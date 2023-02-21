/*

CubicSpline.c   File for cubic spline smoothing

This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by 
the Free Software Foundation.

This file is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public 
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>   

Copyright 2009 
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa

*/

#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <fmt/core.h>

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/cubic_spline_postsmoothing.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class CubicSpline {
  public:
    /*
      Wrapper function for cubic-spline postsmoothing. 
      Assumes input is unsmmothed equivalents (and associated data)
      for random groups design ('R'), single group design ('S'), or CINEG design ('C').
      See Kolen and Brennan (2004, pp. 84-89).  Note in particular that final results
      are the average of cubic spline for x to y and inverse of cubic spline for y to x.

      In a sense, logically, the design type is not needed for cubic spline 
      postsmoothing.  However, knowing the design simplifies obtaining the 
      range [low,high] within which the cubic spline is determined.  This is because
      percentile ranks have already been determined and stored in the USATS or BSTATS
      structure(s). (Recall that RG usus two USTATS structures, SG uses one BSTATS 
      structure, and CG uses two BSTATS structures.) 
      
      Strictly speacking the design parameter could be elminated from the 
      calling sequence since it is avaialble in every PDATA structure.  However,
      requiring the design parameter forces the user to be careful about the 
      functions that are called prior to calling Wrapper_Smooth_CubSpl().  As 
      discussed more fully below, Wrapper_Smooth_CubSpl('R' ...) must be paired with
      two calls to Wrapper_RN(), Wrapper_Smooth_CubSpl('S' ...) must be paired with
      two calls to Wrapper_SN(), and Wrapper_Smooth_CubSpl('C' ...) must be paired 
      with two calls to Wrapper_RN().

      For the 'R' design, two prior calls to ReadRawGet_USTATS() are required, 
      one returning a structure, say, &x, and the other returning a structure,
      say, &y. One call to Wrapper_RN() should use &x followed by &y.  
      The other call should use &y followed by &x.

      For the CINEG design, two calls to ReadRawGet_BSTATS() are required, 
      one returning a structure, say, &xv, and the other returning a structure,
      say, &yv. One call to Wrapper_CN() should use &xv followed by &yv.  
      The other call should use &yv followed by &xv. The method for the 
      x to y and the y to x CINEG designs must be the same. The method
      variable can be E, F, or C; i.e., only one type of equipercentile
      equating is smoothed per call to Wrapper_Smooth_CubSpl(). This
      restriction considerably simplifies the code.  (Complexity
      arises because for frequency estimation, the low and high scores
      for use of the cubic spline need to be determined relative to
      synthestic densities; same is true for modified frequency 
      estimation but the synthetic densities are different.)

      Care must be taken with the 'S' design.  It is assumed here that
      there have been two calls to ReadRawGet_BSTATS(). The first call 
      reads data in order x then y with results stored in, say, &xy.  
      The second call reads data in order y then x, 
      with results stored in a different structure, say, &yx.  Then Wrapper_SN
      needs to be called twice, once using &xy with results stored in 
      structures, say, &pdxy and &rxy, and once using &yx with results stored in 
      structures, say, &pdyx and &ryx.

      NOTE:  As of 6/30/08 this function has not been fully checked for the single
      group design, for the CINEG design, or with the bootstrap

      Input

        design = 'R', 'S', or 'C'
        xtoy = PDATA structure for x to scale of y 
        r_xtoy = ERAW_RESULTS structure for x to scale of y 
        se_xtoy = standard errors for x to scale of y
        cs_xtoy = CS_SMOOTH structure for x to scale of y
        ytox = PDATA structure for y to scale of x 
        r_ytox = ERAW_RESULTS structure for y to scale of x 
        se_ytox = standard errors for y to scale of x
        cs_ytox = CS_SMOOTH structure for y to scale of x
          prlow = percentile rank for low-end truncation  
        prhigh = percentile rank for high-end truncation 
        s = smoothing or "fidelity" parameter
        rep = replication number (set to 0 for actual equating)

      NOTE: cubic spline determined for 
            [lowest score that has a pr >= prlow, highest 
          score that has a pr <= prhigh]
          linear interpolation used outside this range

      Output

        inall = PDATA structure for final cubic-spline equating 
      r = ERAW-RESULTS structure for final cubic-spline equating

      Function calls other than C or NR utilities:
        Smooth_CubSpl()
                                                    
      Author: Robert L. Brennan
      Date of last revision: 3-17-09  
    */
    void runCubicSpline(const EquatingRecipes::Structures::Design& design,
                        EquatingRecipes::Structures::PData& pDataXToY,
                        EquatingRecipes::Structures::EquatedRawScoreResults& equaredRawScoreResultsXToY,
                        Eigen::VectorXd& standardErrorXToY,
                        EquatingRecipes::Structures::CubicSplinePostsmoothing& cubicSplinePostsmoothingXToY,
                        EquatingRecipes::Structures::PData pDataYToX,
                        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResultsYToX,
                        Eigen::VectorXd& standardErrorYToX,
                        EquatingRecipes::Structures::CubicSplinePostsmoothing& cubicSplinePostsmoothingYToX,
                        double& percentileRankLow,
                        double& percentileRankHigh,
                        double& smoothingParameter,
                        size_t& replicationNumber,
                        EquatingRecipes::Structures::PData& pData,
                        EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      std::vector<std::string> methodNames {"     equiv"};
      /* error checking */

      if (design != pDataXToY.design || design != pDataYToX.design) {
        throw std::runtime_error("\nThere is a design mismatch between Wrapper_Smooth_CubSpl() and "
                                 "one or more prior calls to a Wrapper function");
      }

      std::string methodCode = EquatingRecipes::Utilities::getMethodCode(pDataXToY.method);

      if (design == EquatingRecipes::Structures::Design::CommonItenNonEquivalentGroups &&
          methodCode != "E" &&
          methodCode != "F" &&
          methodCode != "C") {
        throw std::runtime_error("\nFor the CINEG design, the method must be F, E, or C.");
      }

      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating */
                                                            /* counting of replications done in Wrapper_Bootstrap() */

      if (pData.bootstrapReplicationNumber == 0) {
        pData.design = design;
        pData.method = pDataXToY.method;
        pData.smoothing = EquatingRecipes::Structures::Smoothing::CUBIC_SPLINE;
        pData.methods = methodNames;

        /* following variables are for x */
        pData.mininumScoreX = pDataXToY.mininumScoreX;
        pData.maximumScoreX = pDataXToY.maximumScoreX;
        pData.scoreIncrementX = pDataXToY.scoreIncrementX;
        pData.scoreFrequenciesX = pDataXToY.scoreFrequenciesX;
        pData.numberOfExaminees = pDataXToY.numberOfExaminees;
      }

      double numberOfScoresX = EquatingRecipes::Utilities::getNumberOfScores(pDataXToY.mininumScoreX,
                                                                             pDataXToY.maximumScoreX,
                                                                             pDataXToY.scoreIncrementX);

      double numberOfScoresY = EquatingRecipes::Utilities::getNumberOfScores(pDataYToX.mininumScoreX,
                                                                             pDataYToX.maximumScoreX,
                                                                             pDataYToX.scoreIncrementX);

      if (pData.bootstrapReplicationNumber <= 1) {
        equatedRawScoreResults.equatedRawScores.resize(1, numberOfScoresX);
        equatedRawScoreResults.equatedRawScoreMoments.resize(1, 4);
      }

      /* populate CS_SMOOTH structures in xtoy and in ytox */
      cubicSplineSmoothing(design,
                           pDataXToY,
                           equaredRawScoreResultsXToY,
                           standardErrorXToY,
                           cubicSplinePostsmoothingXToY,
                           percentileRankLow,
                           percentileRankHigh,
                           smoothingParameter,
                           replicationNumber,
                           numberOfScoresX,
                           numberOfScoresY,
                           false); /* cubic spline for x to y */

      cubicSplineSmoothing(design,
                           pDataYToX,
                           equatedRawScoreResultsYToX,
                           standardErrorYToX,
                           cubicSplinePostsmoothingYToX,
                           percentileRankLow,
                           percentileRankHigh,
                           smoothingParameter,
                           replicationNumber,
                           numberOfScoresY,
                           numberOfScoresX,
                           true); /* inv of cub spl for y to x */

      for (size_t scoreLocation = 0; scoreLocation < numberOfScoresX; scoreLocation++) {
        equatedRawScoreResults.equatedRawScores(0, scoreLocation) = (pDataXToY.cubicSplinePostsmoothing.cubicSplineSmoothedEquivalents(scoreLocation) +
                                                                     pDataYToX.cubicSplinePostsmoothing.cubicSplineSmoothedEquivalents(scoreLocation)) /
                                                                    2.0;
      }

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(
          equatedRawScoreResults.equatedRawScores.row(0),
          pData.scoreFrequenciesX);

      equatedRawScoreResults.equatedRawScoreMoments.row(0) = moments.momentValues;
    }

  private:
    /*
      Populates CS_SMOOTH structure in z (which is either xtoy or ytox).
      If inverse==0, z is xtoy, ns1 = nsx, and ns2 = nsy;
      if inverse==1, x is ytox, ns1 = nsy, ns2 = nsx, and inverse is taken

      Input:

      design = 'R', 'S', or 'C'
      z = PDATA structure for x to y, or for y to x
      r = ERAW_RESULTS structure for x to y, or for y to x
      se = standard errors for raw scores for x to y, or for y to x
      prlow = percentile rank for low-end interpolation
      prhigh = percentile rank for hign-end interpolation
      s = smoothing or "fidelity" parameter
      rep = replication number (0 for actual equating)
        ns1  = number of score categories for first variable
            (e.g. if z = xtoy, x is first variable)
        ns2  = number of score categories for second variable
            (e.g. if z = xtoy, y is second variable)
        inverse = 0 --> no inverse
              = 1 --> get inverse of cubic spline

      Function calls other than C or NR utilities:
        postSmooth()
      inversePostSmooth()
                                                    
      Author: Robert L. Brennan
      Date of last revision: 6/30/08  
    */
    void cubicSplineSmoothing(const EquatingRecipes::Structures::Design& design,
                              EquatingRecipes::Structures::PData& pDataZ,
                              EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults,
                              Eigen::VectorXd& standardErrors,
                              EquatingRecipes::Structures::CubicSplinePostsmoothing& cubicSplinePostsmoothingXToY,
                              double& percentileRankLow,
                              double& percentileRankHigh,
                              double& smoothingParameter,
                              const size_t& replicationNumber,
                              const size_t& numberOfScores1,
                              const size_t& numberOfScores2,
                              const bool& getCubicSplineInverse) {
      size_t numberOfScores;                  /* number of score categories */
      double percentileRankLowestScore;       /* lowest score with pr >= prlow */
      double percentileRankHighestScore;      /* highest score with pr <= prhigh */
      size_t boundedNumberOfScores;           /* number of scores in [low, high] = high-low+1; 'b' --> bounded */
      double scoreIncrement;                  /* score increment    */
      double lowScoresInterpolationConstant;  /* constant for interpolation of low scores (Equation 3.13 in K&B) */
      double highScoresInterpolationConstant; /* constant for interpolation of high scores (Equation 3.13 in K&B) */

      Eigen::VectorXd pseduoRawScores;  /* pointer to vector of pseudo raw scores
      			                      i.e., raw scores in the sense of 0,1,...,ns-1 */
      Eigen::VectorXd pseduoRawScoresX; /* pseudo raw scores always associated with x;
      								                 needed for inversePostSmooth() */
      Eigen::VectorXd cumulativeRelFreqDist;
      Eigen::VectorXd percentileRankDist;
      Eigen::VectorXd smoothingSplineCoefficientMatrix; /* coefficient matrix of the smoothing spline */
      Eigen::VectorXd cubicSplineEquatedScores;         /* vector of equated scores */
      Eigen::VectorXd cubicSplineStandardErrors;        /* vector of standard errors */
      Eigen::VectorXd cubicSplineInverseValues;         /* vector of inverse values  */
      Eigen::VectorXd cubicSplineEquatedEquivalents;    /* vector of equated equivalents smoothed by qubic splne */

      /* get ns as well as low and high scores associated with
           percentile ranks of prlow and prhigh, respectively,
           for various designs (R, S, or G) */

      if (design == EquatingRecipes::Structures::Design::RandomGroups) {
        numberOfScores = pDataZ.summaryRawDataX.numberOfScores;

        size_t scoreLocation;

        for (scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
          if (percentileRankLow <= pDataZ.summaryRawDataX.percentileRankDist(scoreLocation)) {
            break;
          }
        }

        percentileRankLowestScore = static_cast<double>(scoreLocation);

        for (scoreLocation = numberOfScores - 1; scoreLocation >= 0; scoreLocation--) {
          if (percentileRankHigh >= pDataZ.summaryRawDataX.percentileRankDist(scoreLocation)) {
            break;
          }
        }

        percentileRankHighestScore = static_cast<double>(scoreLocation);

      } else if (design == EquatingRecipes::Structures::Design::SingleGroup) {
        numberOfScores = pDataZ.summaryRawDataXY.univariateStatisticsRow.numberOfScores;

        size_t scoreLocation;

        for (scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
          if (percentileRankLow <= pDataZ.summaryRawDataXY.univariateStatisticsRow.percentileRankDist(scoreLocation)) {
            break;
          }
        }

        percentileRankLowestScore = static_cast<double>(scoreLocation);

        for (scoreLocation = numberOfScores - 1; scoreLocation >= 0; scoreLocation--) {
          if (percentileRankHigh >= pDataZ.summaryRawDataXY.univariateStatisticsRow.percentileRankDist(scoreLocation)) {
            break;
          }
        }

        percentileRankHighestScore = static_cast<double>(scoreLocation);

      } else if (design == EquatingRecipes::Structures::Design::CommonItenNonEquivalentGroups) {
        numberOfScores = pDataZ.summaryRawDataXV.univariateStatisticsRow.numberOfScores;

        /* In next section of code, low and high are determined using
      		   actual raw score distribtuion for chained (C), synthetic raw score
      		   distrbution [0] (i.e., fxs[0] or gys[0]) for frequency
      	       estimation (E), and synthetic raw score distribution [1]
      	       (i.e., fxs[1] or gys[1]) for modified frequency estimation (F).
      	       Note that for methods E, F, or C the equivalents are
      	       stored in eraw[0]. See comments in Wrapper_CN() for
      	       further explanation. */

        if (pDataZ.method == EquatingRecipes::Structures::Method::CHAINED) {
          size_t scoreLocation;

          for (scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
            if (percentileRankLow <= pDataZ.summaryRawDataXV.univariateStatisticsRow.percentileRankDist(scoreLocation)) {
              break;
            }
          }

          percentileRankLowestScore = static_cast<double>(scoreLocation);

          for (scoreLocation = numberOfScores - 1; scoreLocation >= 0; scoreLocation--) {
            if (percentileRankHigh >= pDataZ.summaryRawDataXV.univariateStatisticsRow.percentileRankDist(scoreLocation)) {
              break;
            }
          }

          percentileRankHighestScore = static_cast<double>(scoreLocation);
        } else if (pDataZ.method == EquatingRecipes::Structures::Method::FE_BH) {
          cumulativeRelFreqDist = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(0,
                                                                                         numberOfScores - 1,
                                                                                         1,
                                                                                         equatedRawScoreResults.relativeFreqDistsX.row(0));

          percentileRankDist = EquatingRecipes::Utilities::percentileRanks(0,
                                                                           numberOfScores - 1,
                                                                           1,
                                                                           cumulativeRelFreqDist);

          size_t scoreLocation;

          for (scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
            if (percentileRankLow <= percentileRankDist(scoreLocation)) {
              break;
            }
          }

          percentileRankLowestScore = static_cast<double>(scoreLocation);

          for (scoreLocation = numberOfScores - 1; scoreLocation >= 0; scoreLocation++) {
            if (percentileRankHigh >= percentileRankDist(scoreLocation)) {
              break;
            }
          }

          percentileRankHighestScore = static_cast<double>(scoreLocation);
        } else if (pDataZ.method == EquatingRecipes::Structures::Method::MFE_BH) {
          cumulativeRelFreqDist = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(0,
                                                                                         numberOfScores - 1,
                                                                                         1,
                                                                                         equatedRawScoreResults.relativeFreqDistsX.row(1));

          percentileRankDist = EquatingRecipes::Utilities::percentileRanks(0,
                                                                           numberOfScores - 1,
                                                                           1,
                                                                           cumulativeRelFreqDist);

          size_t scoreLocation;

          for (scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
            if (percentileRankLow <= percentileRankDist(scoreLocation)) {
              break;
            }
          }

          percentileRankLowestScore = static_cast<double>(scoreLocation);

          for (scoreLocation = numberOfScores - 1; scoreLocation >= 0; scoreLocation--) {
            if (percentileRankHigh >= percentileRankDist(scoreLocation)) {
              break;
            }
          }

          percentileRankHighestScore = static_cast<double>(scoreLocation);
        }
      } else {
        throw std::runtime_error("Invalid method submitted to cubic spline.");
      }

      boundedNumberOfScores = percentileRankHighestScore - percentileRankLowestScore + 1;

      if (numberOfScores != numberOfScores1) {
        throw std::runtime_error("\nns should be equal to ns1");
      }

      //   inc = z->inc;                                     /* raw score increment*/

      if (pDataZ.bootstrapReplicationNumber <= 1) {
        cubicSplinePostsmoothingXToY.numberOfScores = numberOfScores;
        cubicSplinePostsmoothingXToY.smoothingParameter = smoothingParameter;
        cubicSplinePostsmoothingXToY.percentileRankLowestScore = percentileRankLow;
        cubicSplinePostsmoothingXToY.percentileRankHighestScore = percentileRankHigh;
        cubicSplinePostsmoothingXToY.lowestSmoothedPseudoRawScorePercentileRank = percentileRankLowestScore;
        cubicSplinePostsmoothingXToY.higestSmoothedPseudoRawScorePercentileRank = percentileRankHighestScore;
        cubicSplinePostsmoothingXToY.boundedNumberOfScores = boundedNumberOfScores;
        cubicSplinePostsmoothingXToY.equipercentileEquivalents = equatedRawScoreResults.equatedRawScores.row(0); /* raw-score equipercentile equivalents */
        cubicSplinePostsmoothingXToY.standardErrors = standardErrors;                                            /* standard errors of r->eraw[0] */
        cubicSplinePostsmoothingXToY.coefficients.resize(4 * numberOfScores);                                    /* coeffs; dimensioned at max */
        cubicSplinePostsmoothingXToY.equipercentileEquivalents.resize(numberOfScores);                           /* cub spl results */

        /* inverse of cub spl y to x */
        if (getCubicSplineInverse) {
          cubicSplinePostsmoothingXToY.inverseCubicSplineSmoothedEquivalents.resize(numberOfScores2);
        }

        pDataZ.cubicSplinePostsmoothing = cubicSplinePostsmoothingXToY;
      }

      /* scale input data so that it starts at 0 with unit increment*/
      cubicSplineEquatedScores.resize(numberOfScores);
      cubicSplineStandardErrors.resize(numberOfScores);

      if (getCubicSplineInverse) {
        cubicSplineInverseValues.resize(numberOfScores2);
      }

      cubicSplineEquatedEquivalents.resize(numberOfScores2);
      smoothingSplineCoefficientMatrix.resize(4 * boundedNumberOfScores);

      for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
        cubicSplineEquatedScores(scoreLocation) = ((equatedRawScoreResults.equatedRawScores(0, scoreLocation) / scoreIncrement) -
                                                   pDataZ.mininumScoreX / scoreIncrement);

        cubicSplineStandardErrors(scoreLocation) = standardErrors(scoreLocation) / scoreIncrement;
      }

      Eigen::VectorXd rawScores(numberOfScores);
      for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
        rawScores(scoreLocation) = static_cast<double>(scoreLocation);
      }

      if (!getCubicSplineInverse) {
        postSmooth(rawScores,
                   cubicSplineEquatedScores,
                   cubicSplineStandardErrors,
                   numberOfScores,
                   smoothingParameter,
                   percentileRankLowestScore,
                   percentileRankHighestScore,
                   static_cast<double>(numberOfScores2 - 1),
                   rawScores,
                   numberOfScores,
                   cubicSplineEquatedEquivalents,
                   smoothingSplineCoefficientMatrix);
      } else {
        Eigen::VectorXd rawScoresX(numberOfScores2);

        for (size_t scoreLocation = 0; scoreLocation < numberOfScores2; scoreLocation++) {
          rawScoresX(scoreLocation) = static_cast<double>(scoreLocation);
        }

        inversePostSmooth(rawScoresX,
                          cubicSplineEquatedScores,
                          cubicSplineStandardErrors,
                          numberOfScores,
                          smoothingParameter,
                          percentileRankLowestScore,
                          percentileRankHighestScore,
                          static_cast<double>(numberOfScores2 - 1),
                          rawScoresX,
                          numberOfScores2,
                          cubicSplineInverseValues,
                          smoothingSplineCoefficientMatrix);

        /* The following code is to get the eeqs[] vector; i.e., the
      	   cubic-spline smoothed equivalents for putting y on scale of x.
      	   This code is necessary because inversePostSmooth() does not
      	   return eeqs[], but it can be obtained from cmat[] plus
      	   linear interpolation at the ends */

        for (size_t scoreLocation = 0; scoreLocation < boundedNumberOfScores; scoreIncrement++) {
          cubicSplineEquatedEquivalents(static_cast<size_t>(percentileRankLowestScore) + scoreLocation) =
              smoothingSplineCoefficientMatrix(scoreLocation);
        }

        cubicSplineEquatedEquivalents(static_cast<size_t>(percentileRankHighestScore)) =
            smoothingSplineCoefficientMatrix(boundedNumberOfScores - 2) +
            smoothingSplineCoefficientMatrix(2 * boundedNumberOfScores - 3) +
            smoothingSplineCoefficientMatrix(3 * boundedNumberOfScores - 4) +
            smoothingSplineCoefficientMatrix(4 * boundedNumberOfScores - 5);

        /* linear interpolation for low scores */
        lowScoresInterpolationConstant = cubicSplineEquatedEquivalents(static_cast<size_t>(percentileRankLowestScore) + 0.5) /
                                         cubicSplineEquatedEquivalents(static_cast<size_t>(percentileRankHighestScore) + 0.5);

        for (size_t scoreLocation = 0; scoreLocation < static_cast<size_t>(percentileRankLowestScore); scoreLocation++) {
          cubicSplineEquatedEquivalents(scoreLocation) = lowScoresInterpolationConstant *
                                                             (static_cast<double>(scoreLocation) + 0.5) -
                                                         0.5;
        }

        //     /* linear interpolation for high scores */
        highScoresInterpolationConstant = (cubicSplineEquatedEquivalents(static_cast<size_t>(percentileRankHighestScore)) -
                                           static_cast<double>(numberOfScores2 - 1) + 0.5) /
                                          (percentileRankHighestScore - static_cast<double>(numberOfScores - 1) + 0.5);

        for (size_t scoreLocation = static_cast<size_t>(percentileRankHighestScore) + 1; scoreLocation < numberOfScores2; scoreLocation++) {
          cubicSplineEquatedEquivalents(scoreLocation) = highScoresInterpolationConstant *
                                                             (static_cast<double>(scoreLocation) - percentileRankHighestScore) +
                                                         cubicSplineEquatedEquivalents(static_cast<size_t>(percentileRankHighestScore));
        }
      }

      /* scale input data back to the original scale*/
      for (size_t scoreLocation = 0; scoreLocation < numberOfScores2; scoreLocation++) {
        cubicSplinePostsmoothingXToY.equipercentileEquivalents(scoreLocation) = cubicSplineEquatedEquivalents(scoreLocation) *
                                                                                    scoreIncrement +
                                                                                pDataZ.mininumScoreX;
      }

      for (size_t scoreLocation = 0; scoreLocation < boundedNumberOfScores; scoreLocation++) {
        cubicSplinePostsmoothingXToY.coefficients(scoreLocation) = smoothingSplineCoefficientMatrix(scoreLocation) * scoreIncrement + pDataZ.mininumScoreX;

        cubicSplinePostsmoothingXToY.coefficients(boundedNumberOfScores + scoreLocation) = smoothingSplineCoefficientMatrix(boundedNumberOfScores + scoreLocation);

        cubicSplinePostsmoothingXToY.coefficients(2 * boundedNumberOfScores + scoreLocation) = smoothingSplineCoefficientMatrix(2 * boundedNumberOfScores + scoreLocation) / scoreIncrement;

        cubicSplinePostsmoothingXToY.coefficients(3 * boundedNumberOfScores + scoreLocation) = smoothingSplineCoefficientMatrix(3 * boundedNumberOfScores + scoreLocation) / std::pow(scoreIncrement, 2);
      }

      if (getCubicSplineInverse) {
        for (size_t scoreLocation = 0; scoreLocation < numberOfScores2; scoreLocation++) {
          cubicSplinePostsmoothingXToY.inverseCubicSplineSmoothedEquivalents(scoreLocation) = cubicSplineInverseValues(scoreLocation) * scoreIncrement + pDataZ.mininumScoreX;
        }
      }
    }

    /* functions related to numerical linear algebra */
    Eigen::MatrixXd vectorToMatrix(const Eigen::VectorXd& vect,
                                   const size_t& numberOfRows,
                                   const size_t& numberOfColumns,
                                   const bool& byRow) {
      Eigen::MatrixXd mat(numberOfRows, numberOfColumns);

      size_t rowIndex = 0;
      size_t columnIndex = 0;

      for (size_t vectorIndex = 0; vectorIndex < vect.size(); vectorIndex++) {
        mat(rowIndex, columnIndex) = vect(vectorIndex);

        if (byRow) {
          columnIndex++;

          if (columnIndex == numberOfColumns) {
            rowIndex++;
            columnIndex = 0;
          }
        } else {
          rowIndex++;

          if (rowIndex == numberOfRows) {
            rowIndex = 0;
            columnIndex++;
          }
        }
      }

      return mat;
    }

    Eigen::MatrixXd vectorToSymmetricMatrix(const Eigen::VectorXd& vect,
                                            const size_t& numberOfRowsOrColumns,
                                            const bool& byRow,
                                            const bool& isLowerTriangle,
                                            const bool& includeDiagonal) {
      Eigen::MatrixXd mat(numberOfRowsOrColumns, numberOfRowsOrColumns);

      size_t rowIndex = 0;
      size_t columnIndex = 0;

      for (size_t vectorIndex = 0; vectorIndex < vect.size(); vectorIndex++) {
        mat(rowIndex, columnIndex) = vect(vectorIndex);

        if ()

          if (columnIndex)

            if ((columnIndex == rowIndex && !includeDiagonal) ||
                (columnIndex > rowIndex)) {
              if (byRow) {
              }
            }
      }
    }

    // void	dpdch(double * mx,int dim);
    // void	dtdmt(double *tdm2,double *tdm1, int rdim1, int cdim1);
    // void	dtdsm(double *tdm2,double c, double *tdm1, int rdim1, int cdim1);
    // void	dtdmm(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1, int cdim2);
    // void	dtdma(double *tdm3, double * tdm1, double * tdm2, int rdim1, int cdim1);
    // void	subbk(double* b, double *utm, int dim);
    // void	subfd(double* b, double *utm, int dim);
    // void	chsol(double *vx, double* vb, double *mpd, int dim);
    // void	dtdmv(double *v2, double * mtd, double * v1, int rdim, int cdim);

    /* functions related to numerical linear/cubic polynomial */

    /* 
   Purpose:                                                            
      Linear polynomial function f passing through (x0,y0) and (x1,y1) 
                                                                       
       f(x0,y0,x1,y1,xvalue)                                           
                y1 - y0                                                 
       = y0 + -----------*(xvalue-x0)                                  
                x1 - x0    

   Input:                                                               
      x0,y0,    x & y coordinate of a point p0                         
      x1,y1,    x & y coordinate of a point p1                         
      xvalue,   point at which to evaluate f() 

   Output:
      f(x0,y0,x1,y1,xvalue)

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 7/4/08
*/
    double linearPoly(const double& x0,
                      const double& y0,
                      const double& x1,
                      const double& y1,
                      const double& xvalue) {
      /* evaluate the function         */
      double slope = (y1 - y0) / (x1 - x0);
      /* check that slope != 0         */
      if (std::abs(slope) < std::pow(10.0, 3.0) * std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error("linearPoly: Input Warning, slope ~ 0 \n");
      }

      double rvalue = y0 + slope * (xvalue - x0);

      return rvalue;
    }

    /* 
   Purpose:                                                            
      Cubic polynomial function f defined by ai,bi,ci,di, i.e.,		   
       f(ai,bi,ci,di,xvalue)                                          
       =ai+bi*(xvalue-xleft)+ci*(xvalue-xleft)^2+di*(xvalue-xleft)^3   
      defined over [xleft, xright]      

   Input:                                                              
      xleft,       left end of the domain for f()                      
      xright,      right end of the domain for f()                     
      ai,bi,ci,di, coefficients of f()                                 
      xvalue,      point at which to evaluate f() 

   Output:
      f(ai,bi,ci,di,xvalue)

   Function calls other than C or NR utilities:
      None.

  Author: Jaehoon Seol
  Date of last revision: 7/3/08
*/
    double cubicPoly(const double& xleft,
                     const double& xright,
                     const double& ai,
                     const double& bi,
                     const double& ci,
                     const double& di,
                     const double& xvalue) {
      /* confirm xleft<=xvalue<=xright */
      if (xvalue < xleft || xvalue > xright) {
        std::string msg = "cubicPoly: Input Error, cubic poly. used outside domain\n";
        msg.append(fmt::format("xleft = {:f}, xright= {:f} xvalue = {:f}\n", xleft, xright, xvalue));
        throw std::runtime_error(msg);
      }

      /* evaluate the function         */
      double tempX = xvalue - xleft;
      double rvalue = ai + bi * tempX + ci * std::pow(tempX, 2.0) + di * std::pow(tempX, 3.0);

      return rvalue;
    }

    /* functions related to post-smoothing method using cubic splines */

    /*
   Purpose: 
     This function sets up a positive definite, (n-1)x(n-1) tridiagonal matrix  
	 defined by

	    t     = 2*(h     + h    )/3,         (1)
         i,i        i-1     i
	    t      = t        = h   /3           (2)
		 i,i+1    i+1,i      i
	 Refer to page 179, Reinsch (1967).

   Input  :
     h =  [  h ,h , ... h     ]
              0  1       n-1       
     
	      where  h    = x    -  x   , i=0, 1,2,...,n-1
	              i      i+1     i
 
     n     = number of elements in the vector h.
	 edist = true,      if h is equi-distance vector.
	       = false,     otherwise

   Output :
      mt, The matrix defined by (1) and (2)

   Function calls other than C or NR utilities:
      None.

   References :
     1. C.H. Reinsch, Smoothing by spline functions, 1967, Num. math. 10, p177-183. 
   Comments:
     -The size of memory for mt should be n*sizeof(double) and
	  the memory should be allocated before calling setup_matrixT.

  Author: Jaehoon Seol
  Date of last revision: 2/18/08

*/
    void setupT(const bool& edist,
                Eigen::MatrixXd& mt,
                size_t& n,
                Eigen::VectorXd& h) {
      mt.setZero(n - 1, n - 1);

      double h0 = h(0); /* will be used if h is equi-distance */

      for (size_t index = 0; index < mt.rows(); index++) {
        if (edist) {
          if (index >= 1) {
            mt(index - 1, index) = h0 / 3.0;  
            mt(index, index - 1) = h0 / 3.0;
          }

          mt(index, index) = 4.0 * h0 / 3.0;
        } else {
          if (index >= 1) {
            mt(index - 1, index) = h(index) / 3.0;
            mt(index, index - 1) = h(index) / 3.0;
          }
          
          if (index < mt.rows() - 1) {
            double x = h(index) + h(index + 1);
            mt(index, index) = x * 2.0 / 3.0;
          }
        }


        

        for (size_t columnIndex = left; columnIndex <= right; columnIndex++) {
          if (edist) {
            if (rowIndex == columnIndex) {
              mt((n - 1) * rowIndex + columnIndex) = 4.0 * h0 / 3.0;
            } else {
              mt((n - 1) * rowIndex + columnIndex) = h0 / 3.0;
            }
          } else {
            if (rowIndex == columnIndex) {
              double x = h(rowIndex) + h(rowIndex + 1);
              mt((n - 1) * rowIndex + columnIndex) = x * 2.0 / 3.0;
            } else if (columnIndex == rowIndex + 1) {
              mt((n - 1) * rowIndex + columnIndex) = h(rowIndex + 1) / 3.0;
            } else if (columnIndex == rowIndex - 1) {
              mt((n - 1) * rowIndex + columnIndex) = h(rowIndex) / 3.0;
            }
          }
        }
      }
    }

    /*
   Purpose: 
     This function sets up a positive definite, (n+1)x(n-1) tridiagonal matrix  
	 defined by

	    q     = -1/h     -1/h                 (1)
         i,i        i-1      i
	    q        = 1/h                        (2)
		  i+1,i       i
	    q        = 1/h                        (3)
		 i-1,i        i-1
	  Refer to page 179, Reinsch (1967).

   Input  :
     h =  [  h ,h , ... h     ]
              0  1       n-1       
     
	      where  h    = x    -  x   , i=0, 1,2,...,n-1
	              i      i+1     i
 
     n     = number of elements in the vector h.
	 edist = true,      if h is equi-distance vector.
	       = false,     otherwise

   Output :
      mt, The matrix defined by (1) and (2)

   Function calls other than C or NR utilities:
      None.

   References :
     1. C.H. Reinsch, Smoothing by spline functions, 1967, Num. math. 10, p177-183. 
   Comments:
     - The size of memory for mt should be n*sizeof(double) and
	   the memory should be allocated before calling setup_matrixT.

  Author: Jaehoon Seol
  Date of last revision: 2/18/08
*/
    void setupQt(const bool& edist,
                 Eigen::MatrixXd& mqt,
                 const size_t n,
                 Eigen::MatrixXd& h) {
      mqt.setZero(n + 1, n - 1);

      double h0 = h(0);

      size_t maxIndex = (mqt.rows() <= mqt.cols()) ? mqt.rows() : mqt.cols();

      for (size_t index = 0; index < maxIndex; index++) {
        if (edist) {
          h0 = h(0, 0);

          if (index >= 1) {
            mqt(index - 1, index) = 1.0 / h0;
          }

          mqt(index, index) = -2.0 / h0;

          mqt(index + 1, index) = 1.0 / h0;
        } else {
          if (index >= 1) {
            mqt(index - 1, index) = 1.0 / h(index, 0);
          }

          mqt(index, index) = -1.0 * (1.0 / h(index - 1, 0) + 1.0 / h(index, 0));

          mqt(index + 1, index) = 1.0 / h(index, 0);
        }
      }
    }

    /*
   Purpose:  
     This function calculates the coefficients

	     ai, bi, ci, di, i=0,1,...n 

	 for the cubic spline 
	                                  2        3
	     f(x) = ai + bi(x-xi)+ci(x-xi)+di(x-xi),   x in [xi,xi+1)

	 The cubic spline f(x) minimizes
                         2
	          Int[ f''(x)  ]dx

	 among all functions f(x) satisfying
	                                   
                        f(xi)-yi      2              2
	          Sum( ----------------- )  <= s,  f in C [x0, xn] ---(1)
			               dyi

   Description of the algorithm:     
	 * start with p=0.0                                   * 
	 *                             T       T  2           *  
	 * (1) cholesky decomposition R  R of Q  D  Q + pT    *  
	 * (2) compute u from R^t Ru = Q^t y and 
	 *             v = D Q u, e=v^t v                     * 
	 * (3) if e is greater than S                         *
	 *      compute f = u^t T u and g = w^t w             *
	 *              where R^T w = T u                     *
	 *      replace p by p+(e-(Se)^(.5))/(f-p*g)          *
	 *      restart with step (1).                        *
	 *     Otherwise                     
	 * (4) compute a = y-D v, c= pu                       *
	 *     compute b and d according to (8) and (9)       * 
   
   Input  : 
     x,   describes x-coordinates of input data.
	      (n+1)x1 matrix, i.e.,
		     x = {x0,x1,x2,...,xn}'
     y,   describes y-coordinates of input data.
	      (n+1)x1 matrix, i.e.,
		     y = {y0,y1,y2,...,yn}'
     dyi, estimates the standard deviation of the ordinate yi
	      (n+1)x1 matrix
     num, number of input data 
	 s,   fidelity constant which controls the closeness of f(xi) to yi in (1)
          This constant is scaled by (x[n]-x[0]+1) within the function to 
		  give "s" more consistent meaning across different applications.

   Precondition:  
     mc is a one-dimensional array of size 4*n*sizeof(double).
	 Enough memory for mc should be allocated before calling 
	 ssspline

   Output :
     mc,   coefficient matrix of the smoothing cubic spline returns 
	 number of iteration used to compute p. If the return value is equal
	 to -1, that means the max. iteration number(25) has been used.

   Note :
	 Additional break statement "if (e<=s) break;" at the end of the loop 
	 has been added to keep the matrix
	                  T  2
	                 Q  D  Q + pT  
	 positive definite. Withough this condition, the matrix may not be 
	 positive definite in some cases in which case Cholesky decomposition
	 can not be used. If the user is not comfortable with this, he/she should
	 consider using QR factorization rather than Cholesky decomposition.

   References :
     1. C.H. Reinsch, Smoothing by spline functions, 1967, Num. Math. 10, p177-183. 
	 2. C.H. Reinsch, Smoothing by spline functions. II, 1970, Num. Math. 16, p451-454 

  Author: Jaehoon Seol
  Date of last revision: 8/10/08
*/
    int sspline(const Eigen::VectorXd& x,
                const Eigen::VectorXd& y,
                const Eigen::VectorXd& dyi,
                const size_t& num,
                double& smoothingParameter,
                Eigen::VectorXd& mc) {
      // size_t i;             /* loop index */
      // size_t j;             /* loop index */
      // size_t ii = 0;        /* main loop index */
      // double e = 0.0;       /* e=v^t v */
      // double np = 0.0;      /* new, updated p */
      // double p = 1.0;       /* p<-p+(e-(Se)^(.5))/(f-p*g) */
      // double g = 0.0;       /* g = w^t w */
      // double f = 0.0;       /* f = u^t T u */

      size_t n = num - 1; /* col. dim of cmat */
      double s = smoothingParameter * (x[n] - x[0] + 1);

      Eigen::MatrixXd va(n + 1, 1);       /* cubic spline coeff. a */
      Eigen::MatrixXd vc(n + 1, 1);       /* cubic spline coeff. c */
      Eigen::MatrixXd vu(n - 1, 1);       /* u: R^t R u = Q^t y */
      Eigen::MatrixXd vw(n - 1, 1);       /* w: R^t R w = T u */
      Eigen::MatrixXd vtw(n - 1, 1);      /* tw: T w */
      Eigen::MatrixXd vh(n, 1);           /* h_i = x_(i+1)-x_i */
      Eigen::MatrixXd vv(n + 1, 1);       /* v= D Q u */
      Eigen::MatrixXd vy2(n + 1, 1);      /* y2: Q^t y */
      Eigen::MatrixXd vtu(n - 1, 1);      /* temporary of size n-1 */
      Eigen::MatrixXd mq(n - 1, n + 1);   /* matrix Q */
      Eigen::MatrixXd mt(n - 1, n - 1);   /* matrix Q^t */
      Eigen::MatrixXd mpt(n - 1, n - 1);  /* mpt = p T */
      Eigen::MatrixXd mqt(n - 1, n + 1);  /* matrix T */
      Eigen::MatrixXd mpt1(n - 1, n - 1); /* part 1:(fixed) Q^t D^2 Q */
      Eigen::MatrixXd mpt2(n - 1, n - 1); /* part 2: Q^t D^2 Q + p T */
      Eigen::MatrixXd mtmp(n - 1, n + 1); /* temp matrix for Q */

      // Eigen::VectorXd vh(n);      /* h_i = x_(i+1)-x_i */
      // Eigen::VectorXd va(n + 1);  // va;   /* cubic spline coeff. a */ (double*)malloc((n + 1) * sizeof(double));
      // Eigen::VectorXd vc(n + 1);  // vc;   /* cubic spline coeff. c */ (double*)malloc((n + 1) * sizeof(double));
      // Eigen::VectorXd vu(n - 1);  // vu;   /* u: R^t R u = Q^t y */ (double*)malloc((n - 1) * sizeof(double));
      // Eigen::VectorXd vw(n - 1);  // vw;   /* w: R^t R w = T u */ (double*)malloc((n - 1) * sizeof(double));
      // Eigen::VectorXd vtw(n - 1); // vtw;  /* tw: T w */ (double*)malloc((n - 1) * sizeof(double));
      // Eigen::VectorXd vv(n);      // vv;   /* v= D Q u */ (double*)malloc((n + 1) * sizeof(double));
      // Eigen::VectorXd vy2(n);     // vy2;  /* y2: Q^t y */ (double*)malloc(n * sizeof(double));
      // Eigen::VectorXd vtu(n);     // vtu;  /* temporary of size n-1 */ (double*)malloc((n + 1) * sizeof(double));
      // Eigen::VectorXd mq(n);      // mq;   /* matrix Q */ (double*)malloc((n - 1) * sizeof(double));
      // Eigen::VectorXd mt(n);      // mt;   /* matrix T */ (double*)malloc((n - 1) * (n - 1) * sizeof(double));
      // Eigen::VectorXd mpt(n);     // mpt;  /* mpt = p T */ (double*)malloc((n - 1) * (n - 1) * sizeof(double));
      // Eigen::VectorXd mqt(n);     // mqt;  /* matrix Q^t */ (double*)malloc((n - 1) * (n + 1) * sizeof(double));
      // Eigen::VectorXd mpt1(n);    // mpt1; /* part 1:(fixed) Q^t D^2 Q */ (double*)malloc((n - 1) * (n + 1) * sizeof(double));
      // Eigen::VectorXd mpt2(n);    // mpt2; /* part 2: Q^t D^2 Q + p T */ (double*)malloc((n - 1) * (n - 1) * sizeof(double));
      // Eigen::VectorXd mtmp(n);    // mtmp; /* temp matrix for Q */ (double*)malloc((n - 1) * (n - 1) * sizeof(double));
      //                             //  (double*)malloc((n - 1) * (n + 1) * sizeof(double));

      /* setup vectors and matrices that does not change thru. loop	*/
      vh.col(0) = x(Eigen::seq(1, n)) - x(Eigen::seq(0, n - 1));

      setupQt(false,
              mqt,
              n,
              vh); /* setup matrix Q^t  */

      mq = mqt.transpose(); /* setup matrix Q */

      setupT(FALSE, mt, n, vh); /* setup matrix T */
      memset(mtmp, 0, (n - 1) * (n + 1) * sizeof(double));
      for (i = 0; i < n + 1; i++) /* compute  D^2 Q */
        for (j = 0; j < n - 1; j++)
          mtmp[(i) * (n - 1) + (j)] = dyi[i] * dyi[i] * mq[(i) * (n - 1) + (j)];
      dtdmm(mpt1, mqt, mtmp, n - 1, n + 1, n - 1); /* mpt1 = Q^t D^2 Q */
      free(mtmp);
      dtdmv(vy2, mqt, y, n - 1, n + 1); /* compute y2 = Q^t y */
      ii = 0;
      while (fabs(np - p) > pow(10.0, 3.0) * DBL_EPSILON && ii++ < 35) {
        /* cholesky decomposition R^t R	of	*
		 *      (mpt2=) Q^t D^2 Q + p*T		*/
        p = np;
        dtdsm(mpt, p, mt, n - 1, n - 1);
        dtdma(mpt2, mpt1, mpt, n - 1, n - 1);
        chsol(vu, vy2, mpt2, n - 1);     /* R^t R u = Q^t y	 */
        dtdmv(vv, mq, vu, n + 1, n - 1); /* compute v=D Q u	 */
        for (i = 0; i < n + 1; i++)
          vv[i] = dyi[i] * vv[i];
        e = 0.0;
        for (i = 0; i < n + 1; i++)
          e += vv[i] * vv[i];
        dtdmv(vtu, mt, vu, n - 1, n - 1);
        /* compute w in R^t R w = T u */
        memcpy(vw, vtu, (n - 1) * sizeof(double));
        subfd(vw, mpt2, n - 1); /* vw has w  */
        subbk(vw, mpt2, n - 1); /* vw has w  */
        dtdmv(vtw, mt, vw, n - 1, n - 1);
        /* compute f, g, and p		  */
        f = 0.0;
        for (i = 0; i < n - 1; i++)
          f += vu[i] * vtu[i];
        g = 0.0;
        for (i = 0; i < n - 1; i++)
          g += vu[i] * vtw[i];
        if (s != 0.0)
          np = p + sqrt(e / s) * (e - sqrt(s * e)) / (f - p * g);
        else
          np = p + (e - sqrt(s * e)) / (f - p * g);
        if (e <= s)
          break;
      }
      /* printf("p-value = %f ii = %d\n",np, ii); */
      /* return matrix mc:						* 
	  * mc(0, ) : coefficients a				*
	  * mc(1, ) : coefficients b				*
	  * mc(2, ) : coefficients c				* 
	  * mc(3, ) : coefficients d				*/
      va[0] = y[0] - dyi[0] * vv[0];
      vc[0] = 0.0;
      /* return coefficients b and d				*/
      for (i = 0; i < n; i++) {
        va[i + 1] = y[i + 1] - dyi[i + 1] * vv[i + 1];
        vc[i + 1] = (i == n - 1 ? 0 : np * vu[i]);
        mc[(3) * (n) + (i)] = (vc[i + 1] - vc[i]) / (3.0 * vh[i]);
        mc[(1) * (n) + (i)] = (va[i + 1] - va[i]) / vh[i] - (vc[i] + mc[(3) * (n) + (i)] * vh[i]) * vh[i];
      }
      /* return remaining coefficients			*/
      memcpy(mc, va, n * sizeof(double));
      memcpy(mc + 2 * n, vc, n * sizeof(double));
      /* free all working areas */
      free(va);
      free(vc);
      free(vh);
      free(vu);
      free(vv);
      free(vw);
      free(vtw);
      free(vtu);
      free(vy2);
      free(mq);
      free(mt);
      free(mqt);
      free(mpt);
      free(mpt1);
      free(mpt2);
      if (ii == 25)
        ii = -1;
      return ii;
    }

    void postSmooth(const Eigen::VectorXd& rawScores,
                    const Eigen::VectorXd& cubicSplineRawScores,
                    const Eigen::VectorXd& cubicSplineStandardErrors,
                    const size_t& numberOfScores,
                    const double& smoothingParameter,
                    const double& percentileRankLowestScore,
                    const double& percentileRankHighestScore,
                    const double& ky,
                    const Eigen::VectorXd& vectX,
                    const int& numb2,
                    Eigen::VectorXd& vectY,
                    Eigen::VectorXd& smoothingSplineCoefficientMatrix) {
    }

    /* functions related to the inverse of the post-smoothing method using cubic splines	*/

    /* Purpose:                                                          
      Use bisection method to compute the xvalue satisfying the       
      following equation                                             
                                                                     
            yvalue = f(ai,bi,ci,di,xvalue)  ...... (FI.1)            
                                                                     
        where f(ai,bi,ci,di,*) is a cubic polynomial defined by         
        ai,bi,ci,di.     

      Input:                                                            
        left,        left end of the domain where f( ) is defined       
        right,       right end of the domain where f( ) is defined      
        ai,bi,ci,di, coefficients of the cubic polynomial               
        yvalue,      yvalue at which f^-1 is evaluated

      Output:                                                           
        xvalue that satisfies Equation (FI.1)    

      Precondition:                                                     
        f(left)<yvalue<f(right)                                         
        If this condition is not met, the inverse can not be computed
    
      Function calls other than C or NR utilities: 
        cubicPoly() 

      Author: Jaehoon Seol
      Date of last revision: 7/3/08
    */
    double inverseCubicPoly(double left, double right, double ai, double bi,
                            double ci, double di, double yvalue) {
      return 0;
    }

    /*
    Purpose:
        compute the inverse of the cubic spline defined by ynodes[] and  
        cmat[]. The cubic spline is defined by cmat[i],cmat[i+n2],       
        cmat[i+2*n2],cmat[i+3*n2] on [ynodes[i],ynodes[i+1]].   

    Input:
        ynodes,    nodes of the cubic spline                             
        cmat,      coefficient matrix of the cubic spline                
        n2,        n2x4 is the dimension of the cmat                     
        vectX,     vectX = dx(vectY)                                     
        n1,        dimension of vector vectX and vectY  

    Output:                                                              
        vectY,    inverse value computed 
  
    Function calls other than C or NR utilities: 
      cubicPoly()
      inverseCubicPoly()

    Author: Jaehoon Seol
    Date of last revision: 7/3/08
  */
    void inverseSSpline(double* ynodes, double* cmat, int n2, double* vectX, int n1, double* vectY) {}

    /* Purpose/functionality
      Implementation of the inverse post-smoothing method dx  (x) defined piecewise
      from -0.5 to Ky+0.5 (see Kolen & Brennan, 2004, pp. 84-89).
        Inverts the following function
          
                  dx(ylow)+0.5             0.5*(dx(ylow)+0.5)
                [----------------]*y+[-.5 + -------------------], 
                    ylow + 0.5                 ylow + 0.5
                                            for -0.5 <= y <  ylow    (Eq. 1)
        dx(y) = sspline(.............)       for ylow <= y <= yhigh   (Eq. 2)
                  dx(yhigh)-(Kx+0.5)                    yhigh*[dx(yhigh)-(Kx+0.5)]
                [---------------------]*y + [dx(yhigh) - --------------------------- ] 
                    yhigh - (Ky+0.5)                      yhigh - (Ky+0.5)
                                            for yhigh < y <= Ky+0.5  (Eq. 3)
      
      Input:
        yvalues,	raw Y scores, usually 0,1,2,...,(num1-1)
        xvalues,	equated raw scores on scale of X; i.e., e_X(y)
        dxi,		std. error of equated raw scores; i.e., SE[e_X(y)]
        num1,		array size for yvalues, xvalues, and dxi
        s,			fidelity constant (called smoothing constant in Kolen & Brennan, 2004)
      ylow,		index associated with the lowest percentile rank (prlow)
                  for obtaining cubic spline results (usually, prlow=0.5).
        yhigh,		index associated with the highest percentile rank (prhigh)
                  for obtaining cubic spline results (usually, prhigh=99.5).
        kx,		number of possible score categories minus 1 -- associated
                  with Form X (Kx in Kolen & Brennan, 2004)
        vectX,		x values where d_X^-1(x), the inverse of d_X(y),
                    is evaluated
      num2,		array size for vectX and vectY

      Output:
        vectY,     evaluation of d_x^-1(x) at vectX 
      cmat,      coefficient matrix of the smoothing cubic spline that has dim
                  (yhigh-ylow)x4

      NOTES:  (a) The use of this function in Equating Recipes assumes 
                  the score categories for Form Y are 0,1,2,...,num1-1, and 
              the score categories for Form X are 0,1,2,..., int(kx).
              Transformation to actual raw scores occurs outside this function.
          (b) For equating purposes set num2 = kx+1.  However, if the user
              wanted to get "between-the-nodes" cubic spline equivalents (e.g., to
            create a nearly smooth plot) then num2 might be a much larger
            number
    
      Function calls other than C or NR utilities: 
          sspline()
        cubicPoly()
        linearPoly()
        inverseSSpline()

      Comments:
      - Allocated memory of vectY and cmat should be passed to postSmooth() 
        before calling this function.
      - 0.5 in the code means half a score category 

      Author: Jaehoon Seol
      Date of last revision: 7/8/08 
    */
    void inversePostSmooth(const Eigen::VectorXd& yvalues,
                           const Eigen::VectorXd& xvalues,
                           const Eigen::VectorXd& dxi,
                           const size_t& num1,
                           const double& smoothingParameter,
                           const double& ylow,
                           const double& yhigh,
                           const double& kx,
                           const Eigen::VectorXd& vectX,
                           const size_t& num2,
                           const Eigen::VectorXd& vectY,
                           const Eigen::VectorXd& smoothingSplineCoefficientMatrix) {
      // int i,          /* loop index */
      //     numcoeff;   /* number of cubic spline pieces */
      // double dx_ylow, /* dx(ylow)  */
      //     dx_yhigh,   /* dx(yhigh) */
      //     ky,         /* maximum y value (Ky in Kolen & Brennan, 2004)*/
      //     temp;       /* temp. var. to store vectX[i] */

      /* compute coefficient matrix of the cubic spline */

      // numcoeff = yhigh - ylow;
      /*
	printf("========== inverse post smoothing ===========\n");
	printf("xvalues[xlow]  = %22.15f\n",xvalues[ylow]);
	printf("yvalues[xlow]  = %22.15f\n",yvalues[ylow]);
	printf("dyi[xlow]      = %22.15f\n",dxi[ylow]);
	printf("numcoeff+1     = %d\n",numcoeff+1);
	printf("s              = %22.15f\n",s);
    */
      // sspline(&yvalues[ylow], &xvalues[ylow], &dxi[ylow], numcoeff + 1, s, cmat);
      // dx_ylow = cubicPoly(yvalues[ylow], yvalues[ylow + 1], cmat[0], cmat[numcoeff],
      //                     cmat[2 * numcoeff], cmat[3 * numcoeff], yvalues[ylow]);
      // dx_yhigh = cubicPoly(yvalues[yhigh - 1], yvalues[yhigh], cmat[numcoeff - 1], cmat[2 * numcoeff - 1],
      //                      cmat[3 * numcoeff - 1], cmat[4 * numcoeff - 1], yvalues[yhigh]);
      // ky = yvalues[num1 - 1];

      // for (i = 0; i < num2; i++) {
      //   temp = vectX[i];
      //   if (-0.5 <= temp && temp < dx_ylow) /* inverse of Eq. 1 */
      //     vectY[i] = linearPoly(-0.5, -0.5, dx_ylow, yvalues[ylow], temp);
      //   else if (dx_yhigh < temp && temp <= kx + 0.5) /* inverse of Eq. 3 */
      //     vectY[i] = linearPoly(dx_yhigh, yvalues[yhigh], kx + 0.5, ky + 0.5, temp);
      //   else if (dx_ylow <= temp && temp <= dx_yhigh) /* inverse of Eq. 2 */
      //     inverseSSpline(&yvalues[ylow], cmat, numcoeff, &vectX[i], 1, &vectY[i]);
      //   else {
      //     /* temp=vectX[i] should not be here */
      //     printf("inversePostSmooth: Input Error\n");
      //     printf("xvalue = %8.5f %8.5f %8.5f \n", temp, dx_ylow, dx_yhigh);
      //     exit(-1);
      //   }
      // }
    }
  };
} // namespace EquatingRecipes

#endif