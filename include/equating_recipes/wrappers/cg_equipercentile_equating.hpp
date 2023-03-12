/* 

CG_EquiEquate.c - contains functions for computing equating results for 
          (a) frequency estimation (FE) and Braun-Holland 
              under FE (if requested);
          (b) modified frequency estimation (MFE) (Wang & Brennan, 2006) 
              and Braun-Holland under MFE (if requested); and
          (c) chained equipercentile equating

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

#ifndef CG_EQUIV_EQUATING_HPP
#define CG_EQUIV_EQUATING_HPP

#include <cmath>
#include <string>
#include <vector>
#include <Eigen/Core>

#include <equating_recipes/wrappers/utilities.hpp>
#include <equating_recipes/structures/cg_equipercentile_equating_results.hpp>

namespace EquatingRecipes {
  struct CGEquipercentileEquating {
    /*
      Computes results for common-item (CI) equipercentile equating for EITHER
      (a) frequency estimation (FE) and Braun-Holland under FE (if requested); OR
      (b) modified frequency estimation (MFE) (Wang & Brennan, 2006) and 
          Braun-Holland under MFE (if requested).

      Input
        w1        weight for pop 1
        internal  internal anchor if non-zero; external anchor if zero
      nsv       number of score categories for v
        minv      minimum score for v
        maxv      maximum score for v      
        nsx       number of score categories for x
        minx      minimum score for x
        maxx      maximum score for x
        nsy       number of score categories for y
        miny      minimum score for y
        maxy      maximum score for y
        inc       increment--assumed to be the same for x, y, and v 
        bxvin[][] biv rel freq dist for x (rows) and v -- input 
        byvin[][] biv rel freq dist for y (rows) and v -- input 
        rv1       reliability for common items in pop 1 (used for MFE)
        rv2       reliability for common items in pop 2 (used for MFE)

      Output
        fxs[]      rel freq dist for x in syn pop  
        gys[]      rel freq dist for y in syn pop  
        eraw[]     equipercentile equated raw scores 
        a[]        slope for Braun-Holland Method 
        b[]        intercept for Braun-Holland Method 
        erawBH[]   Braun-Holland linear equated raw scores 

      NOTES: (1) Results are for FE if rv1 == 0 or rv2 == 0;
                  otherwise results are for MFE
              (2) Assumes space allocated for fxs, gys, and eraw
                  before function called
              (3) Braun-Holland results (a, b, erawBH) provided
                  if a != NULL and b !=NULL and erawBH != NULL
              (4) If Braun-Holland results are desired, then before
                  function is called, space must be allocated for erawBH;
                  also, a and b variables must be declared and their
                  addresses passed
              (5) Braun-Holland results are provided along with 
                  equipercentile results; Braun Holland results not
                  provided in isolation (refer to 2, above)

      Function calls other than C or NR utilities:
      SyntheticDensities()
        cum_rel_freqs()
        perc_rank()
        EquiEquate()
        BH_LinEq()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */

    // w1        weight for pop 1
    //   internal  internal anchor if non-zero; external anchor if zero
    // nsv       number of score categories for v
    //   minv      minimum score for v
    //   maxv      maximum score for v
    //   nsx       number of score categories for x
    //   minx      minimum score for x
    //   maxx      maximum score for x
    //   nsy       number of score categories for y
    //   miny      minimum score for y
    //   maxy      maximum score for y
    //   inc       increment--assumed to be the same for x, y, and v
    //   bxvin[][] biv rel freq dist for x (rows) and v -- input
    //   byvin[][] biv rel freq dist for y (rows) and v -- input
    //   rv1       reliability for common items in pop 1 (used for MFE)
    //   rv2       reliability for common items in pop 2 (used for MFE)

    // Output
    //   fxs[]      rel freq dist for x in syn pop
    //   gys[]      rel freq dist for y in syn pop
    //   eraw[]     equipercentile equated raw scores
    //   a[]        slope for Braun-Holland Method
    //   b[]        intercept for Braun-Holland Method
    //   erawBH[]   Braun-Holland linear equated raw scores

    // double w1, int internal, int nsv, 
    //             int nsx, double minx, double maxx, 
    //             int nsy, double miny, double maxy, double inc,
    //             double **bxvin, double **byvin, double rv1, double rv2,
    //             double *fxs, double *gys, double *eraw,
    //             double *a, double *b, double *erawBH
    
    EquatingRecipes::Structures::CGEquipercentileEquatingResults feOrMFEEquipEquating(const double& population1Weight,
                                                                                      const bool& isInternalAnchor,
                                                                                      const size_t& numberOfScoresV,
                                                                                      const size_t& numberOfScoresX,
                                                                                      const double& minimumScoreX,
                                                                                      const double& maximumScoreX,
                                                                                      const size_t& numberOfScoresY,
                                                                                      const double& minimumScoreY,
                                                                                      const double& maximumScoreY,
                                                                                      const double& scoreIncrement,
                                                                                      const bool& doBraunHollandLinearEquating,
                                                                                      Eigen::MatrixXd& bivariateRelativeFreqDistXV,
                                                                                      Eigen::MatrixXd& bivariateRelativeFreqDistYV,
                                                                                      const double& reliabilityCommonItemsPopulation1 = 0.0,
                                                                                      const double& reliabilityCommonItemsPopulation2 = 0.0) {
      EquatingRecipes::Structures::CGEquipercentileEquatingResults results;

      Eigen::VectorXd relFreqDistXSynPop; // (numberOfScoresX);
      Eigen::VectorXd relFreqDistYSynPop; // (numberOfScoresY);

      /* relFreqDistXSynPop and relFreqDistYSynPop: densities for x and y in synthetic population 
      Note: if 0 < reliabilityCommonItemsPopulation1, reliabilityCommonItemsPopulation2 <= 1, results are for MFE;
           otherwise, results are for FE */

      syntheticDensities(population1Weight,
                         isInternalAnchor,
                         numberOfScoresV,
                         numberOfScoresX,
                         bivariateRelativeFreqDistXV,
                         numberOfScoresY,
                         bivariateRelativeFreqDistYV,
                         reliabilityCommonItemsPopulation1,
                         reliabilityCommonItemsPopulation2,
                         relFreqDistXSynPop,
                         relFreqDistYSynPop);

      /* Fxs, Gys, and prxs */
      Eigen::VectorXd cumRelFreqDistXSynPop = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(minimumScoreX,
                                                                                                     maximumScoreX,
                                                                                                     scoreIncrement,
                                                                                                     relFreqDistXSynPop);

      Eigen::VectorXd cumRelFreqDistYSynPop = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(minimumScoreY,
                                                                                                     maximumScoreY,
                                                                                                     scoreIncrement,
                                                                                                     relFreqDistYSynPop);

      Eigen::VectorXd percentileRanksXSynPop(numberOfScoresX);
      for (size_t scoreLocation = 0; scoreLocation < numberOfScoresX; scoreLocation++) {
        double score = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                            minimumScoreX,
                                                            scoreIncrement);

        percentileRanksXSynPop(scoreLocation) = EquatingRecipes::Utilities::getPercentileRank(minimumScoreX,
                                                                                              maximumScoreX,
                                                                                              scoreIncrement,
                                                                                              cumRelFreqDistXSynPop,
                                                                                              score);
      }

      /* Equipercentile equating */
      results.equatedRawScores = EquatingRecipes::Utilities::getEquipercentileEquivalents(numberOfScoresY,
                                                                                          minimumScoreY,
                                                                                          scoreIncrement,
                                                                                          cumRelFreqDistYSynPop,
                                                                                          numberOfScoresX,
                                                                                          percentileRanksXSynPop);

      /* Braun-Holland linear equating */
      if (doBraunHollandLinearEquating) {
        double slope;
        double intercept;

        braunHollandLinearEquating(minimumScoreX,
                                   maximumScoreX,
                                   minimumScoreY,
                                   maximumScoreY,
                                   scoreIncrement,
                                   relFreqDistXSynPop,
                                   relFreqDistYSynPop,
                                   slope,
                                   intercept);

        results.slope = slope;
        results.intercept = intercept;
        Eigen::VectorXd braunHollandEquatedRawScores(numberOfScoresX);

        for (size_t scoreLocation = 0; scoreLocation < numberOfScoresX; scoreLocation++) {
          double score = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                              minimumScoreX,
                                                              scoreIncrement);

          braunHollandEquatedRawScores(scoreLocation) = slope * score + intercept;
        }

        results.braunHollandEquatedRawScores = braunHollandEquatedRawScores;
      }

      return results;
    }

    /*
      Synthetic population densities for x and y

      Input
        w1       weight for population 1 (x)
        internal internal anchor if non-zero; external anchor if zero
        nsv      number of score categories for common items (v)
        nsx      number of score categories for x
        bxv      starts as bivariate rel freq dist for v (rows) and x (cols);
                ends as biv rel freq dist of x given v with v as rows
        nsy      number of score categories for y
        byv      starts as bivariate rel freq dist for y (rows) and v (cols);
                ends as biv rel freq dist of y given v with y as rows
        rv1      rel of common items in pop 1 (0-->FE; non-0-->MFE)
        rv2      rel of common items in pop 2 (0-->FE; non-0-->MFE)

      Output
        fxs    synthetic pop distribution (relative frequencies) for x
        gys    synthetic pop distribution (relative frequencies) for y

      NOTE: space for fxs and gys allocated before call

      Function calls other than C or NR utilities:
        MixSmooth()
        CondBivDist(() 
        ModCondBivDist()
        runerror()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 2/3/09 
    */
    void syntheticDensities(const double& population1Weight,
                            const bool& isInternalAnchor,
                            const size_t& numberOfScoresV,
                            const size_t& numberOfScoresX,
                            Eigen::MatrixXd& bivariateRelativeFreqDistXV,
                            const size_t& numberOfScoresY,
                            Eigen::MatrixXd& bivariateRelativeFreqDistYV,
                            const double& reliabilityCommonItemsPopulation1,
                            const double& reliabilityCommonItemsPopulation2,
                            Eigen::VectorXd& syntheticPopulationRelativeFreqDistX,
                            Eigen::VectorXd& syntheticPopulationRelativeFreqDistY) {
      Eigen::VectorXd postSmoothingRelFreqDistVPop1(numberOfScoresV); /* rel freq dist for v in pop 1 after smoothing */
      Eigen::VectorXd postSmoothingRelFreqDistVPop2(numberOfScoresV); /* rel freq dist for v in pop 2 after smoothing */
      Eigen::VectorXd postSmoothingRelFreqDistXPop1(numberOfScoresX); /* rel freq dist of x in pop 1 after smoothing */
      Eigen::VectorXd postSmoothingRelFreqDistYPop2(numberOfScoresY); /* rel freq dist of y in pop 2 after smoothing */
      double postSmoothingMeanVPop1 = 0.0;                            /* mean of v in pop 1 after smoothing */
      double postSmoothingMeanVPop2 = 0.0;                            /* mean of v in pop 2 after smoothing */
      double uniformMixingProportion = 1.0e-10;

      /* Smooth f1(x,v); get h1(v) & f1(x) */
      mixSmooth(numberOfScoresV,
                numberOfScoresX,
                uniformMixingProportion,
                bivariateRelativeFreqDistXV,
                postSmoothingRelFreqDistXPop1,
                postSmoothingRelFreqDistVPop1);

      uniformMixingProportion = 1.0e-10;

      /* Smooth g2(y,v); get h2(v) & g2(y) */
      mixSmooth(numberOfScoresV,
                numberOfScoresY,
                uniformMixingProportion,
                bivariateRelativeFreqDistYV,
                postSmoothingRelFreqDistYPop2,
                postSmoothingRelFreqDistVPop2);

      /* get f1(x|v) with v as rows */
      conditionalBivariateDistribution(numberOfScoresX,
                                       numberOfScoresV,
                                       bivariateRelativeFreqDistXV,
                                       postSmoothingRelFreqDistVPop1);

      /* get g2(y|v) with v as rows */
      conditionalBivariateDistribution(numberOfScoresY,
                                       numberOfScoresV,
                                       bivariateRelativeFreqDistYV,
                                       postSmoothingRelFreqDistVPop2);

      /* modified conditional bivariate distributions for MFE */
      if (reliabilityCommonItemsPopulation1 != 0.0 &&
          reliabilityCommonItemsPopulation2 != 0.0) {
        if (reliabilityCommonItemsPopulation1 < 0.0 ||
            reliabilityCommonItemsPopulation1 > 1.0 ||
            reliabilityCommonItemsPopulation2 < 0.0 ||
            reliabilityCommonItemsPopulation2 > 1.0) {
          throw std::runtime_error("\nrv1 or rv2 invalid");
        }

        for (size_t scoreLocation = 0; scoreLocation < numberOfScoresV; scoreLocation++) {
          postSmoothingMeanVPop1 += static_cast<double>(scoreLocation) * postSmoothingRelFreqDistVPop1(scoreLocation);
          postSmoothingMeanVPop2 += static_cast<double>(scoreLocation) * postSmoothingRelFreqDistVPop2(scoreLocation);
        }

        /* mod f2(x|v) */
        modifiedConditionalBivariateDistribution(isInternalAnchor,
                                                 numberOfScoresV,
                                                 numberOfScoresX,
                                                 reliabilityCommonItemsPopulation1,
                                                 reliabilityCommonItemsPopulation2,
                                                 postSmoothingMeanVPop1,
                                                 postSmoothingMeanVPop2,
                                                 bivariateRelativeFreqDistXV);

        /* mod g1(y|v) */
        modifiedConditionalBivariateDistribution(isInternalAnchor,
                                                 numberOfScoresV,
                                                 numberOfScoresX,
                                                 reliabilityCommonItemsPopulation2,
                                                 reliabilityCommonItemsPopulation1,
                                                 postSmoothingMeanVPop2,
                                                 postSmoothingMeanVPop1,
                                                 bivariateRelativeFreqDistYV);
      }

      /* synthetic densities using Kolen and Brennan (2004, Eq 5.8) */
      syntheticPopulationRelativeFreqDistX.resize(numberOfScoresX);
      syntheticPopulationRelativeFreqDistY.resize(numberOfScoresY);

      /* syn density for x */
      for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
        syntheticPopulationRelativeFreqDistX(scoreLocationX) = population1Weight * postSmoothingRelFreqDistXPop1(scoreLocationX);

        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          syntheticPopulationRelativeFreqDistX(scoreLocationX) += (1.0 - population1Weight) *
                                                                  bivariateRelativeFreqDistXV(scoreLocationX, scoreLocationV) *
                                                                  postSmoothingRelFreqDistVPop2(scoreLocationV);
        }
      }

      /* syn density for y */
      for (size_t scoreLocationY = 0; scoreLocationY < numberOfScoresY; scoreLocationY++) {
        syntheticPopulationRelativeFreqDistY(scoreLocationY) = (1.0 - population1Weight) * postSmoothingRelFreqDistYPop2(scoreLocationY);

        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          syntheticPopulationRelativeFreqDistY(scoreLocationY) += population1Weight *
                                                                  bivariateRelativeFreqDistYV(scoreLocationY, scoreLocationV) *
                                                                  postSmoothingRelFreqDistVPop1(scoreLocationV);
        }
      }
    }

    /*
      "Smooth" bivariate distribution of v and x by mixing it with a 
      "small" uniform distribution.  The primary purpose of doing so is to 
      replace f(v)=0 with slightly positive relative frequencies.
      Note that descriptions and code are presented in terms of x and v,
      but function obviously applies to y and v, as well.
      
      Input
        nsv     number of categories for v
        nsx	    number of categories for x
        unimix	mixing proportion for uniform distribution; 1-unimix
                is mixing proportion for observed distribution.
        bxv	    observed bivariate rel freqs for x (rows) and v (cols)

      Output
        bxv	    "smoothed" bivariate probabilities of x and v 
        fx		"smoothed" marginal univariate probabilities for x 
        hv		"smoothed" marginal univariate probabilities for v

        NOTE: space for fx and hv must be allocated before function called

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates/revisions by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void mixSmooth(const size_t& numberOfScoresV,
                   const size_t& numberOfScoresX,
                   double& uniformMixingProportion,
                   Eigen::MatrixXd& bivariateProbabilitiesXV,
                   Eigen::VectorXd& marginalProbabilitiesX,
                   Eigen::VectorXd& marginalProbabilitiesV) {
      marginalProbabilitiesV.setZero(numberOfScoresV);
      marginalProbabilitiesX.setZero(numberOfScoresX);

      double observedMixingProportion = 1.0 - uniformMixingProportion;                          /* weight for actual relative freq in a cell */
      double uniformProbability = 1.0 / static_cast<double>(numberOfScoresV * numberOfScoresX); /* cell probability = 1/(nrow*ncol) */

      uniformMixingProportion *= uniformProbability; /* unimix scaled so that sum of "smoothed" bxv[][]=1 */

      for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          double smoothedRelativeFrequency = (observedMixingProportion *
                                              bivariateProbabilitiesXV(scoreLocationX, scoreLocationV)) +
                                             uniformMixingProportion;
          bivariateProbabilitiesXV(scoreLocationX, scoreLocationV) = smoothedRelativeFrequency; /* "smoothed" elements of bxv[][] */
          marginalProbabilitiesV(scoreLocationV) += smoothedRelativeFrequency;                  /* "smoothed" rel freq of i-th element of v */
        }
      }

      for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          marginalProbabilitiesX(scoreLocationX) += bivariateProbabilitiesXV(scoreLocationX, scoreLocationV); /* "smoothed rel freq of j-th element of fx */
        }
      }
    }

    /*
      Conditional bivariate distribution of x given v. (Obviously applies
      as well to conditional bivariate distribution of y given v.)

      Input
        nsv    number of categoies for v
        nsx    number of categories for x 
        bxv    bivariate relative freq distribution of x (rows) and v (cols)
        hv     mariginal rel freq distribution of v
    
      Output
        bxv	   conditional distribution of x given v (or y given v)
                  where v is rows

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void conditionalBivariateDistribution(const size_t& numberOfScoresV,
                                          const size_t& numberOfScoresX,
                                          Eigen::MatrixXd& bivariateRelativeFreqDistXV,
                                          const Eigen::VectorXd& marginalRelativeFreqDistV) {
      for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          bivariateRelativeFreqDistXV(scoreLocationX, scoreLocationV) /= marginalRelativeFreqDistV(scoreLocationV);
        }
      }
    }

    /*
      Braun-Holland linear equating (based on FE or MFE systhetic densities

      Input
        minx  min score for x
        maxx  max score for x
        miny  min score for y
        maxy  max score for y
        inc   increment for both x and y
        fxs   synthetic density for x
        gsy   synthetic density for y

      Output
        a    slope
        b    intercept

      Function calls other than C or NR utilities:
        MomentsFromRFD()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void
    braunHollandLinearEquating(const size_t& minimumScoreX,
                               const size_t& maximumScoreX,
                               const size_t& minimumScoreY,
                               const size_t& maximumScoreY,
                               const double& scoreIncrement,
                               const Eigen::VectorXd& syntheticDensityX,
                               const Eigen::VectorXd& syntheticDensityY,
                               double& slope,
                               double& intercept) {
      /* moments for x in syn pop */
      EquatingRecipes::Structures::Moments momentsXSynPop = EquatingRecipes::Utilities::momentsFromScoreFrequencies(syntheticDensityX,
                                                                                                                          minimumScoreX,
                                                                                                                          maximumScoreX,
                                                                                                                          scoreIncrement);

      /* moments for y in syn pop */
      EquatingRecipes::Structures::Moments momentsYSynPop = EquatingRecipes::Utilities::momentsFromScoreFrequencies(syntheticDensityY,
                                                                                                                          minimumScoreY,
                                                                                                                          maximumScoreY,
                                                                                                                          scoreIncrement);

      slope = momentsYSynPop.momentValues(1) / momentsXSynPop.momentValues(1);
      intercept = momentsYSynPop.momentValues(0) - slope * momentsXSynPop.momentValues(0);
    }

    /*
      For MFE method, modify conditional bivariate distribution of
      x given v with x as rows.  Description here worded 
      in terms of x given v, and code written in same manner, but 
      obviously function applies to y given v, as well.

      Input
        internal internal anchor if non-zero; external anchor if zero
        nsv      number of score categories for v
        nsx      number of score categories for x
        rv1      rel of common items in pop 1
        rv2      rel of common items in pop 2
        muv1     mean of v in pop 1
        muv2     mean of v in pop 2
        bxv      bivariate distribution of x given v in pop 1 with x as rows

      Output
        bxv      modified biv dist of x given v in pop 2 with x as rows

      NOTE:  There is a complexity here for an internal anchor. Let 
              nsu = nsx - nsv = # categories for non-equating items. Recall that
              x = 0,...,nsx-1 and v = 0,..., nsv-1.  In
              the bxv matrix there is a structural 0 whenever
              (x<v || x>nsu-1+v), which creates linear interpolation problems
              that are circumvented here by collapsing the xv matrix to a 
              uv matrix, doing the linear interpolation, and then expanding
              the uv back to an xv matrix. (Note that there is a one-to-one 
              relationship between elements of original xv matrix and 
              collapsed uv matrix.)

      Function calls other than C or NR utilities: 
        interpolate()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void modifiedConditionalBivariateDistribution(const bool& isInternalAnchor,
                                                  const size_t& numberOfScoresV,
                                                  const size_t& numberOfScoresX,
                                                  const double& reliabilityCommonItemsPopulation1,
                                                  const double& reliabilityCommonItemsPopulation2,
                                                  const double& population1MeanV,
                                                  const double& population2MeanV,
                                                  Eigen::MatrixXd& bivariateDistributionXVPopulation) {
      //   double v1,                   /* non integer v in pop 1 associated with integer v in pop 2 */

      /* slope for MFE */
      double slope = std::sqrt(reliabilityCommonItemsPopulation2 /
                               reliabilityCommonItemsPopulation1);

      /* intercept for MFE */
      double intercept = ((1.0 - std::sqrt(reliabilityCommonItemsPopulation2)) / std::sqrt(reliabilityCommonItemsPopulation1)) * population2MeanV -
                         ((1.0 - std::sqrt(reliabilityCommonItemsPopulation1)) / std::sqrt(reliabilityCommonItemsPopulation1)) * population1MeanV;

      Eigen::MatrixXd bXVWorkingValues(numberOfScoresX, numberOfScoresV);
      Eigen::MatrixXd bUVCollapsedValues;
      Eigen::MatrixXd bUVInterpolatedValues;

      //       intercept = ((1 - sqrt(rv2)) / sqrt(rv1)) * muv2 -
      //                   ((1 - sqrt(rv1)) / sqrt(rv1)) * muv1,
      //       **temp,                                           /* working values of bxv */
      //       **temp_collapsed,                                 /* working values of buv */
      //       **temp_interpolated,                              /* interpolated values of buv */
      //       sum;
      //   int v, j, u,
      //       v2,  /* integer v score in pop 2 */
      //       nsu; /* number of non-anchor score categories */

      //   temp = dmatrix(0, nsx - 1, 0, nsv - 1); /* allocate */

      if (!isInternalAnchor) {
        /****************************external anchor ****************************/
        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          double transformedScore = intercept + slope * static_cast<double>(scoreLocationV);

          for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
            bXVWorkingValues(scoreLocationX, scoreLocationV) = EquatingRecipes::Utilities::interpolate(transformedScore,
                                                                                                       numberOfScoresV,
                                                                                                       bivariateDistributionXVPopulation.row(scoreLocationX));
          }
        }
      } else {
        /****************************internal anchor ****************************/
        size_t numberOfScoresU = numberOfScoresX - numberOfScoresV + 1;
        bUVCollapsedValues.resize(numberOfScoresU, numberOfScoresV);
        bUVInterpolatedValues.resize(numberOfScoresU, numberOfScoresV);

        /* store non-structural-zero elements of bxv[][]
        in temp_collapsed[][]; in essense, collapse bxv[][] to buv[][] */
        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          size_t scoreLocationU = 0;

          for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
            if (scoreLocationX < scoreLocationV || scoreLocationX > numberOfScoresU - 1 + scoreLocationV) {
              /* do nothing: skip structural zeros */
            } else {
              bUVCollapsedValues(scoreLocationU++, scoreLocationV) = bivariateDistributionXVPopulation(scoreLocationX,
                                                                                                       scoreLocationV);
            }
          }
        }

        /* get interpolated values under MFE for internal anchor */
        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          double score = intercept + slope * static_cast<double>(scoreLocationV);

          for (size_t scoreLocationU; scoreLocationU < numberOfScoresU; scoreLocationU++) {
            bUVInterpolatedValues(scoreLocationU, scoreLocationV) = EquatingRecipes::Utilities::interpolate(score,
                                                                                                            numberOfScoresV,
                                                                                                            bUVCollapsedValues.row(scoreLocationU));
          }
        }

        /* expand temp_interpolated[0...nsu][v2] to temp[0...nsx][v2] */
        for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
          size_t scoreLocationU = 0;

          for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
            if (scoreLocationX < scoreLocationV || scoreLocationX > numberOfScoresU - 1 + scoreLocationV) {
              bXVWorkingValues(scoreLocationX, scoreLocationV) = 0.0;
            } else {
              bXVWorkingValues(scoreLocationX, scoreLocationV) = bUVInterpolatedValues(scoreLocationU++, scoreLocationV);
            }
          }
        }
      }

      /* for both external and internal anchor
      normalize such that sum-over-x of f(x|v) = 1 */
      for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
        double sum = 0.0;

        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          sum += bXVWorkingValues(scoreLocationX, scoreLocationV);
        }

        /* transfer to bxv[][] */
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          if (sum > 0) {
            bivariateDistributionXVPopulation(scoreLocationX, scoreLocationV) = bXVWorkingValues(scoreLocationX, scoreLocationV) / sum;
          } else {
            bivariateDistributionXVPopulation(scoreLocationX, scoreLocationV) = 0;
          }
        }
      }
    }

    /*
      Chained equipercentile equating using the composed function
      discussed in Kolen and Brennan (2004, pp. 145-147)

      Input

        nsx    = number of score categories for X in pop 1
        prx1[] = PR for X in pop 1
        minv   = min score for V in both pop 1 and pop 2
        maxv   = max score for V in both pop 1 and pop 2
        incv   = increment for V in both pop 1 and pop 2
        nsv    = # of score categories for V in both pops
        Gv1[]  = crfd for V in pop 1

        miny   = min score for Y in pop 2
        incy   = increment for Y in pop 2
        nsy    = number of score categories for Y in pop 2
        Gy2[]  = crfd for Y in pop 2
        Gv2[]  = crfd for V in pop 2

      Output

        eraw[] = chained equipercentile Y-equivalents of X
        
      Notes.  (a) It is assumed that space for eraw[] allocated
                  prior to function call
              (b) minv, maxv, incv, and nsv are 
                  necessarily the same for both pops

      Function calls other than C or NR utilities:
        EquiEquate()
        perc_rank()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    Eigen::VectorXd chainedEquipercentileEquating(const size_t& numberOfScoresX,
                                                  const Eigen::VectorXd& percentileRanksX,
                                                  const double& minimumScoreV,
                                                  const double& maximumScoreV,
                                                  const double& scoreIncrementV,
                                                  const size_t& numberOfScoresV,
                                                  const Eigen::VectorXd& cumlativeRelativeFreqDistVPop1,
                                                  const double& minimumScoreY,
                                                  const double& scoreIncrementY,
                                                  const size_t& numberOfScoresY,
                                                  const Eigen::VectorXd& cumlativeRelativeFreqDistYPop2,
                                                  const Eigen::VectorXd& cumlativeRelativeFreqDistVPop2) {
      Eigen::VectorXd percentileRanksVEquivalents(numberOfScoresX);

      /* Put X on scale of V in pop 1; there are nsx
        (non-integer) V equivalents in extov[] */
      Eigen::VectorXd vEquivalentsXPop1 = EquatingRecipes::Utilities::getEquipercentileEquivalents(numberOfScoresV,
                                                                                                   minimumScoreV,
                                                                                                   scoreIncrementV,
                                                                                                   cumlativeRelativeFreqDistVPop1,
                                                                                                   numberOfScoresX,
                                                                                                   percentileRanksX);

      /* Get PRs (relative to V for pop 2)
        for the non-integer V scores in extov[].
        Note that there are nsx equivalents */
      for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
        percentileRanksVEquivalents(scoreLocationX) = EquatingRecipes::Utilities::getPercentileRank(minimumScoreV,
                                                                                                    maximumScoreV,
                                                                                                    scoreIncrementV,
                                                                                                    cumlativeRelativeFreqDistVPop2,
                                                                                                    vEquivalentsXPop1(scoreLocationX));
      }

      /* Using the PRs in prv2[] get the Y equivalents of X;
      i.e., put V equivalents on scale of Y */

      Eigen::VectorXd equatedRawScores = EquatingRecipes::Utilities::getEquipercentileEquivalents(numberOfScoresY,
                                                                                                  minimumScoreY,
                                                                                                  scoreIncrementY,
                                                                                                  cumlativeRelativeFreqDistYPop2,
                                                                                                  numberOfScoresX,
                                                                                                  percentileRanksVEquivalents);

      return equatedRawScores;
    }

    // void Print_SynDens(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
  };
} // namespace EquatingRecipes

#endif