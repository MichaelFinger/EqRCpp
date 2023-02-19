/* 	
  LogLinear.c  code for log-linear equating

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

  File contains functions used for:
      (a) univariate log-linear fitting under the multinomial model,
          which populates struct ULL_SMOOTH 
      (b) bivariate log-linear fitting under the multinomial model,
          which populates struct BLL_SMOOTH 

  Code follows procedures (and usually the notation) in Holland and
  Thayer (1987), which is sometimes abbreviated H&T in comments.

  A somewhat unique feature of these functions is that the user 
  can select from among various convergence criteria, including
  criteria based on moments for the actual raw scores (e.g., 
  mean, sd, skew, kurt, etc for the scores obtained using
  score() in ERutilities.c).  
*/

#ifndef LOG_LINEAR_EQUATING_HPP
#define LOG_LINEAR_EQUATING_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <fmt/core.h>

#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/score_statistics.hpp>
#include <equating_recipes/utilities.hpp>

#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  class LogLinearEquating {
  public:
    // ctype = comparison type for criterion:
    //         0 --> absolute; 1 --> relative
    // Btype = type of design matrix and, hence, type of moments
    //         for criterion:
    //         0 --> use B (could be scaled or unscaled as indicated
    //               in design_matrix()) --- see note below
    //         1 --> use B_raw and central moments based on it

    enum class CriterionComparisonType {
      ABSOLUTE,
      RELATIVE
    };

    enum class DesignMatrixType {
      SOLUITON,
      RAW_SCORE
    };

    /*
      Wrapper to do univariate log-linear smoothing.

      Input
        x   =  UTATS structure
        c   = number of degrees for polynomial smoothing
        scale = type of scaling:
              0 --> no scaling; 
                1 --> scale such that each column of B has
                      sum (elements) = 0 and sum (elements^2) = 1
        Btype = type of moments for criterion mathching:
                0 --> moments based on B 
                (if scale = 0, design matrix is based on raw scores,
                which means that the moments are based on raw scores;
                if scale = 1, design matrix is based on scaled raw scores,
                which means that the moments are based on scaled raw scores) 
                1 -->  moments based on B_raw, whether scale is 0 or 1
      ctype = comparison type for criterion:
              0 --> means use absolute criterion; 
          1 --> means use relative criterion
      crit = convergence criterion value
        fp  = pointer to output file (if NULL then no output printed;
            in particular, results for each iteration step are not printed)

      Output
        populates s, which is a ULL_SMOOTH structure

      Function calls other than C or NR utilities:  
        Smooth_ULL()

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    void runUnivariateLogLinearSmoothing(const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsX,
                                         const size_t& numberOfDegreesSmoothing,
                                         const bool& useStandardizedScale,
                                         const DesignMatrixType& designMatrixType,
                                         const CriterionComparisonType& criterionComparisonType,
                                         const double& criterion,
                                         EquatingRecipes::Structures::UnivariateLogLinearSmoothing& smoothingResults) {
      smoothUnivaraiteLogLinear(univariateStatisticsX.numberOfExaminees,
                                univariateStatisticsX.numberOfScores,
                                univariateStatisticsX.minimumScore,
                                univariateStatisticsX.scoreIncrement,
                                univariateStatisticsX.freqDistDouble,
                                numberOfDegreesSmoothing,
                                useStandardizedScale,
                                designMatrixType,
                                criterionComparisonType,
                                criterion,
                                smoothingResults);
    }

    /*
      Wrapper for doing equipercentile equating with RG design
        and log-linear smoothing smoothing
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep

      Input
      
        design = 'R'(random groups)
        method = 'E'(equipercentile)
        smoothing = 'L' (log-linear smoothing)  
        *x = pointer to struct USTATS (new form)
        *y = pointer to struct USTATS (old form)
        *ullx = pointer to struct BB_SMOOTH (new form)
        *ully = pointer to struct BB_SMOOTH (old form)
        rep = replication number for bootstrap; should be set to 0
              for actual equating;  
        
      Output
        
        struct PDATA *inall:   populates selected values of inall 
        
        struct ERAW_RESULTS *r: populates            

          **eraw: equated raw scores;          
                  method (rows) by raw score (columns) matrix
                  of equated scores. Here there is only one method.
                  So, memory allocated for eraw[][] is: 
                  eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                                (loc(x->max,x->min,x>-inc)]
                  because we are getting equated raw scores for x 
          **mts:  moments for equated raw scores           
          
      NOTE: If Wrapper_RL() is called in a bootstrap loop,
            then in the calling function struct ERAW_RESULTS must
            be different from struct ERAW_RESULTS for the actual
            equating. 

      Function calls other than C or NR utilities:                   
        EquiEquate()
        MomentsFromFD()  
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08       
    */
    void runRGEquiEquatingWithLoglinearSmoothing(const EquatingRecipes::Structures::Design& design,
                                                 const EquatingRecipes::Structures::Method& method,
                                                 const EquatingRecipes::Structures::Smoothing& smoothing,
                                                 const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsX,
                                                 const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsY,
                                                 const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& univariateLogLinearSmoothingX,
                                                 const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& univariateLogLinearSmoothingY,
                                                 const size_t& replicationNumber,
                                                 EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating */
                                                            /* counting of replications done in Wrapper_Bootstrap() */

      /* allocation and assignments for struct PDATA inall
        Note that for every assignment of the form inall->(var) = x->(var)
        or inall->(var) = y->(var), values vary depending on whether x or y 
        is for actual equating or a bootstrap sample; all other values are 
        the same for the actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) {
        pData.summaryRawDataX = univariateStatisticsX;
        pData.summaryRawDataY = univariateStatisticsY;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        pData.methods.push_back("   y-equiv");
        pData.mininumScoreX = univariateStatisticsX.minimumScore;
        pData.maximumScoreX = univariateStatisticsX.maximumScore;
        pData.scoreIncrementX = univariateStatisticsX.scoreIncrement;
        pData.scoreFrequenciesX = univariateStatisticsX.freqDistDouble;
        pData.numberOfExaminees = univariateStatisticsX.numberOfExaminees;
        pData.univariateLogLinearSmoothingX = univariateLogLinearSmoothingX;
        pData.univariateLogLinearSmoothingY = univariateLogLinearSmoothingY;
      }

      if (pData.bootstrapReplicationNumber <= 1) {
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.mininumScoreX,
                                                                                   pData.scoreIncrementX);
        equatedRawScoreResults.equatedRawScores.resize(pData.methods.size(), maximumScoreLocation + 1);
        equatedRawScoreResults.equatedRawScoreMoments.resize(1, 4);
      }

      /* Compute equating results */
      Eigen::VectorXd equatedRawScores = EquatingRecipes::Utilities::getEquipercentileEquivalents(univariateStatisticsY.numberOfScores,
                                                                                                  univariateStatisticsY.minimumScore,
                                                                                                  univariateStatisticsY.scoreIncrement,
                                                                                                  univariateLogLinearSmoothingY.fittedRawScoreCumulativeRelativeDist,
                                                                                                  univariateStatisticsX.numberOfScores,
                                                                                                  univariateLogLinearSmoothingX.fittedRawScorePercentileRankDist);
      for (size_t index = 0; index < equatedRawScores.size(); index++) {
        equatedRawScoreResults.equatedRawScores(0, index) = equatedRawScores(index);
      }

      /* get moments */
      EquatingRecipes::Structures::Moments moments = EquatingRecipes::ScoreStatistics::momentsFromScoreFrequencies(equatedRawScores,
                                                                                                                   pData.scoreFrequenciesX);

      for (size_t index = 0; index < moments.momentValues.size(); index++) {
        equatedRawScoreResults.equatedRawScoreMoments(0, index) = moments.momentValues(index);
      }
    }

    /*
      Wrapper for doing equipercentile equating with SG design
      and log-linear smoothing.  
      
      NOTE: This is for the SG design in which x and y do not share any items in 
      common, which means that functionally this is the external anchor case.  
      The bivariate log-linear smoothing procedure needs to know this. So, when
      Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
      anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
      (with anchor set to 0) and Wrapper_SL() can still be used, but convergence 
      of the smoothing algorithm may be compromised because of dependencies between
      x and y. (To date my experience does not suggest this is a problem.)
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep

      Input
      
        design = 'S' (single group)
        method = 'E' (equipercentile)
        smoothing = 'L' (log-linear smoothing)  
        xy = struct BSTATS 
      bllxy = struct BLL_SMOOTH
        rep = replication number for bootstrap; should be set to 0
              for actual equating;  
        
        NOTE: it is assumed that the first variable
              in xy is indeed x and the second variable is y.
          For bllxy, data memebers with '_x' are for x,
          and data members with '_v' are for y.  This somewhat
          inconsistent notation arises because BLL_SMOOTH is
          usually used with the CG design in which the variables
          are more naturally designated x (or u) and v.

      Output
        
        struct PDATA inall   populates selected values of inall  

        struct ERAW_RESULTS *r: populates            

          **eraw: equated raw scores;          
                  method (rows) by raw score (columns) matrix
                  of equated scores. Here there is only one method.
                  So, memory allocated for eraw[][] is: 
                  eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                                (loc(x->max,x->min,x>-inc)]
                  because we are getting equated raw scores for x 
          **mts:  moments for equated raw scores            
          
      NOTE: If Wrapper_SL() is called in a bootstrap loop,
            then in the calling function struct ERAW_RESULTS must
            be different from struct ERAW_RESULTS for the actual
            equating. 
                                                
      Function calls other than C or NR utilities:
        EquiEquate()
        MomentsFromFD()  
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08   
    */
    void runSGEquiEquatingWithLoglinearSmoothing(const EquatingRecipes::Structures::Design& design,
                                                 const EquatingRecipes::Structures::Method& method,
                                                 const EquatingRecipes::Structures::Smoothing& smoothing,
                                                 const EquatingRecipes::Structures::BivariateStatistics& xy,
                                                 const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingXY,
                                                 const size_t& replicationNumber,
                                                 EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating */
                                                            /* counting of replications done in Wrapper_Bootstrap() */

      /* Allocation and assignments for struct PDATA inall>
    Note that for every assignment of the form inall->(var) = x->(var)
    or inall->(var) = y->(var), values vary depending on whether x or y 
 	  is for actual equating or a bootstrap sample; all other values are 
	  the same for the actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) {
        pData.summaryRawDataXY = xy;
        pData.bivariateLogLinearSmoothingXY = bivariateLogLinearSmoothingXY;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        pData.isInternalAnchor = false; /* implicitly, anchor is external for biv log-linear
						                        smoothing with the SG design */

        pData.methods.push_back("   y-equiv");

        pData.mininumScoreX = xy.univariateStatisticsRow.minimumScore;
        pData.maximumScoreX = xy.univariateStatisticsRow.maximumScore;
        pData.scoreIncrementX = xy.univariateStatisticsRow.scoreIncrement;
        pData.scoreFrequenciesX = xy.univariateStatisticsRow.freqDistDouble;
        pData.numberOfExaminees = xy.numberOfExaminees;
      }

      if (pData.bootstrapReplicationNumber <= 1) {
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.mininumScoreX,
                                                                                   pData.scoreIncrementX);
        equatedRawScoreResults.equatedRawScores(pData.methods.size(), maximumScoreLocation + 1);
        equatedRawScoreResults.equatedRawScoreMoments.resize(1, 4);
      }

      /* Compute equating results. Put x on scale of y.
        Note that in struct xy, '1' designates x and '2' designates y; 
        in struct bllxy, '_x' designates x and '_v' designates y. So: 
          xy->ns2 = number of score categories for y
        xy->min2 = minimum score for y
        xy->inc2 = increment for y
        bllxy->crfd_v = log-linear smoothed cum rel fd for y
        xy->ns1 = number of score categories for x
          bllxy->prd_x = log-linear smoothed PR dist for x
        r->eraw[0] = y equivalents for x (output) */

      Eigen::VectorXd equatedRawScores = EquatingRecipes::Utilities::getEquipercentileEquivalents(xy.univariateStatisticsColumn.numberOfScores,
                                                                                                  xy.univariateStatisticsColumn.minimumScore,
                                                                                                  xy.univariateStatisticsColumn.scoreIncrement,
                                                                                                  bivariateLogLinearSmoothingXY.fittedRawScoreCumulativeRelativeFreqDistV,
                                                                                                  xy.univariateStatisticsRow.numberOfScores,
                                                                                                  bivariateLogLinearSmoothingXY.fittedRawScorePercentileRankDistX);

      for (size_t index = 0; index < equatedRawScores.size(); index++) {
        equatedRawScoreResults.equatedRawScores(0, index) = equatedRawScores(index);
      }

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::ScoreStatistics::momentsFromScoreFrequencies(equatedRawScoreResults.equatedRawScores,
                                                                                                                   pData.scoreFrequenciesX);

      for (size_t index = 0; index < moments.momentValues.size(); index++) {
        equatedRawScoreResults.equatedRawScoreMoments(0, index) = moments.momentValues(index);
      }
    }

    /*
      Wrapper for equipercentile equating for CG design with log-linear smoothing. 
      Equipercentile equating includes frequency estimation with 
      Braun-Holland (linear) results, modified frequency estimation with 
      Braun-Holland (linear) results, and chained equipercentile equating
        
      Assumes that in xv, score 1 is for x and score 2 is for v
      Assumes that in yv, score 1 is for y and score 2 is for v
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep
      
      Input:
      
        design = 'C' (CINEG)

        method:  'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
                'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
                'G' = FE + BH-FE + MFE + BH-MFE
                'C' = Chained
          'H' = FE + BH-FE + Chained
                'A' = FE + BH-FE + MFE + BH-MFE + Chained
                  
        smoothing = 'L' (log-linear) 

        w1 = weight for pop. 1 (associated with xv)
            [0,1] except that for any number outside this 
            range, proportional weights are used -- i.e.,
            w1 = xv->n/(xv->n + yv->n)
        anchor = 0 --> external; otherwise internal
        rv1 = reliability of common items for population 1 
              (set to 0 for all methods except 'F', 'G, and 'A')
        rv2 = reliability of common items for population 2
              (set to 0 for all methods except 'F', 'G, and 'A')
        xv = struct BSTATS
        yv = struct BSTATS 
        bllxv = struct BLL_SMOOTH; uses brfd, prd_x,  and crfd_v 
      bllyv = struct BLL_SMOOTH; uses brfd, crfd_x, and crfd_v
              (Note that bllyv->crfd_x is really crfd for y in pop 2)
        rep = replication number for bootstrap; should be set to 0
              for actual equating; 
      
        NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
        
      Output:
        
        struct PDATA inall:   populates selected values of inall 
        
        struct ERAW_RESULTS r          
    
          a[] = slopes for Braun-Holland
          b[] = intercepts for Braun-Holland
          eraw[][]:  equated raw scores
          mts[][]:  moments for equated raw scores   
              
          NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
                of equated scores; memory allocated here;
                eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
                because we are getting equated raw scores for x.
                eraw[][] stored in this "row"  manner so that 
                Equated_ss() can be used directly on 
                eraw[k] where k is the method number  
              
      NOTE: Whenever method differs, there must be different structures
            passed as struct PDATA and struct ERAW_RESULTS 
        
      NOTE: If Wrapper_CL() is called in a bootstrap loop, then in
            the calling function struct ERAW_RESULTS must be different
            from struct ERAW_RESULTS for the actual equating. 
                                                
      Function calls other than C or NR utilities:
        FEorMFE_EE()
        Chained_EE()
        runerror()
                                                  
      R. L. Brennan

      Date of last revision: 6/30/08   
    */
    void runCGEquiEquatingWithLoglinearSmoothing(const EquatingRecipes::Structures::Design& design,
                                                 const EquatingRecipes::Structures::Method& method,
                                                 const EquatingRecipes::Structures::Smoothing& smoothing,
                                                 const double& populationWeight1,
                                                 const bool& isInternalAnchor,
                                                 const double& reliabilityCommonItemsPopulation1,
                                                 const double& reliabilityCommonItemsPopulation2,
                                                 const EquatingRecipes::Structures::BivariateStatistics& xv,
                                                 const EquatingRecipes::Structures::BivariateStatistics& yv,
                                                 EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingXV,
                                                 EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingYV,
                                                 const size_t& replicationNumber,
                                                 EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      // char *names[] ={"        FE", "     BH-FE", "       MFE", "    BH-MFE",
      //             "  ChainedE"};
      std::vector<std::string> names {"        FE",
                                      "     BH-FE",
                                      "       MFE",
                                      "    BH-MFE",
                                      "  ChainedE"};

      std::string methodCode = EquatingRecipes::Utilities::getMethodCode(method);

      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating. */
                                                            /* Counting of replications done in Wrapper_Bootstrap(), 
             which is why this statement cannot be in the if statement below */

      /* allocation and assignments for inall
     Note that for every assignment of the form inall->(var) = xv->(var)
     or inall->(var) = yv->(var) values vary depending on whether xv or yv
     is for actual equating or a bootstrap sample; all other values are 
	   the same for the actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) {
        pData.summaryRawDataXV = xv;
        pData.summaryRawDataXV = xv;
        pData.bivariateLogLinearSmoothingXV = bivariateLogLinearSmoothingXV;
        pData.bivariateLogLinearSmoothingYV = bivariateLogLinearSmoothingYV;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        if (populationWeight1 < 0.0 || populationWeight1 > 1.0) {
          /* proportional wts if w1 outside [0,1] */
          pData.weightSyntheticPopulation1 = static_cast<double>(xv.numberOfExaminees) /
                                             static_cast<double>(xv.numberOfExaminees + yv.numberOfExaminees);
        } else {
          pData.weightSyntheticPopulation1 = populationWeight1;
        }

        pData.isInternalAnchor = isInternalAnchor;
        pData.reliabilityCommonItemsPopulation1 = reliabilityCommonItemsPopulation1;
        pData.reliabilityCommonItemsPopulation2 = reliabilityCommonItemsPopulation2;

        if ((methodCode == "F" || methodCode == "G" || methodCode == "A") &&
            (reliabilityCommonItemsPopulation1 == 0.0 || reliabilityCommonItemsPopulation2 == 0)) {
          throw std::runtime_error("\nMFE cannot be conducted since rv1 == 0 or rv2 == 0");
        }

        pData.methods.resize(5, "");

        if (methodCode == "E") {
          pData.methods[0] = names[0];
          pData.methods[1] = names[1];
        } else if (methodCode == "F") {
          pData.methods[0] = names[2];
          pData.methods[1] = names[3];
        } else if (methodCode == "G") {
          pData.methods[0] = names[0];
          pData.methods[1] = names[1];
          pData.methods[0] = names[2];
          pData.methods[1] = names[3];
        } else if (methodCode == "C") {
          pData.methods[0] = names[4];
        } else if (methodCode == "H") {
          pData.methods[0] = names[0];
          pData.methods[1] = names[1];
          pData.methods[2] = names[4];
        } else if (methodCode == "A") {
          pData.methods = names;
        } else {
          throw std::runtime_error("Invalid method in log linear equating.");
        }

        pData.mininumScoreX = xv.univariateStatisticsRow.minimumScore;
        pData.maximumScoreX = xv.univariateStatisticsRow.maximumScore;
        pData.scoreIncrementX = xv.univariateStatisticsRow.scoreIncrement;
        pData.scoreFrequenciesX = xv.univariateStatisticsRow.freqDistDouble;
        pData.numberOfExaminees = xv.numberOfExaminees;
      }

      if (pData.bootstrapReplicationNumber <= 1) {
        size_t maximumScoreLocationXV = EquatingRecipes::Utilities::getScoreLocation(xv.univariateStatisticsRow.maximumScore,
                                                                                     xv.univariateStatisticsRow.minimumScore,
                                                                                     xv.univariateStatisticsRow.scoreIncrement);

        size_t maximumScoreLocationYV = EquatingRecipes::Utilities::getScoreLocation(yv.univariateStatisticsRow.maximumScore,
                                                                                     yv.univariateStatisticsRow.minimumScore,
                                                                                     yv.univariateStatisticsRow.scoreIncrement);

        equatedRawScoreResults.equatedRawScores.resize(pData.methods.size(),
                                                       maximumScoreLocationXV + 1);
        equatedRawScoreResults.equatedRawScoreMoments.resize(pData.methods.size(),
                                                             4);
        equatedRawScoreResults.relativeFreqDistsX.resize(1, maximumScoreLocationXV);
        equatedRawScoreResults.relativeFreqDistsY.resize(1, maximumScoreLocationYV + 1);
      }

      /* Equipercentile results, including Braun-Holland (BH) linear. 
        Note: For FE syn densities are in fxs[0] and gys[0]
        For MFE syn densities are in fxs[1] and gys[1] 
        For BH under FE, slope in a[0] and intercept in b[0]
        For BH under MFE, slope in a[1] and intercept in b[1] */

      /* FE + BH-FE in positions 0 and 1*/
      EquatingRecipes::CGEquipercentileEquating cgEquipercentileEquating;

      if (methodCode == "E" || methodCode == "G" || methodCode == "A" || methodCode == "H") {
        EquatingRecipes::Structures::CGEquipercentileEquatingResults cgResults = cgEquipercentileEquating.feOrMFEEquipEquating(pData.weightSyntheticPopulation1,
                                                                                                                               pData.isInternalAnchor,
                                                                                                                               xv.univariateStatisticsColumn.numberOfScores,
                                                                                                                               xv.univariateStatisticsRow.numberOfScores,
                                                                                                                               xv.univariateStatisticsRow.minimumScore,
                                                                                                                               xv.univariateStatisticsRow.maximumScore,
                                                                                                                               yv.univariateStatisticsRow.numberOfScores,
                                                                                                                               yv.univariateStatisticsRow.minimumScore,
                                                                                                                               yv.univariateStatisticsRow.maximumScore,
                                                                                                                               yv.univariateStatisticsRow.scoreIncrement,
                                                                                                                               true,
                                                                                                                               bivariateLogLinearSmoothingXV.fittedBivariateRelativeFreqDistXV,
                                                                                                                               bivariateLogLinearSmoothingYV.fittedBivariateRelativeFreqDistXV,
                                                                                                                               0,
                                                                                                                               0);
        equatedRawScoreResults.relativeFreqDistsX.row(0) = cgResults.syntheticPopulationRelativeFreqDistX;
        equatedRawScoreResults.relativeFreqDistsY.row(0) = cgResults.syntheticPopulationRelativeFreqDistY;
        equatedRawScoreResults.equatedRawScores.row(0) = cgResults.equatedRawScores;

        if (cgResults.slope.has_value()) {
          equatedRawScoreResults.slope(0) = cgResults.slope.value();
        }

        if (cgResults.intercept.has_value()) {
          equatedRawScoreResults.intercept(0) = cgResults.intercept.value();
        }

        if (cgResults.braunHollandEquatedRawScores.has_value()) {
          equatedRawScoreResults.equatedRawScores.row(1) = cgResults.braunHollandEquatedRawScores.value();
        }
      }

      if (methodCode == "F") {
        EquatingRecipes::Structures::CGEquipercentileEquatingResults cgResults =
            cgEquipercentileEquating.feOrMFEEquipEquating(pData.weightSyntheticPopulation1,
                                                          pData.isInternalAnchor,
                                                          xv.univariateStatisticsColumn.numberOfScores,
                                                          xv.univariateStatisticsRow.numberOfScores,
                                                          xv.univariateStatisticsRow.minimumScore,
                                                          xv.univariateStatisticsRow.maximumScore,
                                                          yv.univariateStatisticsRow.numberOfScores,
                                                          yv.univariateStatisticsRow.minimumScore,
                                                          yv.univariateStatisticsRow.maximumScore,
                                                          yv.univariateStatisticsRow.scoreIncrement,
                                                          true,
                                                          bivariateLogLinearSmoothingXV.fittedBivariateRelativeFreqDistXV,
                                                          bivariateLogLinearSmoothingYV.fittedBivariateRelativeFreqDistXV,
                                                          pData.reliabilityCommonItemsPopulation1,
                                                          pData.reliabilityCommonItemsPopulation2);

        equatedRawScoreResults.relativeFreqDistsX.row(1) = cgResults.syntheticPopulationRelativeFreqDistX;
        equatedRawScoreResults.relativeFreqDistsY.row(1) = cgResults.syntheticPopulationRelativeFreqDistY;
        equatedRawScoreResults.equatedRawScores.row(0) = cgResults.equatedRawScores;
        if (cgResults.slope.has_value()) {
          equatedRawScoreResults.slope(1) = cgResults.slope.value();
        }

        if (cgResults.intercept.has_value()) {
          equatedRawScoreResults.intercept(1) = cgResults.intercept.value();
        }

        if (cgResults.braunHollandEquatedRawScores.has_value()) {
          equatedRawScoreResults.equatedRawScores.row(1) = cgResults.braunHollandEquatedRawScores.value();
        }
      }

      if (methodCode == "G" || methodCode == "A") {
        EquatingRecipes::Structures::CGEquipercentileEquatingResults cgResults =
            cgEquipercentileEquating.feOrMFEEquipEquating(pData.weightSyntheticPopulation1,
                                                          pData.isInternalAnchor,
                                                          xv.univariateStatisticsColumn.numberOfScores,
                                                          xv.univariateStatisticsRow.numberOfScores,
                                                          xv.univariateStatisticsRow.minimumScore,
                                                          xv.univariateStatisticsRow.maximumScore,
                                                          yv.univariateStatisticsRow.numberOfScores,
                                                          yv.univariateStatisticsRow.minimumScore,
                                                          yv.univariateStatisticsRow.maximumScore,
                                                          yv.univariateStatisticsRow.scoreIncrement,
                                                          true,
                                                          bivariateLogLinearSmoothingXV.fittedBivariateRelativeFreqDistXV,
                                                          bivariateLogLinearSmoothingYV.fittedBivariateRelativeFreqDistXV,
                                                          pData.reliabilityCommonItemsPopulation1,
                                                          pData.reliabilityCommonItemsPopulation2);

        // if (method == 'G' || method == 'A') /* MFE + BH-MFE in positions 2 and 3 */
        //   FEorMFE_EE(inall->w1, inall->anchor, xv->ns2, xv->ns1, xv->min1, xv->max1,
        //              yv->ns1, yv->min1, yv->max1, yv->inc1,
        //              bllxv->brfd, bllyv->brfd, inall->rv1, inall->rv2,
        //              r->fxs[1], r->gys[1], r->eraw[2],
        //              &r->a[1], &r->b[1], r->eraw[3]);

        equatedRawScoreResults.relativeFreqDistsX.row(1) = cgResults.syntheticPopulationRelativeFreqDistX;
        equatedRawScoreResults.relativeFreqDistsY.row(1) = cgResults.syntheticPopulationRelativeFreqDistY;
        equatedRawScoreResults.equatedRawScores.row(2) = cgResults.equatedRawScores;
        if (cgResults.slope.has_value()) {
          equatedRawScoreResults.slope(1) = cgResults.slope.value();
        }

        if (cgResults.intercept.has_value()) {
          equatedRawScoreResults.intercept(1) = cgResults.intercept.value();
        }

        if (cgResults.braunHollandEquatedRawScores.has_value()) {
          equatedRawScoreResults.equatedRawScores.row(3) = cgResults.braunHollandEquatedRawScores.value();
        }
      }

      /* Chained equipercentile method. Note that smoothing is bivariate
      log-linear smoothing, not univariate log-linear smoothing.
      if(method == 'C')  Chained in position 0 
      if(method == 'A')  Chained in position 4 
      if(method == 'H')  Chained in position 2 */

      size_t equatedRawScoresRowIndex;
      bool isChainedEquipercentileMethod = true;

      if (methodCode == "C") {
        equatedRawScoresRowIndex = 0;
      } else if (methodCode == "A") {
        equatedRawScoresRowIndex = 4;
      } else if (methodCode == "H") {
        equatedRawScoresRowIndex = 2;
      } else {
        isChainedEquipercentileMethod = false;
      }

      if (isChainedEquipercentileMethod) {
        Eigen::VectorXd equatedRawScores = cgEquipercentileEquating.chainedEquipercentileEquating(xv.univariateStatisticsRow.numberOfScores,
                                                                                                  bivariateLogLinearSmoothingXV.fittedRawScorePercentileRankDistX,
                                                                                                  xv.univariateStatisticsColumn.minimumScore,
                                                                                                  xv.univariateStatisticsColumn.maximumScore,
                                                                                                  xv.univariateStatisticsColumn.scoreIncrement,
                                                                                                  xv.univariateStatisticsColumn.numberOfScores,
                                                                                                  bivariateLogLinearSmoothingXV.fittedRawScoreCumulativeRelativeFreqDistV,
                                                                                                  yv.univariateStatisticsRow.minimumScore,
                                                                                                  yv.univariateStatisticsRow.scoreIncrement,
                                                                                                  yv.univariateStatisticsRow.numberOfScores,
                                                                                                  bivariateLogLinearSmoothingYV.fittedRawScoreCumulativeRelativeFreqDistX,
                                                                                                  bivariateLogLinearSmoothingXV.fittedRawScoreCumulativeRelativeFreqDistV);

        equatedRawScoreResults.equatedRawScores.row(equatedRawScoresRowIndex) = equatedRawScores;
      }

      /* get moments */

      for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
        EquatingRecipes::Structures::Moments moments =
            EquatingRecipes::ScoreStatistics::momentsFromScoreFrequencies(equatedRawScoreResults.equatedRawScores.row(methodIndex),
                                                                          xv.univariateStatisticsRow.freqDistDouble);

        equatedRawScoreResults.equatedRawScoreMoments.row(methodIndex) = moments.momentValues;
      }
    }

    /*    
      void Print_ULL(FILE *fp, char tt[], struct USTATS *x,
                    struct ULL_SMOOTH *s, int print_dm, int print_mts);
      void Print_BLL(FILE *fp, char tt[], struct BSTATS *xv, struct BLL_SMOOTH *s, 
              int print_dm, int print_mts, int print_freq, int print_bfd);
      void Print_RL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
      void Print_SL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
      void Print_CL(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
    */

  private:
    /*
      Performs univariate log-linear smoothing in terms of the 
      multinomial model as described by Holland and Thayer (1987),
      abbreviated here as H&T.  

      Input
      
        n = number of persons
        ns = number of score categories
        min = minimum raw score
        inc = increment in raw scores
        fd[] = frequency distribution
        c = number of degrees for polynomial smoothing
        scale = type of scaling:
              0 --> no scaling; 
                1 --> scale such that each column of B has
                      sum (elements) = 0 and sum (elements^2) = 1
        Btype = type of moments for criterion mathching:
                0 --> moments based on B 
                (if scale = 0, design matrix is based on raw scores,
                which means that the moments are based on raw scores;
                if scale = 1, design matrix is based on scaled raw scores,
                which means that the moments are based on scaled raw scores) 
                1 -->  moments based on B_raw, whether scale is 0 or 1
      ctype = comparison type for criterion:
              0 --> means use absolute criterion; 
          1 --> means use relative criterion
      crit = convergence criterion value
        *fp = pointer to output file
        
      Output: populates struct ULL_SMOOTH s 

      Function calls other than C or NR utilities:
        design_matrix()
        iteration()

      R. L. Brennan

      Date of last revision: 6/30/08

    */
    void smoothUnivaraiteLogLinear(const size_t& numberOfExaminees,
                                   const size_t& numberOfScores,
                                   const double& minimumScore,
                                   const double& scoreIncrement,
                                   const Eigen::VectorXd frequencyDistribution,
                                   const size_t& numberOfDegreesSmoothing,
                                   const bool& useStandardizedScale,
                                   const DesignMatrixType& Btype,
                                   const CriterionComparisonType& ctype,
                                   const double& criterion,
                                   EquatingRecipes::Structures::UnivariateLogLinearSmoothing& univariateLogLinearSmoothing) {
      size_t maximumNumberOfIterations = 40; /* maximum number of iterations */

      univariateLogLinearSmoothing.numberOfExaminees = numberOfExaminees;
      univariateLogLinearSmoothing.numberOfScores = numberOfScores;
      univariateLogLinearSmoothing.mininumRawScore = minimumScore;
      univariateLogLinearSmoothing.rawScoreIncrement = scoreIncrement;
      univariateLogLinearSmoothing.degreesOfSmoothing = numberOfDegreesSmoothing;
      univariateLogLinearSmoothing.observedFrequencies = frequencyDistribution;
      univariateLogLinearSmoothing.useScalingForBDesignMatrix = useStandardizedScale;
      univariateLogLinearSmoothing.useBRawAndCentralMoments = (Btype == DesignMatrixType::RAW_SCORE);
      univariateLogLinearSmoothing.useRelativeCriterionComparison = (ctype == CriterionComparisonType::RELATIVE);
      univariateLogLinearSmoothing.convergenceCriterion = criterion;
      univariateLogLinearSmoothing.rawScoreDesignMatrix.resize(univariateLogLinearSmoothing.numberOfScores,
                                                               univariateLogLinearSmoothing.degreesOfSmoothing);
      univariateLogLinearSmoothing.solutionDesignMatrix.resize(univariateLogLinearSmoothing.numberOfScores,
                                                               univariateLogLinearSmoothing.degreesOfSmoothing);
      univariateLogLinearSmoothing.fittedFrequencies.resize(univariateLogLinearSmoothing.numberOfScores);
      univariateLogLinearSmoothing.betaCoefficients.resize(univariateLogLinearSmoothing.degreesOfSmoothing);
      univariateLogLinearSmoothing.bObservedMoments.resize(univariateLogLinearSmoothing.degreesOfSmoothing);
      univariateLogLinearSmoothing.bFittedMoments.resize(univariateLogLinearSmoothing.degreesOfSmoothing);
      univariateLogLinearSmoothing.observedCentralMoments.resize(univariateLogLinearSmoothing.degreesOfSmoothing);
      univariateLogLinearSmoothing.fittedCentralMoments.resize(univariateLogLinearSmoothing.degreesOfSmoothing);

      univariateLogLinearSmoothing.likelihoodRatioChiSquare = 0.0;
      univariateLogLinearSmoothing.numberOfZeros = 0;

      univariateLogLinearSmoothing.fittedRawScoreDist.resize(univariateLogLinearSmoothing.numberOfScores);
      univariateLogLinearSmoothing.fittedRawScoreCumulativeRelativeDist.resize(univariateLogLinearSmoothing.numberOfScores);
      univariateLogLinearSmoothing.fittedRawScorePercentileRankDist.resize(univariateLogLinearSmoothing.numberOfScores);

      Eigen::MatrixXi crossProductMomentDesignations;

      designMatrix(univariateLogLinearSmoothing.numberOfScores,
                   univariateLogLinearSmoothing.mininumRawScore,
                   univariateLogLinearSmoothing.rawScoreIncrement,
                   0,
                   0,
                   0,
                   univariateLogLinearSmoothing.degreesOfSmoothing,
                   0,
                   0,
                   crossProductMomentDesignations,
                   univariateLogLinearSmoothing.useScalingForBDesignMatrix,
                   univariateLogLinearSmoothing.rawScoreDesignMatrix,
                   univariateLogLinearSmoothing.solutionDesignMatrix);

      /* iteration-step results not printed if first parameter is NULL;
      results are printed if first parameter is fp!=NULL. Note that fp
	    can be set to NULL in Wrapper_Smooth_ULL() */

      std::optional<Eigen::VectorXd> uConstants;

      univariateLogLinearSmoothing.numberOfIterations = iteration(
          univariateLogLinearSmoothing.solutionDesignMatrix,
          univariateLogLinearSmoothing.rawScoreDesignMatrix,
          frequencyDistribution,
          uConstants,
          numberOfDegreesSmoothing,
          0,
          0,
          crossProductMomentDesignations,
          maximumNumberOfIterations,
          ctype,
          Btype,
          criterion,
          univariateLogLinearSmoothing.betaCoefficients,
          univariateLogLinearSmoothing.fittedFrequencies,
          univariateLogLinearSmoothing.bObservedMoments,
          univariateLogLinearSmoothing.bFittedMoments,
          univariateLogLinearSmoothing.observedCentralMoments,
          univariateLogLinearSmoothing.fittedCentralMoments,
          univariateLogLinearSmoothing.likelihoodRatioChiSquare,
          univariateLogLinearSmoothing.numberOfZeros,
          univariateLogLinearSmoothing.cllNormalizingConstant,
          false);

      univariateLogLinearSmoothing.fittedRawScoreDist = univariateLogLinearSmoothing.fittedFrequencies /
                                                        static_cast<double>(univariateLogLinearSmoothing.numberOfExaminees);

      univariateLogLinearSmoothing.fittedRawScoreCumulativeRelativeDist = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(0,
                                                                                                                                 univariateLogLinearSmoothing.numberOfScores - 1,
                                                                                                                                 1,
                                                                                                                                 univariateLogLinearSmoothing.fittedRawScoreDist);

      univariateLogLinearSmoothing.fittedRawScorePercentileRankDist = EquatingRecipes::Utilities::percentileRanks(0,
                                                                                                                  univariateLogLinearSmoothing.numberOfScores - 1,
                                                                                                                  1,
                                                                                                                  univariateLogLinearSmoothing.fittedRawScoreCumulativeRelativeDist);
    }

    /*
      Create design matrix.  Code and variable names are for a
      bivariate u*v distribution, where u is rows and v is columns.
      Note that x = u + v = total score, with v = common-item score 
      and u = non-common-item score.  We never create a design matrix
      with x and v where v is internal to x, because there would be
      a great deal of collinearity. 
      
      Input
        nsu = number of score categories for u
        minu = minimum score for u
        incu = increment for u
        nsv = number of score categories for v
        minv = minimum score for v
        incv = increment for v
        cu = number of degrees of smoothing for u
        cv = number of degrees of smoothing for v
        cuv = number of cross-product moments
        cpm[cuv-1][2] = zero-offset matrix designating cross-
                        product moments.  Example: let cuv = 3,
                        and the desired cross-product moments be 
                        (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                        Then cpm[0] = 1,1; cpm[1] = 1,2; and
                        cpm[2] = 2,1.  
        scale: 0 --> no scaling; 
              1 --> scale such that each column of B has
                    sum (elements) = 0 and sum (elements^2) = 1
                    
      Output
        B_raw[][] = zero-offset design matrix for raw scores; 
                    #rows = (nsu*nsv) and #cols = (cu+cv+cuv);
                    space already allocated for B_raw 
        B[][]     = zero-offset design matrix used for solution;
                    involves scaling if scale==1 
                    #rows = (nsu*nsv) and #cols = (cu+cv+cuv);
                    space already allocated for B 
        
      NOTE: For univariate smoothing, set 
            nsv=0, cv=0, cuv=0, cpm = NULL.  In this case, 
            obviously, u plays a generic role 
            (i.e., any single variable such as x or y)
            
      NOTE: Usually for bivariate smoothing in equating, set cuv=1
            and cpm[0] = 1,1. Otherwise, the conditional 
            distributions are not necessarily stochastically ordered
            (see Rosenbaum & Thayer, 1987, p. 46).

      NOTE: Using minu!=0 and incu!=1 changes both the B and 
            B_raw matrices; similarly for minv and incv.  
            If convergence not achieved, it may be wise to set
            minu = 0 and incu = 1 (same for minv and incv)

      Function calls other than C or NR utilities:  
        score()
      runerror()

      R. L. Brennan

      Date of last revision: 6/30/08

    */
    void designMatrix(const size_t& numberOfScoresU,
                      const double& minimumScoreU,
                      const double& scoreIncrementU,
                      const size_t& numberOfScoresV,
                      const double& minimumScoreV,
                      const double& scoreIncrementV,
                      const size_t& numberOfDegreesOfSmoothingU,
                      const size_t& numberOfDegreesOfSmoothingV,
                      const size_t& numberOfCrossProductMoments,
                      const Eigen::MatrixXi& crossProductMomentDesignations,
                      const bool& useStandardizedScale,
                      Eigen::MatrixXd& rawScoreDesignMatrix,
                      Eigen::MatrixXd& solutionDesignMatrix) {
      size_t numberOfDesignMatrixColumns = numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV + numberOfCrossProductMoments; /* # columns in design matrix */
      size_t numberOfDesigmMatrixRows = (numberOfScoresV > 0) ? numberOfScoresU * numberOfScoresV : numberOfScoresU;                /* # rows in design mat */

      /* univariate LL smoothing */
      if (numberOfScoresV == 0) {
        for (size_t scoreLocationU = 0; scoreLocationU < numberOfScoresU; scoreLocationU++) {
          double scoreU = EquatingRecipes::Utilities::getScore(scoreLocationU, minimumScoreU, scoreIncrementU);

          for (size_t degreesOfSmoothingUIndex = 1; degreesOfSmoothingUIndex <= numberOfDegreesOfSmoothingU; degreesOfSmoothingUIndex++) {
            rawScoreDesignMatrix(scoreLocationU, degreesOfSmoothingUIndex - 1) = std::pow(scoreU, static_cast<double>(degreesOfSmoothingUIndex));
            solutionDesignMatrix(scoreLocationU, degreesOfSmoothingUIndex - 1) = rawScoreDesignMatrix(scoreLocationU, degreesOfSmoothingUIndex - 1);
          }
        }
      } else if (numberOfDegreesOfSmoothingV == 0) {
        throw std::runtime_error("Number of degrees of smoothing for V or number of score categories for V is misspecified.");
      } else {
        /* bivariate LL smoothing */

        /* u polynomials */
        size_t cell = 0;
        for (size_t scoreLocationU; scoreLocationU < numberOfScoresU; scoreLocationU++) {
          double scoreU = EquatingRecipes::Utilities::getScore(scoreLocationU, minimumScoreU, scoreIncrementU);

          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            for (size_t degreesOfSmoothingUIndex = 1; degreesOfSmoothingUIndex <= numberOfDegreesOfSmoothingU; degreesOfSmoothingUIndex++) {
              rawScoreDesignMatrix(cell, degreesOfSmoothingUIndex - 1) = std::pow(scoreU, static_cast<double>(degreesOfSmoothingUIndex));
              cell++;
            }
          }

          /* v polynomials */
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            double scoreV = EquatingRecipes::Utilities::getScore(scoreLocationV,
                                                                 minimumScoreV,
                                                                 scoreIncrementV);

            for (size_t scoreLocationU = 0; scoreLocationU < numberOfScoresU; scoreLocationU++) {
              for (size_t degreesOfSmoothingVIndex = 1; degreesOfSmoothingVIndex <= numberOfDegreesOfSmoothingV; degreesOfSmoothingVIndex++) {
                rawScoreDesignMatrix(scoreLocationV + scoreLocationU * numberOfScoresV,
                                     numberOfDegreesOfSmoothingU + scoreLocationV - 1) = std::pow(scoreV, static_cast<double>(degreesOfSmoothingVIndex));
              }
            }
          }

          if (numberOfDegreesOfSmoothingV > 0 && (crossProductMomentDesignations.rows() == 0 || crossProductMomentDesignations.cols() == 0)) {
            throw std::runtime_error("Number of degrees of smoothing for V or cross-product moment designations misspecified.");
          }
        }

        /* cross-product polynomials */
        cell = 0;

        for (size_t scoreLocationU = 0; scoreLocationU < numberOfScoresU; scoreLocationU++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            for (size_t degreesOfSmoothingVIndex = 0; degreesOfSmoothingVIndex < numberOfDegreesOfSmoothingV; degreesOfSmoothingVIndex++) {
              rawScoreDesignMatrix(cell, numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV + degreesOfSmoothingVIndex) = rawScoreDesignMatrix(cell, crossProductMomentDesignations(degreesOfSmoothingVIndex, 0) - 1) *
                                                                                                                                 rawScoreDesignMatrix(cell, numberOfDegreesOfSmoothingU + crossProductMomentDesignations(degreesOfSmoothingVIndex, 1) - 1);

              cell++;
            }
          }
        }

        /* copy B_raw to B */
        for (size_t rowIndex = 0; rowIndex < numberOfDesigmMatrixRows; rowIndex++) {
          for (size_t columnIndex = 0; columnIndex < numberOfDesignMatrixColumns; columnIndex++) {
            solutionDesignMatrix(rowIndex, columnIndex) = rawScoreDesignMatrix(rowIndex, columnIndex);
          }
        }
      }

      /* scaling such that for each column of B,
        sum (elements) = 0 and sum (elements^2) = 1 */

      if (useStandardizedScale) {
        Eigen::VectorXd columnMeans = rawScoreDesignMatrix.colwise().mean();
        Eigen::VectorXd columnSumSquaredDeviations = (rawScoreDesignMatrix.cwiseProduct(rawScoreDesignMatrix)).colwise().sum() -
                                                     static_cast<double>(numberOfDesignMatrixColumns) * columnMeans.cwiseProduct(columnMeans);

        for (size_t columnIndex = 0; columnIndex < numberOfDesignMatrixColumns; columnIndex++) {
          Eigen::VectorXd vectorOfMean = Eigen::VectorXd::Constant(numberOfDesigmMatrixRows, columnMeans(columnIndex));

          rawScoreDesignMatrix.col(columnIndex) = (rawScoreDesignMatrix.col(columnIndex) - vectorOfMean) / std::sqrt(columnSumSquaredDeviations(columnIndex));
        }
      }
    }

    /*
      Get nct[] from bfd[][]

      Convert bivariate fd (xv->bfd[][]) to a vector nct[] with row j
      elements followed by row j+1 elements.  Note that for an 
      internal anchor the rows are x = u + v and the cols are v,
      which means that there are structural zeros, and we want nct[]
      to contain only the u*v elements.  See example below in which
      - indicates a structural 0 and the within-matrix numbers are 
      the cell locations in nct[].

                    v                                   v
                0  1  2  3                          0  1  2  3
          ---------------                      --------------
          0 |  0  -  -  -                      0|  0  1  2  3
          1 |  4  1  -  -    nsx = 9           1|  4  5  6  7
          2 |  8  5  2  -    nsv = 4 -->     u 2|  8  9 10 11  
          3 | 12  9  6  3    nsu = 6           3| 12 13 14 15
        x 4 | 16 13 10  7                      4| 16 17 18 19
          5 | 20 17 14 11                      5| 20 21 22 23
          6 |  - 21 18 15
          7 |  -  - 22 19
          8 |  -  -  - 23
    

      Input
        anchor : 0 --> external; 1 --> internal
        nsx = number of score categories for total scores (x)
        nsv = number of score categories for common-item scores (v)
        bfd[][] = bivariate freq dist for x and v

      Output
        nct[] = vector version of bfd[][], where nct[] is "collaped",
                as discussed above, if anchor is internal.
                Assumes space already allocated for nct[]
    
      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08      
    */
    void getNctBfd(const bool& isInternalAnchor,
                   const size_t& numberOfScoresX,
                   const size_t& numberOfScoresV,
                   const Eigen::MatrixXd& bivariateFreqDist,
                   Eigen::VectorXd& nct) {
      if (!isInternalAnchor) {
        /* external anchor */

        size_t cells = 0;
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            nct(cells++) = bivariateFreqDist(scoreLocationX, scoreLocationV);
          }
        }
      } else {
        /* internal anchor */

        size_t numberOfScoresU = numberOfScoresX - numberOfScoresV + 1; /* number of score categories for non-common items */

        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            if (scoreLocationX < scoreLocationV || scoreLocationX > numberOfScoresU - 1 + scoreLocationV) {
              continue;
            } else {
              nct((scoreLocationX - scoreLocationV) * numberOfScoresV + scoreLocationV) = bivariateFreqDist(scoreLocationX, scoreLocationV);
            }
          }
        }
      }
    }

    /*
      Get bfd[][] from mct[]

      For an external anchor, directly convert the row major 
      vector mct[nsu*nsv] to bfd[nsu][nsv]. For an 
      internal anchor, convert mct[nsu*nsv] to
      bfd[nsu+nsv][nsv] be adding structural zeros.
      See comments for get_nct_bfd().

      Input
        anchor : 0 --> external; 1 --> internal
        nsx = number of score categories for total scores (x)
        nsv = number of score categories for common-item scores (v)
        mct[] = row-major vector for fitted bivariate frequencies
                for non-common items by common items

      Output
        bfd[][] = fitted bivariate frequencies for total scores by
                  common-item scores. That is mct[] is mapped into 
                  bfd[][] with structural zeros added if anchor is 
                  internal (see comments for get_nct_bfd()). 
                  Assumes space already allocated for bfd[][]
    
      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void getBfdMct(const bool& isInternalAnchor,
                   const size_t& numberOfScoresX,
                   const size_t& numberOfScoresV,
                   const Eigen::VectorXd& mct,
                   Eigen::MatrixXd& bivariateFreqDist) {
      if (!isInternalAnchor) {
        /* external anchor */

        size_t cells = 0;
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            bivariateFreqDist(scoreLocationX, scoreLocationV) = mct(cells++);
          }
        }
      } else {
        /* internal anchor */

        size_t numberOfScoresU = numberOfScoresX - numberOfScoresV + 1; /* number of score categories for non-common items */
        for (size_t scoreLocationX = 0; scoreLocationX < numberOfScoresX; scoreLocationX++) {
          for (size_t scoreLocationV = 0; scoreLocationV < numberOfScoresV; scoreLocationV++) {
            if (scoreLocationX < scoreLocationV || scoreLocationX > numberOfScoresU - 1 + scoreLocationV) {
              bivariateFreqDist(scoreLocationX, scoreLocationV) = 0.0;
            } else {
              bivariateFreqDist(scoreLocationX, scoreLocationV) = mct((scoreLocationX - scoreLocationV) * numberOfScoresV + scoreLocationV);
            }
          }
        }
      }
    }

    /*
      Input:
        B[][] = design matrix (ns x nc)
        m[]   = fitted frequencies (ns x 1)

      Output:
        BtSmB[][] = Bt x Sm x B (nc x nc)
                  = minus the 2nd derivative of log-likelihood
                  (Eq. 22 and 32 in Holland & Thayer, 1987)

      Removed From Input:
        ns    = number of score categories (rows in design matrix)
        nc    = number of columns in design matrix
        N     = total of all frequencies

      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    Eigen::MatrixXd getBtSmB(const Eigen::MatrixXd& designMatrix,
                             const Eigen::VectorXd& fittedFrequencies) {
      Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(designMatrix.cols(),
                                                      designMatrix.cols());

      Eigen::VectorXd weightedSumDesignMatrixColumns =
          designMatrix.cwiseProduct(fittedFrequencies.replicate(designMatrix.cols(), 1)).colwise().sum();

      double fittedFrequencySum = fittedFrequencies.sum();

      for (size_t rowIndex = 0; rowIndex < designMatrix.cols(); rowIndex++) {
        for (size_t columnIndex = rowIndex; columnIndex < designMatrix.cols(); columnIndex++) {
          hessian(rowIndex, columnIndex) =
              (designMatrix.col(rowIndex) * designMatrix.col(columnIndex) * fittedFrequencies).sum();

          hessian(rowIndex, columnIndex) /= fittedFrequencySum;
          hessian(columnIndex, rowIndex) = hessian(rowIndex, columnIndex);
        }
      }

      return hessian;
    }

    /*
      Input:
        B[][] = design matrix (ns x nc)
        n[]   = actual frequencies (ns x 1)
        m[]   = fitted frequencies (ns x 1)

      Output:
        Btnm[] = 1st derivative of log-likelihood (nc x 1)
                (Eq. 19 and 33 in Holland & Thayer, 1987)

      Removed From Input:
        ns    = number of score categories (rows in design matrix)
        nc    = number of columns in design matrix

      Function calls other than C or NR utilities: None 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    Eigen::VectorXd getBtnm(const Eigen::MatrixXd& designMatrix,
                            const Eigen::VectorXd& observedFrequencies,
                            const Eigen::VectorXd& fittedFrequencies) {
      Eigen::VectorXd gradient(designMatrix.cols());

      for (size_t columnIndex = 0; columnIndex < designMatrix.cols(); columnIndex++) {
        gradient(columnIndex) = (designMatrix.col(columnIndex) * (observedFrequencies - fittedFrequencies)).sum();
      }

      return gradient;
    }

    /*
      Input:
        B[][] = zero-offset design matrix (ns x nc) 
        n[]   = zero-offset actual frequencies (ns x 1)
        
        fp    = output file pointer for debugging
                (NULL --> no output)

      Output:
        Beta0[] = zero-offset initial values of Beta (nc x 1);
                  space already allocated;
                  based on Eq 49 (and next line) of Holland and
                    Thayer (1987) -- abbreviated H&T below

      Removed From Input:
        N     = total of frequencies
        ns    = number of score categories (rows in design matrix)
        nc    = number of columns in design matrix

      Function calls other than C or NR utilities: 
        Print_vector()
        Print_matrix() 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    Eigen::VectorXd getBeta0(const Eigen::MatrixXd& designMatrix,
                             const Eigen::VectorXd& observedFrequencies,
                             const bool& debug) {
      size_t numberOfRows = designMatrix.rows();
      size_t numberOfColumns = designMatrix.cols();
      double observedFrequenciesSum = observedFrequencies.sum();

      /* get a; 0.8 can be changed to any value in (0,1) */
      double rho = 0.8;

      Eigen::VectorXd a = (rho * observedFrequencies) + Eigen::VectorXd::Constant(numberOfRows, observedFrequenciesSum / static_cast<double>(numberOfRows));

      if (debug) {
        std::cout << fmt::format("a debug:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(a, false));
      }

      /* get BtSaB -- first term on left side of Eq 49 in H&T */
      Eigen::MatrixXd BtSaB = getBtSmB(designMatrix, a);

      if (debug) {
        std::cout << fmt::format("BtSaB debug:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(BtSaB));
      }

      /*  get right side of Eq 49 in H&T, which is computed using
      Equation 38 in Holland and Thayer (2000) with all mu = 0 */
      double aLogA = a.array().log().sum();

      Eigen::VectorXd BtSaloga(numberOfColumns);

      for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
        double bALogA = (designMatrix.col(columnIndex).array() * a.array() * a.array().log()).sum();
        double bA = (designMatrix.col(columnIndex) * a).sum();

        BtSaloga(columnIndex) = bALogA - (bA * aLogA / observedFrequenciesSum);
      }

      /* get Beta0 using NR ludcmp() and lubksb() */

      /* right[] is one-offset B in Ax=B */
      /* left[][] is one-offset A in Ax=B */
      Eigen::MatrixXd left = Eigen::MatrixXd::Zero(numberOfColumns + 1, numberOfColumns + 1);
      Eigen::VectorXd right = Eigen::VectorXd::Zero(numberOfColumns + 1);

      right(Eigen::seq(1, numberOfColumns)) = BtSaloga;

      left.block(1, 1, numberOfColumns, numberOfColumns) = BtSaB;

      Eigen::VectorXd solution = left.partialPivLu().solve(right);

      /* Beta0 is zero-offset solution */
      Eigen::VectorXd beta0(numberOfColumns + 1);
      beta0(Eigen::seq(1, numberOfColumns + 1)) = solution;

      return beta0;
    }

    /*
      Get m using Equation 8 in Holland and Thayer (1987)

      Input
        B[][]  = design matrix (ns x nc)
        Beta[] = parameter estimates (nc x 1)
        uin[]  = u constants (ns x 1); if NULL, set elements to 0 
        N      = total of frequencies
        ns     = number of rows of design matrix
        nc     = number of columns of design matrix
        fp     = output file pointer for debugging
                (NULL --> no output)

      Output
        m[]    = estimated frequencies (space already allocated)

      NOTE: DBL_MIN is defined in <float.h>.  It is the minimum normalized 
          floating point number. For Visual Studio DBL_MIN = 2.225074 E-308
        log(base e) of DBL_MIN in Visual Studio is -708.3964

      Return ap = alpha' --- see top of p. 3 in H&T

      Function calls other than C or NR utilities: 
        Print_vector() 

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    double getEstimatedFrequences(const Eigen::MatrixXd& designMatrix,
                                  const Eigen::VectorXd& betaParameterEstimates,
                                  const std::optional<Eigen::VectorXd>& uConstants,
                                  const double& sumOfFrequencies,
                                  const bool& debug,
                                  Eigen::VectorXd& estimatedFrequencies) {
      size_t numberOfRows = designMatrix.rows();

      Eigen::VectorXd u = uConstants.value_or(Eigen::VectorXd::Zero(numberOfRows));
      Eigen::VectorXd BBeta = designMatrix * betaParameterEstimates;

      if (debug) {
        std::cout << fmt::format("BBeta debug:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(BBeta, false));
      }

      double normalizingConstant = 0.0; // = ap
      for (size_t rowIndex = 0; rowIndex < numberOfRows; rowIndex++) {
        double uPlusBBeta = u(rowIndex) + BBeta(rowIndex);

        normalizingConstant += uPlusBBeta < std::numeric_limits<double>::min() ? 0.0 : std::exp(uPlusBBeta);
      }

      if (normalizingConstant < std::numeric_limits<double>::min()) {
        normalizingConstant = std::log(sumOfFrequencies);
      } else {
        normalizingConstant = std::log(sumOfFrequencies) - std::log(normalizingConstant);
      }

      if (debug) {
        std::cout << fmt::format("\n\nap = {:12.5f}\n", normalizingConstant);
      }

      estimatedFrequencies.resize(numberOfRows);

      for (size_t rowIndex = 0; rowIndex < numberOfRows; rowIndex++) {
        double z = normalizingConstant + u(rowIndex) + BBeta(rowIndex);

        if (z < std::numeric_limits<double>::min()) {
          estimatedFrequencies[rowIndex] = 0.0;
        } else {
          estimatedFrequencies[rowIndex] = std::exp(z);
        }
      }

      return normalizingConstant; /* ap is the ; */
    }

    /* Iterate to a solution for log-linear models under a mutinomial 
      distribution as discussed in Holland and Thayer (1987), abbreviated
      H&T here.

      This function is called from Smooth_ULL() and Smooth_BLL()

      Input
        fp      = file pointer for output
                  (if fp==NULL) no output written
        B[][]   = design matrix (ns x nc) where nc = cu + cv + cuv
                  used for solution
        B_raw[][]   = design matrix for raw scores;
                      used to get central moments               
        nct[]   = frequencies (ns x 1)
        N       = total of frequencies
        uin[]   = constant vector (see Eq 8 in H&T);
                  if NULL, all elements set to 0
        ns      = number of frequencies (rows in design matrix)
        cu      = number of degrees of smoothing for u
        cv      = number of degrees of smoothing for v
        cuv     = number of cross-product moments
        cpm[cuv-1][2] = zero-offset matrix designating cross-
                        product moments.  Example: let cuv = 3,
                        and the desired cross-product moments be 
                        (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                        Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                        cpm[2] = {2,1}. 
        max_nit = maximum number of iterations
        ctype = comparison type for criterion:
                0 --> absolute; 1 --> relative
        Btype = type of design matrix and, hence, type of moments 
                for criterion:
                0 --> use B (could be scaled or unscaled as indicated
                      in design_matrix()) --- see note below
                1 --> use B_raw and central moments based on it
        crit = criterion.  See crit_mts() for discussion 
        
      Output
        Beta[]  = coefficients for variables (columns) of B (nc x 1)
        mct[]   = fitted frequencies (ns x 1) 
        n_mts[] = actual moments based on B (nc x 1)
        m_mts[] = fitted moments based on B (nc x 1)
        n_mts_raw[] = actual central moments based on B_raw (nc x 1)
        m_mts_raw[] = fitted central moments based on B_raw (nc x 1)  
        *lrc        = likelihood-ratio chi-square; see Agresti, 2007, p. 37 
        *nzero      = number of times that nct[i]/mct[i] involves 
                      at least one zero frequency; df correction for *lrc  
      *ap     = alpha (used as normalizing constant in CLL---added by TW)

        NOTE: space already allocated for 
              Beta[], mct[], n_mts[], m_mts[], n_mts_raw, m_mts_raw;
              lrc and nzero declared in calling function

        NOTE: if scale==0 for design_matrix(),
              then n_mts and m_mts are non-central moments;
              if scale==1 for design_matrix(),
              then n_mts and m_mts are analogous to central moments

      Return: number of iterations to convergence 

      Function calls other than C or NR utilities: 
        get_LLmoments()
        Print_iteration_heading() 
        get_Beta0()
        get_mct()
        crit_mts()
        get_BtSmB()
        ludcmp() from NR.c
      lubksb() from NR.c
        get_Btnm()
        mmult_a_v()
        runerror()

      R. L. Brennan

      Date of last revision: 6/30/08

    */
    size_t iteration(const Eigen::MatrixXd& solutionDesignMatrix,
                     const Eigen::MatrixXd& rawScoreDesigmMatrix,
                     const Eigen::VectorXd& frequencies,
                     const std::optional<Eigen::VectorXd>& uConstants,
                     const size_t& numberOfDegreesOfSmoothingU,
                     const size_t& numberOfDegreesOfSmoothingV,
                     const size_t& numberOfCrossProductMoments,
                     const Eigen::MatrixXi& crossProductMomentDesignations,
                     const size_t& maximumIterations,
                     const CriterionComparisonType& criterionComparisonType,
                     const DesignMatrixType& designMatrixType,
                     const double& criterion,
                     Eigen::VectorXd& betaParameterEstimates,
                     Eigen::VectorXd& fittedFrequencies,
                     Eigen::VectorXd& observedMoments,
                     Eigen::VectorXd& fittedMoments,
                     Eigen::VectorXd& observedCentralMoments,
                     Eigen::VectorXd& fittedCentralMoments,
                     double& likelihoodRatioChiSquare,
                     size_t& numberOfZeroFrequencies,
                     double& normalizingConstant,
                     const bool& debug) {
      size_t numberOfRows = solutionDesignMatrix.rows();

      /* number of columns of design matrix */
      size_t numberOfColumns = numberOfDegreesOfSmoothingU +
                               numberOfDegreesOfSmoothingV +
                               numberOfCrossProductMoments;

      size_t iterationNumber;
      size_t numberOfCellsToPrint = std::min(numberOfRows, 250UL);

      getLogLinearMoments(solutionDesignMatrix,
                          rawScoreDesigmMatrix,
                          frequencies,
                          numberOfDegreesOfSmoothingU,
                          numberOfDegreesOfSmoothingV,
                          numberOfCrossProductMoments,
                          crossProductMomentDesignations,
                          observedMoments,
                          observedCentralMoments);

      if (debug) {
        PrintIterationHeading(numberOfRows,
                              numberOfColumns,
                              frequencies,
                              observedMoments,
                              observedCentralMoments,
                              criterionComparisonType,
                              designMatrixType,
                              criterion);
      }

      betaParameterEstimates = getBeta0(solutionDesignMatrix,
                                        frequencies,
                                        debug);

      /***** begin iteration loop *****/
      for (iterationNumber = 0; iterationNumber <= maximumIterations; iterationNumber++) {
        likelihoodRatioChiSquare = 0.0;
        numberOfZeroFrequencies = 0;

        /* get ap (the normalizing constant which is needed in CLL) and fitted frequencies (mct) */
        normalizingConstant = getEstimatedFrequences(solutionDesignMatrix,
                                                     betaParameterEstimates,
                                                     uConstants,
                                                     frequencies.sum(),
                                                     debug,
                                                     fittedFrequencies);

        for (size_t rowIndex = 0; rowIndex < numberOfRows; rowIndex++) {
          if (frequencies(rowIndex) != 0.0 && fittedFrequencies(rowIndex) != 0.0) {
            likelihoodRatioChiSquare += frequencies(rowIndex) * std::log(frequencies(rowIndex) / fittedFrequencies(rowIndex));
          } else {
            numberOfZeroFrequencies++;
          }
        }

        likelihoodRatioChiSquare *= 2.0;

        getLogLinearMoments(solutionDesignMatrix,
                            rawScoreDesigmMatrix,
                            fittedFrequencies,
                            numberOfDegreesOfSmoothingU,
                            numberOfDegreesOfSmoothingV,
                            numberOfCrossProductMoments,
                            crossProductMomentDesignations,
                            fittedMoments,
                            fittedCentralMoments);

        if (debug) {
          std::cout << fmt::format("\n        {:2d}    ", iterationNumber);
          std::cout << fmt::format("  {10.5f}", normalizingConstant);
          for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
            std::cout << fmt::format("  {10.5f}\n", uConstants.has_value() ? (uConstants.value())(columnIndex) : 0.0);
          }

          for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
            std::cout << fmt::format("{:12.3f}\n", betaParameterEstimates(columnIndex));
          }

          for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
            std::cout << fmt::format("{:12.5f}\n", fittedMoments(columnIndex));
          }

          for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
            std::cout << fmt::format("{:12.5f}\n", fittedCentralMoments(columnIndex));
          }

          std::cout << fmt::format("{:12.5f}\n", likelihoodRatioChiSquare);

          std::cout << fmt::format("{:7d}\n", numberOfZeroFrequencies);

          for (size_t cellIndex = 0; cellIndex < numberOfCellsToPrint; cellIndex++) {
            std::cout << fmt::format("{:12.5f}\n", fittedFrequencies(cellIndex));
          }

          if (numberOfCellsToPrint > 250) {
            std::cout << "  ...\n";
          }
        }

        if (momentsCriterion(numberOfColumns,
                             numberOfDegreesOfSmoothingU,
                             criterionComparisonType,
                             designMatrixType,
                             (designMatrixType == DesignMatrixType::SOLUITON) ? observedMoments : observedCentralMoments,
                             (designMatrixType == DesignMatrixType::SOLUITON) ? fittedMoments : fittedCentralMoments,
                             criterion)) {
          break;
        }

        /* Bt x Sm x B --- see Eq 30 and 32 in H&T */
        Eigen::MatrixXd BtSmB = getBtSmB(solutionDesignMatrix, fittedFrequencies);

        /* Bt x (n - m) --- see Eq 30 and 33 in H&T */
        Eigen::VectorXd Btnm = getBtnm(solutionDesignMatrix,
                                       frequencies,
                                       fittedFrequencies);

        Eigen::MatrixXd left(numberOfColumns + 1, numberOfColumns + 1); /* left[1...nc][1...nc] set to BtSmB[0...nc-1][0...nc-1] */
        Eigen::VectorXd right(numberOfColumns + 1);                     /* right[1..nc] set to Btnm[0...nc-1] */

        /* solve for delta^r in ER Equation 10.14 == H&T Equation 30*/
        left(Eigen::seq(1, numberOfColumns), Eigen::seq(1, numberOfColumns)) = BtSmB; /* left[][] is one-offset A in Ax=B */
        right(Eigen::seq(1, numberOfColumns)) = Btnm;                                 /* right[] is one-offset B in Ax=B */

        Eigen::VectorXd solution = left.partialPivLu().solve(right);

        Eigen::VectorXd delta(numberOfColumns + 1); /* delta[] is zero-offset version of right[], see Eq 31 in H&T */

        delta = right(Eigen::seq(1, numberOfColumns));

        betaParameterEstimates += delta;
      }

      if (iterationNumber > maximumIterations) {
        throw std::runtime_error("Criterion not satisfied after max_nit iterations");
      }

      return iterationNumber;
    }

    /*
      Get moments based on B and B_raw design matrices for log-linear model

      Input
        B[][] = design matrix (ns x nc)
                could be based on scale = 0 or scale = 1 as specified
                in design_matrix()
        B_raw[][] = design matrix (ns x nc) for raw scores
                    (e.g., first column is min, min+inc, ...)
        f[]   = frequencies (ns x 1) (could be nct[] or mct[])
        ns    = number of frequencies (rows in design matrix)
        cu    = number of degrees of smoothing for u
        cv    = number of degrees of smoothing for v
        cuv   = number of cross-product moments
        cpm[cuv-1][2] = zero-offset matrix designating cross-
                        product moments.  Example: let cuv = 3,
                        and the desired cross-product moments be 
                        (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                        Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                        cpm[2] = {2,1}.  
      Output
        mts[] = moments (nc x 1) based on B (space already allocated)
        mts_raw[] = central moments (nc x 1) based on B_raw 
                    (space already allocated);
                    e.g., for j = 2,...(cu-1)), the (j+1)th central moment is 
                          mts_raw[j] = 
                          sum_i(B_raw[i][j] - xbar)^(j+1)*f[i]/sd^((j+1)),
                          where xbar = sum_i(B_raw[i][0]*f[i]).
                          Note that mts_raw[0] = 0 by construction, 
                          and we replace it with xbar 
                          (i.e., we set mts_raw[0] = xbar);
                          simlarly, we set mts_raw[1] = sd

      Function calls other than C or NR utilities: none

      R. L. Brennan

      Date of last revision: 6/30/08
    */
    void getLogLinearMoments(const Eigen::MatrixXd& solutionDesignMatrix,
                             const Eigen::MatrixXd& rawScoreDesigmMatrix,
                             const Eigen::VectorXd& frequencies,
                             const size_t& numberOfDegreesOfSmoothingU,
                             const size_t& numberOfDegreesOfSmoothingV,
                             const size_t& numberOfCrossProductMoments,
                             const Eigen::MatrixXi& crossProductMomentDesignations,
                             Eigen::VectorXd& fittedMoments,
                             Eigen::VectorXd& fittedCentralMoments) {
      // int i, j,
      //     nc = cu + cv + cuv;
      // double mnu, sdu, mnv, sdv,
      //     *rf; /* relative frequency */

      size_t numberOfColumns = numberOfDegreesOfSmoothingU +
                               numberOfDegreesOfSmoothingV +
                               numberOfCrossProductMoments; /* number of columns in design matrix */

      Eigen::VectorXd relativeFrequencies = frequencies / frequencies.sum();

      fittedMoments.setZero(numberOfColumns);
      fittedCentralMoments.setZero(numberOfColumns);

      /*** moments based on B ***/
      fittedMoments = solutionDesignMatrix * relativeFrequencies;

      /*** "typical" central moments based on B_raw ***/
      /* for cu columns associated with u;  scores are in column 0 */
      double meanU = rawScoreDesigmMatrix.col(0).cwiseProduct(relativeFrequencies).sum();

      double sdU = (rawScoreDesigmMatrix.col(0).array().square()).cwiseProduct(relativeFrequencies.array()).sum();
      sdU = std::sqrt(sdU - std::pow(meanU, 2));

      if (numberOfDegreesOfSmoothingU >= 1) {
        fittedCentralMoments(0) = meanU;
      }

      if (numberOfDegreesOfSmoothingU >= 2) {
        fittedCentralMoments(1) = sdU;
      }

      Eigen::VectorXd deviation = rawScoreDesigmMatrix.col(0) -
                                  Eigen::VectorXd::Constant(rawScoreDesigmMatrix.rows(), meanU);

      for (size_t momentIndex = 2; momentIndex < numberOfDegreesOfSmoothingU; momentIndex++) {
        fittedCentralMoments(momentIndex) += ((deviation.array().pow(momentIndex + 1)) * relativeFrequencies.array()).sum();
        fittedCentralMoments(momentIndex) /= std::pow(sdU, static_cast<double>(momentIndex + 1));
      }

      // /* for cv columns associated with v; scores are in column cu*/
      double meanV = 0.0;
      double sdV = 0.0;
      if (numberOfDegreesOfSmoothingV > 0) {
        meanV = rawScoreDesigmMatrix.col(numberOfDegreesOfSmoothingU).cwiseProduct(relativeFrequencies).sum();
        sdV = (rawScoreDesigmMatrix.col(numberOfDegreesOfSmoothingU).array().pow(2)).cwiseProduct(relativeFrequencies.array()).sum();

        sdV = std::sqrt(sdV - std::pow(meanV, static_cast<double>(2)));

        if (numberOfDegreesOfSmoothingV >= 1) {
          fittedCentralMoments(numberOfDegreesOfSmoothingU) = meanV;
        }

        if (numberOfDegreesOfSmoothingV >= 2) {
          fittedCentralMoments(numberOfDegreesOfSmoothingU + 1) = sdV;
        }

        Eigen::VectorXd deviations = rawScoreDesigmMatrix.col(numberOfDegreesOfSmoothingU) - Eigen::VectorXd::Constant(rawScoreDesigmMatrix.rows(), meanV);

        for (size_t momentIndex = numberOfDegreesOfSmoothingU + 2; momentIndex < numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV; momentIndex++) {
          fittedCentralMoments(momentIndex) = ((deviations.array().pow(momentIndex - numberOfDegreesOfSmoothingU + 1)).cwiseProduct(relativeFrequencies.array())).sum();
          fittedCentralMoments(momentIndex) /= std::pow(sdV, static_cast<double>(momentIndex - numberOfDegreesOfSmoothingU + 1));
        }
      }

      /* for cuv columns associated with cross products;
          scores are in columns 0 and cu */
      if (numberOfCrossProductMoments > 0 &&
          numberOfDegreesOfSmoothingU >= 2 &&
          numberOfDegreesOfSmoothingV >= 2) {
        for (size_t columnIndex = numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV;
             columnIndex < numberOfColumns;
             columnIndex++) {
          for (size_t rowIndex = 0; rowIndex < rawScoreDesigmMatrix.rows(); rowIndex++) {
            fittedCentralMoments(columnIndex) += std::pow(rawScoreDesigmMatrix(rowIndex, 0) - meanU,
                                                          static_cast<double>(crossProductMomentDesignations(columnIndex - numberOfDegreesOfSmoothingU - numberOfDegreesOfSmoothingV, 0))) *
                                                 std::pow(rawScoreDesigmMatrix(rowIndex, numberOfDegreesOfSmoothingU),
                                                          static_cast<double>(crossProductMomentDesignations(columnIndex - numberOfDegreesOfSmoothingU - numberOfDegreesOfSmoothingV, 1))) *
                                                 relativeFrequencies(rowIndex);
          }

          fittedCentralMoments(columnIndex) /= std::pow(sdU,
                                                        static_cast<double>(crossProductMomentDesignations(columnIndex - numberOfDegreesOfSmoothingU - numberOfDegreesOfSmoothingV, 0))) *
                                               std::pow(sdV,
                                                        static_cast<double>(crossProductMomentDesignations(columnIndex - numberOfDegreesOfSmoothingU - numberOfDegreesOfSmoothingV, 1)));
        }
      }
    }

    /*
      Moments criterion for iteration loop

      Input
        nc = number of cmoments
        cu = number of moments for u
        ctype = comparison type for criterion:
                0 --> absolute; 1 --> relative
        Btype = type of moments for criterion mathching:
                0 --> moments based on B 
                (if scale = 0, design matrix is based on raw scores,
                which means that the moments are based on raw scores;
                if scale = 1, design matrix is based on scaled raw scores,
                which means that the moments are based on scaled raw scores) 
                1 -->  moments based on B_raw, whether scale is 0 or 1
          NOTE: see design_matrix() for comments about scale
        n_mts[] = moments based on actual frequencies
        m_mts[] = moments based on fitted frequencies
        crit = criterion  (see notes below)

      return: 0 --> criterion not met ---- continue iteration
              1 --> criterion met --- end iteration

      NOTES.

      If ctype==0 (absolute comparison), criterion is met if
      |n_mts[i] - m_mts[i]| <= crit for all nc moments.

      If ctype = 1 (relative comparison), criterion is met if
      |(n_mts[i] - m_mts[i])/n_mts[i]| <= crit for all nc moments.  
      An error occurs if n_mts[i] = 0.

      If Btype==1, the first central moment for both u and v
      is 0 whether the frequencies are actual (n) or fitted (m).
      These are the moments associated with column 0 and cu in B[][].
      Further, when Btype==1, these moments are  replaced by u and v 
      means, respectively, in get LLmoments().Therefore, when Btype==1,
      no test is made for |n_mts[0] - m_mts[0]| or for
      |n_mts[cu] - m_mts[cu]|. 

      Function calls other than C or NR utilities:
        runerror()

      R. L. Brennan

      Date of last revision: 6/30/08

    */
    bool momentsCriterion(const size_t& numberOfCMoments,
                          const size_t& numberOfMomentsForU,
                          const CriterionComparisonType& criterionComparisonType,
                          const DesignMatrixType& designMatrixType,
                          const Eigen::VectorXd& observedMoments,
                          const Eigen::VectorXd& fittedMoments,
                          const double& criterion) {
      // int i;

      if (designMatrixType == DesignMatrixType::SOLUITON) {
        for (size_t momentIndex = 0; momentIndex < numberOfCMoments; momentIndex++) {
          if (criterionComparisonType == CriterionComparisonType::ABSOLUTE) {
            if (std::abs(observedMoments(momentIndex) - fittedMoments(momentIndex)) > criterion) {
              return false;
            }
          } else {
            if (observedMoments(momentIndex) == 0.0) {
              throw std::runtime_error("\n\nRelative criterion and n_mts[] = 0");
            }

            if (std::abs((observedMoments(momentIndex) - fittedMoments(momentIndex)) / observedMoments(momentIndex)) > criterion) {
              return false;
            }
          }
        }
      } else {
        for (size_t momentIndex = 0; momentIndex < numberOfCMoments; momentIndex++) {
          if (momentIndex == 0 || momentIndex == numberOfMomentsForU) {
            if (criterionComparisonType == CriterionComparisonType::ABSOLUTE) {
              if (std::abs(observedMoments(momentIndex) - fittedMoments(momentIndex)) > criterion) {
                return false;
              }
            } else {
              if (observedMoments(momentIndex) == 0.0) {
                throw std::runtime_error("\n\nRelative criterion and n_mts[] = 0");
              }

              if (std::abs((observedMoments(momentIndex) - fittedMoments(momentIndex)) / observedMoments(momentIndex)) > criterion) {
                return false;
              }
            }
          }
        }
      }

      return true;
    }

    /*
      For log-linear iterations print (a) actual results based on nct[]
      and (b) heading for fitted results

      Input
        fp          = file pointer for output
        ns          = number of score categories 
                      (number of rows in design matrix)
        nc          = number of columns in design matrix
        nct[]       = actual frequencies
        n_mts[]     = moments based on nct[] and design matrix B
        n_mts_raw[] = moments based on nct[] and design matrix B_raw
        ctype = comparison type for criterion:
                0 --> absolute; 1 --> relative
        Btype = type of design matrix and, hence, type of moments 
                for criterion:
                0 --> use B (could be scaled or unscaled as indicated
                      in design_matrix()) 
                1 --> use B_raw and central moments based on it
        crit = criterion.  See crit_mts() for discussion 

      Function calls other than C or NR utilities: none

      R. L. Brennan

      Date of last revision: 6/30/08


    */
    void PrintIterationHeading(const size_t& numberOfDesignMarixRows,
                               const size_t& numberOfDesignMatrixColumns,
                               const Eigen::VectorXd& observedFrequencies,
                               const Eigen::VectorXd& observedMomentsSolutionDesignMatrix,
                               const Eigen::VectorXd& observedMomentsRawScoreDesignMatrix,
                               const CriterionComparisonType& criterionComparisonType,
                               const DesignMatrixType& designMatrixType,
                               const double& criterion) {
      // int i,
      //     pcells = (ns > 250) ? 250 : ns;

      // /* print headings */

      // fprintf(fp, "\n\n\nRESULTS FOR ITERATIONS");

      // fprintf(fp, "\n\n              ");
      // for (i = 1; i <= (12 * (nc + 2)); i++)
      //   fprintf(fp, " ");

      // for (i = 1; i <= ((12 * nc) - 18) / 2; i++)
      //   fprintf(fp, " ");
      // fprintf(fp, "Moments Based on B");
      // for (i = 1; i <= ((12 * nc) - 18) / 2; i++)
      //   fprintf(fp, " ");

      // for (i = 1; i <= ((12 * nc) - 30) / 2; i++)
      //   fprintf(fp, " ");
      // fprintf(fp, "Central Moments Based on B_raw");
      // for (i = 1; i <= ((12 * nc) - 30) / 2; i++)
      //   fprintf(fp, " ");

      // fprintf(fp, "\n              ");
      // for (i = 1; i <= (12 * (nc + 2)); i++)
      //   fprintf(fp, " ");
      // fprintf(fp, " ");
      // for (i = 1; i <= 12 * nc - 1; i++)
      //   fprintf(fp, "-");
      // fprintf(fp, " ");
      // for (i = 1; i <= 12 * nc - 1; i++)
      //   fprintf(fp, "-");

      // /* print actual results based on nct[] */

      // fprintf(fp, "\n              ");
      // for (i = 1; i <= (12 * (nc + 2)); i++)
      //   fprintf(fp, " ");
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "       n[%2d]", i + 1);
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "   n_raw[%2d]", i + 1);
      // for (i = 0; i < 19; i++)
      //   fprintf(fp, " ");
      // for (i = 0; i < pcells; i++)
      //   fprintf(fp, "    nct[%3d]", i);
      // if (ns > 250)
      //   fprintf(fp, "  ...");
      // fprintf(fp, "\n\n        Actual");
      // for (i = 1; i <= (12 * (nc + 2)); i++)
      //   fprintf(fp, " ");
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "%12.5f", n_mts[i]);
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "%12.5f", n_mts_raw[i]);
      // for (i = 0; i < 19; i++)
      //   fprintf(fp, " ");
      // for (i = 0; i < pcells; i++)
      //   fprintf(fp, "%12.5f", nct[i]);
      // if (ns > 250)
      //   fprintf(fp, "  ...");

      // /* print criterion */

      // fprintf(fp, "\n\n    Criterion:  ");
      // if (Btype == 0 && ctype == 0)
      //   fprintf(fp, "|(n[i] - m[i]| <= %15.10f", crit);
      // else if (Btype == 0 && ctype == 1)
      //   fprintf(fp, "|n[i] - m[i])/n[i]| <= %15.10f", crit);
      // else if (Btype == 1 && ctype == 0)
      //   fprintf(fp, "|(n_raw[i] - m_raw[i]| <= %15.10f", crit);
      // else
      //   fprintf(fp, "|n_raw[i] - m_raw[i])/n_raw[i]|",
      //           " <= %15.10f", crit);
      // fprintf(fp, "  for i = 1 to %d", nc);

      // /* print fitted-results heading for iterations */

      // fprintf(fp, "\n\n     Iteration");
      // fprintf(fp, "          ap");
      // fprintf(fp, "           u");
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "    Beta[%2d]", i + 1);
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "       m[%2d]", i + 1);
      // for (i = 0; i < nc; i++)
      //   fprintf(fp, "   m_raw[%2d]", i + 1);
      // fprintf(fp, "     LR Chi2");
      // fprintf(fp, "  nzero");
      // for (i = 0; i < pcells; i++)
      //   fprintf(fp, "    mct[%3d]", i);
      // if (ns > 250)
      //   fprintf(fp, "  ...");
      // fprintf(fp, "\n");
    }
  };
} // namespace EquatingRecipes

#endif