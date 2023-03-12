/*
  RGandSG_NoSmooth.c 

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


  File contains random groups and single group drivers and print function, along 
  with linear equating functions.  Equipercentile functions are in ERutilities.c                                             
*/

#ifndef RG_AND_SG_NO_SMOOTH_HPP
#define RG_AND_SG_NO_SMOOTH_HPP

#include <string>
#include <vector>
#include <Eigen/Core>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/wrappers/utilities.hpp>

namespace EquatingRecipes {
  struct RandomAndSingleGroupEquating {
    /*
      Wrapper for getting mean, linear, or equipercentile equating for RG design
        with no smoothing

      Wrapper_RN() does the following:
          (i) allocate space for r->eraw[0][] and r->mts[0][] 
        (ii) populate elements of inall 
        (iii) compute r->a[0] = slope (if mean or linear), 
              r->b[0] = intercept (if mean or linear), 
              r->eraw[0][] and r->mts[0][] 
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep

      Input
      
        design = 'R' (random groups)
        method = 'M', 'L', or 'E'(mean, linear, or equipercentile)
        smoothing = 'N' (none)  
        x = struct USTATS
        y = struct USTATS 
        rep = replication number for bootstrap; should be set to 0
              for actual equating;  
        
      Output
        
        struct PDATA inall   populates selected values of inall 
        
        struct ERAW_RESULTS r          

          a[0] = slope (if method = 'M', 'L', or 'E')
          b[0] = intercept (if method = 'M', 'L', or 'E')
          eraw[][]: equated raw scores;          
                    method (rows) by raw score (columns) matrix
                    of equated scores. Here there is only one method.
                    So, memory allocated for eraw[][] is: 
                    eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                                  (loc(x->max,x->min,x>-inc)]
                    because we are getting equated raw scores for x   

          mts[][]:  moments for equated raw scores         
          
      NOTE: If Wrapper_RN() is called in a bootstrap loop, then in the 
            calling function struct ERAW_RESULTS must be different 
            from struct ERAW_RESULTS for the actual equating. 
                                                
      Function calls other than C or NR utilities:
        RGandSG_LinEq()
        EquiEquate()
        MomentsFromFD()  
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void randomGroupEquating(const EquatingRecipes::Structures::Design& design,
                             const EquatingRecipes::Structures::Method& method,
                             const EquatingRecipes::Structures::Smoothing& smoothing,
                             const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsX,
                             const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsY,
                             const size_t& replication,
                             EquatingRecipes::Structures::PData& pData,
                             EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      /* method names --- 10 characters; right justified */
      std::vector<std::string> methodNames {"   y-equiv"};

      /* should be set to 0 for actual equating */
      /* counting of replications done in Wrapper_Bootstrap() */
      pData.bootstrapReplicationNumber = replication;

      /* allocation and assignments for in
      Note that for every assignment of the form inall->(var) = xv->(var)
      values vary depending on whether x (or y) is for actual equating or
      a bootstrap sample; all other values are the same for the 
      actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) {
        pData.summaryRawDataX = univariateStatisticsX;
        pData.summaryRawDataY = univariateStatisticsY;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        pData.methods = methodNames;

        pData.minimumScoreX = univariateStatisticsX.minimumScore;
        pData.maximumScoreX = univariateStatisticsX.maximumScore;
        pData.scoreIncrementX = univariateStatisticsX.scoreIncrement;
        pData.scoreFrequenciesX = univariateStatisticsX.freqDistDouble;
        pData.numberOfExaminees = univariateStatisticsX.numberOfExaminees;
      }

      /* allocation and assignments for r */

      if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.minimumScoreX,
                                                                                   pData.scoreIncrementX);

        equatedRawScoreResults.equatedRawScores.setZero(1, maximumScoreLocation + 1);
        equatedRawScoreResults.equatedRawScoreMoments.setZero(1, 4);
      }

      /* Compute equating results */
      double slope;
      double intercept;
      Eigen::VectorXd equatedRawScores;

      if (method != EquatingRecipes::Structures::Method::EQUIPERCENTILE) {
        equatedRawScoreResults.slope.resize(1);
        equatedRawScoreResults.intercept.resize(1);

        equatedRawScores = linearEquating(univariateStatisticsX.momentValues(0),
                                          univariateStatisticsX.momentValues(1),
                                          univariateStatisticsY.momentValues(0),
                                          univariateStatisticsY.momentValues(1),
                                          method,
                                          pData.minimumScoreX,
                                          pData.maximumScoreX,
                                          pData.scoreIncrementX,
                                          equatedRawScoreResults.slope(0),
                                          equatedRawScoreResults.intercept(0));

      } else {
        equatedRawScores = EquatingRecipes::Utilities::getEquipercentileEquivalents(univariateStatisticsY.numberOfScores,
                                                                                    univariateStatisticsY.minimumScore,
                                                                                    univariateStatisticsY.scoreIncrement,
                                                                                    univariateStatisticsY.cumulativeRelativeFreqDist,
                                                                                    univariateStatisticsX.numberOfScores,
                                                                                    univariateStatisticsX.percentileRankDist);
      }

      equatedRawScoreResults.equatedRawScores.row(0) = equatedRawScores;

      /* get moments */
      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(equatedRawScoreResults.equatedRawScores.row(0),
                                                                                                                   pData.scoreFrequenciesX);

      equatedRawScoreResults.equatedRawScoreMoments.row(0) = moments.momentValues;
    }

    /*
      Wrapper for getting mean, linear, or equipercentile equating for SG design
        with no smoothing

      Wrapper_SN() does the following:
          (i) allocate space for r->eraw[0][] and r->mts[0][] 
        (ii) populate elements of inall 
        (iii) compute r->a[0] = slope (if mean or linear),
              r->b[0] = intercept (if mean or linear), and
              r->eraw[0][] and r->mts[0][] 
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep
      NOTE: required input is one struct BSTATS

      Input
      
        design = 'S' (single group)
        method = 'M', 'L', or 'E'(mean, linear, or equipercentile)
        smoothing = 'N' (none)  
        xy = struct BSTATS 
        rep = replication number for bootstrap; should be set to 0
              for actual equating;  
        
        NOTE: it is assumed that the first variable
              in xy is indeed x and the second variable is y

      Output
        
        struct PDATA inall   populates selected values of inall 
        
        struct ERAW_RESULTS r          

          a[0] = slope (if mean or linear)
          b[0] = intercept (if mean or linear)
          eraw[][]: equated raw scores;          
                    method (rows) by raw score (columns) matrix
                    of equated scores. Here there is only one method.
                    So, memory allocated for eraw[][] is: 
                    eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                                  (loc(x->max,x->min,x>-inc)]
                    because we are getting equated raw scores for x   

          mts[][]:  moments for equated raw scores          
          
      NOTE: If Wrapper_SN() is called in a bootstrap loop,
            then in the calling function struct ERAW_RESULTS must
            be different from struct ERAW_RESULTS for the actual
            equating. 
                                                
      Function calls other than C or NR utilities:
        RGandSG_LinEq()
        EquiEquate()
        MomentsFromFD()  
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08   
    */
    void singleGroupEquating(const EquatingRecipes::Structures::Design& design,
                             const EquatingRecipes::Structures::Method& method,
                             const EquatingRecipes::Structures::Smoothing& smoothing,
                             const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXY,
                             const size_t& replication,
                             EquatingRecipes::Structures::PData& pData,
                             EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      // void Wrapper_SN(char design, char method, char smoothing,  struct BSTATS *xy,
      //                int rep, struct PDATA *inall, struct ERAW_RESULTS *r)

      pData.methods.push_back("   y-equiv");

      /* should be set to 0 for actual equating */
      /* counting of replications done in Wrapper_Bootstrap() */
      pData.bootstrapReplicationNumber = replication;

      /* allocation and assignments for in
        Note that for every assignment of the form inall->(var) = xv->(var)
        values vary depending on whether xy is for actual equating or
        a bootstrap sample; all other values are the same for the 
        actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) { /* no assignment or stor alloc for bootstrap reps */
        pData.summaryRawDataXY = bivariateStatisticsXY;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;

        pData.minimumScoreX = bivariateStatisticsXY.univariateStatisticsRow.minimumScore;
        pData.maximumScoreX = bivariateStatisticsXY.univariateStatisticsRow.maximumScore;
        pData.scoreIncrementX = bivariateStatisticsXY.univariateStatisticsRow.scoreIncrement;
        pData.scoreFrequenciesX = bivariateStatisticsXY.univariateStatisticsRow.freqDistDouble;
        pData.numberOfExaminees = bivariateStatisticsXY.univariateStatisticsRow.numberOfExaminees;
      }

      /* allocation and assignments for r */
      /* no storage allocation for bootstrap reps >1 */
      if (pData.bootstrapReplicationNumber <= 1) {
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.minimumScoreX,
                                                                                   pData.scoreIncrementX);

        equatedRawScoreResults.equatedRawScores.setZero(pData.methods.size(), maximumScoreLocation + 1);
        equatedRawScoreResults.equatedRawScoreMoments.setZero(pData.methods.size(), 4);
      }

      /* Compute equating results */

      //   if(method != 'E')
      if (method == EquatingRecipes::Structures::Method::EQUIPERCENTILE) {
        equatedRawScoreResults.equatedRawScores = EquatingRecipes::Utilities::getEquipercentileEquivalents(bivariateStatisticsXY.univariateStatisticsColumn.numberOfScores,
                                                                                                           bivariateStatisticsXY.univariateStatisticsColumn.minimumScore,
                                                                                                           bivariateStatisticsXY.univariateStatisticsColumn.scoreIncrement,
                                                                                                           bivariateStatisticsXY.univariateStatisticsColumn.cumulativeRelativeFreqDist,
                                                                                                           bivariateStatisticsXY.univariateStatisticsColumn.numberOfScores,
                                                                                                           bivariateStatisticsXY.univariateStatisticsColumn.percentileRankDist);
      } else {
        equatedRawScoreResults.equatedRawScores = linearEquating(bivariateStatisticsXY.univariateStatisticsColumn.momentValues(0),
                                                                 bivariateStatisticsXY.univariateStatisticsColumn.momentValues(1),
                                                                 bivariateStatisticsXY.univariateStatisticsColumn.momentValues(0),
                                                                 bivariateStatisticsXY.univariateStatisticsColumn.momentValues(1),
                                                                 method,
                                                                 pData.minimumScoreX,
                                                                 pData.maximumScoreX,
                                                                 pData.scoreIncrementX,
                                                                 equatedRawScoreResults.slope(0),
                                                                 equatedRawScoreResults.intercept(0));
      }

      /* get moments */

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(equatedRawScoreResults.equatedRawScores,
                                                                                                                   pData.scoreFrequenciesX);

      equatedRawScoreResults.equatedRawScoreMoments = moments.momentValues;
    }

    /*  
      Linear equating for RG and SG designs

      Input
        mnx = mean for x
        sdx = sd for x
        mny = mean for y
        sdy = sd for y
        method = 'M' or 'L'
        min = minimum score for x
        max = maximum score for x
        inc = increment for x

      Output
        a = slope
        b = intercept
        eraw[] = equated raw scores

      Function calls other than C or NR utilities:
        score() 
                                                  
      R. L. Brennan

      Date of last revision: 6/30/08   
    */
    Eigen::VectorXd linearEquating(const double& meanX,
                                   const double& sdX,
                                   const double& meanY,
                                   const double& sdY,
                                   const EquatingRecipes::Structures::Method& method,
                                   const double& minimumScoreX,
                                   const double& maximumScoreX,
                                   const double& scoreIncrementX,
                                   double& slope,
                                   double& intercept) {
      if (method == EquatingRecipes::Structures::Method::MEAN) {
        slope = 1.0;
      } else {
        slope = sdY / sdX;
      }

      intercept = meanY - slope * meanX;

      /* get equated raw scores */
      size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(maximumScoreX,
                                                                                 minimumScoreX,
                                                                                 scoreIncrementX);

      Eigen::VectorXd equatedRawScores(maximumScoreLocation + 1);

      for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
        equatedRawScores(scoreLocation) = intercept + slope * EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                                   minimumScoreX,
                                                                                                   scoreIncrementX);
      }

      return equatedRawScores;
    }
  };
} // namespace EquatingRecipes

#endif