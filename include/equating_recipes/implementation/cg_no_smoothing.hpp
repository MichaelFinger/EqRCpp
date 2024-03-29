/*
  CG_NoSmooth.c 

  Wrapper, print, and linear functions for common-item non-equivalent groups 
  design with no smoothing.

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

#ifndef CG_NO_SMOOTHING_HPP
#define CG_NO_SMOOTHING_HPP

#include <cmath>
#include <string>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>
#include <equating_recipes/implementation/cg_equipercentile_equating.hpp>

namespace EquatingRecipes {
  namespace Implementation {
    class CGEquatingNoSmoothing {
    public:
      /*
      Wrapper for getting mean, linear, or equipercentile equating for CINEG design
      with no smoothing. Equipercentile equating includes frequency estimation with 
      Braun-Holland (linear) results, modified frequency estimation with 
      Braun-Holland (linear) results, and chained equipercentile equating

      Purposes of Wrapper_CN() are to:
          (i) allocate space for r->eraw[][] and r->mts[0][] 
        (ii) populate elements of inall, including the 
              determination of proportional weights, if requested
        (iii) get equating results and store them in r
        (iv) get moments
        
      Assumes that in xv, score 1 is for x and score 2 is for v
      Assumes that in yv, score 1 is for y and score 2 is for v
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep
      
      Input:
      
        design = 'C' (CINEG)
        method:  'M' = mean, 
                'L' = linear (Tucker, Levine-obs, Levine-true, chained)
                'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
                'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
                'G' = FE + BH-FE + MFE + BH-MFE
                'C' = Chained
                'H' = FE + BH-FE + Chained
                'A' = FE + BH-FE + MFE + BH-MFE + Chained
                  
        smoothing = 'N' (none)  
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
        rep = replication number for bootstrap; should be set to 0
              for actual equating; 
      
        NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
        
      Output:
        
        struct PDATA inall:   populates selected values of inall 
        
        struct ERAW_RESULTS r          

          msx[] = means for x for synthetic pop 
          msy[] = means for y for synthetic pop 
          dssx[] = sd's for x for synthetic pop 
          ssy[] = sd's for y for synthetic pop 
          gamma1[] = gamma's for pop 1 
          gamma2[] = gamma's for pop 2 
          a[] = slopes
          b[] = intercepts
          eraw[][]:  equated raw scores
          mts[][]:  moments for equated raw scores  
          
          NOTE: vectors in ERAW_RESULTS for nm = 4 linear methods 
                give results for following methods: 0 - Tucker,
                1 - Levine Observed, 2 - Levine True, 3 - Chained  
              
          NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
                of equated scores; memory allocated here;
                eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
                because we are getting equated raw scores for x.
                eraw[][] stored in this "row"  manner so that 
                Equated_ss() can be used directly on 
                eraw[k] where k is the method number  
              
      NOTE: Whenever method differs, there must be different structures
            passed as struct PDATA and struct ERAW_RESULTS 
        
      NOTE: If Wrapper_CN() is called in a bootstrap loop, then in
            the calling function struct ERAW_RESULTS must be different
            from struct ERAW_RESULTS for the actual equating. 
                                                
      Function calls other than C or NR utilities:
        CI_LinEq() 
        FEorMFE_EE()
        Chained_EE()
        runerror()
                                                  
      R. L. Brennan

      Date of last revision: 6/30/08   
    */
      void runWithNoSmoothing(const EquatingRecipes::Structures::Design& design,
                              const EquatingRecipes::Structures::Method& method,
                              const EquatingRecipes::Structures::Smoothing& smoothing,
                              const double& population1Weight,
                              const bool& isInternalAnchor,
                              const double& reliabilityCommonItemsPopulation1,
                              const double& reliabilityCommonItemsPopulation2,
                              EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXV,
                              EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsYV,
                              const size_t& bootstrapReplicationNumber,
                              EquatingRecipes::Structures::PData& pData,
                              EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
        std::vector<std::string> methodNames = {"    Tucker",
                                                "   Lev Obs",
                                                "  Lev True",
                                                "  ChainedL",
                                                "        FE",
                                                "     BH-FE",
                                                "       MFE",
                                                "    BH-MFE",
                                                "  ChainedE"};

        /* should be set to 0 for actual equating. */
        /* Counting of replications done in Wrapper_Bootstrap(),
      which is why this statement cannot be in the if statement below */
        pData.bootstrapReplicationNumber = bootstrapReplicationNumber;

        std::string methodCode = EquatingRecipes::Implementation::Utilities::getMethodCode(method);

        /* allocation and assignments for in
          Note that for every assignment of the form inall->(var) = xv->(var)
          values vary depending on whether xv is for actual equating or
          a bootstrap sample; all other values are the same for the
          actual equating and a bootstrap sample */
        if (pData.bootstrapReplicationNumber == 0) { /* no assignment or stor alloc for bootstrap reps */
          pData.summaryRawDataXV = bivariateStatisticsXV;
          pData.summaryRawDataYV = bivariateStatisticsYV;
          pData.design = design;
          pData.method = method;
          pData.smoothing = smoothing;

          if (population1Weight < 0.0 || population1Weight > 1.0) {
            /* proportional wts if w1 outside [0,1] */
            pData.weightSyntheticPopulation1 = static_cast<double>(bivariateStatisticsXV.numberOfExaminees) /
                                               static_cast<double>(bivariateStatisticsXV.numberOfExaminees + bivariateStatisticsYV.numberOfExaminees);
          } else {
            pData.weightSyntheticPopulation1 = population1Weight;
          }

          pData.isInternalAnchor = isInternalAnchor;
          pData.reliabilityCommonItemsPopulation1 = reliabilityCommonItemsPopulation1;
          pData.reliabilityCommonItemsPopulation2 = reliabilityCommonItemsPopulation2;

          pData.methodCode = methodCode;

          if ((methodCode == "F" ||
               methodCode == "G" ||
               methodCode == "A") &&
              (reliabilityCommonItemsPopulation1 == 0.0 ||
               reliabilityCommonItemsPopulation2 == 0.0)) {
            throw std::runtime_error("\nMFE cannot be conducted since reliabilityCommonItemsPopulation1 == 0 or reliabilityCommonItemsPopulation2 == 0");
          }

          size_t methodStartIndex;
          size_t methodEndIndex;

          if (methodCode == "M" ||
              methodCode == "L") {
            methodStartIndex = 0;
            methodEndIndex = 3;

          } else if (methodCode == "E") {
            methodStartIndex = 4;
            methodEndIndex = 5;

          } else if (methodCode == "F") {
            methodStartIndex = 6;
            methodEndIndex = 7;

          } else if (methodCode == "G") {
            methodStartIndex = 4;
            methodEndIndex = 7;

          } else if (methodCode == "C") {
            methodStartIndex = 8;
            methodEndIndex = 8;

          } else if (methodCode == "H") {
            methodStartIndex = 4;
            methodEndIndex = 5;
          } else if (methodCode == "A") {
            methodStartIndex = 4;
            methodEndIndex = 8;
          } else {
            throw std::runtime_error("Invalid method code.");
          }

          for (size_t methodIndex = methodStartIndex; methodIndex <= methodEndIndex; methodIndex++) {
            pData.methods.push_back(methodNames[methodIndex]);
          }

          if (methodCode == "H") {
            pData.methods.push_back(methodNames[8]);
          }

          pData.minimumScoreX = bivariateStatisticsXV.univariateStatisticsRow.minimumScore;
          pData.maximumScoreX = bivariateStatisticsXV.univariateStatisticsRow.maximumScore;
          pData.scoreIncrementX = bivariateStatisticsXV.univariateStatisticsRow.scoreIncrement;
          pData.scoreFrequenciesX = bivariateStatisticsXV.univariateStatisticsRow.freqDistDouble;
          pData.numberOfExaminees = bivariateStatisticsXV.numberOfExaminees;
        }

        /* allocation and assignments for r */

        if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
          size_t maxScoreLocationX = EquatingRecipes::Implementation::Utilities::getScoreLocation(bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                                                                                  bivariateStatisticsXV.univariateStatisticsRow.minimumScore,
                                                                                                  bivariateStatisticsXV.univariateStatisticsRow.scoreIncrement);

          size_t maxScoreLocationY = EquatingRecipes::Implementation::Utilities::getScoreLocation(bivariateStatisticsYV.univariateStatisticsRow.maximumScore,
                                                                                                  bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                                                  bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement);

          equatedRawScoreResults.equatedRawScores.resize(pData.methods.size(), maxScoreLocationX + 1);
          equatedRawScoreResults.equatedRawScoreMoments.resize(pData.methods.size(), 4);
          equatedRawScoreResults.relativeFreqDistsX.resize(2, maxScoreLocationX + 1);
          equatedRawScoreResults.relativeFreqDistsY.resize(2, maxScoreLocationY + 1);
          equatedRawScoreResults.slope.resize(pData.methods.size());
          equatedRawScoreResults.intercept.resize(pData.methods.size());
        }

        /* Mean and Linear equating results:

         xv and yv variables relative to variables in CI_LinEq()

         Function  Function Call
          mnx1      xv->mts1[0]
          sdx1      xv->mts1[1]
          mnv1      xv->mts2[0]
          sdv1      xv->mts2[1]
          covxv1    xv->cov

          mny2      yv->mts1[0]
          sdy2      yv->mts1[1]
          mnv2      yv->mts2[0]
          sdv2      yv->mts2[1]
          covyv2    yv->cov

          In function, 1 and 2 are populations; in xv and yv, 1 and 2 are variables
      */

        if (methodCode == "M" || methodCode == "L") {
          commonItemLinearEquating(bivariateStatisticsXV.univariateStatisticsRow.momentValues(0),
                                   bivariateStatisticsXV.univariateStatisticsRow.momentValues(1),
                                   bivariateStatisticsXV.univariateStatisticsColumn.momentValues(0),
                                   bivariateStatisticsXV.univariateStatisticsColumn.momentValues(1),
                                   bivariateStatisticsXV.covariance,
                                   bivariateStatisticsYV.univariateStatisticsRow.momentValues(0),
                                   bivariateStatisticsYV.univariateStatisticsRow.momentValues(1),
                                   bivariateStatisticsYV.univariateStatisticsColumn.momentValues(0),
                                   bivariateStatisticsYV.univariateStatisticsColumn.momentValues(1),
                                   bivariateStatisticsYV.covariance,
                                   pData.weightSyntheticPopulation1,
                                   isInternalAnchor,
                                   pData.method,
                                   bivariateStatisticsXV.univariateStatisticsRow.minimumScore,
                                   bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                   bivariateStatisticsXV.univariateStatisticsRow.scoreIncrement,
                                   pData.methods.size(),
                                   equatedRawScoreResults.xSyntheticPopulationMean,
                                   equatedRawScoreResults.ySyntheticPopulationMean,
                                   equatedRawScoreResults.xSyntheticPopulationSD,
                                   equatedRawScoreResults.ySyntheticPopulationSD,
                                   equatedRawScoreResults.gammaPopulation1,
                                   equatedRawScoreResults.gammaPopulation2,
                                   equatedRawScoreResults.slope,
                                   equatedRawScoreResults.intercept,
                                   equatedRawScoreResults.equatedRawScores);
        }

        /* Equipercentile results, including Braun-Holland (BH) linear.
         Note: For FE syn densities are in fxs[0] and gys[0]
               For MFE syn densities are in fxs[1] and gys[1]
               For BH under FE, slope in a[0] and intercept in b[0]
               For BH under MFE, slope in a[1] and intercept in b[1]
        */

        EquatingRecipes::Implementation::CGEquipercentileEquating cgEquipercentileEquating;

        /* FE + BH-FE in positions 0 and 1 */
        if (methodCode == "E" || methodCode == "G" || methodCode == "A" || methodCode == "H") {
          EquatingRecipes::Structures::CGEquipercentileEquatingResults cgResults =
              cgEquipercentileEquating.runCGEquipercentileEquating(pData.weightSyntheticPopulation1,
                                                                   pData.isInternalAnchor,
                                                                   bivariateStatisticsXV.univariateStatisticsColumn.numberOfScores,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.numberOfScores,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.minimumScore,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.numberOfScores,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.maximumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement,
                                                                   true,
                                                                   bivariateStatisticsXV.bivariateProportions,
                                                                   bivariateStatisticsYV.bivariateProportions,
                                                                   0.0,
                                                                   0.0);

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
          /* MFE + BH-MFE in positions 0 and 1 */
          EquatingRecipes::Structures::CGEquipercentileEquatingResults cgResults =
              cgEquipercentileEquating.runCGEquipercentileEquating(pData.weightSyntheticPopulation1,
                                                                   pData.isInternalAnchor,
                                                                   bivariateStatisticsXV.univariateStatisticsColumn.numberOfScores,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.numberOfScores,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.minimumScore,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.numberOfScores,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.maximumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement,
                                                                   true,
                                                                   bivariateStatisticsXV.bivariateProportions,
                                                                   bivariateStatisticsYV.bivariateProportions,
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
          /* MFE + BH-MFE in positions 2 and 3 */
          EquatingRecipes::Structures::CGEquipercentileEquatingResults cgResults =
              cgEquipercentileEquating.runCGEquipercentileEquating(pData.weightSyntheticPopulation1,
                                                                   pData.isInternalAnchor,
                                                                   bivariateStatisticsXV.univariateStatisticsColumn.numberOfScores,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.numberOfScores,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.minimumScore,
                                                                   bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.numberOfScores,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.maximumScore,
                                                                   bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement,
                                                                   true,
                                                                   bivariateStatisticsXV.bivariateProportions,
                                                                   bivariateStatisticsYV.bivariateProportions,
                                                                   pData.reliabilityCommonItemsPopulation1,
                                                                   pData.reliabilityCommonItemsPopulation2);

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

        /* Chained equipercentile method */
        if (methodCode == "C") {
          /* Chained in position 0 */
          equatedRawScoreResults.equatedRawScores.row(0) =
              cgEquipercentileEquating.runChainedEquipercentileEquating(bivariateStatisticsXV.univariateStatisticsRow.numberOfScores,
                                                                        bivariateStatisticsXV.univariateStatisticsRow.percentileRankDist,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.minimumScore,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.maximumScore,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.scoreIncrement,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.numberOfScores,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.cumulativeRelativeFreqDist,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.numberOfScores,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.cumulativeRelativeFreqDist,
                                                                        bivariateStatisticsYV.univariateStatisticsColumn.cumulativeRelativeFreqDist);
        }

        /* All methods: FE, BF under FE, MFE, BH under MFE, Chained */
        if (methodCode == "A") {
          /* Chained in position 4 */
          equatedRawScoreResults.equatedRawScores.row(4) =
              cgEquipercentileEquating.runChainedEquipercentileEquating(bivariateStatisticsXV.univariateStatisticsRow.numberOfScores,
                                                                        bivariateStatisticsXV.univariateStatisticsRow.percentileRankDist,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.minimumScore,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.maximumScore,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.scoreIncrement,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.numberOfScores,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.cumulativeRelativeFreqDist,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.numberOfScores,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.cumulativeRelativeFreqDist,
                                                                        bivariateStatisticsYV.univariateStatisticsColumn.cumulativeFreqDist);
        }

        //  /* FE, BF under FE, Chained */
        if (methodCode == "H") {
          /* Chained in position 2 */
          equatedRawScoreResults.equatedRawScores.row(2) =
              cgEquipercentileEquating.runChainedEquipercentileEquating(bivariateStatisticsXV.univariateStatisticsRow.numberOfScores,
                                                                        bivariateStatisticsXV.univariateStatisticsRow.percentileRankDist,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.minimumScore,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.maximumScore,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.scoreIncrement,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.numberOfScores,
                                                                        bivariateStatisticsXV.univariateStatisticsColumn.cumulativeRelativeFreqDist,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.numberOfScores,
                                                                        bivariateStatisticsYV.univariateStatisticsRow.cumulativeRelativeFreqDist,
                                                                        bivariateStatisticsYV.univariateStatisticsColumn.cumulativeFreqDist);
        }

        /* get moments */
        for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
          EquatingRecipes::Structures::Moments moments = EquatingRecipes::Implementation::Utilities::momentsFromScoreFrequencies(equatedRawScoreResults.equatedRawScores.row(methodIndex),
                                                                                                                                 bivariateStatisticsXV.univariateStatisticsRow.freqDistDouble);

          equatedRawScoreResults.equatedRawScoreMoments.row(methodIndex) = moments.momentValues;
        }
      }

    private:
      /*
      computes results for common-item linear equating
      of x to scale of y for methods:
        0. Tucker
        1. Levine Observed
        2. Levine True
        3. Chained
                  
      Kolen and Brennan (2004) notation used here
      
      Input:
          
          mnx1 = mean for x for pop 1    
          sdx1 = sd for x for pop 1    
          mnv1 = mean for v for pop 1    
          sdv1 = sd for v for pop 1     
          covxv1 = cov for x and v in pop 1  

          mny2 = mean for y in pop 2     
          sdy2 = sd for y in pop 2     
          mnv2 = mean for v in pop 2     
          sdv2 = sd for v in pop 2    
          covyv2 = cov for y and v in pop 2
          
          w1 = weight for pop 1
          anchor = 0 --> external; otherwise internal
          method = 'M' --> mean equating (i.e., slope = 1)
          min = min score for x 
          max = max score for x
          inc = increment for x
          fdx[] = fd for x 
      
      Output:
      
          msx[] = vector of means for x for synthetic pop (for the methods)
          msy[] = vector of means for y for synthetic pop (for the methods)
          ssx[] = vector of sd's for x for synthetic pop (for the methods) 
          ssy[] = vector of sd's for y for synthetic pop (for the methods) 
          gamma1[] = vector of gamma's for pop 1 (for the methods)
          gamma2[] = vector of gamma's for pop 2 (for the methods)      
          a[] = vector of slopes (for the methods) 
          b[] = vector of intercepts (for the methods)
          eraw[][] =  method (rows) by raw score (columns) matrix
                      of equated scores; see Wrapper_CLN()for more details
                    
          NOTE:  it is assumed that storage already allocated for
                all vectors and eraw[][]

        Function calls other than C or NR utilities:
          CI_LinObsEq()
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08                                   
    */
      void commonItemLinearEquating(const double& meanXPop1,
                                    const double& sdXPop1,
                                    const double& meanVPop1,
                                    const double& sdVPop1,
                                    const double& covXVPop1,
                                    const double& meanYPop2,
                                    const double& sdYPop2,
                                    const double& meanVPop2,
                                    const double& sdVPop2,
                                    const double& covYVPop2,
                                    const double& population1Weight,
                                    const bool& isInternalAnchor,
                                    const EquatingRecipes::Structures::Method& method,
                                    const double& minimumScore,
                                    const double& maximumScore,
                                    const double& scoreIncrement,
                                    const size_t& numberOfMethods,
                                    Eigen::VectorXd& meanVectorXSynPop,
                                    Eigen::VectorXd& meanVectorYSynPop,
                                    Eigen::VectorXd& sdVectorXSynPop,
                                    Eigen::VectorXd& sdVectorYSynPop,
                                    Eigen::VectorXd& gammaVectorPop1,
                                    Eigen::VectorXd& gammaVectorPop2,
                                    Eigen::VectorXd& slopes,
                                    Eigen::VectorXd& intercepts,
                                    Eigen::MatrixXd& methodByEquatedRawScores) {
        if (method == EquatingRecipes::Structures::Method::MEAN) {
          /* slopes are unity */
          slopes.setOnes(4);

          /* no synthetic sd's for mean equating */
          sdVectorXSynPop.setZero(4);
          sdVectorYSynPop.setZero(4);
        }

        /* For Levine true, synthetic population is irrelevant */
        meanVectorXSynPop(2) = 0.0;
        meanVectorYSynPop(2) = 0.0;
        sdVectorXSynPop(2) = 0.0;
        sdVectorYSynPop(2) = 0.0;

        /* Tucker */
        gammaVectorPop1(0) = covXVPop1 / std::pow(sdVPop1, 2);
        gammaVectorPop2(0) = covYVPop2 / std::pow(sdVPop2, 2);

        commonItemLinearObservedEquating(meanXPop1,
                                         sdXPop1,
                                         meanVPop1,
                                         sdVPop1,
                                         meanYPop2,
                                         sdYPop2,
                                         meanVPop2,
                                         sdVPop2,
                                         population1Weight,
                                         method,
                                         gammaVectorPop1(0),
                                         gammaVectorPop2(0),
                                         meanVectorXSynPop(0),
                                         meanVectorYSynPop(0),
                                         sdVectorXSynPop(0),
                                         sdVectorYSynPop(0),
                                         slopes(0),
                                         intercepts(0));

        /* Levine */

        if (isInternalAnchor) {
          gammaVectorPop1(1) = std::pow(sdXPop1, 2) / covXVPop1;
          gammaVectorPop2(1) = std::pow(sdYPop2, 2) / covYVPop2;
        } else {
          gammaVectorPop1(1) = (std::pow(sdXPop1, 2) + covXVPop1) / (std::pow(sdVPop1, 2) + covXVPop1);
          gammaVectorPop2(1) = (std::pow(sdYPop2, 2) + covYVPop2) / (std::pow(sdVPop2, 2) + covYVPop2);
        }

        /* Levine lin obs */
        commonItemLinearObservedEquating(meanXPop1,
                                         sdXPop1,
                                         meanVPop1,
                                         sdVPop1,
                                         meanYPop2,
                                         sdYPop2,
                                         meanVPop2,
                                         sdVPop2,
                                         population1Weight,
                                         method,
                                         gammaVectorPop1(1),
                                         gammaVectorPop2(1),
                                         meanVectorXSynPop(1),
                                         meanVectorYSynPop(1),
                                         sdVectorXSynPop(1),
                                         sdVectorYSynPop(1),
                                         slopes(1),
                                         intercepts(1));

        gammaVectorPop1(2) = gammaVectorPop1(1);
        gammaVectorPop2(2) = gammaVectorPop2(1);

        if (method != EquatingRecipes::Structures::Method::MEAN) {
          slopes(2) = gammaVectorPop2(2) / gammaVectorPop1(2); /* Levine lin true slope */
        }

        /* Levine lin true inter */
        intercepts(2) = (meanYPop2 - slopes(2) * meanXPop1) + gammaVectorPop2(2) * (meanVPop1 - meanVPop2);

        /* Chained (Note:  Chained true = Levine true) */
        gammaVectorPop1(3) = sdXPop1 / sdVPop1;
        gammaVectorPop2(3) = sdYPop2 / sdVPop2;

        commonItemLinearObservedEquating(meanXPop1,
                                         sdXPop1,
                                         meanVPop1,
                                         sdVPop1,
                                         meanYPop2,
                                         sdYPop2,
                                         meanVPop2,
                                         sdVPop2,
                                         population1Weight,
                                         method,
                                         gammaVectorPop1(3),
                                         gammaVectorPop2(3),
                                         meanVectorXSynPop(3),
                                         meanVectorYSynPop(3),
                                         sdVectorXSynPop(3),
                                         sdVectorYSynPop(3),
                                         slopes(3),
                                         intercepts(3));

        /* NOTE: The synthetic group means and sd's don't really apply to
      chained linear equating in the sense that the slope and intercept are
      invariant with respect to choice of weights.  However, 
      CI_LinObsEq() is a convenient way to get the slope and intercept */

        /* get equated raw scores */
        size_t maximumScoreLocation = EquatingRecipes::Implementation::Utilities::getScoreLocation(maximumScore,
                                                                                                   minimumScore,
                                                                                                   scoreIncrement);

        for (size_t methodIndex = 0; methodIndex < numberOfMethods; methodIndex++) {
          for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
            double score = EquatingRecipes::Implementation::Utilities::getScore(scoreLocation,
                                                                                minimumScore,
                                                                                scoreIncrement);
            methodByEquatedRawScores(methodIndex, scoreLocation) = intercepts(methodIndex) + slopes(methodIndex) * score;
          }
        }
      }

      /* 
      For a linear observed score equating method and a CINEG design,
      get synthetic pop means and sd's, and get slope (*a) and intercept (*b)  

      Kolen and Brennan (2004) notation used here
    
      Input:
        
        mnx1 = mean for x for pop 1    
        sdx1 = sd for x for pop 1    
        mnv1 = mean for v for pop 1    
        sdv1 = sd for v for pop 1     

        mny2 = mean for y in pop 2     
        sdy2 = sd for y in pop 2     
        mnv2 = mean for v in pop 2     
        sdv2 = sd for v in pop 2    
        
        w1 = weight for pop 1
        method = 'M' --> mean equating (i.e., slope = 1)
        gamma1 = gamma for pop 1 
        gamma2 = gamma for pop 2
    
    Output:

        msx = mean for x for synthetic pop 
        msy = mean for y for synthetic pop 
        ssx = sd for x for synthetic pop  
        ssy = sd for y for synthetic pop       
        a = slope  
        b = intercept 
      
      Function calls other than C or NR utilities: None
                                                  
      R. L. Brennan

      Date of last revision: 6/30/08         
    */
      void commonItemLinearObservedEquating(const double& meanXPop1,
                                            const double& sdXPop1,
                                            const double& meanVPop1,
                                            const double& sdVPop1,
                                            const double& meanYPop2,
                                            const double& sdYPop2,
                                            const double& meanVPop2,
                                            const double& sdVPop2,
                                            const double& population1Weight,
                                            const EquatingRecipes::Structures::Method& method,
                                            const double& gammaPop1,
                                            const double& gammaPop2,
                                            double& meanXSynPop,
                                            double& meanYSynPop,
                                            double& sdXSynPop,
                                            double& sdYSynPop,
                                            double& slope,
                                            double& intercept) {
        double population2Weight = 1.0 - population1Weight;
        double meanUSynPopX = meanXPop1 - population2Weight * gammaPop1 * (meanVPop1 - meanVPop2);
        double meanUSynPopY = meanYPop2 - population1Weight * gammaPop2 * (meanVPop1 - meanVPop2);

        if (method != EquatingRecipes::Structures::Method::MEAN) {
          double varXSynPop = std::pow(sdXPop1, 2) -
                              population2Weight * std::pow(gammaPop1, 2) * (std::pow(sdVPop1, 2) - std::pow(sdVPop2, 2)) +
                              population1Weight * population2Weight * std::pow(gammaPop1, 2) * std::pow(meanVPop1 - meanVPop2, 2);

          double varYSynPop = std::pow(sdYPop2, 2) -
                              population1Weight * std::pow(gammaPop2, 2) * (std::pow(sdVPop1, 2) - std::pow(sdVPop2, 2)) +
                              population1Weight * population2Weight * std::pow(gammaPop2, 2) * std::pow(meanVPop1 - meanVPop2, 2);

          slope = std::sqrt(varYSynPop / varXSynPop);
          sdXSynPop = std::sqrt(varXSynPop);
          sdYSynPop = std::sqrt(varYSynPop);
        }

        intercept = meanUSynPopY - slope * meanUSynPopX;
        meanXSynPop = meanUSynPopX;
        meanYSynPop = meanUSynPopY;
      }
    };
  } // namespace Implementation
} // namespace EquatingRecipes

#endif