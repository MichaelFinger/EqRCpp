/* Bootstrap.c

   For different equating design/method/smoothing the only statements 
   that need to be changed are those between the 'start' and 'end' comments
   in Wrapper_Bootstrap
   
   NOTES:
   
     Boot_BSTATS and Boot_USTATS() use the following functions from NR:
       ran2() and sort()

     Parametric_boot_univ_BB(), Parametric_boot_univ_ULL(), 
	 and Parametric_boot_biv(),uses ran2()
       
     See comments for Equated_ss() for conventions used to find 
       (a) raw scores associated with locations in vectors
       (b) locations in vectors associated with raw scores
 
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

#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <algorithm>
#include <optional>
#include <random>
#include <vector>
#include <Eigen/Core>

#include <equating_recipes/beta_binomial.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/cg_no_smoothing.hpp>
#include <equating_recipes/log_linear_equating.hpp>
#include <equating_recipes/rg_and_sg_equating.hpp>

#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/bootstrap_equated_raw_score_results.hpp>
#include <equating_recipes/structures/bootstrap_equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class Bootstrap {
  public:
    /*
      Wrapper_Bootstrap

      Input
      
        struct PDATA inall = input for actual equating
        nrep = number of bootstrap replications
        idum = seed                         
          
      Output
      
        struct BOOT_ERAW_RESULTS t = boot results for raw scores
        struct BOOT_ESS_RESULTS u = boot results for scale scores 
                              (NULL indicates no bootstrap for scale scores)
                              
      NOTES:
        (a) struct PDATA inall is the input associated with the actual equating
            (starts out with inall->rep = 0); 
        (b) for any bootstrap replication, everything in struct PDATA inall
            remains the same except the USTATS and BSTATS structures 
            and the replication number
              (1) SG design requires one bootstrap struct BSTATS
              (2) RG design requires two bootstrap struct USTATS
              (3) CI design requires two bootstrap struct BSTATS 
        (c) For different equating design/method/smoothing the only statements 
            that need to be changed are those between 
            the '*****start...*****' and '*****end...******' comments   
          
      Function calls other than C or NR utilities:
        Boot_USTATS()
        Boot_BSTATS()  
        Wrapper_RN()
        Wrapper_SN()
        Wrapper_CN()
        Wrapper_RB()
      Wrapper_RL()    
      Wrapper_SL()
        Parametric_boot_univ_BB()
        Parametric_boot_univ_ULL()
      Parametric_boot_biv()
      Boot_initialize_eraw()
        Boot_accumulate_eraw()
        Boot_accumulate_ess()
        Boot_se_eraw()
        Boot_se_ess()
      Wrapper_ESS()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08           
    */
    void runBootstrap(EquatingRecipes::Structures::PData& pData,
                      const size_t& numberOfReplications,
                      EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& bootstrapEquatedRawScoreResults,
                      std::optional<EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults>& bootstrapEquatedScaledScoresResults) {
      this->gen = EquatingRecipes::Utilities::getSeedEngine();
      this->distrib = EquatingRecipes::Utilities::getUniformRealDistribution();

      EquatingRecipes::Structures::BivariateStatistics bootstrapXV;                                     /* for bootstrap reps and CI design */
      EquatingRecipes::Structures::BivariateStatistics bootstrapYV;                                     /* for bootstrap reps and CI design */
      EquatingRecipes::Structures::UnivariateStatistics bootstrapX;                                     /* for bootstrap reps and RG design */
      EquatingRecipes::Structures::UnivariateStatistics bootstrapY;                                     /* for bootstrap reps and RG design */
      EquatingRecipes::Structures::BivariateStatistics bootstrapXY;                                     /* for bootstrap reps and SG design */
      EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;                       /* for bootstrap reps; not dependent on design */
      EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;                /* for bootstrap reps; not dependent on design */
      EquatingRecipes::Structures::BetaBinomialSmoothing bootstrapBetaBinomialSmoothingX;               /* for boot reps with REB */
      EquatingRecipes::Structures::BetaBinomialSmoothing bootstrapBetaBinomialSmoothingY;               /* for boot reps with REB */
      EquatingRecipes::Structures::UnivariateLogLinearSmoothing bootstrapUnivariateLogLinearSmoothingX; /* for boot reps with REL */
      EquatingRecipes::Structures::UnivariateLogLinearSmoothing bootstrapUnivariateLogLinearSmoothingY; /* for boot reps with REL */
      EquatingRecipes::Structures::BivariateLogLinearSmoothing bootstrapBivariateLogLinearSmoothingXY;  /* for boot reps with SEL */
      EquatingRecipes::Structures::BivariateLogLinearSmoothing bootstrapBivariateLogLinearSmoothingXV;  /* for boot reps with CEL */
      EquatingRecipes::Structures::BivariateLogLinearSmoothing bootstrapBivariateLogLinearSmoothingYV;  /* for boot reps with CEL */

      pData.numberOfBootstrapReplications = numberOfReplications;

      bootstrapInitializeEquatedRawScoreResults(pData, bootstrapEquatedRawScoreResults);

      if (bootstrapEquatedScaledScoresResults.has_value()) {
        EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults value = bootstrapEquatedScaledScoresResults.value();

        bootstraptInitializeEquatedScaledScores(pData, value);

        bootstrapEquatedScaledScoresResults = value;
      }

      for (size_t replicationNumber = 1; replicationNumber <= numberOfReplications; replicationNumber++) {
        pData.bootstrapReplicationNumber = replicationNumber; /* rep always stored in struct PDATA inall */

        /***** start for different designs/methods/smoothing *****/

        /* Notes: (a) Since struct PDATA has the name "in" here,
                  the name of struct PDATA must be "in" everywhere; 
              (b) when smoothing=='N', storage allocation done in
                  boot_USTATS() and boot_BSTATS(); otherwise,
                  storage allocation done in Parametric_boot_univ_BB(),
				  Parametric_boot_univ_ULL(), or Parametric_boot_biv().
              (c) Drivers here are the same as those for actual
                  equating; here the last three parameters for all of them 
                  must be: rep,inall,&br */

        /* nonparametric bootstrap procedures (no smoothing involved) */
        EquatingRecipes::RandomAndSingleGroupEquating rgSgEquating;
        EquatingRecipes::CGEquatingNoSmoothing cgEquatingNoSmoothing;
        EquatingRecipes::BetaBinomial betaBinomial;
        EquatingRecipes::LogLinearEquating logLinearEquating;

        if (pData.design == EquatingRecipes::Structures::Design::RANDOM_GROUPS &&
            pData.smoothing == EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED) {
          bootstrapUnivariateStatistics(pData.summaryRawDataX.value(),
                                        replicationNumber,
                                        bootstrapX);

          bootstrapUnivariateStatistics(pData.summaryRawDataY.value(),
                                        replicationNumber,
                                        bootstrapY);

          rgSgEquating.randomGroupEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                           pData.method,
                                           EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                           bootstrapX,
                                           bootstrapY,
                                           replicationNumber,
                                           pData,
                                           equatedRawScoreResults);
        } else if (pData.design == EquatingRecipes::Structures::Design::SINGLE_GROUP &&
                   pData.smoothing == EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED) {
          bootstrapBivariateStatistics(pData.summaryRawDataXY.value(), replicationNumber, bootstrapXY);

          rgSgEquating.singleGroupEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                           pData.method,
                                           EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                           bootstrapXY,
                                           replicationNumber,
                                           pData,
                                           equatedRawScoreResults);
        } else if (pData.design == EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS &&
                   pData.smoothing == EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED) {
          bootstrapBivariateStatistics(pData.summaryRawDataXV.value(), replicationNumber, bootstrapXV);
          bootstrapBivariateStatistics(pData.summaryRawDataYV.value(), replicationNumber, bootstrapYV);

          cgEquatingNoSmoothing.run(EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS,
                                    pData.method,
                                    EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                    pData.weightSyntheticPopulation1,
                                    pData.isInternalAnchor,
                                    pData.reliabilityCommonItemsPopulation1,
                                    pData.reliabilityCommonItemsPopulation2,
                                    bootstrapXV,
                                    bootstrapYV,
                                    replicationNumber,
                                    pData,
                                    equatedRawScoreResults);
        }

        /* parametric bootstrap procedures (involve smoothing) */

        else if (pData.design == EquatingRecipes::Structures::Design::RANDOM_GROUPS &&
                 pData.smoothing == EquatingRecipes::Structures::Smoothing::BETA_BINOMIAL) {
          EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomalSmoothingX;

          parametricBootstraptUnivariateBetaBinomial(betaBinomalSmoothingX,
                                                     replicationNumber,
                                                     bootstrapBetaBinomialSmoothingX);

          pData.betaBinomalSmoothingX = bootstrapBetaBinomialSmoothingX;

          EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomalSmoothingY;

          parametricBootstraptUnivariateBetaBinomial(betaBinomalSmoothingY,
                                                     replicationNumber,
                                                     bootstrapBetaBinomialSmoothingY);

          pData.betaBinomalSmoothingY = bootstrapBetaBinomialSmoothingY;

          equatedRawScoreResults = betaBinomial.randomGroupsEquipercentileEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                                                                   EquatingRecipes::Structures::Method::EQUIPERCENTILE,
                                                                                   EquatingRecipes::Structures::Smoothing::BETA_BINOMIAL,
                                                                                   pData.summaryRawDataX.value(),
                                                                                   pData.summaryRawDataY.value(),
                                                                                   bootstrapBetaBinomialSmoothingX,
                                                                                   bootstrapBetaBinomialSmoothingY,
                                                                                   replicationNumber,
                                                                                   pData);
        } else if (pData.design == EquatingRecipes::Structures::Design::RANDOM_GROUPS &&
                   pData.smoothing == EquatingRecipes::Structures::Smoothing::LOG_LINEAR) {
          parametricBootstraptUnivariateLogLinear(pData.univariateLogLinearSmoothingX.value(),
                                                  replicationNumber,
                                                  bootstrapUnivariateLogLinearSmoothingX);

          parametricBootstraptUnivariateLogLinear(pData.univariateLogLinearSmoothingY.value(),
                                                  replicationNumber,
                                                  bootstrapUnivariateLogLinearSmoothingY);

          logLinearEquating.runRGEquiEquatingWithLoglinearSmoothing(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                                                    EquatingRecipes::Structures::Method::EQUIPERCENTILE,
                                                                    EquatingRecipes::Structures::Smoothing::LOG_LINEAR,
                                                                    pData.summaryRawDataX.value(),
                                                                    pData.summaryRawDataY.value(),
                                                                    bootstrapUnivariateLogLinearSmoothingX,
                                                                    bootstrapUnivariateLogLinearSmoothingY,
                                                                    replicationNumber,
                                                                    pData,
                                                                    equatedRawScoreResults);

        } else if (pData.design == EquatingRecipes::Structures::Design::SINGLE_GROUP &&
                   pData.smoothing == EquatingRecipes::Structures::Smoothing::LOG_LINEAR) {
          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingXY;

          parametricBootstraptBivariate(bivariateLogLinearSmoothingXY,
                                        replicationNumber,
                                        bootstrapBivariateLogLinearSmoothingXY);

          pData.bivariateLogLinearSmoothingXY = bivariateLogLinearSmoothingXY;

          logLinearEquating.runSGEquiEquatingWithLoglinearSmoothing(EquatingRecipes::Structures::Design::SINGLE_GROUP,
                                                                    pData.method,
                                                                    EquatingRecipes::Structures::Smoothing::LOG_LINEAR,
                                                                    pData.summaryRawDataXY.value(),
                                                                    bootstrapBivariateLogLinearSmoothingXY,
                                                                    replicationNumber,
                                                                    pData,
                                                                    equatedRawScoreResults);
        } else if (pData.design == EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS &&
                   pData.smoothing == EquatingRecipes::Structures::Smoothing::LOG_LINEAR) {
          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingXV;

          parametricBootstraptBivariate(bivariateLogLinearSmoothingXV,
                                        replicationNumber,
                                        bootstrapBivariateLogLinearSmoothingXV);

          pData.bivariateLogLinearSmoothingXV = bivariateLogLinearSmoothingXV;

          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingYV;

          parametricBootstraptBivariate(bivariateLogLinearSmoothingYV,
                                        replicationNumber,
                                        bootstrapBivariateLogLinearSmoothingYV);

          pData.bivariateLogLinearSmoothingYV = bivariateLogLinearSmoothingYV;

          logLinearEquating.runCGEquiEquatingWithLoglinearSmoothing(EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS,
                                                                    pData.method,
                                                                    EquatingRecipes::Structures::Smoothing::LOG_LINEAR,
                                                                    pData.weightSyntheticPopulation1,
                                                                    pData.isInternalAnchor,
                                                                    pData.reliabilityCommonItemsPopulation1,
                                                                    pData.reliabilityCommonItemsPopulation2,
                                                                    pData.summaryRawDataXV.value(),
                                                                    pData.summaryRawDataYV.value(),
                                                                    bootstrapBivariateLogLinearSmoothingXV,
                                                                    bootstrapBivariateLogLinearSmoothingYV,
                                                                    replicationNumber,
                                                                    pData,
                                                                    equatedRawScoreResults);
        }

        /***** end for different designs/methods/smoothing *****/

        bootstrapAccumulateEquatedRawScores(pData,
                                            equatedRawScoreResults,
                                            bootstrapEquatedRawScoreResults);

        if (bootstrapEquatedScaledScoresResults.has_value()) {
          EquatingRecipes::Utilities::runEquatedScaledScores(pData,
                                                             equatedRawScoreResults,
                                                             pData.minimumRawScoreYct,
                                                             pData.maximumRawScoreYct,
                                                             pData.scoreIncrementYct,
                                                             pData.rawToScaledScoreTable.value(),
                                                             pData.roundToNumberOfDecimalPlaces,
                                                             pData.lowestObservableRoundedScaledScore,
                                                             pData.highestObservableRoundedScaledScore,
                                                             equatedScaledScoreResults);

          EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults btEqScaledScoreResults;
          bootstrapAccumulateEquatedScaledScores(pData,
                                                 equatedScaledScoreResults,
                                                 btEqScaledScoreResults);

          bootstrapEquatedScaledScoresResults = btEqScaledScoreResults;
        }

        bootstrapStandardErrorsEquatedRawScores(pData,
                                                bootstrapEquatedRawScoreResults);

        if (bootstrapEquatedScaledScoresResults.has_value()) {
          EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults btEquatedScaledScoresResults;
          bootstrapStandardErrorsEquatedScaledScores(pData,
                                                     btEquatedScaledScoresResults);

          bootstrapEquatedScaledScoresResults = btEquatedScaledScoresResults;
        }

        pData.bootstrapReplicationNumber = 0; /* bootstrap concluded */
      }
    }

  private:
    std::mt19937_64 gen;
    std::uniform_real_distribution<> distrib;

    /*
      Input
        struct PDATA inall
          
      Output  
        struct BOOT_ERAW_RESULTS t
          double mn[][] = matrix of sum and mean for scores
          double sd[][] = matrix of sum2 and sd for scores
          double bse[] = overall bootstrap se's 

      Function calls other than C or NR utilities:
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08          
    */
    void bootstrapInitializeEquatedRawScoreResults(EquatingRecipes::Structures::PData& pData,
                                                   EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& bootstrapEquatedRawScoreResults) {
      size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                 pData.minimumScoreX,
                                                                                 pData.scoreIncrementX);

      bootstrapEquatedRawScoreResults.sumAndMeanScores.setZero(pData.methods.size(), maximumScoreLocation + 1);
      bootstrapEquatedRawScoreResults.sumSquareAndSDScores.setZero(pData.methods.size(), maximumScoreLocation + 1);
      bootstrapEquatedRawScoreResults.bootstrapStandardErrors.setZero(pData.methods.size(), maximumScoreLocation + 1);
    }

    /*
      Based on USTATS x for actual equating, get a bootstrap sample 
      and store results in USTATS xb
      
      NOTE:  some initialization done here rather than in Boot_initialize-eraw()
            because otherwise there would have to be several versions of 
            Boot_initialize_eraw() for different designs. 
      
      Input: 
        struct USTATS x
        idum = seed
        rep = replication number
        
      Output
        struct USTATS xb:  xb must be different from x

      Function calls other than C or NR utilities:
        ran2()
        sort()
        cum_rel_freqs()
        perc_rank()
        score()
        MomentsFromFD()
        runerror()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08      
    */
    void bootstrapUnivariateStatistics(EquatingRecipes::Structures::UnivariateStatistics& x,
                                       const size_t& replicationNumber,
                                       EquatingRecipes::Structures::UnivariateStatistics& bootstrapX) {
      size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(x.maximumScore,
                                                                                 x.minimumScore,
                                                                                 x.scoreIncrement);

      if (replicationNumber == 1) {
        bootstrapX.numberOfExaminees = x.numberOfExaminees;
        bootstrapX.minimumScore = x.minimumScore;
        bootstrapX.maximumScore = x.maximumScore;
        bootstrapX.scoreIncrement = x.scoreIncrement;
        bootstrapX.numberOfScores = x.numberOfScores;
        bootstrapX.freqDist.resize(maximumScoreLocation + 1);
        bootstrapX.freqDistDouble.resize(maximumScoreLocation + 1);
        bootstrapX.cumulativeFreqDist.resize(maximumScoreLocation + 1);
        bootstrapX.relativeFreqDist.resize(maximumScoreLocation + 1);
        bootstrapX.cumulativeRelativeFreqDist.resize(maximumScoreLocation + 1);
        bootstrapX.percentileRankDist.resize(maximumScoreLocation + 1);
      }

      std::vector<int> pvec;

      bootstrapX.freqDist.setZero();

      for (size_t examineeIndex = 0; examineeIndex < bootstrapX.numberOfExaminees; examineeIndex++) {
        double randomValue = EquatingRecipes::Utilities::getUniformDoubleRandomNumber(this->distrib,
                                                                                      this->gen);

        int randomExamineeNumber = static_cast<int>(static_cast<double>(bootstrapX.numberOfExaminees) * randomValue + 1.0);
      }

      std::sort(pvec.begin(),
                pvec.end());

      /* get bootstrap fd */
      int sumFrequencies = 0;

      size_t randomExamineeIndex = 0;
      bool exitLoop = false;

      for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
        sumFrequencies += x.freqDist(scoreLocation);

        while (static_cast<int>(static_cast<double>(pvec[randomExamineeIndex++] + 0.00001) <= sumFrequencies)) {
          bootstrapX.freqDist(scoreLocation) += 1;

          if (randomExamineeIndex == bootstrapX.numberOfExaminees) {
            exitLoop = true;
            break;
          }
        }

        if (exitLoop) {
          break;
        } else {
          randomExamineeIndex--; /* gone 1 too far */
        }
      }

      /* get actual minimum and maximum scores in bootstrap data */
      for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
        if (bootstrapX.freqDist(scoreLocation) != 0) {
          bootstrapX.freqDistMinimumScore = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                 bootstrapX.minimumScore,
                                                                                 bootstrapX.scoreIncrement);
          break;
        }
      }

      for (size_t scoreLocation = maximumScoreLocation; scoreLocation >= 0; scoreLocation++) {
        if (bootstrapX.freqDist(scoreLocation) != 0) {
          bootstrapX.freqDistMaximumScore = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                 bootstrapX.minimumScore,
                                                                                 bootstrapX.scoreIncrement);

          break;
        }
      }

      /* get cfd, rfd, crfd, prd in bootstrap data */

      for (size_t scoreLocation = 1; scoreLocation <= maximumScoreLocation; scoreLocation++) {
        bootstrapX.cumulativeFreqDist(scoreLocation) = bootstrapX.cumulativeFreqDist(scoreLocation - 1) +
                                                       bootstrapX.freqDist(scoreLocation);
      }

      bootstrapX.freqDistDouble = bootstrapX.freqDist.cast<double>();
      bootstrapX.relativeFreqDist = bootstrapX.freqDistDouble / static_cast<double>(bootstrapX.numberOfExaminees);
      bootstrapX.cumulativeRelativeFreqDist = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(bootstrapX.minimumScore,
                                                                                                     bootstrapX.maximumScore,
                                                                                                     bootstrapX.scoreIncrement,
                                                                                                     bootstrapX.relativeFreqDist);
      bootstrapX.percentileRankDist = EquatingRecipes::Utilities::percentileRanks(bootstrapX.minimumScore,
                                                                                  bootstrapX.maximumScore,
                                                                                  bootstrapX.scoreIncrement,
                                                                                  bootstrapX.cumulativeRelativeFreqDist);

      /*  get moments */

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(bootstrapX.freqDistDouble,
                                                                                                                   bootstrapX.minimumScore,
                                                                                                                   bootstrapX.maximumScore,
                                                                                                                   bootstrapX.scoreIncrement);

      bootstrapX.momentValues = moments.momentValues;
    }

    void buildBootstrapMarginalAndBivariateFreqDists(const size_t& maximumScoreRow,
                                                     const size_t& maximumScoreColumn,
                                                     const std::vector<int>& pvec,
                                                     const EquatingRecipes::Structures::BivariateStatistics& xv,
                                                     EquatingRecipes::Structures::BivariateStatistics& bootstrapXV) {
      double sumBivarateFrequencies = 0.0;
      size_t selectedExamineeIndex = 0;

      for (size_t scoreLocationRow = 0; scoreLocationRow <= maximumScoreRow; scoreLocationRow++) {
        for (size_t scoreLocationColumn = 0; scoreLocationColumn <= maximumScoreColumn; scoreLocationColumn++) {
          sumBivarateFrequencies += xv.bivariateFreqDist(scoreLocationRow, scoreLocationColumn);

          while (static_cast<int>(static_cast<double>(pvec[selectedExamineeIndex]) + 0.00001) <= sumBivarateFrequencies) {
            bootstrapXV.bivariateFreqDist(scoreLocationRow, scoreLocationColumn) += 1;
            bootstrapXV.univariateStatisticsRow.freqDist(scoreLocationRow) += 1;
            bootstrapXV.univariateStatisticsColumn.freqDist(scoreLocationColumn) += 1;

            double scoreRow = EquatingRecipes::Utilities::getScore(scoreLocationRow,
                                                                   bootstrapXV.univariateStatisticsRow.minimumScore,
                                                                   bootstrapXV.univariateStatisticsRow.scoreIncrement);

            double scoreColumn = EquatingRecipes::Utilities::getScore(scoreLocationColumn,
                                                                      bootstrapXV.univariateStatisticsColumn.minimumScore,
                                                                      bootstrapXV.univariateStatisticsColumn.scoreIncrement);

            bootstrapXV.covariance += scoreRow * scoreColumn;

            if (selectedExamineeIndex + 1 == bootstrapXV.numberOfExaminees) {
              return;
            }

            selectedExamineeIndex++;
          }

          selectedExamineeIndex--; /* gone 1 too far */
        }
      }
    }

    /*
      Based on BSTATS xv for actual equating, get a bootstrap sample 
      and store results in BSTATS xvb
      
      NOTE:  some initialization done here rather than in Boot_initialize-eraw()
            because otherwise there would have to be several versions of 
            Boot_initialize_eraw() for different designs. 
      
      Input: 
        sort()
        struct BSTATS xv
        idum = seed
        rep = replication number
        
      Output
        struct BSTATS xvb:  xvb must be different from xv

      Function calls other than C or NR utilities:
        ran2()
        sort()
        cum_rel_freqs()
        perc_rank()
        score()
        MomentsFromFD()
        runerror()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08      
    */
    void bootstrapBivariateStatistics(EquatingRecipes::Structures::BivariateStatistics& xv,
                                      const size_t& replicationNumber,
                                      EquatingRecipes::Structures::BivariateStatistics& bootstrapXV) {
      size_t maximumScoreRow = EquatingRecipes::Utilities::getScoreLocation(xv.univariateStatisticsRow.maximumScore,
                                                                            xv.univariateStatisticsRow.minimumScore,
                                                                            xv.univariateStatisticsRow.scoreIncrement);

      size_t maximumScoreColumn = EquatingRecipes::Utilities::getScoreLocation(xv.univariateStatisticsColumn.maximumScore,
                                                                               xv.univariateStatisticsColumn.minimumScore,
                                                                               xv.univariateStatisticsColumn.scoreIncrement);

      if (replicationNumber == 1) {
        bootstrapXV.numberOfExaminees = xv.numberOfExaminees;
        bootstrapXV.univariateStatisticsRow.minimumScore = xv.univariateStatisticsRow.minimumScore;
        bootstrapXV.univariateStatisticsRow.maximumScore = xv.univariateStatisticsRow.maximumScore;
        bootstrapXV.univariateStatisticsRow.scoreIncrement = xv.univariateStatisticsRow.scoreIncrement;
        bootstrapXV.univariateStatisticsRow.numberOfScores = xv.univariateStatisticsRow.numberOfScores;
        bootstrapXV.univariateStatisticsColumn.minimumScore = xv.univariateStatisticsColumn.minimumScore;
        bootstrapXV.univariateStatisticsColumn.maximumScore = xv.univariateStatisticsColumn.maximumScore;
        bootstrapXV.univariateStatisticsColumn.scoreIncrement = xv.univariateStatisticsColumn.scoreIncrement;
        bootstrapXV.univariateStatisticsColumn.numberOfScores = xv.univariateStatisticsColumn.numberOfScores;

        bootstrapXV.univariateStatisticsRow.freqDist.resize(maximumScoreRow + 1);
        bootstrapXV.univariateStatisticsColumn.freqDist.resize(maximumScoreColumn + 1);
        bootstrapXV.bivariateFreqDist.resize(maximumScoreRow + 1, maximumScoreColumn + 1);

        bootstrapXV.univariateStatisticsRow.freqDistDouble.resize(maximumScoreRow + 1);
        bootstrapXV.univariateStatisticsColumn.freqDistDouble.resize(maximumScoreColumn + 1);
        bootstrapXV.bivariateFreqDistDouble.resize(maximumScoreRow + 1, maximumScoreColumn + 1);

        bootstrapXV.univariateStatisticsRow.cumulativeFreqDist.resize(maximumScoreRow + 1);
        bootstrapXV.univariateStatisticsRow.relativeFreqDist.resize(maximumScoreRow + 1);
        bootstrapXV.univariateStatisticsRow.cumulativeRelativeFreqDist.resize(maximumScoreRow + 1);
        bootstrapXV.univariateStatisticsRow.percentileRankDist.resize(maximumScoreRow + 1);

        bootstrapXV.univariateStatisticsColumn.cumulativeFreqDist.resize(maximumScoreColumn + 1);
        bootstrapXV.univariateStatisticsColumn.relativeFreqDist.resize(maximumScoreColumn + 1);
        bootstrapXV.univariateStatisticsColumn.cumulativeRelativeFreqDist.resize(maximumScoreColumn + 1);
        bootstrapXV.univariateStatisticsColumn.percentileRankDist.resize(maximumScoreColumn + 1);

        bootstrapXV.bivariateProportions.resize(maximumScoreRow + 1, maximumScoreColumn + 1);
      }

      std::vector<int> pvec;

      bootstrapXV.covariance = 0.0;
      bootstrapXV.univariateStatisticsRow.freqDist.setZero();
      bootstrapXV.univariateStatisticsColumn.freqDist.setZero();
      bootstrapXV.bivariateFreqDist.setZero();

      for (size_t examineeIndex = 1; examineeIndex <= bootstrapXV.numberOfExaminees; examineeIndex++) {
        double randomValue = EquatingRecipes::Utilities::getUniformDoubleRandomNumber(this->distrib,
                                                                                      this->gen);

        int selectedExaminee = static_cast<int>(static_cast<double>(bootstrapXV.numberOfExaminees) * randomValue + 1.0);

        pvec.push_back(selectedExaminee);
      }

      std::sort(pvec.begin(),
                pvec.end());

      buildBootstrapMarginalAndBivariateFreqDists(maximumScoreRow,
                                                  maximumScoreColumn,
                                                  pvec,
                                                  xv,
                                                  bootstrapXV);

      /*  get observed (actual) minimum and maximum scores in bootstrap data */
      for (size_t scoreLocation = 0; scoreLocation <= maximumScoreRow; scoreLocation++) {
        if (bootstrapXV.univariateStatisticsRow.freqDist(scoreLocation) != 0) {
          bootstrapXV.univariateStatisticsRow.freqDistMinimumScore = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                                          bootstrapXV.univariateStatisticsRow.minimumScore,
                                                                                                          bootstrapXV.univariateStatisticsRow.scoreIncrement);

          break;
        }
      }

      for (size_t scoreLocation = maximumScoreRow; scoreLocation >= 0; scoreLocation++) {
        if (bootstrapXV.univariateStatisticsRow.freqDist(scoreLocation) != 0) {
          bootstrapXV.univariateStatisticsRow.freqDistMaximumScore = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                                          bootstrapXV.univariateStatisticsRow.minimumScore,
                                                                                                          bootstrapXV.univariateStatisticsRow.scoreIncrement);
        }
      }

      for (size_t scoreLocation = 0; scoreLocation <= maximumScoreColumn; scoreLocation++) {
        if (bootstrapXV.univariateStatisticsColumn.freqDist(scoreLocation) != 0) {
          bootstrapXV.univariateStatisticsColumn.freqDistMinimumScore = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                                             bootstrapXV.univariateStatisticsColumn.minimumScore,
                                                                                                             bootstrapXV.univariateStatisticsColumn.scoreIncrement);

          break;
        }
      }

      for (size_t scoreLocation = maximumScoreColumn; scoreLocation >= 0; scoreLocation++) {
        if (bootstrapXV.univariateStatisticsColumn.freqDist(scoreLocation) != 0) {
          bootstrapXV.univariateStatisticsColumn.freqDistMaximumScore = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                                                             bootstrapXV.univariateStatisticsColumn.minimumScore,
                                                                                                             bootstrapXV.univariateStatisticsColumn.scoreIncrement);
        }
      }

      /* get cfd1,rfd1,crfd1,prd1, and cfd2,rfd2,crfd2,prd2 in boot data */

      bootstrapXV.univariateStatisticsRow.cumulativeFreqDist(0) = bootstrapXV.univariateStatisticsRow.freqDist(0);
      for (size_t scoreLocation = maximumScoreRow; scoreLocation >= 0; scoreLocation++) {
        bootstrapXV.univariateStatisticsRow.cumulativeFreqDist(scoreLocation) = bootstrapXV.univariateStatisticsRow.cumulativeFreqDist(scoreLocation - 1) +
                                                                                bootstrapXV.univariateStatisticsRow.freqDist(scoreLocation);
      }

      bootstrapXV.univariateStatisticsRow.freqDistDouble = bootstrapXV.univariateStatisticsRow.freqDist.cast<double>();
      bootstrapXV.univariateStatisticsRow.relativeFreqDist = bootstrapXV.univariateStatisticsRow.freqDistDouble / static_cast<double>(bootstrapXV.numberOfExaminees);
      bootstrapXV.univariateStatisticsRow.cumulativeRelativeFreqDist = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(bootstrapXV.univariateStatisticsRow.minimumScore,
                                                                                                                              bootstrapXV.univariateStatisticsRow.maximumScore,
                                                                                                                              bootstrapXV.univariateStatisticsRow.scoreIncrement,
                                                                                                                              bootstrapXV.univariateStatisticsRow.relativeFreqDist);
      bootstrapXV.univariateStatisticsRow.percentileRankDist = EquatingRecipes::Utilities::percentileRanks(bootstrapXV.univariateStatisticsRow.minimumScore,
                                                                                                           bootstrapXV.univariateStatisticsRow.maximumScore,
                                                                                                           bootstrapXV.univariateStatisticsRow.scoreIncrement,
                                                                                                           bootstrapXV.univariateStatisticsRow.cumulativeRelativeFreqDist);

      bootstrapXV.univariateStatisticsColumn.cumulativeFreqDist(0) = bootstrapXV.univariateStatisticsColumn.freqDist(0);
      for (size_t scoreLocation = maximumScoreColumn; scoreLocation >= 0; scoreLocation++) {
        bootstrapXV.univariateStatisticsColumn.cumulativeFreqDist(scoreLocation) = bootstrapXV.univariateStatisticsColumn.cumulativeFreqDist(scoreLocation - 1) +
                                                                                   bootstrapXV.univariateStatisticsColumn.freqDist(scoreLocation);
      }

      bootstrapXV.univariateStatisticsColumn.freqDistDouble = bootstrapXV.univariateStatisticsColumn.freqDist.cast<double>();
      bootstrapXV.univariateStatisticsColumn.relativeFreqDist = bootstrapXV.univariateStatisticsColumn.freqDistDouble / static_cast<double>(bootstrapXV.numberOfExaminees);
      bootstrapXV.univariateStatisticsColumn.cumulativeRelativeFreqDist = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(bootstrapXV.univariateStatisticsColumn.minimumScore,
                                                                                                                                 bootstrapXV.univariateStatisticsColumn.maximumScore,
                                                                                                                                 bootstrapXV.univariateStatisticsColumn.scoreIncrement,
                                                                                                                                 bootstrapXV.univariateStatisticsColumn.relativeFreqDist);
      bootstrapXV.univariateStatisticsColumn.percentileRankDist = EquatingRecipes::Utilities::percentileRanks(bootstrapXV.univariateStatisticsColumn.minimumScore,
                                                                                                              bootstrapXV.univariateStatisticsColumn.maximumScore,
                                                                                                              bootstrapXV.univariateStatisticsColumn.scoreIncrement,
                                                                                                              bootstrapXV.univariateStatisticsColumn.cumulativeRelativeFreqDist);

      /*  get moments, cov, and corr in bootstrap data */
      EquatingRecipes::Structures::Moments rowMoments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(bootstrapXV.univariateStatisticsRow.freqDistDouble,
                                                                                                                      bootstrapXV.univariateStatisticsRow.minimumScore,
                                                                                                                      bootstrapXV.univariateStatisticsRow.maximumScore,
                                                                                                                      bootstrapXV.univariateStatisticsRow.scoreIncrement);

      bootstrapXV.univariateStatisticsRow.momentValues = rowMoments.momentValues;

      EquatingRecipes::Structures::Moments columnMoments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(bootstrapXV.univariateStatisticsColumn.freqDistDouble,
                                                                                                                         bootstrapXV.univariateStatisticsColumn.minimumScore,
                                                                                                                         bootstrapXV.univariateStatisticsColumn.maximumScore,
                                                                                                                         bootstrapXV.univariateStatisticsColumn.scoreIncrement);

      bootstrapXV.univariateStatisticsColumn.momentValues = columnMoments.momentValues;

      bootstrapXV.covariance = (bootstrapXV.covariance / static_cast<double>(bootstrapXV.numberOfExaminees)) -
                               (bootstrapXV.univariateStatisticsRow.momentValues(0) *
                                bootstrapXV.univariateStatisticsColumn.momentValues(0));

      bootstrapXV.correlation = bootstrapXV.covariance / (bootstrapXV.univariateStatisticsRow.momentValues(1) *
                                                          bootstrapXV.univariateStatisticsColumn.momentValues(1));

      /* bivariate proportions for frequency estimation, as well as
      double version of xvb->bfd[][] */
      for (size_t scoreLocationColumn = 0; scoreLocationColumn < bootstrapXV.univariateStatisticsColumn.numberOfScores; scoreLocationColumn++) {
        for (size_t scoreLocationRow = 0; scoreLocationRow < bootstrapXV.univariateStatisticsRow.numberOfScores; scoreLocationRow++) {
          bootstrapXV.bivariateProportions(scoreLocationRow, scoreLocationColumn) = static_cast<double>(bootstrapXV.bivariateFreqDist(scoreLocationRow, scoreLocationColumn)) /
                                                                                    static_cast<double>(bootstrapXV.numberOfExaminees);
        }
      }

      bootstrapXV.bivariateFreqDistDouble = bootstrapXV.bivariateFreqDist.cast<double>();
    }

    /*
      accumulate b->eraw[][] in t->mn[][]
      accumulate b->eraw[][]*b->eraw[][] in t->sd[][]
      
      
      Input
        struct PDATA inall: following variables are used
          nm = number of methods 
          min = min raw score 
          max = max raw score 
          inc = raw score increment
          
        struct ERAW_RESULTS b 
          double eraw[][] --- for a boot rep, not actual equating
      
      Output   
        struct BOOT_ERAW_RESULTS s
          double **mn: for matrix of sum and mean for scores 
          double **sd; for matrix of sum2 and sd for scores

      Function calls other than C or NR utilities: 
        loc()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08       
    */
    void bootstrapAccumulateEquatedRawScores(EquatingRecipes::Structures::PData& pData,
                                             EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults,
                                             EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& bootstrapEquatedRawScoreResults) {
      size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                 pData.minimumScoreX,
                                                                                 pData.scoreIncrementX);

      for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
        for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
          bootstrapEquatedRawScoreResults.sumAndMeanScores(methodIndex, scoreLocation) +=
              equatedRawScoreResults.equatedRawScores(methodIndex, scoreLocation);

          bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) +=
              std::pow(equatedRawScoreResults.equatedRawScores(methodIndex, scoreLocation), 2);
        }
      }
    }

    /*
      final computations for s->mn[][] and s->sd[][] == bootstrap se's
      
      Input
        struct PDATA inall: following variables are used
          nm = number of methods 
          min = min raw score
          max = max raw score
          inc = raw score increment
          fdx[] = fd for x
          nrep = number of replications
          n = sample size associated with fdx[] 
          
        struct BOOT_ERAW_RESULTS t
          mn[][] = for matrix of sum's for scores 
          sd[][] = for matrix of sum2's for scores 
      
      Output   
        struct BOOT_ERAW_RESULTS t
          mn[][] = for matrix of mean for scores 
          sd[][] = for matrix of sd for scores 
          bse[] = overall bootstrap se's (one for each method)

      NOTE: bootstrap se's require a divisor of nrep-1, not nrep;
            see "divisor" statement in code

      Function calls other than C or NR utilities: 
        loc()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08  
    */
    void bootstrapStandardErrorsEquatedRawScores(EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& bootstrapEquatedRawScoreResults) {
      size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                 pData.minimumScoreX,
                                                                                 pData.scoreIncrementX);

      for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
        for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
          bootstrapEquatedRawScoreResults.sumAndMeanScores(methodIndex, scoreLocation) /= static_cast<double>(pData.numberOfBootstrapReplications);
          bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) =
              bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) / static_cast<double>(pData.numberOfBootstrapReplications);

          bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) -=
              std::pow(bootstrapEquatedRawScoreResults.sumAndMeanScores(methodIndex, scoreLocation), 2);

          bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) *=
              static_cast<double>(pData.numberOfBootstrapReplications) / static_cast<double>(pData.numberOfBootstrapReplications - 1);

          bootstrapEquatedRawScoreResults.bootstrapStandardErrors(methodIndex) += pData.scoreFrequenciesX(scoreLocation) *
                                                                                  bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation);

          if (bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) > 0.00001) {
            bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) =
                std::sqrt(bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation));
          } else {
            bootstrapEquatedRawScoreResults.sumSquareAndSDScores(methodIndex, scoreLocation) = 0.0;
          }
        }

        bootstrapEquatedRawScoreResults.bootstrapStandardErrors(methodIndex) = std::sqrt(bootstrapEquatedRawScoreResults.bootstrapStandardErrors(methodIndex)) /
                                                                               static_cast<double>(pData.numberOfExaminees);
      }
    }

    // void Print_Boot_se_eraw(FILE *fp, char tt[], struct PDATA *inall,
    //                         struct ERAW_RESULTS *r,
    //                         struct BOOT_ERAW_RESULTS *t, int mdiff);

    /*
      Input
        
        struct PDATA inall: variables used:
          nm = number of methods
          min = minimum raw score for x 
          max = maximum raw score for x
        inc = increment in raw scores for x
      
      Output             
        
        struct BOOT_ESS_RESULTS u: variables used:
          int rep = replication number
          mnu[][] = matrix of sum and mean for unrounded scale scores
          sdu[][] = matrix of sum2 and sd for unrounded scale scores
          bseu[] = overall bootstrap se's for unrounded scale scales 
          mnr[][] = matrix of sum and mean for rounded scale scores
          sdr[[]] = matrix of sum2 and sd for rounded scale scores
          bser[] = overall bootstrap se's for rounded scale scales

      Function calls other than C or NR utilities: 
        loc()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08    
    */
    void bootstraptInitializeEquatedScaledScores(EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& bootstrapEquatedScaledScoresResults) {
      size_t numberOfMethods = pData.methods.size();
      size_t numberOfScores = EquatingRecipes::Utilities::getNumberOfScores(pData.maximumScoreX,
                                                                            pData.minimumScoreX,
                                                                            pData.scoreIncrementX);

      bootstrapEquatedScaledScoresResults.unroundedScaledScoresSumsAndMeans.setZero(numberOfMethods, numberOfScores);
      bootstrapEquatedScaledScoresResults.unroundedScaledScoresSumSquaresAndSDs.setZero(numberOfMethods, numberOfScores);
      bootstrapEquatedScaledScoresResults.unroundedScaledScoresBoostrapStandardErrors.setZero(numberOfScores);

      bootstrapEquatedScaledScoresResults.roundedScaledScoresSumsAndMeans.setZero(numberOfMethods, numberOfScores);
      bootstrapEquatedScaledScoresResults.roundedScaledScoresSumSquaresAndSDs.setZero(numberOfMethods, numberOfScores);
      bootstrapEquatedScaledScoresResults.roundedScaledScoresBoostrapStandardErrors.setZero(numberOfScores);
    }

    /*
      accumulate s->essu[][] and store in u->mnu[][]
      accumulate s->essu[][]*s->essu[][] and store in u->sdu[][]
      accumulate s->essr[][] and store in s->mnr[][]
      accumulate s->essr[][]*s->essr[][] and store in u->sdr[][]
      
      Input
        
        struct PDATA inall: variables used:
          nm = number of methods
          min = min raw score for x 
          max = max raw score for x 
          inc = raw score increment for x          
          round: if round = i, then round to i-th place
      
        struct ESS_RESULTS s (must be same struct as in Wrapper_ESS()
                              within the bootstrap loop)
    
          essu[][] = unrounded equated scale scores for a replication
          essr[][] = rounded equated scale scores for a replication
      
      Output
        
        struct Boot_ess_results u (must be same struct as in 
                                  Boot_initialize_ess())
                                
        rep = replication number for bootstrap
        mnu[][] = sum and mean for unrounded scale scores
        sdu[][] = sum2 and sd for unrounded scale scores 
        mnr[][] = sum and mean for rounded scale scores
        sdr[][] = sum2 and sd for rounded scale scores

      Function calls other than C or NR utilities: None
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08    
    */
    void bootstrapAccumulateEquatedScaledScores(EquatingRecipes::Structures::PData& pData,
                                                EquatingRecipes::Structures::EquatedScaledScoresResults& equatedScaledScoresResults,
                                                EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& bootstrapEquatedScaledScoresResults) {
      size_t numberOfScores = EquatingRecipes::Utilities::getNumberOfScores(pData.maximumScoreX,
                                                                            pData.minimumScoreX,
                                                                            pData.scoreIncrementX);

      for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          bootstrapEquatedScaledScoresResults.unroundedScaledScoresSumsAndMeans(methodIndex, scoreIndex) +=
              equatedScaledScoresResults.unroundedEquatedScaledScores(methodIndex, scoreIndex);

          bootstrapEquatedScaledScoresResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) +=
              std::pow(equatedScaledScoresResults.unroundedEquatedScaledScores(methodIndex, scoreIndex), 2);

          if (pData.roundToNumberOfDecimalPlaces > 0) {
            bootstrapEquatedScaledScoresResults.roundedScaledScoresSumsAndMeans(methodIndex, scoreIndex) +=
                equatedScaledScoresResults.roundedEquatedScaledScores(methodIndex, scoreIndex);

            bootstrapEquatedScaledScoresResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) +=
                std::pow(equatedScaledScoresResults.roundedEquatedScaledScores(methodIndex, scoreIndex), 2);
          }
        }
      }
    }

    /*
    final computations for
      u->muu[][] and u->sdu[][] == u->bseu[][] = bootstrap unrounded se's 
      u->mur[][] and u->sdr[][] == u->bser[][] = bootstrap rounded se's
        
    Input
      
      struct PDATA inall: variables used:
        nm = number of methods
        min = min raw score for x 
        max = max raw score for x
        inc = raw score increment for x
        nrep = number of replications
        fdx[] = FD for old form x
        n = sample size associated with fdx[]
        round: if round = i, then round to i-th place
      
      struct BOOT_ESS_RESULTS u (must be same struct as in 
                                  Boot_initialize_ess() and
                                  Boot_accumulate_ess())
    
        mnu[][] = sums for unrounded scale scores
        sdu[][] = sum2's for unrounded scale scores 
        mnr[][] = sums for rounded scale scores
        sdr[][] = sum2's for rounded scale scores
        bser[] = overall bootstrap se's for rounded scale scales
        
    Output
      struct BOOT_ESS_RESULTS u (must be same struct as in 
                                  Boot_initialize_ess() and
                                  Boot_accumulate_ess())
    
        mnu[][] = means for unrounded scale scores
        sdu[][] = sd's for unrounded scale scores
        bseu[] = overall bootstrap se's for unrounded scale scales 
        mnr[][] = means for rounded scale scores
        sdr[][] = sd's for rounded scale scores
        bser[] = overall bootstrap se's for rounded scale scales

    NOTE: bootstrap se's require a divisor of nrep-1, not nrep;
          see "divisor" statements in code

    Function calls other than C or NR utilities:
      loc()
      score()
                                                  
    R. L. Brennan

    Date of last revision: 6/30/08    
  */
    void bootstrapStandardErrorsEquatedScaledScores(EquatingRecipes::Structures::PData& pData,
                                                    EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& bootstrapEqautedScaledScoreResults) {
      size_t numberOfScores = EquatingRecipes::Utilities::getNumberOfScores(pData.maximumScoreX,
                                                                            pData.minimumScoreX,
                                                                            pData.scoreIncrementX);
      for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
        /* for unrounded scale scores */
        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumsAndMeans(methodIndex, scoreIndex) /= static_cast<double>(pData.numberOfBootstrapReplications);

          bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) =
              bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) / static_cast<double>(pData.numberOfBootstrapReplications);

          bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) -=
              std::pow(bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumsAndMeans(methodIndex, scoreIndex), 2);

          bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) *=
              (static_cast<double>(pData.numberOfBootstrapReplications)) / (static_cast<double>(pData.numberOfBootstrapReplications - 1));

          if (bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) > 0.00001) {
            bootstrapEqautedScaledScoreResults.unroundedScaledScoresBoostrapStandardErrors(methodIndex) +=
                pData.scoreFrequenciesX(scoreIndex) * bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex);
          }

          if (bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) > 0.00001) {
            bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) =
                std::sqrt(bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex));
          } else {
            bootstrapEqautedScaledScoreResults.unroundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) = 0.0;
          }
        }

        bootstrapEqautedScaledScoreResults.unroundedScaledScoresBoostrapStandardErrors(methodIndex) =
            std::sqrt(bootstrapEqautedScaledScoreResults.unroundedScaledScoresBoostrapStandardErrors(methodIndex) /
                      static_cast<double>(pData.numberOfExaminees));

        /* for rounded scale scores */
        if (pData.roundToNumberOfDecimalPlaces > 0) {
          for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
            bootstrapEqautedScaledScoreResults.roundedScaledScoresSumsAndMeans(methodIndex, scoreIndex) /= static_cast<double>(pData.numberOfBootstrapReplications);

            bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) =
                bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) /
                static_cast<double>(pData.numberOfBootstrapReplications);

            bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) -=
                std::pow(bootstrapEqautedScaledScoreResults.roundedScaledScoresSumsAndMeans(methodIndex, scoreIndex), 2);

            bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) *=
                static_cast<double>(pData.numberOfBootstrapReplications) / static_cast<double>(pData.numberOfBootstrapReplications - 1);

            if (bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) > 0.00001) {
              bootstrapEqautedScaledScoreResults.roundedScaledScoresBoostrapStandardErrors(methodIndex) += pData.scoreFrequenciesX(scoreIndex) *
                                                                                                           bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex);
            }

            if (bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) > 0.00001) {
              bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) =
                  std::sqrt(bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex));
            } else {
              bootstrapEqautedScaledScoreResults.roundedScaledScoresSumSquaresAndSDs(methodIndex, scoreIndex) = 0.0;
            }
          }

          bootstrapEqautedScaledScoreResults.roundedScaledScoresBoostrapStandardErrors(methodIndex) =
              std::sqrt(bootstrapEqautedScaledScoreResults.roundedScaledScoresBoostrapStandardErrors(methodIndex) /
                        static_cast<double>(pData.numberOfExaminees));
        }
      }
    }

    // void Print_Boot_se_ess(FILE *fp, char tt[], struct PDATA *inall,
    //                        struct ESS_RESULTS *s,
    //                        struct BOOT_ESS_RESULTS *u, int mdiff);

    /*
      Parametric bootstrap sample for a univariate distribution

      Calling sequence from Wrapper_Bootstrap():
        Parametric_boot_univ_BB(inall->bbx, idum, rep, bt_bbx)

      Input
        x    = struct BB_SMOOTH for actual fitted distribution;
              taken as population distribution
        idum = seed
      rep  = replication number

      Output
        btx = struct BB_SMOOTH for bootstrap sample from x; only elements
            needed are btx->crfd and btx->prd

      Function calls other than C or NR utilities:
        ran2()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08    
    */
    void parametricBootstraptUnivariateBetaBinomial(EquatingRecipes::Structures::BetaBinomialSmoothing& x,
                                                    size_t& replicationNumber,
                                                    EquatingRecipes::Structures::BetaBinomialSmoothing& bootstrapX) {
      size_t numberOfScores = x.numberOfItems + 1;

      if (replicationNumber == 1) {
        bootstrapX.fittedRawScoreCumulativeRelativeFreqDist.resize(numberOfScores);
        bootstrapX.fittedRawScorePercentileRankDist.resize(numberOfScores);
      }

      bootstrapX.fittedRawScoreCumulativeRelativeFreqDist.setZero();

      for (size_t examineeIndex = 0; examineeIndex < x.numberOfExaminees; examineeIndex++) {
        double randomValue = EquatingRecipes::Utilities::getUniformDoubleRandomNumber(this->distrib,
                                                                                      this->gen);
        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          if (randomValue < x.fittedRawScoreCumulativeRelativeFreqDist(scoreIndex)) {
            bootstrapX.fittedRawScoreCumulativeRelativeFreqDist(scoreIndex) += 1;
            break;
          }
        }
      }

      /* here btx->crfd[] is frequencies; 
       next two statemenst makes btx->crfd[] cum rel freqs */

      bootstrapX.fittedRawScoreCumulativeRelativeFreqDist /= static_cast<double>(x.numberOfExaminees);

      for (size_t scoreIndex = 1; scoreIndex < numberOfScores; scoreIndex++) {
        bootstrapX.fittedRawScoreCumulativeRelativeFreqDist(scoreIndex) += bootstrapX.fittedRawScoreCumulativeRelativeFreqDist(scoreIndex - 1);
      }

      bootstrapX.fittedRawScorePercentileRankDist = EquatingRecipes::Utilities::percentileRanks(0,
                                                                                                numberOfScores - 1,
                                                                                                1,
                                                                                                bootstrapX.fittedRawScoreCumulativeRelativeFreqDist);
    }

    /*
      Parametric bootstrap sample for a univariate distribution

      Calling sequence from Wrapper_Bootstrap():
        Parametric_boot_univ_ULL(inall->ullx, idum, rep, bt_ullx)

      Input
        x    = struct ULL_SMOOTH for actual fitted distribution;
              taken as population distribution
        idum = seed
      rep  = replication number

      Output
        btx = struct ULL_SMOOTH for bootstrap sample from x; only elements
            needed are btx->crfd and btx->prd

      Function calls other than C or NR utilities:
        ran2()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08    
    */
    void
    parametricBootstraptUnivariateLogLinear(EquatingRecipes::Structures::UnivariateLogLinearSmoothing& x,
                                            size_t& replicationNumber,
                                            EquatingRecipes::Structures::UnivariateLogLinearSmoothing& bootstrapX) {
      if (replicationNumber == 1) {
        bootstrapX.fittedRawScoreCumulativeRelativeDist.resize(x.numberOfScores);
        bootstrapX.fittedRawScorePercentileRankDist.resize(x.numberOfScores);
      }

      bootstrapX.fittedRawScoreCumulativeRelativeDist.setZero();

      for (size_t examineeIndex = 1; examineeIndex <= x.numberOfExaminees; examineeIndex++) {
        double randomValue = EquatingRecipes::Utilities::getUniformDoubleRandomNumber(this->distrib,
                                                                                      this->gen);

        for (size_t scoreIndex = 0; scoreIndex < x.numberOfScores; scoreIndex++) {
          if (randomValue < x.fittedRawScoreCumulativeRelativeDist(scoreIndex)) {
            bootstrapX.fittedRawScoreCumulativeRelativeDist(scoreIndex) += 1;
            break;
          }
        }
      }

      /* here btx->crfd[] is frequencies; 
      next two statemenst makes btx->crfd[] cum rel freqs */
      bootstrapX.fittedRawScoreCumulativeRelativeDist /= static_cast<double>(x.numberOfExaminees);

      for (size_t scoreIndex = 1; scoreIndex < x.numberOfScores; scoreIndex++) {
        bootstrapX.fittedRawScoreCumulativeRelativeDist(scoreIndex) += bootstrapX.fittedRawScoreCumulativeRelativeDist(scoreIndex - 1);
      }

      bootstrapX.fittedRawScorePercentileRankDist = EquatingRecipes::Utilities::percentileRanks(0,
                                                                                                x.numberOfScores - 1,
                                                                                                1,
                                                                                                bootstrapX.fittedRawScoreCumulativeRelativeDist);
    }

    /*
      Parametric bootstrap sample for a bivariate distribution

      Input
        xv      = BLL_SMOOTH for original data (smoothed);  
                xv->brd[][] is taken as population distribution
        idum    = seed
      rep     = replication number

      Output
        btxv   = struct BLL_SMOOTH for bootstrap sample; not
                all data elements are needed

      Function calls other than C or NR utilities:
        ran2()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08    
    */
    void parametricBootstraptBivariate(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& xv,
                                       size_t& replicationNumber,
                                       EquatingRecipes::Structures::BivariateLogLinearSmoothing& bootstrapXV) {
      size_t totalNumberOfScores = xv.numberOfScoresX * xv.numberOfScoresV;

      if (replicationNumber == 1) {
        bootstrapXV.isInternalAnchor = xv.isInternalAnchor;
        bootstrapXV.numberOfExaminees = xv.numberOfExaminees;
        bootstrapXV.numberOfScoresX = xv.numberOfScoresX;
        bootstrapXV.minimumRawScoreX = xv.minimumRawScoreX;
        bootstrapXV.scoreIncrementX = xv.scoreIncrementX;
        bootstrapXV.numberOfScoresV = xv.numberOfScoresV;
        bootstrapXV.minimumRawScoreV = xv.minimumRawScoreV;
        bootstrapXV.scoreIncrementV = xv.scoreIncrementV;
        bootstrapXV.totalNumberOfScores = xv.totalNumberOfScores;

        bootstrapXV.fittedBivariateFreqDist.resize(xv.numberOfScoresX, xv.numberOfScoresV);

        bootstrapXV.fittedFrequencesX.resize(xv.numberOfScoresX);
        bootstrapXV.fittedRawScoreDensityX.resize(xv.numberOfScoresX);
        bootstrapXV.fittedRawScoreCumulativeRelativeFreqDistX.resize(xv.numberOfScoresX);
        bootstrapXV.fittedRawScorePercentileRankDistX.resize(xv.numberOfScoresX);

        bootstrapXV.fittedFrequencesV.resize(xv.numberOfScoresV);
        bootstrapXV.fittedRawScoreDensityV.resize(xv.numberOfScoresV);
        bootstrapXV.fittedRawScoreCumulativeRelativeFreqDistV.resize(xv.numberOfScoresV);
        bootstrapXV.fittedRawScorePercentileRankDistV.resize(xv.numberOfScoresV);

        bootstrapXV.cumulativeRelativeFreqDistRowMajorVector.resize(totalNumberOfScores);
        bootstrapXV.fittedBivariateRelativeFreqDistXV.resize(xv.numberOfScoresX, xv.numberOfScoresV);
      }

      bootstrapXV.cumulativeRelativeFreqDistRowMajorVector.setZero();

      /* In the following code, xv->crfd_vector_bfd[] is the crfd  for the 
        fitted actual data, while btxv->crfd_vector_bfd[] is the bootstrap fd 
      at end of for loop */
      for (size_t examineeIndex = 1; examineeIndex <= xv.numberOfExaminees; examineeIndex++) {
        double randomValue = EquatingRecipes::Utilities::getUniformDoubleRandomNumber(this->distrib,
                                                                                      this->gen);

        for (size_t index = 0; index < totalNumberOfScores; index++) {
          if (randomValue < xv.cumulativeRelativeFreqDistRowMajorVector(index)) {
            bootstrapXV.cumulativeRelativeFreqDistRowMajorVector(index)++;
            break;
          }
        }
      }

      /* Here btxv->crfd_vector_bfd[] is bootstrap frequencies in vector format; 
      next statements put it in matrix format and store it in btxv->bfd[][] */
      size_t index = 0;
      for (size_t rowIndex = 0; rowIndex < xv.numberOfScoresX; rowIndex++) {
        for (size_t columnIndex = 0; columnIndex < xv.numberOfScoresV; columnIndex++) {
          bootstrapXV.fittedBivariateFreqDist(rowIndex, columnIndex) = bootstrapXV.cumulativeRelativeFreqDistRowMajorVector(index);

          index++;
        }
      }

      /* get row and col marginal fd, rfd, crfd, and PR */
      bootstrapXV.fittedFrequencesV.setZero();
      bootstrapXV.fittedFrequencesX.setZero();
      index = 0;
      for (size_t rowIndex = 0; rowIndex < xv.numberOfScoresX; rowIndex++) {
        for (size_t columnIndex = 0; columnIndex < xv.numberOfScoresV; columnIndex++) {
          bootstrapXV.fittedFrequencesX(rowIndex) += bootstrapXV.fittedBivariateFreqDist(rowIndex, columnIndex);
          bootstrapXV.fittedFrequencesV(columnIndex) += bootstrapXV.fittedBivariateFreqDist(rowIndex, columnIndex);
        }

        bootstrapXV.fittedRawScoreDensityX(rowIndex) = bootstrapXV.fittedFrequencesX(rowIndex) / static_cast<double>(xv.numberOfExaminees);
      }

      bootstrapXV.fittedRawScoreDensityV = bootstrapXV.fittedFrequencesV / static_cast<double>(xv.numberOfExaminees);

      bootstrapXV.fittedRawScoreCumulativeRelativeFreqDistX = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(0,
                                                                                                                     xv.numberOfScoresX - 1,
                                                                                                                     1,
                                                                                                                     bootstrapXV.fittedRawScoreDensityX);

      bootstrapXV.fittedRawScorePercentileRankDistX = EquatingRecipes::Utilities::percentileRanks(0,
                                                                                                  xv.numberOfScoresX - 1,
                                                                                                  1,
                                                                                                  bootstrapXV.fittedRawScoreDensityX);

      bootstrapXV.fittedRawScoreCumulativeRelativeFreqDistV = EquatingRecipes::Utilities::cumulativeRelativeFreqDist(0,
                                                                                                                     xv.numberOfScoresV - 1,
                                                                                                                     1,
                                                                                                                     bootstrapXV.fittedRawScoreDensityV);

      bootstrapXV.fittedRawScorePercentileRankDistV = EquatingRecipes::Utilities::percentileRanks(0,
                                                                                                  xv.numberOfScoresV - 1,
                                                                                                  1,
                                                                                                  bootstrapXV.fittedRawScoreDensityV);

      /* following code gets rel fd version of btxv->bfd[]; i.e., 
      btxv->brfd[][] is the bootstrap rel freq biv dist for x by v; 
	    needed for FEorMFE_EE() */

      bootstrapXV.fittedBivariateRelativeFreqDistXV = bootstrapXV.fittedBivariateFreqDist / static_cast<double>(xv.numberOfExaminees);
    }
  };
} // namespace EquatingRecipes

#endif