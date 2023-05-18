/*
  struct USTATS x,y;  
  struct PDATA pdREN;
  struct ERAW_RESULTS rREN;
  struct ESS_RESULTS sREN;
  struct BOOT_ERAW_RESULTS tREN;
  struct BOOT_ESS_RESULTS uREN;

  long idum = -15L;            
  FILE *outf;                                         
                   
  outf = fopen("Chap 8 out","w");
   
  Random Groups Design: 
     Kolen and Brennan (2004): Chapter 2 example: 
     Bootstrap standard errors (see also p. 238)

  ReadFdGet_USTATS("actmathfreq.dat",1,2,0,40,1,'X',&x);
  ReadFdGet_USTATS("actmathfreq.dat",1,3,0,40,1,'Y',&y);

  Wrapper_RN('R','E','N',&x,&y,0,&pdREN,&rREN);
  Wrapper_ESS(&pdREN,&rREN,0,40,1,"yctmath.TXT",1,1,36,&sREN); 
  
  Wrapper_Bootstrap(&pdREN,1000,&idum,&tREN,&uREN);
  Print_Boot_se_eraw(outf,"ACT Math---Equipercentile",
                     &pdREN,&rREN,&tREN,0);
  Print_Boot_se_ess(outf,"ACT Math---Equipercentile",
                    &pdREN,&sREN,&uREN,0);
*/

#ifndef TESTS_EXAMPLES_CHAPTER_8_HPP
#define TESTS_EXAMPLES_CHAPTER_8_HPP

#include <equating_recipes/analyses/bootstrap.hpp>
#include <equating_recipes/analyses/equated_scaled_scores.hpp>
#include <equating_recipes/analyses/random_groups_equating.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/json/json_document.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/rg_and_sg_equating.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>

#include "datasets/actmathfreq.hpp"
#include "datasets/yctmath.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter8 {
        void operator()() {
          /* Random Groups Design: 
          Kolen and Brennan (2004): Chapter 2 example: 
          Bootstrap standard errors (see also p. 238) */

          EquatingRecipes::Tests::Examples::Datasets::ACTMathFreq actMathFreq;
          EquatingRecipes::Tests::Examples::Datasets::YctMath yctMath;

          EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataX;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataY;

          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;

          inputDataX.title = "ACT Math X";
          inputDataX.datasetName = "ACT Math";
          inputDataX.variableName = "X";
          inputDataX.scoreFrequencies = actMathFreq.freqX;
          inputDataX.minimumScore = 0;
          inputDataX.maximumScore = 40;
          inputDataX.scoreIncrement = 1;
          inputDataX.id = "X";

          nlohmann::json univariateStatisticsXJson = univariateStatistics(inputDataX,
                                                                          univariateStatisticsX);

          inputDataY.title = "ACT Math Y";
          inputDataY.datasetName = "ACT Math";
          inputDataY.variableName = "Y";
          inputDataY.scoreFrequencies = actMathFreq.freqY;
          inputDataY.minimumScore = 0;
          inputDataY.maximumScore = 40;
          inputDataY.scoreIncrement = 1;
          inputDataY.id = "Y";

          nlohmann::json univariateStatisticsYJson = univariateStatistics(inputDataY,
                                                                          univariateStatisticsY);

          EquatingRecipes::Analyses::RandomGroupsEquating linearEquatingRandomGroups;
          EquatingRecipes::Analyses::RandomGroupsEquating::InputData inputData;
          EquatingRecipes::Analyses::RandomGroupsEquating::OutputData outputData;

          inputData.title = "ACT Math---Equipercentile";
          inputData.datasetName = actMathFreq.datasetName;
          inputData.design = EquatingRecipes::Structures::Design::RANDOM_GROUPS;
          inputData.method = EquatingRecipes::Structures::Method::EQUIPERCENTILE;
          inputData.smoothing = EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED;
          inputData.univariateStatisticsX = univariateStatisticsX;
          inputData.univariateStatisticsY = univariateStatisticsY;

          EquatingRecipes::Analyses::RandomGroupsEquating::OutputData linearEquatingRandomGroupsOutputData;

          nlohmann::json linearEquatingRandomGroupsJson = linearEquatingRandomGroups(inputData,
                                                                                     linearEquatingRandomGroupsOutputData);

          EquatingRecipes::Analyses::EquatedScaledScores equatedScaledScores;
          EquatingRecipes::Analyses::EquatedScaledScores::InputData inputDataScaledScores;
          EquatingRecipes::Analyses::EquatedScaledScores::OutputData outputDataScaledScores;

          inputDataScaledScores.datasetName = actMathFreq.datasetName;
          inputDataScaledScores.equatedRawScoreResults = linearEquatingRandomGroupsOutputData.equatedRawScoreResults;
          inputDataScaledScores.pData = linearEquatingRandomGroupsOutputData.pData;
          inputDataScaledScores.title = "ACT Math---Equipercentile";
          inputDataScaledScores.lowestObservableEquatedRawScore = 0;
          inputDataScaledScores.highestObservableEquatedRawScore = 40;
          inputDataScaledScores.scoreIncrementEquatedRawScore = 1;
          inputDataScaledScores.lowestObservableScaledScore = 1;
          inputDataScaledScores.highestObservableScaledScore = 36;
          inputDataScaledScores.rawToScaledScoreTable = yctMath.rawToScaledScoreTable;

          nlohmann::json equatedScaledScoresJson = equatedScaledScores(inputDataScaledScores,
                                                                       outputDataScaledScores);
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif