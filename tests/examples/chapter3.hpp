#ifndef TESTS_EXAMPLES_CHAPTER_3_HPP
#define TESTS_EXAMPLES_CHAPTER_3_HPP

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>

#include "datasets/actmathfreq.hpp"
#include "datasets/yctmath.hpp"

#include <equating_recipes/analyses/random_groups_equating.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/analyses/equated_scaled_scores.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter3 {
        void operator()() {
          EquatingRecipes::Tests::Examples::Datasets::ACTMathFreq actMathFreq;
          EquatingRecipes::Tests::Examples::Datasets::YctMath yctMath;

          EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;

          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataX;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataY;

          inputDataX.datasetName = "ACT Math Freq";
          inputDataX.variableName = "X";
          inputDataX.scoreFrequencies = actMathFreq.freqX;
          inputDataX.minimumScore = 0;
          inputDataX.maximumScore = 40;
          inputDataX.scoreIncrement = 1;
          inputDataX.id = "X";

          inputDataY.datasetName = "ACT Math Freq";
          inputDataY.variableName = "Y";
          inputDataY.scoreFrequencies = actMathFreq.freqY;
          inputDataY.minimumScore = 0;
          inputDataY.maximumScore = 40;
          inputDataY.scoreIncrement = 1;
          inputDataY.id = "Y";

          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;

          nlohmann::json univariateStatisticsXJson = univariateStatistics(inputDataX, univariateStatisticsX);
          nlohmann::json univariateStatisticsYJson = univariateStatistics(inputDataY, univariateStatisticsY);

          EquatingRecipes::Analyses::RandomGroupsEquating::InputData inputData;
          inputData.title = actMathFreq.datasetName + "---Linear";
          inputData.datasetName = actMathFreq.datasetName;
          inputData.design = EquatingRecipes::Structures::Design::RANDOM_GROUPS;
          inputData.method = EquatingRecipes::Structures::Method::LINEAR;
          inputData.smoothing = EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED;

          inputData.univariateStatisticsX = univariateStatisticsX;
          inputData.univariateStatisticsY = univariateStatisticsY;
          
          EquatingRecipes::Analyses::RandomGroupsEquating linearEquatingRandomGroups;

          EquatingRecipes::Analyses::RandomGroupsEquating::OutputData linearEquatingRandomGroupsOutputData;

          nlohmann::json linearEquatingRandomGroupsJson = linearEquatingRandomGroups(inputData,
                                                                                     linearEquatingRandomGroupsOutputData);

          EquatingRecipes::Analyses::EquatedScaledScores equatedScaledScores;
          EquatingRecipes::Analyses::EquatedScaledScores::InputData inputDataScaledScores;
          inputDataScaledScores.datasetName = "ACT Math";
          inputDataScaledScores.equatedRawScoreResults = linearEquatingRandomGroupsOutputData.equatedRawScoreResults;
          inputDataScaledScores.pData = linearEquatingRandomGroupsOutputData.pData;
          inputDataScaledScores.lowestObservableEquatedRawScore = 0;
          inputDataScaledScores.highestObservableEquatedRawScore = 40;
          inputDataScaledScores.scoreIncrementEquatedRawScore = 1;
          inputDataScaledScores.lowestObservableScaledScore = 1;
          inputDataScaledScores.highestObservableScaledScore = 36;
          inputDataScaledScores.rawToScaledScoreTable = yctMath.rawToScaledScoreTable;

          nlohmann::json j = {univariateStatisticsXJson,
                              univariateStatisticsYJson,
                              linearEquatingRandomGroupsJson};

          EquatingRecipes::JSON::JsonDocument jsonDoc;
          jsonDoc.setJson(j);
          jsonDoc.toTextFile("chapter3.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif