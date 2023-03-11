#ifndef TESTS_EXAMPLES_CHAPTER_3_HPP
#define TESTS_EXAMPLES_CHAPTER_3_HPP

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/wrappers/utilities.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>
#include <equating_recipes/wrappers/rg_and_sg_equating.hpp>

#include "datasets/actmathfreq.hpp"
#include "datasets/yctmath.hpp"

#include <equating_recipes/analyses/linear_equating_random_groups.hpp>

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter3 {
        void operator()() {
          EquatingRecipes::Tests::Examples::Datasets::ACTMathFreq actMathFreq;
          EquatingRecipes::Tests::Examples::Datasets::YctMath yctMath;

          // EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX = EquatingRecipes::Utilities::univariateFromScoreFrequencies(actMathFreq.freqX,
          //                                                                                                                                      0,
          //                                                                                                                                      40,
          //                                                                                                                                      1,
          //                                                                                                                                      "X");

          // EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY = EquatingRecipes::Utilities::univariateFromScoreFrequencies(actMathFreq.freqY,
          //                                                                                                                                      0,
          //                                                                                                                                      40,
          //                                                                                                                                      1,
          //                                                                                                                                      "Y");

          // EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;
          // EquatingRecipes::Structures::PData pData;
          // EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          // EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;

          // randomAndSingleGroupEquating.randomGroupEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
          //                                                  EquatingRecipes::Structures::Method::LINEAR,
          //                                                  EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
          //                                                  univariateStatisticsX,
          //                                                  univariateStatisticsY,
          //                                                  0,
          //                                                  pData,
          //                                                  equatedRawScoreResults);

          // EquatingRecipes::Utilities::runEquatedScaledScores(pData,
          //                                                    equatedRawScoreResults,
          //                                                    0,
          //                                                    40,
          //                                                    1,
          //                                                    yctMath.rawToScaledScoreTable,
          //                                                    1,
          //                                                    1,
          //                                                    36,
          //                                                    equatedScaledScoreResults);

          // EquatingRecipes::JSON::JsonDocument jsonDoc;

          // nlohmann::json j = nlohmann::json::object();
          // j["EquatingAnalysisType"] = "Linear Equating With Random Groups Design";
          // j["UnivariateStatisticsX"] = univariateStatisticsX;
          // j["UnivariateStatisticsY"] = univariateStatisticsY;
          // j["PData"] = pData;
          // j["EquatedRawScoreResults"] = equatedRawScoreResults;
          // j["EquatedScaledScoreResults"] = equatedScaledScoreResults;

          // jsonDoc.setJson(j);
          // jsonDoc.toTextFile("chapter3.json");

          EquatingRecipes::Analyses::LinearEquatingRandomGroups::InputData inputData;

          inputData.scoreFrequenciesX = actMathFreq.freqX;
          inputData.minimumScoreX = 0;
          inputData.maximumScoreX = 40;
          inputData.scoreIncrementX = 1;
          inputData.scoreFrequenciesY = actMathFreq.freqY;
          inputData.minimumScoreY = 0;
          inputData.maximumScoreY = 40;
          inputData.scoreIncrementY = 1;
          inputData.lowestObservableEquatedRawScore = 0;
          inputData.highestObservableEquatedRawScore = 40;
          inputData.scoreIncrementEquatedRawScore = 1;
          inputData.lowestObservableScaledScore = 1;
          inputData.highestObservableScaledScore = 36;
          inputData.rawToScaledScoreTable = yctMath.rawToScaledScoreTable;

          EquatingRecipes::Analyses::LinearEquatingRandomGroups linearEquatingRandomGroups;

          nlohmann::json j = linearEquatingRandomGroups(inputData);

          EquatingRecipes::JSON::JsonDocument jsonDoc;
          jsonDoc.setJson(j);
          jsonDoc.toTextFile("chapter3.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif