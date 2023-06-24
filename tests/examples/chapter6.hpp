#ifndef TESTS_EXAMPLES_CHAPTER_6_HPP
#define TESTS_EXAMPLES_CHAPTER_6_HPP

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/json/json_document.hpp>
#include <equating_recipes/analyses/bivariate_statistics.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/analyses/common_item_nonequivalent_groups.hpp>

#include "datasets/mondatx.hpp"
#include "datasets/mondaty.hpp"

namespace EquatingRecipes {
  namespace Tests {
    namespace Examples {
      struct Chapter6 {
        void operator()() {
          /* Common-item Nonequivalent Groups Design: 
            Kolen and Brennan (2004) Chapter 5 example:
            Equipercentile equating (see p. 151) */

          EquatingRecipes::Tests::Examples::Datasets::MondatX mondatX;
          EquatingRecipes::Tests::Examples::Datasets::MondatY mondatY;

          EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;
          EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataXV;
          EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;

          inputDataXV.title = "Mondat X";
          inputDataXV.datasetName = mondatX.datasetName;
          inputDataXV.rowVariableName = "X";
          inputDataXV.columnVariableName = "V";
          inputDataXV.rowScores = mondatX.rawScores.col(0);
          inputDataXV.rowMinimumScore = 0;
          inputDataXV.rowMaximumScore = 36;
          inputDataXV.rowScoreIncrement = 1;
          inputDataXV.columnScores = mondatX.rawScores.col(1);
          inputDataXV.columnMinimumScore = 0;
          inputDataXV.columnMaximumScore = 12;
          inputDataXV.columnScoreIncrement = 1;
          inputDataXV.rowScoreId = "X";
          inputDataXV.columnScoreId = "V";

          nlohmann::json bivariateStatisticsXVJson = bivariateStatistics(inputDataXV,
                                                                         bivariateStatisticsXV);

          EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataYV;
          EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;

          inputDataYV.title = "Mondat Y";
          inputDataYV.datasetName = mondatY.datasetName;
          inputDataYV.rowVariableName = "Y";
          inputDataYV.columnVariableName = "V";
          inputDataYV.rowScores = mondatY.rawScores.col(0);
          inputDataYV.rowMinimumScore = 0;
          inputDataYV.rowMaximumScore = 36;
          inputDataYV.rowScoreIncrement = 1;
          inputDataYV.columnScores = mondatY.rawScores.col(1);
          inputDataYV.columnMinimumScore = 0;
          inputDataYV.columnMaximumScore = 12;
          inputDataYV.columnScoreIncrement = 1;
          inputDataYV.rowScoreId = "Y";
          inputDataYV.columnScoreId = "V";

          nlohmann::json bivariateStatisticsYVJson = bivariateStatistics(inputDataYV,
                                                                         bivariateStatisticsYV);

          EquatingRecipes::Analyses::CommonItemNonequivalentGroups commonItemNonequivalentGroups;
          EquatingRecipes::Analyses::CommonItemNonequivalentGroups::InputData inputData;
          EquatingRecipes::Analyses::CommonItemNonequivalentGroups::OutputData outputData;

          inputData.bootstrapReplicationNumber = 0;
          inputData.datasetName = "Mondat X, Y";
          inputData.design = EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS;
          inputData.isInternalAnchor = true;
          inputData.method = EquatingRecipes::Structures::Method::FE_BH_MFE_BH_CHAINED;
          inputData.populationWeight = 1;
          inputData.reliabilityCommonItemsPopulation1 = 0.5584431;
          inputData.reliabilityCommonItemsPopulation2 = 0.5735077;
          inputData.smoothing = EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED;
          inputData.title = "Chapter 5 in K&B: w1 = 1";
          inputData.xVariableName = inputDataXV.rowVariableName;
          inputData.yVariableName = inputDataYV.rowVariableName;
          inputData.bivariateStatisticsXV = bivariateStatisticsXV;
          inputData.bivariateStatisticsYV = bivariateStatisticsYV;

          nlohmann::json commonItemNonequivalentGroupsJson = commonItemNonequivalentGroups(inputData,
                                                                                           outputData);

          nlohmann::json j = commonItemNonequivalentGroupsJson;

          EquatingRecipes::JSON::JsonDocument jsonDoc;
          jsonDoc.setJson(j);
          jsonDoc.toTextFile("chapter6.json");
        }
      };
    } // namespace Examples
  }   // namespace Tests
} // namespace EquatingRecipes

#endif