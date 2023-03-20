#ifndef COMMON_ITEM_NONEQUIV_GROUPS_HPP
#define COMMON_ITEM_NONEQUIV_GROUPS_HPP

#include <string>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/cg_no_smoothing.hpp>
#include <equating_recipes/utilities.hpp>

#include <equating_recipes/analyses/bivariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct CommonItemNonequivalentGroups {
      struct InputData {
        std::string datasetName;
        std::string xVariableName;
        std::string yVariableName;
        EquatingRecipes::Structures::Design design = EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS;
        EquatingRecipes::Structures::Method method = EquatingRecipes::Structures::Method::LINEAR;
        EquatingRecipes::Structures::Smoothing smoothing = EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED;
        double populationWeight = -1;
        bool isInternalAnchor = true;
        double reliabilityCommonItemsPopulation1 = 0;
        double reliabilityCommonItemsPopulation2 = 0;
        size_t bootstrapReplicationNumber = 0;
      };

      struct OutputData {
        EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
        EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::CommonItemNonequivalentGroups::InputData& inputData,
                                EquatingRecipes::Analyses::CommonItemNonequivalentGroups::OutputData& outputData) {
        EquatingRecipes::CGEquatingNoSmoothing cgEquatingNoSmoothing;

        cgEquatingNoSmoothing.run(inputData.design,
                                  inputData.method,
                                  inputData.smoothing,
                                  inputData.populationWeight,
                                  inputData.isInternalAnchor,
                                  inputData.reliabilityCommonItemsPopulation1,
                                  inputData.reliabilityCommonItemsPopulation2,
                                  outputData.bivariateStatisticsXV,
                                  outputData.bivariateStatisticsYV,
                                  inputData.bootstrapReplicationNumber,
                                  outputData.pData,
                                  outputData.equatedRawScoreResults);

        nlohmann::json results = nlohmann::json::object();
        results["DatasetName"] = inputData.datasetName;
        results["RowVariableName"] = inputData.xVariableName;
        results["ColumnwVariableName"] = inputData.yVariableName;
        results["PData"] = outputData.pData;
        results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;

        nlohmann::json j = {{"common_item_nonequivalent_groups", results}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif