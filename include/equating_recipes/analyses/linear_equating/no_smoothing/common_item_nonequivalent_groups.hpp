#ifndef ANALYSES_LINEAR_EQ_CINE_NO_SMOOTHING_HPP
#define ANALYSES_LINEAR_EQ_CINE_NO_SMOOTHING_HPP

#include <string>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/cg_equipercentile_equating.hpp>
#include <equating_recipes/implementation/cg_no_smoothing.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace LinearEquating {
      namespace NoSmoothing {
        struct CommonItemNonequivalentGroups {
          struct InputData {
            std::string title;
            std::string datasetName;
            std::string xVariableName;
            std::string yVariableName;
            EquatingRecipes::Structures::Design design;
            EquatingRecipes::Structures::Method method;
            EquatingRecipes::Structures::Smoothing smoothing;
            double populationWeight = -1;
            bool isInternalAnchor = true;
            double reliabilityCommonItemsPopulation1 = 0;
            double reliabilityCommonItemsPopulation2 = 0;
            size_t bootstrapReplicationNumber = 0;
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;
          };

          struct OutputData {
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          };

          nlohmann::json operator()(const EquatingRecipes::Analyses::CommonItemNonequivalentGroups::InputData& inputData,
                                    EquatingRecipes::Analyses::CommonItemNonequivalentGroups::OutputData& outputData) {
            EquatingRecipes::Implementation::CGEquatingNoSmoothing cgEquatingNoSmoothing;

            outputData.bivariateStatisticsXV = inputData.bivariateStatisticsXV;
            outputData.bivariateStatisticsYV = inputData.bivariateStatisticsYV;

            cgEquatingNoSmoothing.runWithNoSmoothing(inputData.design,
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

            nlohmann::json j = nlohmann::json::array();

            j.push_back({{"analysis_title", inputData.title},
                         {"analysis_type", "bivariate_statistics"},
                         {"analysis_results", outputData.bivariateStatisticsXV}});

            j.push_back({{"analysis_title", inputData.title},
                         {"analysis_type", "bivariate_statistics"},
                         {"analysis_results", outputData.bivariateStatisticsYV}});

            j.push_back({{"analysis_title", inputData.title},
                         {"analysis_type", "common_item_nonequivalent_groups"},
                         {"analysis_results", results}});

            return j;
          }
        };
      } // namespace NoSmoothing
    }   // namespace LinearEquating
  }     // namespace Analyses
} // namespace EquatingRecipes

#endif