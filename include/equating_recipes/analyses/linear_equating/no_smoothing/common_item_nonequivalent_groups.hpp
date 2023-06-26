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
#include <equating_recipes/analyses/bivariate_statistics.hpp>

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
            EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataXV;
            EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataYV;
          };

          struct OutputData {
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          };

          nlohmann::json operator()(const InputData& inputData,
                                    OutputData& outputData) {
            EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;
            EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataXV;
            EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataYV;
            EquatingRecipes::Analyses::BivariateStatistics::OutputData outputDataXV;
            EquatingRecipes::Analyses::BivariateStatistics::OutputData outputDataYV;

            nlohmann::json bivariateStatisticsXVResults = bivariateStatistics(inputData.bivariateStatisticsInputDataXV,
                                                                              outputDataXV);

            nlohmann::json bivariateStatisticsYVResults = bivariateStatistics(inputData.bivariateStatisticsInputDataYV,
                                                                              outputDataYV);

            EquatingRecipes::Implementation::CGEquatingNoSmoothing cgEquatingNoSmoothing;

            outputData.bivariateStatisticsXV = outputDataXV.bivariateStatistics;
            outputData.bivariateStatisticsYV = outputDataYV.bivariateStatistics;

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

            nlohmann::json cineGroupsEquatingResults = nlohmann::json {{"analysis_type", "common_item_nonequivalent_groups"},
                                                                       {"analysis_results", results}};

            nlohmann::json j = nlohmann::json::array();

            j.push_back(bivariateStatisticsXVResults);
            j.push_back(bivariateStatisticsYVResults);
            j.push_back(cineGroupsEquatingResults);

            return j;
          }
        };
      } // namespace NoSmoothing
    }   // namespace LinearEquating
  }     // namespace Analyses
} // namespace EquatingRecipes

#endif