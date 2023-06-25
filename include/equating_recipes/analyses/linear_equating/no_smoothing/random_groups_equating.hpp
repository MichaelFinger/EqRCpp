#ifndef ANALYSES_LINEAR_EQ_RANDOM_GROUPS_NO_SMOOTHING_HPP
#define ANALYSES_LINEAR_EQ_RANDOM_GROUPS_NO_SMOOTHING_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/univariate_statistics_input_data.hpp>
#include <equating_recipes/implementation/rg_and_sg_equating.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace LinearEquating {
      namespace NoSmoothing {

        struct RandomGroupsEquating {
          struct InputData {
            EquatingRecipes::Structures::Design design;
            EquatingRecipes::Structures::Method method;
            EquatingRecipes::Structures::Smoothing smoothing;
            EquatingRecipes::Structures::UnivariateStatisticsInputData univariateStatisticsInputDataX;
            EquatingRecipes::Structures::UnivariateStatisticsInputData univariateStatisticsInputDataY;
          };

          struct OutputData {
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
            EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
          };

          nlohmann::json operator()(const std::string& title,
                                    const std::string& datasetName,
                                    const EquatingRecipes::Analyses::RandomGroupsEquating::InputData& inputData,
                                    EquatingRecipes::Analyses::RandomGroupsEquating::OutputData& outputData) {
            EquatingRecipes::Analyses::UnivariateStatistics univariateStatisticsAnalysis;
            univariateStatisticsAnalysis(title,
                                         datasetName,
                                         inputData.univariateStatisticsInputDataX,
                                         outputData.univariateStatisticsX);

            EquatingRecipes::Implementation::RandomAndSingleGroupEquating randomAndSingleGroupEquating;

            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;

            randomAndSingleGroupEquating.runRandomGroupEquating(inputData.design,
                                                                inputData.method,
                                                                inputData.smoothing,
                                                                inputData.univariateStatisticsX,
                                                                inputData.univariateStatisticsY,
                                                                0,
                                                                outputData.pData,
                                                                outputData.equatedRawScoreResults);

            nlohmann::json results = nlohmann::json::object();
            results["DatasetName"] = inputData.datasetName;
            results["RowVariableName"] = inputData.univariateStatisticsX.variableName;
            results["ColumnwVariableName"] = inputData.univariateStatisticsY.variableName;
            results["PData"] = outputData.pData;
            results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;

            nlohmann::json j = {{"analysis_title", inputData.title},
                                {"analysis_type", "random_groups_equating"},
                                {"analysis_results", results}};

            return j;
          }
        };
      } // namespace NoSmoothing
    }   // namespace LinearEquating
  }     // namespace Analyses
} // namespace EquatingRecipes

#endif