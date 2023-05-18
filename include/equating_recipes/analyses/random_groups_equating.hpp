#ifndef ANALYSES_LINEAR_EQUATING_RANDOM_GROUPS_HPP
#define ANALYSES_LINEAR_EQUATING_RANDOM_GROUPS_HPP

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
#include <equating_recipes/rg_and_sg_equating.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct RandomGroupsEquating {
      struct InputData {
        std::string title;
        std::string datasetName;
        EquatingRecipes::Structures::Design design;
        EquatingRecipes::Structures::Method method;
        EquatingRecipes::Structures::Smoothing smoothing;
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
        size_t roundToNumberOfDecimalPlaces = 1;
      };

      struct OutputData {
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::RandomGroupsEquating::InputData& inputData,
                                EquatingRecipes::Analyses::RandomGroupsEquating::OutputData& outputData) {
        EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;

        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;

        randomAndSingleGroupEquating.randomGroupEquating(inputData.design,
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
  } // namespace Analyses
} // namespace EquatingRecipes

#endif