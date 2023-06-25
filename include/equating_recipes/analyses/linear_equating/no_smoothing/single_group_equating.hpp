#ifndef ANALYSES_LINEAR_EQ_SINGLE_GROUP_NO_SMOOTHING_HPP
#define ANALYSES_LINEAR_EQ_SINGLE_GROUP_NO_SMOOTHING_HPP

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
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/implementation/rg_and_sg_equating.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace LinearEquating {
      namespace NoSmoothing {
        struct SingleGroupEquating {
          struct InputData {
            std::string title;
            std::string datasetName;
            EquatingRecipes::Structures::Design design;
            EquatingRecipes::Structures::Method method;
            EquatingRecipes::Structures::Smoothing smoothing;
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXY;
            size_t roundToNumberOfDecimalPlaces = 1;
          };

          struct OutputData {
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          };

          nlohmann::json operator()(const EquatingRecipes::Analyses::SingleGroupEquating::InputData& inputData,
                                    EquatingRecipes::Analyses::SingleGroupEquating::OutputData& outputData) {
            EquatingRecipes::Implementation::RandomAndSingleGroupEquating randomAndSingleGroupEquating;

            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;

            randomAndSingleGroupEquating.runSingleGroupEquating(inputData.design,
                                                                inputData.method,
                                                                inputData.smoothing,
                                                                inputData.bivariateStatisticsXY,
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
                                {"analysis_type", "single_group_equating"},
                                {"analysis_results", results}};

            return j;
          }
        };
      } // namespace NoSmoothing
    }   // namespace LinearEquating
  }     // namespace Analyses
} // namespace EquatingRecipes

#endif