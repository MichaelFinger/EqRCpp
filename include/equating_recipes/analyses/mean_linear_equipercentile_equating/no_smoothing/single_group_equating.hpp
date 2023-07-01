#ifndef ANALYSES_SG_MEAN_LIN_EQUI_EQ_NO_SMOOTHING_HPP
#define ANALYSES_SG_MEAN_LIN_EQUI_EQ_NO_SMOOTHING_HPP

#include <optional>
#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <fmt/core.h>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/rg_and_sg_equating.hpp>
#include <equating_recipes/analyses/bivariate_statistics.hpp>
#include <equating_recipes/analyses/equated_scaled_scores.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace MeanLinearEquipercentileEquating {
      namespace NoSmoothing {
        struct SingleGroupEquating {
          struct InputData {
            std::string title;
            std::string datasetName;
            EquatingRecipes::Structures::Design design;
            EquatingRecipes::Structures::Method method;
            EquatingRecipes::Structures::Smoothing smoothing;
            EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataXY;
            std::optional<EquatingRecipes::Structures::RawToScaledScoreTable> rawToScaledScoreTable;
          };

          struct OutputData {
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            std::optional<EquatingRecipes::Structures::EquatedScaledScoresResults> equatedScaledScoreResults;
            EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXY;
          };

          nlohmann::json operator()(const InputData& inputData,
                                    OutputData& outputData) {
            EquatingRecipes::Analyses::BivariateStatistics bivariateStatisticsAnalysis;

            EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsXYOutputData;

            nlohmann::json bivariateStatisticsXYResults = bivariateStatisticsAnalysis(inputData.bivariateStatisticsInputDataXY,
                                                                                      bivariateStatisticsXYOutputData);

            outputData.bivariateStatisticsXY = bivariateStatisticsXYOutputData.bivariateStatistics;

            EquatingRecipes::Implementation::RandomAndSingleGroupEquating randomAndSingleGroupEquating;

            randomAndSingleGroupEquating.runSingleGroupEquating(inputData.design,
                                                                inputData.method,
                                                                inputData.smoothing,
                                                                outputData.bivariateStatisticsXY,
                                                                0,
                                                                outputData.pData,
                                                                outputData.equatedRawScoreResults);

            nlohmann::json equatedScaledScoresResults;

            if (inputData.rawToScaledScoreTable.has_value()) {
              EquatingRecipes::Analyses::EquatedScaledScores equatedScaledScores;
              EquatingRecipes::Analyses::EquatedScaledScores::InputData inputDataScaledScores;
              EquatingRecipes::Analyses::EquatedScaledScores::OutputData outputDataScaledScores;

              inputDataScaledScores.datasetName = "ACT Math";
              inputDataScaledScores.equatedRawScoreResults = outputData.equatedRawScoreResults;
              inputDataScaledScores.pData = outputData.pData;
              inputDataScaledScores.rawToScaledScoreTable = inputData.rawToScaledScoreTable.value();

              equatedScaledScoresResults = equatedScaledScores(inputDataScaledScores,
                                                               outputDataScaledScores);

              outputData.equatedScaledScoreResults = outputDataScaledScores.equatedScaledScoreResults;

              outputData.pData = outputDataScaledScores.pData;
            }            

            nlohmann::json results = nlohmann::json::object();
            results["DatasetName"] = inputData.datasetName;
            results["RowVariableName"] = outputData.bivariateStatisticsXY.univariateStatisticsRow.id;
            results["ColumnVariableName"] = outputData.bivariateStatisticsXY.univariateStatisticsRow.id;
            results["PData"] = outputData.pData;
            results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;

            if (inputData.rawToScaledScoreTable.has_value()) {
              results["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults.value();
            }

            std::string analysisType = fmt::format("{}_mean_linear_equipercentile_equating",
                                                   EquatingRecipes::Implementation::Utilities::getDesignName(inputData.design));

            nlohmann::json singleGroupEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                        {"analysis_results", results}};

            nlohmann::json j = nlohmann::json::array();
            j.push_back(bivariateStatisticsXYResults);
            j.push_back(singleGroupEquatingResults);

            if (inputData.rawToScaledScoreTable.has_value()) {
              j.push_back(equatedScaledScoresResults);
            }

            return j;
          }
        };
      } // namespace NoSmoothing
    }   // namespace MeanLinearEquipercentileEquating
  }     // namespace Analyses
} // namespace EquatingRecipes

#endif