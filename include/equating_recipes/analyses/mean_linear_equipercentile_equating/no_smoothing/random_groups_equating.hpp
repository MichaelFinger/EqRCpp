#ifndef ANALYSES_RG_MEAN_LIN_EQUI_EQ_NO_SMOOTHING_HPP
#define ANALYSES_RG_MEAN_LIN_EQUI_EQ_NO_SMOOTHING_HPP

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
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/analyses/equated_scaled_scores.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace MeanLinearEquipercentileEquating {
      namespace NoSmoothing {
        struct RandomGroupsEquating {
          struct InputData {
            std::string datasetName;
            EquatingRecipes::Structures::Design design;
            EquatingRecipes::Structures::Method method;
            EquatingRecipes::Structures::Smoothing smoothing;
            EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataX;
            EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataY;
            std::optional<EquatingRecipes::Structures::RawToScaledScoreTable> rawToScaledScoreTable;
          };

          struct OutputData {
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            std::optional<EquatingRecipes::Structures::EquatedScaledScoresResults> equatedScaledScoreResults;
            EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
            EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
          };

          nlohmann::json operator()(const InputData& inputData,
                                    OutputData& outputData) {
            EquatingRecipes::Analyses::UnivariateStatistics univariateStatisticsAnalysis;

            EquatingRecipes::Analyses::UnivariateStatistics::OutputData univariateStatisticsXOutputData;
            EquatingRecipes::Analyses::UnivariateStatistics::OutputData univariateStatisticsYOutputData;

            nlohmann::json univariateStatisticsXResults = univariateStatisticsAnalysis(inputData.univariateStatisticsInputDataX,
                                                                                       univariateStatisticsXOutputData);

            nlohmann::json univariateStatisticsYResults = univariateStatisticsAnalysis(inputData.univariateStatisticsInputDataY,
                                                                                       univariateStatisticsYOutputData);

            outputData.univariateStatisticsX = univariateStatisticsXOutputData.univariateStatistics;
            outputData.univariateStatisticsY = univariateStatisticsYOutputData.univariateStatistics;

            EquatingRecipes::Implementation::RandomAndSingleGroupEquating randomAndSingleGroupEquating;

            randomAndSingleGroupEquating.runRandomGroupEquating(inputData.design,
                                                                inputData.method,
                                                                inputData.smoothing,
                                                                outputData.univariateStatisticsX,
                                                                outputData.univariateStatisticsY,
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
            results["RowVariableName"] = outputData.univariateStatisticsX.variableName;
            results["ColumnVariableName"] = outputData.univariateStatisticsY.variableName;
            results["PData"] = outputData.pData;
            results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;

            if (inputData.rawToScaledScoreTable.has_value()) {
              results["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults.value();
            }

            std::string analysisType = fmt::format("{}_mean_linear_equipercentile_equating",
                                                   EquatingRecipes::Implementation::Utilities::getDesignName(inputData.design));
            nlohmann::json randomGroupsEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                         {"analysis_results", results}};

            nlohmann::json j = nlohmann::json::array();
            j.push_back(univariateStatisticsXResults);
            j.push_back(univariateStatisticsYResults);
            j.push_back(randomGroupsEquatingResults);

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