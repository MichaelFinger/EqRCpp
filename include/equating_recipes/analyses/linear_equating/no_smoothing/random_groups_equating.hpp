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
#include <equating_recipes/implementation/rg_and_sg_equating.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/analyses/equated_scaled_scores.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace LinearEquating {
      namespace NoSmoothing {
        struct RandomGroupsEquating {
          struct InputData {
            std::string title;
            std::string datasetName;
            EquatingRecipes::Structures::Design design;
            EquatingRecipes::Structures::Method method;
            EquatingRecipes::Structures::Smoothing smoothing;
            EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataX;
            EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataY;
            EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;
          };

          struct OutputData {
            EquatingRecipes::Structures::PData pData;
            EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;
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

            EquatingRecipes::Analyses::EquatedScaledScores equatedScaledScores;
            EquatingRecipes::Analyses::EquatedScaledScores::InputData inputDataScaledScores;
            EquatingRecipes::Analyses::EquatedScaledScores::OutputData outputDataScaledScores;
            
            inputDataScaledScores.datasetName = "ACT Math";
            inputDataScaledScores.equatedRawScoreResults = outputData.equatedRawScoreResults;
            inputDataScaledScores.pData = outputData.pData;
            inputDataScaledScores.rawToScaledScoreTable = inputData.rawToScaledScoreTable;

            nlohmann::json equatedScaledScoresResults = equatedScaledScores(inputDataScaledScores,
                                                                            outputDataScaledScores);

            outputData.pData = outputDataScaledScores.pData;
            outputData.equatedScaledScoreResults = outputDataScaledScores.equatedScaledScoreResults;

            nlohmann::json results = nlohmann::json::object();
            results["DatasetName"] = inputData.datasetName;
            results["RowVariableName"] = outputData.univariateStatisticsX.variableName;
            results["ColumnwVariableName"] = outputData.univariateStatisticsY.variableName;
            results["PData"] = outputData.pData;
            results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
            results["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults;

            nlohmann::json randomGroupsEquatingResults = nlohmann::json {{"analysis_type", "random_groups_equating"},
                                                                         {"analysis_results", results}};

            nlohmann::json j = nlohmann::json::array();
            j.push_back(univariateStatisticsXResults);
            j.push_back(univariateStatisticsYResults);
            j.push_back(randomGroupsEquatingResults);
            j.push_back(equatedScaledScoresResults);

            return j;
          }
        };
      } // namespace NoSmoothing
    }   // namespace LinearEquating
  }     // namespace Analyses
} // namespace EquatingRecipes

#endif