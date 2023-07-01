#ifndef ANALYSES_RG_POST_SMOOTHING_HPP
#define ANALYSES_RG_POST_SMOOTHING_HPP

#include <optional>
#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <fmt/core.h>

#include <equating_recipes/analyses/equated_scaled_scores.hpp>
#include <equating_recipes/analyses/mean_linear_equipercentile_equating/no_smoothing/random_groups_equating.hpp>
#include <equating_recipes/implementation/analytic_standard_errors.hpp>
#include <equating_recipes/implementation/cubic_spline.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace MeanLinearEquipercentileEquating {
      namespace PostSmoothing {
        namespace CubicSpline {
          struct RandomGroupsEquating {
            struct InputData {
              std::string title;
              std::string datasetName;
              EquatingRecipes::Structures::Design design;
              EquatingRecipes::Structures::Method method;
              EquatingRecipes::Structures::Smoothing smoothing;
              EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataX;
              EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataY;
              std::optional<EquatingRecipes::Structures::RawToScaledScoreTable> rawToScaledScoreTable;
              double percentileRankLow = 0.5;
              double percentileRankHigh = 99.5;
              double smoothingParameter = 0.2;
              size_t replicationNumber = 0;
            };

            struct OutputData {
              EquatingRecipes::Structures::PData pData;
              EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
              std::optional<EquatingRecipes::Structures::EquatedScaledScoresResults> equatedScaledScoreResults;
              EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
              EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
              EquatingRecipes::Analyses::MeanLinearEquipercentileEquating::NoSmoothing::RandomGroupsEquating::OutputData rgOutputDataXY;
              EquatingRecipes::Analyses::MeanLinearEquipercentileEquating::NoSmoothing::RandomGroupsEquating::OutputData rgOutputDataYX;
              Eigen::VectorXd seXY;
              Eigen::VectorXd seYX;
              EquatingRecipes::Structures::CubicSplinePostsmoothing cubicSplinePostsmoothingXY;
              EquatingRecipes::Structures::CubicSplinePostsmoothing cubicSplinePostsmoothingYX;
            };

            nlohmann::json operator()(const InputData& inputData,
                                      OutputData& outputData) {
              EquatingRecipes::Analyses::MeanLinearEquipercentileEquating::NoSmoothing::RandomGroupsEquating randomGroupsEquatingNoSmoothing;
              EquatingRecipes::Analyses::MeanLinearEquipercentileEquating::NoSmoothing::RandomGroupsEquating::InputData rgInputDataXY;
              EquatingRecipes::Analyses::MeanLinearEquipercentileEquating::NoSmoothing::RandomGroupsEquating::InputData rgInputDataYX;

              rgInputDataXY.datasetName = inputData.datasetName;
              rgInputDataXY.design = inputData.design;
              rgInputDataXY.method = inputData.method;
              rgInputDataXY.smoothing = inputData.smoothing;
              rgInputDataXY.univariateStatisticsInputDataX = inputData.univariateStatisticsInputDataX;
              rgInputDataXY.univariateStatisticsInputDataY = inputData.univariateStatisticsInputDataY;

              nlohmann::json randomGroupsEquatingResultsXY = randomGroupsEquatingNoSmoothing(rgInputDataXY, outputData.rgOutputDataXY);

              rgInputDataYX.datasetName = inputData.datasetName;
              rgInputDataYX.design = inputData.design;
              rgInputDataYX.method = inputData.method;
              rgInputDataYX.smoothing = inputData.smoothing;
              rgInputDataYX.univariateStatisticsInputDataX = inputData.univariateStatisticsInputDataY;
              rgInputDataYX.univariateStatisticsInputDataY = inputData.univariateStatisticsInputDataX;

              nlohmann::json randomGroupsEquatingResultsYX = randomGroupsEquatingNoSmoothing(rgInputDataYX, outputData.rgOutputDataYX);

              EquatingRecipes::Implementation::AnalyticStandardErrors analyticStandardErrors;
              outputData.seXY = analyticStandardErrors.calculate(outputData.rgOutputDataXY.univariateStatisticsY.numberOfScores,
                                                                 outputData.rgOutputDataXY.univariateStatisticsY.numberOfExaminees,
                                                                 outputData.rgOutputDataXY.univariateStatisticsY.cumulativeRelativeFreqDist,
                                                                 outputData.rgOutputDataXY.univariateStatisticsX.numberOfScores,
                                                                 outputData.rgOutputDataXY.univariateStatisticsX.scoreIncrement,
                                                                 outputData.rgOutputDataXY.univariateStatisticsX.numberOfExaminees,
                                                                 outputData.rgOutputDataXY.univariateStatisticsX.percentileRankDist);

              outputData.seYX = analyticStandardErrors.calculate(outputData.rgOutputDataXY.univariateStatisticsX.numberOfScores,
                                                                 outputData.rgOutputDataXY.univariateStatisticsX.numberOfExaminees,
                                                                 outputData.rgOutputDataXY.univariateStatisticsX.cumulativeRelativeFreqDist,
                                                                 outputData.rgOutputDataXY.univariateStatisticsY.numberOfScores,
                                                                 outputData.rgOutputDataXY.univariateStatisticsY.scoreIncrement,
                                                                 outputData.rgOutputDataXY.univariateStatisticsY.numberOfExaminees,
                                                                 outputData.rgOutputDataXY.univariateStatisticsY.percentileRankDist);

              EquatingRecipes::Implementation::CubicSpline cubicSpline;

              cubicSpline.runCubicSpline(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                         outputData.rgOutputDataXY.pData,
                                         outputData.rgOutputDataXY.equatedRawScoreResults,
                                         outputData.seXY,
                                         outputData.cubicSplinePostsmoothingXY,
                                         outputData.rgOutputDataYX.pData,
                                         outputData.rgOutputDataYX.equatedRawScoreResults,
                                         outputData.seYX,
                                         outputData.cubicSplinePostsmoothingYX,
                                         inputData.percentileRankLow,
                                         inputData.percentileRankHigh,
                                         inputData.smoothingParameter,
                                         inputData.replicationNumber,
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

              nlohmann::json randomGroupsCubicSplineEquatingResults = nlohmann::json::object();
              randomGroupsCubicSplineEquatingResults["DatasetName"] = inputData.datasetName;
              randomGroupsCubicSplineEquatingResults["RowVariableName"] = outputData.univariateStatisticsX.variableName;
              randomGroupsCubicSplineEquatingResults["ColumnVariableName"] = outputData.univariateStatisticsY.variableName;
              randomGroupsCubicSplineEquatingResults["PData"] = outputData.pData;
              randomGroupsCubicSplineEquatingResults["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;

              if (inputData.rawToScaledScoreTable.has_value()) {
                randomGroupsCubicSplineEquatingResults["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults.value();
              }

              std::string analysisType = fmt::format("{}_random_groups_cubic_spline",
                                                     EquatingRecipes::Implementation::Utilities::getDesignName(inputData.design));
              nlohmann::json randomGroupsEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                           {"analysis_results", randomGroupsCubicSplineEquatingResults}};

              nlohmann::json j = nlohmann::json::array();
              j.push_back(randomGroupsEquatingResultsXY);
              j.push_back(randomGroupsEquatingResultsYX);
              j.push_back(randomGroupsEquatingResults);

              // if (inputData.rawToScaledScoreTable.has_value()) {
              //   j.push_back(equatedScaledScoresResults);
              // }

              return j;
            }
          };
        } // namespace CubicSpline
      }   // namespace PostSmoothing
    }     // namespace MeanLinearEquipercentileEquating
  }       // namespace Analyses
} // namespace EquatingRecipes

#endif