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
#include <equating_recipes/wrappers/rg_and_sg_equating.hpp>
#include <equating_recipes/wrappers/utilities.hpp>
#include <equating_recipes/analyses/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct LinearEquatingRandomGroups {
      struct InputData {
        std::string datasetName;
        std::string xVariableName;
        std::string yVariableName;

        Eigen::VectorXd scoreFrequenciesX;
        double minimumScoreX;
        double maximumScoreX;
        double scoreIncrementX;
        std::string idX = "X";

        Eigen::VectorXd scoreFrequenciesY;
        double minimumScoreY;
        double maximumScoreY;
        double scoreIncrementY;
        std::string idY = "Y";

        double lowestObservableEquatedRawScore;
        double highestObservableEquatedRawScore;
        double scoreIncrementEquatedRawScore;
        double lowestObservableScaledScore;
        double highestObservableScaledScore;
        EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;

        size_t roundToNumberOfDecimalPlaces = 1;
      };

      struct OutputData {
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
        EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
        EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::LinearEquatingRandomGroups::InputData& inputData,
                                EquatingRecipes::Analyses::LinearEquatingRandomGroups::OutputData& outputData) {
        EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;
        EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataX;

        inputDataX.datasetName = inputData.datasetName;
        inputDataX.id = inputData.idX;
        inputDataX.variableName = inputData.xVariableName;
        inputDataX.minimumScore = inputData.minimumScoreX;
        inputDataX.maximumScore = inputData.maximumScoreX;
        inputDataX.scoreIncrement = inputData.scoreIncrementX;
        inputDataX.scoreFrequencies = inputData.scoreFrequenciesX;

        EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataY;
        inputDataY.datasetName = inputData.datasetName;
        inputDataY.id = inputData.idY;
        inputDataY.variableName = inputData.yVariableName;
        inputDataY.minimumScore = inputData.minimumScoreY;
        inputDataY.maximumScore = inputData.maximumScoreY;
        inputDataY.scoreIncrement = inputData.scoreIncrementY;
        inputDataY.scoreFrequencies = inputData.scoreFrequenciesY;

        nlohmann::json univariateStatisticsXResults = univariateStatistics.operator()(inputDataX, outputData.univariateStatisticsX);
        nlohmann::json univariateStatisticsYResults = univariateStatistics.operator()(inputDataY, outputData.univariateStatisticsY);

        EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
        EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;

        randomAndSingleGroupEquating.randomGroupEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                                         EquatingRecipes::Structures::Method::LINEAR,
                                                         EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                                         outputData.univariateStatisticsX,
                                                         outputData.univariateStatisticsY,
                                                         0,
                                                         pData,
                                                         equatedRawScoreResults);

        EquatingRecipes::Utilities::runEquatedScaledScores(pData,
                                                           equatedRawScoreResults,
                                                           inputData.lowestObservableEquatedRawScore,
                                                           inputData.highestObservableEquatedRawScore,
                                                           inputData.scoreIncrementEquatedRawScore,
                                                           inputData.rawToScaledScoreTable,
                                                           inputData.roundToNumberOfDecimalPlaces,
                                                           inputData.lowestObservableScaledScore,
                                                           inputData.highestObservableScaledScore,
                                                           equatedScaledScoreResults);

        outputData.pData = pData;
        outputData.equatedRawScoreResults = equatedRawScoreResults;
        outputData.equatedScaledScoreResults = equatedScaledScoreResults;

        nlohmann::json j = nlohmann::json::object();

        j["AnalysisTitle"] = "Linear Equating With Random Groups Design";
        j["AnalysisType"] = "linear_equating_random_groups";
        j["DatasetName"] = inputData.datasetName;
        j["RowVariableName"] = inputData.xVariableName;
        j["ColumnwVariableName"] = inputData.yVariableName;
        j["UnivariateStatisticsX"] = outputData.univariateStatisticsX;
        j["UnivariateStatisticsY"] = outputData.univariateStatisticsY;
        j["PData"] = outputData.pData;
        j["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
        j["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults;

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif