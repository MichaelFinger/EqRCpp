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
    struct LinearEquatingRandomGroups {
      struct InputData {
        std::string datasetName;
        std::string xVariableName;
        std::string yVariableName;

        // Eigen::VectorXd scoreFrequenciesX;
        // double minimumScoreX;
        // double maximumScoreX;
        // double scoreIncrementX;
        // std::string idX = "X";

        // Eigen::VectorXd scoreFrequenciesY;
        // double minimumScoreY;
        // double maximumScoreY;
        // double scoreIncrementY;
        // std::string idY = "Y";

        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;

        double lowestObservableEquatedRawScore;
        double highestObservableEquatedRawScore;
        double scoreIncrementEquatedRawScore;
        double lowestObservableScaledScore;
        double highestObservableScaledScore;
        EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;

        size_t roundToNumberOfDecimalPlaces = 1;
      };

      struct OutputData {
        EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
        EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::LinearEquatingRandomGroups::InputData& inputData,
                                EquatingRecipes::Analyses::LinearEquatingRandomGroups::OutputData& outputData) {
        EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
        EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;

        randomAndSingleGroupEquating.randomGroupEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                                         EquatingRecipes::Structures::Method::LINEAR,
                                                         EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                                         inputData.univariateStatisticsX,
                                                         inputData.univariateStatisticsY,
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

        nlohmann::json results = nlohmann::json::object();
        results["DatasetName"] = inputData.datasetName;
        results["RowVariableName"] = inputData.xVariableName;
        results["ColumnwVariableName"] = inputData.yVariableName;
        results["PData"] = outputData.pData;
        results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
        results["EquatedScaledScoreResults"] = outputData.equatedScaledScoreResults;

        nlohmann::json j = {{"linear_equating_random_groups", results}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif