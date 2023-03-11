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

namespace EquatingRecipes {
  namespace Analyses {
    struct LinearEquatingRandomGroups {
      struct InputData {
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

      nlohmann::json operator()(const EquatingRecipes::Analyses::LinearEquatingRandomGroups::InputData& inputData) {
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX = EquatingRecipes::Utilities::univariateFromScoreFrequencies(inputData.scoreFrequenciesX,
                                                                                                                                             inputData.minimumScoreX,
                                                                                                                                             inputData.maximumScoreX,
                                                                                                                                             inputData.scoreIncrementX,
                                                                                                                                             inputData.idX);

        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY = EquatingRecipes::Utilities::univariateFromScoreFrequencies(inputData.scoreFrequenciesY,
                                                                                                                                             inputData.minimumScoreY,
                                                                                                                                             inputData.maximumScoreY,
                                                                                                                                             inputData.scoreIncrementY,
                                                                                                                                             inputData.idY);

        EquatingRecipes::RandomAndSingleGroupEquating randomAndSingleGroupEquating;
        EquatingRecipes::Structures::PData pData;
        EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
        EquatingRecipes::Structures::EquatedScaledScoresResults equatedScaledScoreResults;

        randomAndSingleGroupEquating.randomGroupEquating(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                                         EquatingRecipes::Structures::Method::LINEAR,
                                                         EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED,
                                                         univariateStatisticsX,
                                                         univariateStatisticsY,
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

        nlohmann::json j = nlohmann::json::object();
        nlohmann::json equatingResults = nlohmann::json::object();

        equatingResults["Title"] = "Linear Equating With Random Groups Design";
        equatingResults["UnivariateStatisticsX"] = univariateStatisticsX;
        equatingResults["UnivariateStatisticsY"] = univariateStatisticsY;
        equatingResults["PData"] = pData;
        equatingResults["EquatedRawScoreResults"] = equatedRawScoreResults;
        equatingResults["EquatedScaledScoreResults"] = equatedScaledScoreResults;

        j["linear_equating_random_groups"] = equatingResults;

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif