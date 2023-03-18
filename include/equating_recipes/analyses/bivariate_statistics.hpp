#ifndef ANALYSES_BIVARIATE_STATISTICS_HPP
#define ANALYSES_BIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/wrappers/utilities.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct BivariateStatistics {
      struct InputData {
        std::string datasetName;
        std::string rowVariableName;
        std::string columnVariableName;

        Eigen::VectorXd rowScores;
        double rowMinimumScore;
        double rowMaximumScore;
        double rowScoreIncrement;

        Eigen::VectorXd columnScores;
        double columnMinimumScore;
        double columnMaximumScore;
        double columnScoreIncrement;

        std::string rowScoreId = "X";
        std::string columnScoreId = "Y";
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::BivariateStatistics::InputData& inputData,
                                EquatingRecipes::Structures::BivariateStatistics& bivariateStatistics) {
        Eigen::MatrixXd scores(inputData.rowScores.size(), 2);
        scores.col(0) = inputData.rowScores;
        scores.col(1) = inputData.columnScores;

        bivariateStatistics = EquatingRecipes::Utilities::bivariateFromScores(scores,
                                                                              inputData.rowMinimumScore,
                                                                              inputData.rowMaximumScore,
                                                                              inputData.rowScoreIncrement,
                                                                              inputData.columnMinimumScore,
                                                                              inputData.columnMaximumScore,
                                                                              inputData.columnScoreIncrement,
                                                                              inputData.rowScoreId,
                                                                              inputData.columnScoreId,
                                                                              inputData.datasetName,
                                                                              inputData.rowVariableName,
                                                                              inputData.columnVariableName);

        nlohmann::json j = {{"bivariate_statistics", bivariateStatistics}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif