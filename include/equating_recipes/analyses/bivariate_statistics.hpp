#ifndef ANALYSES_BIVARIATE_STATISTICS_HPP
#define ANALYSES_BIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>

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

      struct OutputData {
        EquatingRecipes::Structures::BivariateStatistics bivariateStatistics;
      };

      nlohmann::json operator()(const InputData& inputData,
                                OutputData& outputData) {
        Eigen::MatrixXd scores(inputData.rowScores.size(), 2);
        scores.col(0) = inputData.rowScores;
        scores.col(1) = inputData.columnScores;

        outputData.bivariateStatistics = EquatingRecipes::Implementation::Utilities::bivariateFromScores(scores,
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

        nlohmann::json j = {{"analysis_type", "bivariate_statistics"},
                            {"analysis_results", outputData.bivariateStatistics}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif