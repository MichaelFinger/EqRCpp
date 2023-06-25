#ifndef ANALYSES_BIVARIATE_STATISTICS_HPP
#define ANALYSES_BIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>
#include <equating_recipes/structures/bivariate_statistics_input_data.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct BivariateStatistics {
      nlohmann::json operator()(const std::string& title,
                                const std::string& datasetName,
                                const EquatingRecipes::Structures::BivariateStatisticsInputData& inputData,
                                EquatingRecipes::Structures::BivariateStatistics& bivariateStatistics) {
        Eigen::MatrixXd scores(inputData.rowScores.size(), 2);
        scores.col(0) = inputData.rowScores;
        scores.col(1) = inputData.columnScores;

        bivariateStatistics = EquatingRecipes::Implementation::Utilities::bivariateFromScores(scores,
                                                                                              inputData.rowMinimumScore,
                                                                                              inputData.rowMaximumScore,
                                                                                              inputData.rowScoreIncrement,
                                                                                              inputData.columnMinimumScore,
                                                                                              inputData.columnMaximumScore,
                                                                                              inputData.columnScoreIncrement,
                                                                                              inputData.rowScoreId,
                                                                                              inputData.columnScoreId,
                                                                                              datasetName,
                                                                                              inputData.rowVariableName,
                                                                                              inputData.columnVariableName);

        nlohmann::json j = {{"analysis_title", title},
                            {"analysis_type", "bivariate_statistics"},
                            {"analysis_results", bivariateStatistics}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif