#ifndef ANALYSES_ANALYTIC_STANDARD_ERRORS_HPP
#define ANALYSES_ANALYTIC_STANDARD_ERRORS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>
#include <equating_recipes/implementation/analytic_standard_errors.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct AnalyticStandardErrors {
      struct InputData {
        std::string title;
        std::string datasetName;
        std::string idX = "X";
        std::string idY = "Y";
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
        EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::AnalyticStandardErrors::InputData& inputData,
                                Eigen::VectorXd& standardErrors) {
        EquatingRecipes::Implementation::AnalyticStandardErrors analyticStandardErrors;

        standardErrors = analyticStandardErrors.calculate(inputData.univariateStatisticsY.numberOfScores,
                                                          inputData.univariateStatisticsY.numberOfExaminees,
                                                          inputData.univariateStatisticsY.cumulativeRelativeFreqDist,
                                                          inputData.univariateStatisticsX.numberOfScores,
                                                          inputData.univariateStatisticsX.scoreIncrement,
                                                          inputData.univariateStatisticsX.numberOfExaminees,
                                                          inputData.univariateStatisticsX.percentileRankDist);

        nlohmann::json j = {{"analysis_title", inputData.title},
                            {"analysis_type", "analytic_standard_errors"},
                            {"analysis_results", {{"standard_errors", standardErrors}}}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif