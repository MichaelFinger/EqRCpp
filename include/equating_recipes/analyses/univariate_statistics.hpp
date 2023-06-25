#ifndef ANALYSES_UNIVARIATE_STATISTICS_HPP
#define ANALYSES_UNIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>

#include <equating_recipes/structures/univariate_statistics_input_data.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct UnivariateStatistics {
      nlohmann::json operator()(const std::string& title,
                                const std::string& datasetName,
                                const EquatingRecipes::Structures::UnivariateStatisticsInputData& inputData,
                                EquatingRecipes::Structures::UnivariateStatistics& univariateStatistics) {
        univariateStatistics = EquatingRecipes::Implementation::Utilities::univariateFromScoreFrequencies(inputData.scoreFrequencies,
                                                                                                          inputData.minimumScore,
                                                                                                          inputData.maximumScore,
                                                                                                          inputData.scoreIncrement,
                                                                                                          inputData.id,
                                                                                                          datasetName,
                                                                                                          inputData.variableName);

        nlohmann::json j = {{"analysis_title", title},
                            {"analysis_type", "univariate_statistics"},
                            {"analysis_results", univariateStatistics}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif