#ifndef ANALYSES_UNIVARIATE_STATISTICS_HPP
#define ANALYSES_UNIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/implementation/utilities.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    struct UnivariateStatistics {
      struct InputData {
        std::string datasetName;
        std::string variableName;
        Eigen::VectorXd scoreFrequencies;
        double minimumScore;
        double maximumScore;
        double scoreIncrement;
        std::string id;
      };

      struct OutputData {
        EquatingRecipes::Structures::UnivariateStatistics univariateStatistics;
      };

      nlohmann::json operator()(const InputData& inputData,
                                OutputData& outputData) {
        outputData.univariateStatistics = EquatingRecipes::Implementation::Utilities::univariateFromScoreFrequencies(inputData.scoreFrequencies,
                                                                                                          inputData.minimumScore,
                                                                                                          inputData.maximumScore,
                                                                                                          inputData.scoreIncrement,
                                                                                                          inputData.id,
                                                                                                          inputData.datasetName,
                                                                                                          inputData.variableName);

        nlohmann::json j = {{"analysis_type", "univariate_statistics"},
                            {"analysis_results", outputData.univariateStatistics}};

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif