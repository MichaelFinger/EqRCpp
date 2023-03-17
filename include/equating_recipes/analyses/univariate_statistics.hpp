#ifndef ANALYSES_UNIVARIATE_STATISTICS_HPP
#define ANALYSES_UNIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/wrappers/utilities.hpp>

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
        std::string id = "X";
      };

      nlohmann::json operator()(const EquatingRecipes::Analyses::UnivariateStatistics::InputData& inputData,
                                EquatingRecipes::Structures::UnivariateStatistics& univariateStatistics) {
        univariateStatistics = EquatingRecipes::Utilities::univariateFromScoreFrequencies(inputData.scoreFrequencies,
                                                                                          inputData.minimumScore,
                                                                                          inputData.maximumScore,
                                                                                          inputData.scoreIncrement,
                                                                                          inputData.id,
                                                                                          inputData.datasetName,
                                                                                          inputData.variableName);

        nlohmann::json j = nlohmann::json::object();

        j["Analysis"] = "univariate_statistics";
        j["Title"] = "Univariate Statistics";
        j["DatasetName"] = inputData.datasetName;
        j["VariableName"] = inputData.variableName;
        j["Results"] = univariateStatistics;

        return j;
      }
    };
  } // namespace Analyses
} // namespace EquatingRecipes

#endif