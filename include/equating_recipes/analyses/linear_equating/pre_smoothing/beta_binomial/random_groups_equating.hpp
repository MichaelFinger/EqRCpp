#ifndef ANALYSES_LINEAR_EQ_RANDOM_GROUPS_BETA_BIN_SMOOTHING_HPP
#define ANALYSES_LINEAR_EQ_RANDOM_GROUPS_BETA_BIN_SMOOTHING_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/implementation/beta_binomial.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace LinearEquating {
      namespace PreSmoothing {
        namespace BetaBinomial {
          struct RandomGroupsEquating {
            struct InputData {
              std::string title;
              std::string datasetName;
              std::string xVariableName;
              std::string yVariableName;
              EquatingRecipes::Structures::Design design;
              EquatingRecipes::Structures::Method method;
              EquatingRecipes::Structures::Smoothing smoothing;
              EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataX;
              EquatingRecipes::Analyses::UnivariateStatistics::InputData inputDataY;
              size_t numberOfBetaParametersX;
              size_t numberOfBetaParametersY;
              double reliabilityX;
              double reliabilityY;
            };

            struct OutputData {
              EquatingRecipes::Structures::PData pData;
              EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
              EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
              EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomialSmoothingX;
              EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomialSmoothingY;
              EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            };

            nlohmann::json operator()(const InputData& inputData,
                                      OutputData& outputData) {
              EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;

              EquatingRecipes::Analyses::UnivariateStatistics::OutputData outputDataX;
              nlohmann::json univariateStatisticsXResults = univariateStatistics(inputData.inputDataX,
                                                                                 outputDataX);
              outputData.univariateStatisticsX = outputDataX.univariateStatistics;

              EquatingRecipes::Analyses::UnivariateStatistics::OutputData outputDataY;
              nlohmann::json univariateStatisticsYResults = univariateStatistics(inputData.inputDataY,
                                                                                 outputDataY);
              outputData.univariateStatisticsY = outputDataY.univariateStatistics;

              EquatingRecipes::Implementation::BetaBinomial betaBinomial;

              outputData.betaBinomialSmoothingX = betaBinomial.runBetaBinomialSmoothing(outputData.univariateStatisticsX,
                                                                                        inputData.numberOfBetaParametersX,
                                                                                        inputData.reliabilityX);

              outputData.betaBinomialSmoothingY = betaBinomial.runBetaBinomialSmoothing(outputData.univariateStatisticsY,
                                                                                        inputData.numberOfBetaParametersY,
                                                                                        inputData.reliabilityY);

              outputData.equatedRawScoreResults = betaBinomial.runRandomGroupsEquipercentileEquating(inputData.design,
                                                                                                     inputData.method,
                                                                                                     inputData.smoothing,
                                                                                                     outputData.univariateStatisticsX,
                                                                                                     outputData.univariateStatisticsY,
                                                                                                     outputData.betaBinomialSmoothingX,
                                                                                                     outputData.betaBinomialSmoothingY,
                                                                                                     0,
                                                                                                     outputData.pData);

              nlohmann::json results = nlohmann::json::object();
              results["DatasetName"] = inputData.datasetName;
              results["RowVariableName"] = inputData.xVariableName;
              results["ColumnwVariableName"] = inputData.yVariableName;
              results["PData"] = outputData.pData;
              results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
              results["BetaBinomialSmoothingX"] = outputData.betaBinomialSmoothingX;
              results["BetaBinomialSmoothingY"] = outputData.betaBinomialSmoothingY;

              nlohmann::json betaBinomialEquatingResults = nlohmann::json {{"analysis_type", "random_groups_beta_binomial_smoothing"},
                                                                           {"analysis_results", results}};

              nlohmann::json j = nlohmann::json::array();

              j.push_back(univariateStatisticsXResults);
              j.push_back(univariateStatisticsYResults);
              j.push_back(betaBinomialEquatingResults);

              return j;
            }
          };
        } // namespace BetaBinomial
      }   // namespace PreSmoothing
    }     // namespace LinearEquating
  }       // namespace Analyses
} // namespace EquatingRecipes

#endif