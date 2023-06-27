#ifndef ANALYSES_LINEAR_EQ_RANDOM_GROUPS_LOG_LIN_SMOOTHING_HPP
#define ANALYSES_LINEAR_EQ_RANDOM_GROUPS_LOG_LIN_SMOOTHING_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/implementation/beta_binomial.hpp>
#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/criterion_comparison_type.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace LinearEquating {
      namespace PreSmoothing {
        namespace LogLinear {
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
              size_t numberOfDegreesSmoothng;
              bool useStandardizedScale;
              EquatingRecipes::Structures::DesignMatrixType designMatrixType;
              EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType;
              double criterion;
            };

            struct OutputData {
              EquatingRecipes::Structures::PData pData;
              EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
              EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
              EquatingRecipes::Structures::UnivariateLogLinearSmoothing univariateLogLinearSmoothingX;
              EquatingRecipes::Structures::UnivariateLogLinearSmoothing univariateLogLinearSmoothingY;
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

              EquatingRecipes::Implementation::LogLinearEquating logLinearEquating;

              outputData.univariateLogLinearSmoothingX = logLinearEquating.runUnivariateLogLinearSmoothing(outputData.univariateStatisticsX,
                                                                                                           inputData.numberOfDegreesSmoothng,
                                                                                                           inputData.useStandardizedScale,
                                                                                                           inputData.designMatrixType,
                                                                                                           inputData.criterionComparisonType,
                                                                                                           inputData.criterion);

              outputData.univariateLogLinearSmoothingY = logLinearEquating.runUnivariateLogLinearSmoothing(outputData.univariateStatisticsY,
                                                                                                           inputData.numberOfDegreesSmoothng,
                                                                                                           inputData.useStandardizedScale,
                                                                                                           inputData.designMatrixType,
                                                                                                           inputData.criterionComparisonType,
                                                                                                           inputData.criterion);
              
              outputData.equatedRawScoreResults = logLinearEquating.runRGEquiEquatingWithLoglinearSmoothing(inputData.design,
                                                                                                            inputData.method,
                                                                                                            inputData.smoothing,
                                                                                                            outputData.univariateStatisticsX,
                                                                                                            outputData.univariateStatisticsY,
                                                                                                            outputData.univariateLogLinearSmoothingX,
                                                                                                            outputData.univariateLogLinearSmoothingY,
                                                                                                            0,
                                                                                                            outputData.pData);

              nlohmann::json results = nlohmann::json::object();
              results["DatasetName"] = inputData.datasetName;
              results["RowVariableName"] = inputData.xVariableName;
              results["ColumnwVariableName"] = inputData.yVariableName;
              results["PData"] = outputData.pData;
              results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
              results["LogLinearSmoothingX"] = outputData.univariateLogLinearSmoothingX;
              results["LogLinearSmoothingY"] = outputData.univariateLogLinearSmoothingY;

              nlohmann::json loglinearEquatingResults = nlohmann::json {{"analysis_type", "random_groups_log_linear_smoothing"},
                                                                           {"analysis_results", results}};

              nlohmann::json j = nlohmann::json::array();

              j.push_back(univariateStatisticsXResults);
              j.push_back(univariateStatisticsYResults);
              j.push_back(loglinearEquatingResults);

              return j;
            }
          };
        } // namespace LogLinear
      }   // namespace PreSmoothing
    }     // namespace LinearEquating
  }       // namespace Analyses
} // namespace EquatingRecipes

#endif