#ifndef ANALYSES_CINE_MEAN_LIN_EQUI_EQ_LOG_LIN_SMOOTHING_HPP
#define ANALYSES_CINE_MEAN_LIN_EQUI_EQ_LOG_LIN_SMOOTHING_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace MeanLinearEquipercentileEquating {
      namespace PreSmoothing {
        namespace LogLinear {
          struct CINEGroupsEquating {
            struct InputData {
              std::string datasetName;
              std::string xVariableName;
              std::string yVariableName;
              EquatingRecipes::Structures::Design design;
              EquatingRecipes::Structures::Method method;
              EquatingRecipes::Structures::Smoothing smoothing;
              double populationWeight1;
              bool isInternalAnchor;
              double reliabilityCommonItemsPopulation1;
              double reliabilityCommonItemsPopulation2;
              EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataXV;
              EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataYV;
              size_t numberOfDegreesSmoothingU;
              size_t numberOfDegreesSmoothingV;
              size_t numberOfCrossProductMoments;
              Eigen::MatrixXi crossProductMatrix;
              bool useStandardizedScale;
              EquatingRecipes::Structures::DesignMatrixType designMatrixType;
              EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType;
              double criterion;
            };

            struct OutputData {
              EquatingRecipes::Structures::PData pData;
              EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
              EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;
              EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingXV;
              EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingYV;
              EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            };

            nlohmann::json operator()(const InputData& inputData,
                                      OutputData& outputData) {
              EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;

              EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsOuputDataXV;
              nlohmann::json bivariateStatisticsResultsXV = bivariateStatistics(inputData.bivariateStatisticsInputDataXV,
                                                                                bivariateStatisticsOuputDataXV);

              EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsOuputDataYV;
              nlohmann::json bivariateStatisticsResultsYV = bivariateStatistics(inputData.bivariateStatisticsInputDataYV,
                                                                                bivariateStatisticsOuputDataYV);

              outputData.bivariateStatisticsXV = bivariateStatisticsOuputDataXV.bivariateStatistics;
              outputData.bivariateStatisticsYV = bivariateStatisticsOuputDataYV.bivariateStatistics;

              EquatingRecipes::Implementation::LogLinearEquating logLinearEquating;

              outputData.equatedRawScoreResults =
                  logLinearEquating.runCGEquiEquatingWithLoglinearSmoothing(inputData.design,
                                                                            inputData.method,
                                                                            inputData.smoothing,
                                                                            inputData.populationWeight1,
                                                                            inputData.isInternalAnchor,
                                                                            inputData.reliabilityCommonItemsPopulation1,
                                                                            inputData.reliabilityCommonItemsPopulation2,
                                                                            outputData.bivariateStatisticsXV,
                                                                            outputData.bivariateStatisticsYV,
                                                                            outputData.bivariateLogLinearSmoothingXV,
                                                                            outputData.bivariateLogLinearSmoothingYV,
                                                                            0,
                                                                            outputData.pData);

              nlohmann::json results = nlohmann::json::object();

              results["DatasetName"] = inputData.datasetName;
              results["RowVariableName"] = inputData.xVariableName;
              results["ColumnwVariableName"] = inputData.yVariableName;
              results["PData"] = outputData.pData;
              results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
              results["LogLinearSmoothingXV"] = outputData.bivariateLogLinearSmoothingXV;
              results["LogLinearSmoothingYV"] = outputData.bivariateLogLinearSmoothingYV;

              nlohmann::json loglinearEquatingResults = nlohmann::json {{"analysis_type", "cine_group_log_linear_smoothing"},
                                                                        {"analysis_results", results}};

              nlohmann::json j = nlohmann::json::array();

              j.push_back(bivariateStatisticsResultsXV);
              j.push_back(bivariateStatisticsResultsYV);
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