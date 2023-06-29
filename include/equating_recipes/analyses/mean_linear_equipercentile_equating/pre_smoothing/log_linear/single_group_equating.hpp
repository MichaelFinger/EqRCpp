#ifndef ANALYSES_SG_MEAN_LIN_EQUI_EQ_LOG_LIN_SMOOTHING_HPP
#define ANALYSES_SG_MEAN_LIN_EQUI_EQ_LOG_LIN_SMOOTHING_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>
#include <fmt/core.h>

#include <equating_recipes/analyses/bivariate_statistics.hpp>
#include <equating_recipes/implementation/beta_binomial.hpp>
#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/json/structures.hpp>
#include <equating_recipes/structures/criterion_comparison_type.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace MeanLinearEquipercentileEquating {
      namespace PreSmoothing {
        namespace LogLinear {
          struct SingleGroupEquating {
            struct InputData {
              std::string datasetName;
              std::string xVariableName;
              std::string yVariableName;
              EquatingRecipes::Structures::Design design;
              EquatingRecipes::Structures::Method method;
              EquatingRecipes::Structures::Smoothing smoothing;
              EquatingRecipes::Analyses::BivariateStatistics::InputData inputDataUV;
              size_t numberOfDegreesSmoothngU;
              size_t numberOfDegreesSmoothngV;
              size_t numberOfCrossProductMoments;
              Eigen::MatrixXi crossProductMatrix;
              bool useStandardizedScale;
              EquatingRecipes::Structures::DesignMatrixType designMatrixType;
              EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType;
              double criterion;
              bool isInternalAnchor;
            };

            struct OutputData {
              EquatingRecipes::Structures::PData pData;
              EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsUV;
              EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingUV;
              EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
            };

            nlohmann::json operator()(const InputData& inputData,
                                      OutputData& outputData) {
              EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;
              EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsOuputData;
              nlohmann::json bivariateStatisticsResults = bivariateStatistics(inputData.inputDataUV,
                                                                              bivariateStatisticsOuputData);

              EquatingRecipes::Implementation::LogLinearEquating logLinearEquating;

              outputData.bivariateLogLinearSmoothingUV = logLinearEquating.runBivariateLogLinearSmoothing(outputData.bivariateStatisticsUV,
                                                                                                          inputData.isInternalAnchor,
                                                                                                          inputData.numberOfDegreesSmoothngU,
                                                                                                          inputData.numberOfDegreesSmoothngV,
                                                                                                          inputData.numberOfCrossProductMoments,
                                                                                                          inputData.crossProductMatrix,
                                                                                                          inputData.useStandardizedScale,
                                                                                                          inputData.designMatrixType,
                                                                                                          inputData.criterionComparisonType,
                                                                                                          inputData.criterion);

              outputData.equatedRawScoreResults = logLinearEquating.runSGEquiEquatingWithLoglinearSmoothing(inputData.design,
                                                                                                            inputData.method,
                                                                                                            inputData.smoothing,
                                                                                                            outputData.bivariateStatisticsUV,
                                                                                                            outputData.bivariateLogLinearSmoothingUV,
                                                                                                            0,
                                                                                                            outputData.pData);

              nlohmann::json results = nlohmann::json::object();
              
              results["DatasetName"] = inputData.datasetName;
              results["RowVariableName"] = inputData.xVariableName;
              results["ColumnVariableName"] = inputData.yVariableName;
              results["PData"] = outputData.pData;
              results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
              results["LogLinearSmoothingUV"] = outputData.bivariateLogLinearSmoothingUV;

              std::string analysisType = fmt::format("{}_mean_linear_equipercentile_log_linear_equating",
                                                   EquatingRecipes::Implementation::Utilities::getDesignName(inputData.design));

              nlohmann::json loglinearEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                        {"analysis_results", results}};

              nlohmann::json j = nlohmann::json::array();

              j.push_back(bivariateStatisticsResults);
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