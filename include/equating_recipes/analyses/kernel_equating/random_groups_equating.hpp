#ifndef ANALYSES_KERNEL_EQUATING_RG_HPP
#define ANALYSES_KERNEL_EQUATING_RG_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/analyses/univariate_statistics.hpp>
#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/implementation/kernel_equating.hpp>
#include <equating_recipes/implementation/utilities.hpp>
#include <equating_recipes/implementation/kernel_equating.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>
#include <equating_recipes/structures/criterion_comparison_type.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace KernelEquating {
      struct RandomGroupsEquating {
        struct InputData {
          std::string title;
          std::string datasetName;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataX;
          EquatingRecipes::Analyses::UnivariateStatistics::InputData univariateStatisticsInputDataY;
          EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;
          size_t numberOfDegreesSmoothing = 6;
          bool useStandardizedScale = true;
          EquatingRecipes::Structures::DesignMatrixType designMatrixType = EquatingRecipes::Structures::DesignMatrixType::RAW_SCORE;
          EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType = EquatingRecipes::Structures::CriterionComparisonType::ABSOLUTE;
          double criterion = 0.000001;
        };

        struct OutputData {
          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX;
          EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY;
          EquatingRecipes::Structures::UnivariateLogLinearSmoothing univariateLogLinearSmoothingX;
          EquatingRecipes::Structures::UnivariateLogLinearSmoothing univariateLogLinearSmoothingY;
          EquatingRecipes::Structures::PData pData;
          EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          Eigen::VectorXd unroundedEquatedScaledScores;
          Eigen::VectorXd roundedEquatedScaledScores;
        };

        nlohmann::json operator()(const InputData& inputData,
                                  OutputData& outputData) {
          EquatingRecipes::Analyses::UnivariateStatistics univariateStatistics;
          EquatingRecipes::Analyses::UnivariateStatistics::OutputData univariateStatisticsOutputDataX;
          EquatingRecipes::Analyses::UnivariateStatistics::OutputData univariateStatisticsOutputDataY;

          nlohmann::json univariateStatisticsXResults = univariateStatistics(inputData.univariateStatisticsInputDataX,
                                                                             univariateStatisticsOutputDataX);

          nlohmann::json univariateStatisticsYResults = univariateStatistics(inputData.univariateStatisticsInputDataY,
                                                                             univariateStatisticsOutputDataY);

          outputData.univariateStatisticsX = univariateStatisticsOutputDataX.univariateStatistics;
          outputData.univariateStatisticsY = univariateStatisticsOutputDataY.univariateStatistics;

          EquatingRecipes::Implementation::LogLinearEquating logLinearEquating;
          outputData.univariateLogLinearSmoothingX = logLinearEquating.runUnivariateLogLinearSmoothing(outputData.univariateStatisticsX,
                                                                                                       inputData.numberOfDegreesSmoothing,
                                                                                                       inputData.useStandardizedScale,
                                                                                                       inputData.designMatrixType,
                                                                                                       inputData.criterionComparisonType,
                                                                                                       inputData.criterion);

          outputData.univariateLogLinearSmoothingY = logLinearEquating.runUnivariateLogLinearSmoothing(outputData.univariateStatisticsY,
                                                                                                       inputData.numberOfDegreesSmoothing,
                                                                                                       inputData.useStandardizedScale,
                                                                                                       inputData.designMatrixType,
                                                                                                       inputData.criterionComparisonType,
                                                                                                       inputData.criterion);

          /*
            Kernel equating for Random Groups Design

  Wrapper_RK('R','E','K', &x, &y, &ullx, &ully, 0, &pdREK,&rREK);
  Print_RK(outf,"ACT Math---Kernel---Log Linear Smoothing",&pdREK, &rREK);

  Wrapper_ESS(&pdREK,&rREK,0,40,1,"yctmath.TXT",1,1,36,&sREK);
  Print_ESS(outf,"ACT Math---Equipercentile",&pdREK,&sREK);

          */

          EquatingRecipes::Implementation::KernelEquating kernelEquating;

          kernelEquating.runWithRGDesign(EquatingRecipes::Structures::Design::RANDOM_GROUPS,
                                         EquatingRecipes::Structures::Method::EQUIPERCENTILE,
                                         EquatingRecipes::Structures::Smoothing::KERNEL,
                                         outputData.univariateStatisticsX,
                                         outputData.univariateStatisticsY,
                                         outputData.univariateLogLinearSmoothingX,
                                         outputData.univariateLogLinearSmoothingY,
                                         0,
                                         outputData.pData,
                                         outputData.equatedRawScoreResults);

          outputData.unroundedEquatedScaledScores(outputData.equatedRawScoreResults.equatedRawScores.size());
          outputData.roundedEquatedScaledScores(outputData.equatedRawScoreResults.equatedRawScores.size());

          EquatingRecipes::Implementation::Utilities::getEquatedScaledScores(inputData.univariateStatisticsInputDataX.minimumScore,
                                                                             inputData.univariateStatisticsInputDataX.maximumScore,
                                                                             inputData.univariateStatisticsInputDataX.scoreIncrement,
                                                                             outputData.pData.minimumRawScoreYct,
                                                                             outputData.pData.maximumRawScoreYct,
                                                                             outputData.pData.scoreIncrementYct,
                                                                             outputData.equatedRawScoreResults.equatedRawScores,
                                                                             inputData.rawToScaledScoreTable,
                                                                             0,
                                                                             inputData.rawToScaledScoreTable.lowestObservableScaledScore,
                                                                             inputData.rawToScaledScoreTable.highestObservableScaledScore,
                                                                             outputData.unroundedEquatedScaledScores,
                                                                             outputData.roundedEquatedScaledScores);

          nlohmann::json results = nlohmann::json::object();

          results["DatasetName"] = inputData.datasetName;
          results["RowVariableName"] = outputData.univariateStatisticsX.variableName;
          results["ColumnVariableName"] = outputData.univariateStatisticsY.variableName;
          results["PData"] = outputData.pData;
          results["UnivariateStatisticsX"] = outputData.univariateStatisticsX;
          results["UnivariateStatisticsY"] = outputData.univariateStatisticsY;
          results["UnivariateLogLinearSmoothingX"] = outputData.univariateLogLinearSmoothingX;
          results["UnivariateLogLinearSmoothingY"] = outputData.univariateLogLinearSmoothingY;
          results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
          results["UnroundedEquatedScaledScores"] = outputData.unroundedEquatedScaledScores;
          results["RoundedEquatedScaledScores"] = outputData.roundedEquatedScaledScores;

          std::string analysisType = fmt::format("{}_kernel_log_linear_equating",
                                                 EquatingRecipes::Implementation::Utilities::getDesignName(EquatingRecipes::Structures::Design::RANDOM_GROUPS));

          nlohmann::json kernelEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                 {"analysis_results", results}};

          nlohmann::json j = nlohmann::json::array();

          j.push_back(univariateStatisticsXResults);
          j.push_back(univariateStatisticsYResults);
          j.push_back(kernelEquatingResults);

          return j;
        }
      };
    } // namespace KernelEquating
  }   // namespace Analyses
} // namespace EquatingRecipes

#endif