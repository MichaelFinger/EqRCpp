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
          EquatingRecipes::Structures::Design design;
          EquatingRecipes::Structures::Method method;
          EquatingRecipes::Structures::Smoothing smoothing;
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

          size_t numberOfScoresX = EquatingRecipes::Implementation::Utilities::getNumberOfScores(inputData.univariateStatisticsInputDataX.minimumScore,
                                                                                                 inputData.univariateStatisticsInputDataX.maximumScore,
                                                                                                 inputData.univariateStatisticsInputDataX.scoreIncrement);

          size_t numberOfScoresY = EquatingRecipes::Implementation::Utilities::getNumberOfScores(inputData.univariateStatisticsInputDataY.minimumScore,
                                                                                                 inputData.univariateStatisticsInputDataY.maximumScore,
                                                                                                 inputData.univariateStatisticsInputDataY.scoreIncrement);
          Eigen::VectorXd scoresX(numberOfScoresX);
          for (size_t scoreIndex = 0; scoreIndex < numberOfScoresX; scoreIndex) {
            scoresX(scoreIndex) = inputData.univariateStatisticsInputDataX.minimumScore +
                                  static_cast<double>(scoreIndex) * inputData.univariateStatisticsInputDataX.scoreIncrement;
          }

          Eigen::VectorXd scoresY(numberOfScoresY);
          for (size_t scoreIndex = 0; scoreIndex < numberOfScoresY; scoreIndex) {
            scoresY(scoreIndex) = inputData.univariateStatisticsInputDataY.minimumScore +
                                  static_cast<double>(scoreIndex) * inputData.univariateStatisticsInputDataY.scoreIncrement;
          }

          EquatingRecipes::Implementation::KernelEquating kernelEquating;
          kernelEquating.runWithRGDesign(inputData.design,
                                         inputData.method,
                                         inputData.smoothing,
                                         outputData.univariateStatisticsX,
                                         outputData.univariateStatisticsY,
                                         outputData.univariateLogLinearSmoothingX,
                                         outputData.univariateLogLinearSmoothingY,
                                         0,
                                         outputData.pData,
                                         outputData.equatedRawScoreResults);

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

          std::string analysisType = fmt::format("{}_kernel_log_linear_equating",
                                                 EquatingRecipes::Implementation::Utilities::getDesignName(inputData.design));

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