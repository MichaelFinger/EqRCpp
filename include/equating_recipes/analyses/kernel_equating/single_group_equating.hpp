#ifndef ANALYSES_KERNEL_EQUATING_SG_HPP
#define ANALYSES_KERNEL_EQUATING_SG_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/analyses/bivariate_statistics.hpp>
#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/implementation/kernel_equating.hpp>
#include <equating_recipes/implementation/utilities.hpp>
#include <equating_recipes/implementation/kernel_equating.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>
#include <equating_recipes/structures/criterion_comparison_type.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace KernelEquating {
      struct SingleGroupsEquating {
        struct InputData {
          std::string title;
          std::string datasetName;
          EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputData;
          EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;
          EquatingRecipes::Structures::DesignMatrixType designMatrixType = EquatingRecipes::Structures::DesignMatrixType::RAW_SCORE;
          EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType = EquatingRecipes::Structures::CriterionComparisonType::ABSOLUTE;
          double criterion = 0.000001;
          bool isInternalAnchor;
          size_t numberOfDegreesSmoothingX = 6;
          size_t numberOfDegreesSmoothingY = 6;
          size_t numberOfCrossProductMoments;
          Eigen::MatrixXi crossProductMatrix;
          bool useStandardizedScale = true;
        };

        struct OutputData {
          EquatingRecipes::Structures::BivariateStatistics bivariateStatistics;
          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothing;
          EquatingRecipes::Structures::PData pData;
          EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          Eigen::VectorXd unroundedEquatedScaledScores;
          Eigen::VectorXd roundedEquatedScaledScores;
        };

        nlohmann::json operator()(const InputData& inputData,
                                  OutputData& outputData) {
          // Bivariate stats
          EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;
          EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsOutputData;

          nlohmann::json bivariateStatisticsResults = bivariateStatistics(inputData.bivariateStatisticsInputData,
                                                                          bivariateStatisticsOutputData);

          outputData.bivariateStatistics = bivariateStatisticsOutputData.bivariateStatistics;

          //log-linear smoothing of the bivariate distribution
          EquatingRecipes::Implementation::LogLinearEquating logLinearEquating;
          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothing =
              logLinearEquating.runBivariateLogLinearSmoothing(outputData.bivariateStatistics,
                                                               inputData.isInternalAnchor,
                                                               inputData.numberOfDegreesSmoothingX,
                                                               inputData.numberOfDegreesSmoothingY,
                                                               inputData.numberOfCrossProductMoments,
                                                               inputData.crossProductMatrix,
                                                               inputData.useStandardizedScale,
                                                               inputData.designMatrixType,
                                                               inputData.criterionComparisonType,
                                                               inputData.criterion);

          EquatingRecipes::Implementation::KernelEquating kernelEquating;
          kernelEquating.runWithSGDesign(EquatingRecipes::Structures::Design::SINGLE_GROUP,
                                         EquatingRecipes::Structures::Method::EQUIPERCENTILE,
                                         EquatingRecipes::Structures::Smoothing::KERNEL,
                                         outputData.bivariateStatistics,
                                         outputData.bivariateLogLinearSmoothing,
                                         0,
                                         outputData.pData,
                                         outputData.equatedRawScoreResults);

          outputData.unroundedEquatedScaledScores(outputData.equatedRawScoreResults.equatedRawScores.size());
          outputData.roundedEquatedScaledScores(outputData.equatedRawScoreResults.equatedRawScores.size());

          EquatingRecipes::Implementation::Utilities::getEquatedScaledScores(outputData.bivariateStatistics.univariateStatisticsRow.minimumScore,
                                                                             outputData.bivariateStatistics.univariateStatisticsRow.maximumScore,
                                                                             outputData.bivariateStatistics.univariateStatisticsRow.scoreIncrement,
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
          results["RowVariableName"] = outputData.bivariateStatistics.rowVariableName;
          results["ColumnVariableName"] = outputData.bivariateStatistics.columnVariableName;
          results["PData"] = outputData.pData;
          results["BivariateStatisticsXY"] = outputData.bivariateStatistics;
          results["BivariateLogLinearSmoothingXY"] = outputData.bivariateLogLinearSmoothing;
          results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
          results["UnroundedEquatedScaledScores"] = outputData.unroundedEquatedScaledScores;
          results["RoundedEquatedScaledScores"] = outputData.roundedEquatedScaledScores;

          std::string analysisType = fmt::format("{}_kernel_log_linear_equating",
                                                 EquatingRecipes::Implementation::Utilities::getDesignName(EquatingRecipes::Structures::Design::RANDOM_GROUPS));

          nlohmann::json kernelEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                 {"analysis_results", results}};

          nlohmann::json j = nlohmann::json::array();

          j.push_back(bivariateStatisticsResults);
            j.push_back(kernelEquatingResults);

          return j;
        }
      };
    } // namespace KernelEquating
  }   // namespace Analyses
} // namespace EquatingRecipes

#endif