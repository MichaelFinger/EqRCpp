#ifndef ANALYSES_KERNEL_EQUATING_CG_HPP
#define ANALYSES_KERNEL_EQUATING_CG_HPP

#include <string>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/analyses/bivariate_statistics.hpp>
#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/implementation/continuized_log_linear_equating.hpp>
#include <equating_recipes/implementation/utilities.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>
#include <equating_recipes/structures/criterion_comparison_type.hpp>

namespace EquatingRecipes {
  namespace Analyses {
    namespace KernelEquating {
      struct CINEGroupsEquating {
        struct InputData {
          std::string title;
          std::string datasetName;
          EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataXV;
          EquatingRecipes::Analyses::BivariateStatistics::InputData bivariateStatisticsInputDataYV;

          EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;
          EquatingRecipes::Structures::DesignMatrixType designMatrixType = EquatingRecipes::Structures::DesignMatrixType::RAW_SCORE;
          EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType = EquatingRecipes::Structures::CriterionComparisonType::ABSOLUTE;
          double criterion = 0.000001;
          bool isInternalAnchor;
          size_t numberOfDegreesSmoothingX = 6;
          size_t numberOfDegreesSmoothingY = 6;
          size_t numberOfDegreesSmoothingV = 6;
          size_t numberOfCrossProductMoments = 3;
          Eigen::MatrixXi crossProductMatrix;
          bool useStandardizedScale = true;
          double population1Weight;
          double reliabilityCommonItemsPopulation1;
          double reliabilityCommonItemsPopulation2;
        };

        struct OutputData {
          EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsXV;
          EquatingRecipes::Structures::BivariateStatistics bivariateStatisticsYV;
          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingXV;
          EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingYV;
          EquatingRecipes::Structures::PData pData;
          EquatingRecipes::Structures::EquatedRawScoreResults equatedRawScoreResults;
          Eigen::VectorXd unroundedEquatedScaledScores;
          Eigen::VectorXd roundedEquatedScaledScores;
        };

        nlohmann::json operator()(const InputData& inputData,
                                  OutputData& outputData) {
          EquatingRecipes::Analyses::BivariateStatistics bivariateStatistics;

          EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsOutputDataXV;
          EquatingRecipes::Analyses::BivariateStatistics::OutputData bivariateStatisticsOutputDataYV;

          nlohmann::json bivariateStatisticsResultsXV = bivariateStatistics(inputData.bivariateStatisticsInputDataXV,
                                                                            bivariateStatisticsOutputDataXV);

          nlohmann::json bivariateStatisticsResultsYV = bivariateStatistics(inputData.bivariateStatisticsInputDataYV,
                                                                            bivariateStatisticsOutputDataYV);

          outputData.bivariateStatisticsXV = bivariateStatisticsOutputDataXV.bivariateStatistics;
          outputData.bivariateStatisticsYV = bivariateStatisticsOutputDataYV.bivariateStatistics;

          EquatingRecipes::Implementation::LogLinearEquating logLinearEquating;

          outputData.bivariateLogLinearSmoothingXV =
              logLinearEquating.runBivariateLogLinearSmoothing(outputData.bivariateStatisticsXV,
                                                               inputData.isInternalAnchor,
                                                               inputData.numberOfDegreesSmoothingX,
                                                               inputData.numberOfDegreesSmoothingV,
                                                               inputData.numberOfCrossProductMoments,
                                                               inputData.crossProductMatrix,
                                                               inputData.useStandardizedScale,
                                                               inputData.designMatrixType,
                                                               inputData.criterionComparisonType,
                                                               inputData.criterion);

          outputData.bivariateLogLinearSmoothingXV =
              logLinearEquating.runBivariateLogLinearSmoothing(outputData.bivariateStatisticsXV,
                                                               inputData.isInternalAnchor,
                                                               inputData.numberOfDegreesSmoothingY,
                                                               inputData.numberOfDegreesSmoothingV,
                                                               inputData.numberOfCrossProductMoments,
                                                               inputData.crossProductMatrix,
                                                               inputData.useStandardizedScale,
                                                               inputData.designMatrixType,
                                                               inputData.criterionComparisonType,
                                                               inputData.criterion);

          EquatingRecipes::Implementation::ContinuizedLogLinearEquating continuizedLogLinearEquating;
          continuizedLogLinearEquating.runWithCGDesign(EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS,
                                                       EquatingRecipes::Structures::Method::EQUIPERCENTILE,
                                                       EquatingRecipes::Structures::Smoothing::CONTINUIZED_LOG_LINEAR_EQUATING,
                                                       inputData.population1Weight,
                                                       inputData.isInternalAnchor,
                                                       inputData.reliabilityCommonItemsPopulation1,
                                                       inputData.reliabilityCommonItemsPopulation2,
                                                       outputData.bivariateStatisticsXV,
                                                       outputData.bivariateStatisticsYV,
                                                       outputData.bivariateLogLinearSmoothingXV,
                                                       outputData.bivariateLogLinearSmoothingYV,
                                                       0,
                                                       outputData.pData,
                                                       outputData.equatedRawScoreResults);

          outputData.unroundedEquatedScaledScores(outputData.equatedRawScoreResults.equatedRawScores.size());
          outputData.roundedEquatedScaledScores(outputData.equatedRawScoreResults.equatedRawScores.size());

          EquatingRecipes::Implementation::Utilities::getEquatedScaledScores(inputData.bivariateStatisticsInputDataXV.rowMinimumScore,
                                                                             inputData.bivariateStatisticsInputDataXV.rowMaximumScore,
                                                                             inputData.bivariateStatisticsInputDataXV.rowScoreIncrement,
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

          results["PData"] = outputData.pData;
          results["BivariateStatisticsXV"] = outputData.bivariateStatisticsXV;
          results["BivariateStatisticsYV"] = outputData.bivariateStatisticsYV;
          results["BivariateLogLinearSmoothingXV"] = outputData.bivariateLogLinearSmoothingXV;
          results["BivariateLogLinearSmoothingYV"] = outputData.bivariateLogLinearSmoothingYV;
          results["EquatedRawScoreResults"] = outputData.equatedRawScoreResults;
          results["UnroundedEquatedScaledScores"] = outputData.unroundedEquatedScaledScores;
          results["RoundedEquatedScaledScores"] = outputData.roundedEquatedScaledScores;

          std::string analysisType = fmt::format("{}_continuized_log_linear_equating",
                                                 EquatingRecipes::Implementation::Utilities::getDesignName(EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS));

          nlohmann::json kernelEquatingResults = nlohmann::json {{"analysis_type", analysisType},
                                                                 {"analysis_results", results}};

          nlohmann::json j = nlohmann::json::array();

          j.push_back(bivariateStatisticsResultsXV);
          j.push_back(bivariateStatisticsResultsYV);
          j.push_back(kernelEquatingResults);

          return j;
        }
      };
    } // namespace KernelEquating
  }   // namespace Analyses
} // namespace EquatingRecipes

#endif