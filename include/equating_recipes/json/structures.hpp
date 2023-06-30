#ifndef JSON_IMPLEMENTATION_HPP
#define JSON_IMPLEMENTATION_HPP

#include <map>
#include <set>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/bootstrap_equated_raw_score_results.hpp>
#include <equating_recipes/structures/bootstrap_equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/cg_equipercentile_equating_results.hpp>
#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/cubic_spline_postsmoothing.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_bootstrap_results.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_score_bootstrapt_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/form_type.hpp>
#include <equating_recipes/structures/irt_equating_results.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_input.hpp>
#include <equating_recipes/structures/irt_method.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/irt_scale_transformation_data.hpp>
#include <equating_recipes/structures/irt_scale_transformation_item_results.hpp>
#include <equating_recipes/structures/irt_scale_transformation_method.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/loss_spec.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/optimization_results.hpp>
#include <equating_recipes/structures/criterion_comparison_type.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>


namespace Eigen {
  template<typename Derived>
  void to_json(nlohmann::json& j, const Eigen::DenseBase<Derived>& matrix) {
    j = nlohmann::json::array();

    for (size_t rowIndex = 0; rowIndex < matrix.rows(); rowIndex++) {
      nlohmann::json row = nlohmann::json::array();

      for (size_t columnIndex = 0; columnIndex < matrix.cols(); columnIndex++) {
        row.push_back(matrix(rowIndex, columnIndex));
      }

      j.push_back(row);
    }
  }

  void from_json(const nlohmann::json& j, Eigen::MatrixXd& matrix) {
    size_t numberOfRows = j.size();
    size_t numberOfColumns = j[0].size();

    matrix.resize(numberOfRows, numberOfColumns);

    for (size_t rowIndex = 0; rowIndex < numberOfRows; rowIndex++) {
      nlohmann::json jsonRow = j[rowIndex];

      for (size_t columnIndex = 0; columnIndex < numberOfColumns; columnIndex++) {
        matrix(rowIndex, columnIndex) = jsonRow[columnIndex];
      }
    }
  }
} // namespace Eigen

namespace EquatingRecipes {
  namespace Structures {
    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::Design, {{EquatingRecipes::Structures::Design::NOT_SPECIFIED, ""},
                                                                       {EquatingRecipes::Structures::Design::RANDOM_GROUPS, "random groups"},
                                                                       {EquatingRecipes::Structures::Design::SINGLE_GROUP, "single group"},
                                                                       {EquatingRecipes::Structures::Design::COMMON_ITEN_NON_EQUIVALENT_GROUPS, "common iten non-equivalent groups"},
                                                                       {EquatingRecipes::Structures::Design::SINGLE_GROUP_COUNTER_BALANCE, "single group counter-balance"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::FormType, {{EquatingRecipes::Structures::FormType::NEW, "new"},
                                                                         {EquatingRecipes::Structures::FormType::OLD, "old"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::IRTMethod, {{EquatingRecipes::Structures::IRTMethod::NOT_SPECIFIED, ""},
                                                                          {EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE, "observed score"},
                                                                          {EquatingRecipes::Structures::IRTMethod::TRUE_SCORE, "true score"},
                                                                          {EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE, "true and observed score"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::IRTModel, {{EquatingRecipes::Structures::IRTModel::NOT_SPECIFIED, ""},
                                                                         {EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC, "3-parameter logistic"},
                                                                         {EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE, "graded response"},
                                                                         {EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT, "generaized partial credit"},
                                                                         {EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE, "nominal response"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::IRTScaleTransformationMethod, {{EquatingRecipes::Structures::IRTScaleTransformationMethod::NOT_SPECIFIED, ""},
                                                                                             {EquatingRecipes::Structures::IRTScaleTransformationMethod::HAEBARA, "Haebara"},
                                                                                             {EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_MEAN, "mean mean"},
                                                                                             {EquatingRecipes::Structures::IRTScaleTransformationMethod::MEAN_SIGMA, "mean sigma"},
                                                                                             {EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD, "Stocking-Lord"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::LossSpec, {{EquatingRecipes::Structures::LossSpec::NOT_SPECIFIED, ""},
                                                                         {EquatingRecipes::Structures::LossSpec::MIX_HA, "Haebara"},
                                                                         {EquatingRecipes::Structures::LossSpec::MIX_SL, "Stocking-Lord"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::Method, {{EquatingRecipes::Structures::Method::NOT_SPECIFIED, ""},
                                                                       {EquatingRecipes::Structures::Method::MEAN, "mean"},
                                                                       {EquatingRecipes::Structures::Method::LINEAR, "linear"},
                                                                       {EquatingRecipes::Structures::Method::EQUIPERCENTILE, "equipercentile"},
                                                                       {EquatingRecipes::Structures::Method::FE_BH, "frequency estimation with Braun-Holland"},
                                                                       {EquatingRecipes::Structures::Method::MFE_BH, "modified frequency estimation with Braun-Holland"},
                                                                       {EquatingRecipes::Structures::Method::FE_BH_MFE_BH, "frequency estimation with Braun-Holland & modified estimation with Braun-Holland"},
                                                                       {EquatingRecipes::Structures::Method::CHAINED, "chained"},
                                                                       {EquatingRecipes::Structures::Method::FE_BH_CHAINED, "chained frequency estimation with Braun-Holland"},
                                                                       {EquatingRecipes::Structures::Method::FE_BH_MFE_BH_CHAINED, "chained frequency estimation with Braun-Holland & modified estimation with Braun-Holland"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::OptimizationResultCode, {{EquatingRecipes::Structures::OptimizationResultCode::SUCCESS, "success"},
                                                                                       {EquatingRecipes::Structures::OptimizationResultCode::FAILED, "failed"},
                                                                                       {EquatingRecipes::Structures::OptimizationResultCode::FUNCTION_STOP_VALUE, "function_stop_value"},
                                                                                       {EquatingRecipes::Structures::OptimizationResultCode::FUNCTION_TOLERANCE, "function_tolerance"},
                                                                                       {EquatingRecipes::Structures::OptimizationResultCode::PARAMETER_TOLERANCE, "parameter_tolerance"},
                                                                                       {EquatingRecipes::Structures::OptimizationResultCode::MAXIMUM_NUMBER_OF_ITERATIONS, "maximum_number_of_iterations"},
                                                                                       {EquatingRecipes::Structures::OptimizationResultCode::MAXIMUM_TIME, "maximum_time"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::Smoothing, {{EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED, ""},
                                                                          {EquatingRecipes::Structures::Smoothing::LOG_LINEAR, "log linear"},
                                                                          {EquatingRecipes::Structures::Smoothing::BETA_BINOMIAL, "beta binomial"},
                                                                          {EquatingRecipes::Structures::Smoothing::CUBIC_SPLINE, "cubic spline"},
                                                                          {EquatingRecipes::Structures::Smoothing::KERNEL, "kernel"},
                                                                          {EquatingRecipes::Structures::Smoothing::CONTINUIZED_LOG_LINEAR_EQUATING, "continuized log linear equating"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::Symmetry, {{EquatingRecipes::Structures::Symmetry::NOT_SPECIFIED, ""},
                                                                         {EquatingRecipes::Structures::Symmetry::NEW_SCALE, "new scale"},
                                                                         {EquatingRecipes::Structures::Symmetry::OLD_SCALE, "old scale"},
                                                                         {EquatingRecipes::Structures::Symmetry::SYMMETRIC, "symmetric"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::CriterionComparisonType, {{EquatingRecipes::Structures::CriterionComparisonType::ABSOLUTE, "absolute"},
                                                                                               {EquatingRecipes::Structures::CriterionComparisonType::RELATIVE, "relative"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::DesignMatrixType, {{EquatingRecipes::Structures::DesignMatrixType::SOLUTION, "solution"},
                                                                                        {EquatingRecipes::Structures::DesignMatrixType::RAW_SCORE, "raw score"}})

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::UnivariateStatistics& rec) {
      j = nlohmann::json {{"id", rec.id},
                          {"datasetName", rec.datasetName},
                          {"variableName", rec.variableName},
                          {"numberOfExaminees", rec.numberOfExaminees},
                          {"minimumScore", rec.minimumScore},
                          {"maximumScore", rec.maximumScore},
                          {"frequencyDisributiontMinimumScoreObserved", rec.freqDistMinimumScore},
                          {"frequencyDisributiontMaximumScoreObserved", rec.freqDistMaximumScore},
                          {"scoreIncrement", rec.scoreIncrement},
                          {"numberOfScores", rec.numberOfScores},
                          {"frequencyDistribution", rec.freqDist},
                          {"frequencyDistributionDouble", rec.freqDistDouble},
                          {"cumulativeFrequencyDistribution", rec.cumulativeFreqDist},
                          {"relativeFrequencyDistribution", rec.relativeFreqDist},
                          {"cumulativeRelativeFrequencyDistribution", rec.cumulativeRelativeFreqDist},
                          {"percentileRankDistribution", rec.percentileRankDist},
                          {"moments", rec.momentValues},
                          {"rawScores", rec.rawScores}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BivariateLogLinearSmoothing& rec) {
      j = nlohmann::json {{"isInternalAnchor", rec.isInternalAnchor},
                          {"numberOfExaminees", rec.numberOfExaminees},
                          {"numberOfScoresU", rec.numberOfScoresU},
                          {"minimumRawScoreU", rec.minimumRawScoreU},
                          {"scoreIncrementU", rec.scoreIncrementU},
                          {"numberOfScoresV", rec.numberOfScoresV},
                          {"minimumRawScoreV", rec.minimumRawScoreV},
                          {"scoreIncrementV", rec.scoreIncrementV},
                          {"totalNumberOfScores", rec.totalNumberOfScores},
                          {"numberOfDegreesOfSmoothingU", rec.numberOfDegreesOfSmoothingU},
                          {"numberOfDegreesOfSmoothingV", rec.numberOfDegreesOfSmoothingV},
                          {"numberOfCrossProductMoments", rec.numberOfCrossProductMoments},
                          {"crossProductMoments", rec.crossProductMoments},
                          {"numberOfColumnsInDesignMatrix", rec.numberOfColumnsInDesignMatrix},
                          {"rawScoreDesignMatrix", rec.rawScoreDesignMatrix},
                          {"solutionDesignMatrix", rec.solutionDesignMatrix},
                          {"observedFrequencies", rec.observedFrequencies},
                          {"fittedFrequencies", rec.fittedFrequencies},
                          {"betaCoefficients", rec.betaCoefficients},
                          {"cllNormalizingConstant", rec.cllNormalizingConstant},
                          {"solutionDesignMatrixObservedMoments", rec.solutionDesignMatrixObservedMoments},
                          {"solutionDesignMatrixFittedMoments", rec.solutionDesignMatrixFittedMoments},
                          {"rawScoreDesignMatrixObservedCentralMoments", rec.rawScoreDesignMatrixObservedCentralMoments},
                          {"rawScoreDesignMatrixFittedCentralMoments", rec.rawScoreDesignMatrixFittedCentralMoments},
                          {"numberOfIterations", rec.numberOfIterations},
                          {"maximumNumberOfIterations", rec.maximumNumberOfIterations},
                          {"criterionComparisonType", rec.criterionComparisonType},
                          {"designMatrixType", rec.designMatrixType},
                          {"useScalingForBDeisngMatrix", rec.useScalingForBDeisngMatrix},
                          {"convergenceCriterion", rec.convergenceCriterion},
                          {"likelihoodRatioChiSquare", rec.likelihoodRatioChiSquare},
                          {"numberOfZeros", rec.numberOfZeros},
                          {"numberOfScoresX", rec.numberOfScoresX},
                          {"minimumRawScoreX", rec.minimumRawScoreX},
                          {"scoreIncrementX", rec.scoreIncrementX},
                          {"fittedBivariateFreqDist", rec.fittedBivariateFreqDist},
                          {"fittedFrequencesX", rec.fittedFrequencesX},
                          {"fittedRawScoreDensityX", rec.fittedRawScoreDensityX},
                          {"fittedRawScoreCumulativeRelativeFreqDistX", rec.fittedRawScoreCumulativeRelativeFreqDistX},
                          {"fittedRawScorePercentileRankDistX", rec.fittedRawScorePercentileRankDistX},
                          {"fittedFrequencesV", rec.fittedFrequencesV},
                          {"fittedRawScoreDensityV", rec.fittedRawScoreDensityV},
                          {"fittedRawScoreCumulativeRelativeFreqDistV", rec.fittedRawScoreCumulativeRelativeFreqDistV},
                          {"fittedRawScorePercentileRankDistV", rec.fittedRawScorePercentileRankDistV},
                          {"fittedBivariateRelativeFreqDistXV", rec.fittedBivariateRelativeFreqDistXV},
                          {"cumulativeRelativeFreqDistRowMajorVector", rec.cumulativeRelativeFreqDistRowMajorVector}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& rec) {
      j = nlohmann::json {{"sumAndMeanScores", rec.sumAndMeanScores},
                          {"sumSquareAndSDScores", rec.sumSquareAndSDScores},
                          {"bootstrapStandardErrors", rec.bootstrapStandardErrors}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& rec) {
      j = nlohmann::json {{"unroundedScaledScoresSumsAndMeans", rec.unroundedScaledScoresSumsAndMeans},
                          {"unroundedScaledScoresSumSquaresAndSDs", rec.unroundedScaledScoresSumSquaresAndSDs},
                          {"unroundedScaledScoresBoostrapStandardErrors", rec.unroundedScaledScoresBoostrapStandardErrors},
                          {"roundedScaledScoresSumsAndMeans", rec.roundedScaledScoresSumsAndMeans},
                          {"roundedScaledScoresSumSquaresAndSDs", rec.roundedScaledScoresSumSquaresAndSDs},
                          {"roundedScaledScoresBoostrapStandardErrors", rec.roundedScaledScoresBoostrapStandardErrors}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::CGEquipercentileEquatingResults& rec) {
      j = nlohmann::json {{"syntheticPopulationRelativeFreqDistX", rec.syntheticPopulationRelativeFreqDistX},
                          {"syntheticPopulationRelativeFreqDistY", rec.syntheticPopulationRelativeFreqDistY},
                          {"equatedRawScores", rec.equatedRawScores}};

      if (rec.slope.has_value()) {
        j["slope"] = rec.slope.value();
      }

      if (rec.intercept.has_value()) {
        j["intercept"] = rec.intercept.value();
      }

      if (rec.braunHollandEquatedRawScores.has_value()) {
        j["braunHollandEquatedRawScores"] = rec.braunHollandEquatedRawScores.value();
      }
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::CubicSplinePostsmoothing& rec) {
      j = nlohmann::json {{"numberOfScores", rec.numberOfScores},
                          {"smoothingParameter", rec.smoothingParameter},
                          {"percentileRankLowestScore", rec.percentileRankLowestScore},
                          {"percentileRankHighestScore", rec.percentileRankHighestScore},
                          {"lowestSmoothedPseudoRawScorePercentileRank", rec.lowestSmoothedPseudoRawScorePercentileRank},
                          {"higestSmoothedPseudoRawScorePercentileRank", rec.higestSmoothedPseudoRawScorePercentileRank},
                          {"boundedNumberOfScores", rec.boundedNumberOfScores},
                          {"equipercentileEquivalents", rec.equipercentileEquivalents},
                          {"standardErrors", rec.standardErrors},
                          {"coefficients", rec.coefficients},
                          {"cubicSplineSmoothedEquivalents", rec.cubicSplineSmoothedEquivalents},
                          {"inverseCubicSplineSmoothedEquivalents", rec.inverseCubicSplineSmoothedEquivalents}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::EquatedRawScoreBootstrapResults& rec) {
      j = nlohmann::json {{"sumAndMeanRawScores", rec.sumAndMeanRawScores},
                          {"sumSquaredAndSDRawScores", rec.sumSquaredAndSDRawScores},
                          {"rawScoresBootstrapStandardErrors", rec.rawScoresBootstrapStandardErrors}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::EquatedRawScoreResults& rec) {
      j = nlohmann::json {{"xSyntheticPopulationMean", rec.xSyntheticPopulationMean},
                          {"ySyntheticPopulationMean", rec.ySyntheticPopulationMean},
                          {"xSyntheticPopulationSD", rec.xSyntheticPopulationSD},
                          {"ySyntheticPopulationSD", rec.ySyntheticPopulationSD},
                          {"gammaPopulation1", rec.gammaPopulation1},
                          {"gammaPopulation2", rec.gammaPopulation2},
                          {"slope", rec.slope},
                          {"intercept", rec.intercept},
                          {"equatedRawScores", rec.equatedRawScores},
                          {"equatedRawScoreMoments", rec.equatedRawScoreMoments},
                          {"relativeFreqDistsX", rec.relativeFreqDistsX},
                          {"relativeFreqDistsY", rec.relativeFreqDistsY}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::EquatedScaledScoreBootstrapResults& rec) {
      j = nlohmann::json {{"sumAndMeanUnroundedScaledScores", rec.sumAndMeanUnroundedScaledScores},
                          {"sumSquaredAndSDUnroundedScaledScores", rec.sumSquaredAndSDUnroundedScaledScores},
                          {"unroundedScaledScoresBootstrapStandardErrors", rec.unroundedScaledScoresBootstrapStandardErrors},
                          {"sumAndMeanRoundedScaledScores", rec.sumAndMeanRoundedScaledScores},
                          {"sumSquaredAndSDRoundedScaledScores", rec.sumSquaredAndSDRoundedScaledScores},
                          {"roundedScaledScoresBootstrapStandardErrors", rec.roundedScaledScoresBootstrapStandardErrors}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::EquatedScaledScoresResults& rec) {
      j = nlohmann::json {{"unroundedEquatedScaledScores", rec.unroundedEquatedScaledScores},
                          {"roundedEquatedScaledScores", rec.roundedEquatedScaledScores},
                          {"unroundedEquatedScaledScoreMoments", rec.unroundedEquatedScaledScoreMoments},
                          {"roundedEquatedScaledScoreMoments", rec.roundedEquatedScaledScoreMoments}};
    };

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTEquatingResults& rec) {
      j = nlohmann::json {{"numberOfRawScoreCategoriesNewForm", rec.numberOfRawScoreCategoriesNewForm},
                          {"minimumTrueScoreNewForm", rec.minimumTrueScoreNewForm},
                          {"minimumTrueScoreOldForm", rec.minimumTrueScoreOldForm},
                          {"thetaEquivalentFormXScore", rec.thetaEquivalentFormXScore},
                          {"unroundedEquatedTrueScore", rec.unroundedEquatedTrueScore},
                          {"roundedEquatedTrueScore", rec.roundedEquatedTrueScore},
                          {"unroundedEquatedObservedScore", rec.unroundedEquatedObservedScore},
                          {"roundedEquatedObservedScore", rec.roundedEquatedObservedScore},
                          {"momentsEquatedTrueScores", rec.momentsEquatedTrueScores},
                          {"momentsEquatedObservedScores", rec.momentsEquatedObservedScores}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTFittedDistribution& rec) {
      j = nlohmann::json {{"numberOfRawScoreCategories", rec.numberOfRawScoreCategories},
                          {"rawScores", rec.rawScores},
                          {"fittedDistributionNewGroup", rec.fittedDistributionNewGroup},
                          {"fittedDistributionOldGroup", rec.fittedDistributionOldGroup},
                          {"fittedDistributionSyntheticGroup", rec.fittedDistributionSyntheticGroup},
                          {"momentsFittedDistributionNewGroup", rec.momentsFittedDistributionNewGroup},
                          {"momentsFittedDistributionOldGroup", rec.momentsFittedDistributionOldGroup},
                          {"momentsFittedDistributionSyntheticGroup", rec.momentsFittedDistributionSyntheticGroup}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Moments& rec) {
      j = nlohmann::json {{"moments", rec.momentValues}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Quadrature& rec) {
      j = nlohmann::json {{"thetaValues", rec.thetaValues},
                          {"thetaWeights", rec.thetaWeights}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& rec) {
      j = nlohmann::json {{"numberOfExaminees", rec.numberOfExaminees},
                          {"numberOfScores", rec.numberOfScores},
                          {"minimumRawScore", rec.minimumRawScore},
                          {"rawScoreIncrement", rec.rawScoreIncrement},
                          {"degreesOfSmoothing", rec.degreesOfSmoothing},
                          {"rawScoreDesignMatrix", rec.rawScoreDesignMatrix},
                          {"solutionDesignMatrix", rec.solutionDesignMatrix},
                          {"observedFrequencies", rec.observedFrequencies},
                          {"fittedFrequencies", rec.fittedFrequencies},
                          {"betaCoefficients", rec.betaCoefficients},
                          {"cllNormalizingConstant", rec.cllNormalizingConstant},
                          {"bObservedMoments", rec.bObservedMoments},
                          {"bFittedMoments", rec.bFittedMoments},
                          {"observedCentralMoments", rec.observedCentralMoments},
                          {"fittedCentralMoments", rec.fittedCentralMoments},
                          {"numberOfIterations", rec.numberOfIterations},
                          {"maximumNumberOfIterators", rec.maximumNumberOfIterators},
                          {"useRelativeCriterionComparison", rec.useRelativeCriterionComparison},
                          {"useBRawAndCentralMoments", rec.useBRawAndCentralMoments},
                          {"useScalingForBDesignMatrix", rec.useScalingForBDesignMatrix},
                          {"convergenceCriterion", rec.convergenceCriterion},
                          {"likelihoodRatioChiSquare", rec.likelihoodRatioChiSquare},
                          {"numberOfZeros", rec.numberOfZeros},
                          {"fittedRawScoreDist", rec.fittedRawScoreDist},
                          {"fittedRawScoreCumulativeRelativeDist", rec.fittedRawScoreCumulativeRelativeDist},
                          {"fittedRawScorePercentileRankDist", rec.fittedRawScorePercentileRankDist}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::RawToScaledScoreTable::Entry& rec) {
      j = nlohmann::json {{"scoreLocation", rec.scoreLocation},
                          {"rawScore", rec.rawScore},
                          {"scaledScore", rec.scaledScore}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::RawToScaledScoreTable& rec) {
      j = nlohmann::json {{"entries", rec.lookup}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BetaBinomialSmoothing& rec) {
      j = nlohmann::json {{"numberOfItems", rec.numberOfItems},
                          {"numberOfExaminees", rec.numberOfExaminees},
                          {"numberOfParameters", rec.numberOfParameters},
                          {"reliablilty", rec.reliablilty},
                          {"lordK", rec.lordK},
                          {"betaParameters;", rec.betaParameters},
                          {"rawScoreMoments;", rec.rawScoreMoments},
                          {"trueScoreMoments;", rec.trueScoreMoments},
                          {"likelihoodRatioChiSq", rec.likelihoodRatioChiSq},
                          {"pearsonChiSq", rec.pearsonChiSq},
                          {"numberOfMomentsFit", rec.numberOfMomentsFit},
                          {"fittedRawScoreDensity", rec.fittedRawScoreDensity},
                          {"fittedRawScoreCumulativeRelativeFreqDist", rec.fittedRawScoreCumulativeRelativeFreqDist},
                          {"fittedRawScorePercentileRankDist", rec.fittedRawScorePercentileRankDist},
                          {"rawScoreUnivariateStatistics", rec.rawScoreUnivariateStatistics}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BivariateStatistics& rec) {
      j = nlohmann::json {{"datasetName", rec.datasetName},
                          {"rowVariableName", rec.rowVariableName},
                          {"columnVariableName", rec.columnVariableName},
                          {"univariateStatisticsRow", rec.univariateStatisticsRow},
                          {"univariateStatisticsColumn", rec.univariateStatisticsColumn},
                          {"numberOfExaminees", rec.numberOfExaminees},
                          {"bivariateFreqDist", rec.bivariateFreqDist},
                          {"bivariateFreqDistDouble", rec.bivariateFreqDistDouble},
                          {"covariance", rec.covariance},
                          {"correlation", rec.correlation},
                          {"bivariateProportions", rec.bivariateProportions}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::CommonItemSpecification& rec) {
      j = nlohmann::json {{"newID", rec.newID},
                          {"oldID", rec.oldID},
                          {"numberOfCategories", rec.numberOfCategories},
                          {"scalingConstant", rec.scalingConstant},
                          {"scoringFunctionValues", rec.scoringFunctionValues},
                          {"newA", rec.newA},
                          {"newB", rec.newB},
                          {"newC", rec.newC},
                          {"oldA", rec.oldA},
                          {"oldB", rec.oldB},
                          {"oldC", rec.oldC},
                          {"irtModel", rec.irtModel}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::ItemSpecification& rec) {
      j = nlohmann::json {{"itemID", rec.itemID},
                          {"numberOfCategories", rec.numberOfCategories},
                          {"scalingConstant", rec.scalingConstant},
                          {"scoringFunctionValues", rec.scoringFunctionValues},
                          {"a", rec.a},
                          {"b", rec.b},
                          {"c", rec.c},
                          {"d", rec.d},
                          {"irtModel", rec.irtModel}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTScaleTransformationItemResults& rec) {
      j = nlohmann::json {{"itemID", rec.itemID},
                          {"irtModel", rec.irtModel},
                          {"numberOfCategories", rec.numberOfCategories},
                          {"scalingConstant", rec.scalingConstant},
                          {"scoringFunctionValues", rec.scoringFunctionValues},
                          {"transformedA", rec.transformedA},
                          {"transformedB", rec.transformedB},
                          {"transformedC", rec.transformedC}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::OptimizationResults& rec) {
      j = nlohmann::json {{"functionValue", rec.functionValue},
                          {"gradientVector", rec.gradientVector},
                          {"parameterEstimates", rec.parameterEstimates},
                          {"parameterStartingValues", rec.parameterStartingValues},
                          {"resultCode", rec.resultCode},
                          {"numberOfIterations", rec.numberOfIterations},
                          {"maximumAbsoluteChangeInFunctionValue", rec.maximumAbsoluteChangeInFunctionValue},
                          {"maximumRelativeChangeInFunctionValue", rec.maximumRelativeChangeInFunctionValue},
                          {"maximumAbsoluteChangeInParameterValues", rec.maximumAbsoluteChangeInParameterValues},
                          {"maximumRelativeChangeInParameterValues", rec.maximumRelativeChangeInParameterValues}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTScaleTransformationData& rec) {
      j = nlohmann::json {{"maximumNumberOfIterations", rec.maximumNumberOfIterations},
                          {"maximumAbsoluteChangeInFunctionValue", rec.maximumAbsoluteChangeInFunctionValue},
                          {"maximumRelativeChangeInFunctionValue", rec.maximumRelativeChangeInFunctionValue},
                          {"maximumAbsoluteChangeInParameterValues", rec.maximumAbsoluteChangeInParameterValues},
                          {"maximumRelativeChangeInParameterValues", rec.maximumRelativeChangeInParameterValues},
                          {"quadratureNewForm", rec.quadratureNewForm},
                          {"quadratureOldForm", rec.quadratureOldForm},
                          {"newItems", rec.newItems},
                          {"oldItems", rec.oldItems},
                          {"commonItems", rec.commonItems},
                          {"symmetryOptions", rec.symmetryOptions},
                          {"standardizations", rec.standardizations},
                          {"slopeStartingValue", rec.slopeStartingValue},
                          {"interceptStartingValue", rec.interceptStartingValue},
                          {"irtScaleTranformationMethods", rec.irtScaleTranformationMethods},
                          {"slopeEstimate", rec.slopeEstimate},
                          {"interceptEstimate", rec.interceptEstimate},
                          {"transformedQuadratureNewForm", rec.transformedQuadratureNewForm},
                          {"itemResultsNewForm", rec.itemResultsNewForm},
                          {"optimization_results", rec.optimizationResults}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTInput& rec) {
      j = nlohmann::json {{"method", rec.method},
                          {"newFormFrequencyDistribution", rec.newFormFrequencyDistribution},
                          {"newItems", rec.newItems},
                          {"oldItems", rec.oldItems},
                          {"irtScaleTransformationData", rec.irtScaleTransformationData},
                          {"newFormIRTFittedDistribution", rec.newFormIRTFittedDistribution},
                          {"oldFormIRTFittedDistribution", rec.oldFormIRTFittedDistribution},
                          {"irtEquatingResults", rec.irtEquatingResults}};
    }

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::PData& rec) {
      j = nlohmann::json {{"design", rec.design},
                          {"method", rec.method},
                          {"smoothing", rec.smoothing},
                          {"methodCode", rec.methodCode},
                          {"minimumRawScoreYct", rec.minimumRawScoreYct},
                          {"maximumRawScoreYct", rec.maximumRawScoreYct},
                          {"scoreIncrementYct", rec.scoreIncrementYct},
                          {"lowestObservableRoundedScaledScore", rec.lowestObservableRoundedScaledScore},
                          {"highestObservableRoundedScaledScore", rec.highestObservableRoundedScaledScore},
                          {"weightSyntheticPopulation1", rec.weightSyntheticPopulation1},
                          {"isInternalAnchor", rec.isInternalAnchor},
                          {"reliabilityCommonItemsPopulation1", rec.reliabilityCommonItemsPopulation1},
                          {"reliabilityCommonItemsPopulation2", rec.reliabilityCommonItemsPopulation2},
                          {"methods", rec.methods},
                          {"minimumScoreX", rec.minimumScoreX},
                          {"maximumScoreX", rec.maximumScoreX},
                          {"scoreIncrementX", rec.scoreIncrementX},
                          {"scoreFrequenciesX", rec.scoreFrequenciesX},
                          {"numberOfExaminees", rec.numberOfExaminees},
                          {"numberOfBootstrapReplications", rec.numberOfBootstrapReplications},
                          {"bootstrapReplicationNumber", rec.bootstrapReplicationNumber},
                          {"roundToNumberOfDecimalPlaces", rec.roundToNumberOfDecimalPlaces}};

      if (rec.summaryRawDataX.has_value()) {
        j["summaryRawDataX"] = rec.summaryRawDataX.value();
      }
      if (rec.summaryRawDataY.has_value()) {
        j["summaryRawDataY"] = rec.summaryRawDataY.value();
      }
      if (rec.summaryRawDataXV.has_value()) {
        j["summaryRawDataXV"] = rec.summaryRawDataXV.value();
      }
      if (rec.summaryRawDataYV.has_value()) {
        j["summaryRawDataYV"] = rec.summaryRawDataYV.value();
      }
      if (rec.summaryRawDataXY.has_value()) {
        j["summaryRawDataXY"] = rec.summaryRawDataXY.value();
      }
      if (rec.rawToScaledScoreTable.has_value()) {
        j["rawToScaledScoreTable"] = rec.rawToScaledScoreTable.value();
      }
      if (rec.betaBinomalSmoothingX.has_value()) {
        j["betaBinomalSmoothingX"] = rec.betaBinomalSmoothingX.value();
      }
      if (rec.betaBinomalSmoothingY.has_value()) {
        j["betaBinomalSmoothingY"] = rec.betaBinomalSmoothingY.value();
      }
      if (rec.univariateLogLinearSmoothingX.has_value()) {
        j["univariateLogLinearSmoothingX"] = rec.univariateLogLinearSmoothingX.value();
      }
      if (rec.univariateLogLinearSmoothingY.has_value()) {
        j["univariateLogLinearSmoothingY"] = rec.univariateLogLinearSmoothingY.value();
      }
      if (rec.bivariateLogLinearSmoothingXV.has_value()) {
        j["bivariateLogLinearSmoothingXV"] = rec.bivariateLogLinearSmoothingXV.value();
      }
      if (rec.bivariateLogLinearSmoothingYV.has_value()) {
        j["bivariateLogLinearSmoothingYV"] = rec.bivariateLogLinearSmoothingYV.value();
      }
      if (rec.bivariateLogLinearSmoothingXY.has_value()) {
        j["bivariateLogLinearSmoothingXY"] = rec.bivariateLogLinearSmoothingXY.value();
      }
      if (rec.cubicSplinePostsmoothing.has_value()) {
        j["cubicSplinePostsmoothing"] = rec.cubicSplinePostsmoothing.value();
      }
      if (rec.irtInput.has_value()) {
        j["irtInput"] = rec.irtInput.value();
      }
    }
  } // namespace Structures
} // namespace EquatingRecipes

#endif