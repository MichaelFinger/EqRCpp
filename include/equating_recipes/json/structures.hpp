#ifndef JSON_STRUCTURES_HPP
#define JSON_STRUCTURES_HPP

#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

#include <nlohmann/json.hpp>

#include <equating_recipes/structures/all_structures.hpp>
#include <equating_recipes/log_linear_equating.hpp>

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
                                                                                             {EquatingRecipes::Structures::IRTScaleTransformationMethod::STOCKING_LORD, "Stocking Lord"}})

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

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::Smoothing, {{EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED, ""},
                                                                          {EquatingRecipes::Structures::Smoothing::LOG_LINEAR, "log linear"},
                                                                          {EquatingRecipes::Structures::Smoothing::BETA_BINOMIAL, "beta binomial"},
                                                                          {EquatingRecipes::Structures::Smoothing::CUBIC_SPLINE, "cubic spline"},
                                                                          {EquatingRecipes::Structures::Smoothing::KERNAL, "kernal"},
                                                                          {EquatingRecipes::Structures::Smoothing::CONTINUIZED_LOG_LINEAR_EQUATING, "continuized log linear equating"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::Structures::Symmetry, {{EquatingRecipes::Structures::Symmetry::NOT_SPECIFIED, ""},
                                                                         {EquatingRecipes::Structures::Symmetry::NEW_SCALE, "new scale"},
                                                                         {EquatingRecipes::Structures::Symmetry::OLD_SCALE, "old scale"},
                                                                         {EquatingRecipes::Structures::Symmetry::SYMMETRIC, "symmetric"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::LogLinearEquating::CriterionComparisonType, {{EquatingRecipes::LogLinearEquating::CriterionComparisonType::ABSOLUTE, "absolute"},
                                                                                               {EquatingRecipes::LogLinearEquating::CriterionComparisonType::RELATIVE, "relative"}})

    NLOHMANN_JSON_SERIALIZE_ENUM(EquatingRecipes::LogLinearEquating::DesignMatrixType, {{EquatingRecipes::LogLinearEquating::DesignMatrixType::SOLUTION, "solution"},
    {EquatingRecipes::LogLinearEquating::DesignMatrixType::RAW_SCORE, "raw score"}})

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::UnivariateStatistics& rec) {
      j = nlohmann::json {{"id", rec.id},
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
                          {"moments", rec.momentValues}};
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
                          {"useRelativeCriterionComparison", rec.useRelativeCriterionComparison},
                          {"useRawAndCentralMoments", rec.useRawAndCentralMoments},
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
                          {"mininumRawScore", rec.mininumRawScore},
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
      j = nlohmann::json {{"univariateStatisticsRow", rec.univariateStatisticsRow},
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

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTScaleTransformationData& rec) {
      j = nlohmann::json {{"minimumRawScoreNewForm", rec.minimumRawScoreNewForm},
                          {"maximumRawScoreNewForm", rec.maximumRawScoreNewForm},
                          {"rawScoreIncrementNewForm", rec.rawScoreIncrementNewForm},
                          {"minimumRawScoreOldForm", rec.minimumRawScoreOldForm},
                          {"maximumRawScoreOldForm", rec.maximumRawScoreOldForm},
                          {"rawScoreIncrementOldForm", rec.rawScoreIncrementOldForm},
                          {"quadratureNewForm", rec.quadratureNewForm},
                          {"quadratureOldForm", rec.quadratureOldForm},
                          {"newItems", rec.newItems},
                          {"oldItems", rec.oldItems},
                          {"commonItems", rec.commonItems},
                          {"runHaebara", rec.runHaebara},
                          {"haebaraSymmetryOption", rec.haebaraSymmetryOption},
                          {"haebaraFunctionStandardization", rec.haebaraFunctionStandardization},
                          {"runStockingLord", rec.runStockingLord},
                          {"stockingLordSymmetryOption", rec.stockingLordSymmetryOption},
                          {"stockingLordFunctionStandardization", rec.stockingLordFunctionStandardization},
                          {"meanMeanSlope", rec.meanMeanSlope},
                          {"meanMeanIntercept", rec.meanMeanIntercept},
                          {"meanSigmaSlope", rec.meanSigmaSlope},
                          {"meanSigmaIntercept", rec.meanSigmaIntercept},
                          {"transformedQuadratureNewForm", rec.transformedQuadratureNewForm},
                          {"itemResultsNewForm", rec.itemResultsNewForm},
                          {"irtScaleTranformationMethod", rec.irtScaleTranformationMethod}};

      if (rec.haebaraSlopeStartingValue.has_value()) {
        j["haebaraSlopeStartingValue"] = rec.haebaraSlopeStartingValue.value();
      }

      if (rec.haebaraInterceptStartingValue.has_value()) {
        j["haebaraInterceptStartingValue"] = rec.haebaraInterceptStartingValue.value();
      }

      if (rec.stockingLordSlopeStartingValue.has_value()) {
        j["stockingLordSlopeStartingValue"] = rec.stockingLordSlopeStartingValue.value();
      }

      if (rec.stockingLordInterceptStartingValue.has_value()) {
        j["stockingLordInterceptStartingValue"] = rec.stockingLordInterceptStartingValue.value();
      }

      if (rec.haebaraSlope.has_value()) {
        j["haebaraSlope"] = rec.haebaraSlope.value();
      }

      if (rec.haebaraIntercept.has_value()) {
        j["haebaraIntercept"] = rec.haebaraIntercept.value();
      }

      if (rec.stockingLordSlope.has_value()) {
        j["stockingLordSlope"] = rec.stockingLordSlope.value();
      }

      if (rec.stockingLordIntercept.has_value()) {
        j["stockingLordIntercept"] = rec.stockingLordIntercept.value();
      }
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
                          {"mininumScoreX", rec.mininumScoreX},
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