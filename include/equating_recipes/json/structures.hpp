#ifndef JSON_STRUCTURES_HPP
#define JSON_STRUCTURES_HPP

#include <string>
#include <vector>
#include <Eigen/Core>

#include <nlohmann/json.hpp>

#include <equating_recipes/structures/all_structures.hpp>

namespace Eigen {
    template<typename Derived>
    void to_json(nlohmann::json& j, const MatrixBase<Derived>& matrix) {
      j = const_cast<Eigen::MatrixBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>&>(matrix);
    }
}

namespace EquatingRecipes {
  namespace Structures {
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

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::CGEquipercentileEquatingResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::CubicSplinePostsmoothing& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Design& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::EquatedRawScoreBootstrapResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::EquatedScaledScoresResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::FormType& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTEquatingResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTFittedDistribution& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTMethod& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTModel& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTScaleTransformationMethod& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTTrueScoreEquating& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::LossSpec& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Method& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Moments& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Quadrature& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Smoothing& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::Symmetry& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::RawToScaledScoreTable& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BetaBinomialSmoothing& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::BivariateStatistics& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::CommonItemSpecification& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::ItemSpecification& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTScaleTransformationItemResults& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTScaleTransformationData& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::IRTInput& rec) {}

    void to_json(nlohmann::json& j, const EquatingRecipes::Structures::PData& rec) {}
  } // namespace Structures
} // namespace EquatingRecipes

#endif