#include <cmath>
#include <limits>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  namespace Structures {
    Moments Moments::getScoreMoments(const Eigen::VectorXd& scores) {
      Moments scoreMoments;

      scoreMoments.numberOfExaminees = scores.size();

      scoreMoments.minimumObservedScore = scores.minCoeff();
      scoreMoments.maximumObservedScore = scores.maxCoeff();

      scoreMoments.momentValues.setZero(4);
      scoreMoments.momentValues(0) = scores.mean();

      Eigen::VectorXd meanVector = Eigen::VectorXd::Constant(scores.size(), scoreMoments.momentValues(0));
      Eigen::VectorXd deviations = scores - meanVector;

      double variance = deviations.array().pow(2).mean();
      double skewness = deviations.array().pow(3).mean();
      double kurtosis = deviations.array().pow(4).mean();

      scoreMoments.momentValues(1) = std::sqrt(variance);
      scoreMoments.momentValues(2) = skewness / std::pow(scoreMoments.momentValues(1), 3);
      scoreMoments.momentValues(3) = kurtosis / std::pow(scoreMoments.momentValues(1), 4);

      return scoreMoments;
    }

    Moments Moments::getScoreMoments(const Eigen::VectorXi& scoreFrequencies,
                                     const double& minimumScore,
                                     const double& maximumScore,
                                     const double& scoreIncrement) {
      Moments scoreMoments;

      size_t numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore,
                                                                         maximumScore,
                                                                         scoreIncrement);

      size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(maximumScore,
                                                                                 minimumScore,
                                                                                 scoreIncrement);

      scoreMoments.minimumObservedScore = std::numeric_limits<double>::max();
      scoreMoments.maximumObservedScore = std::numeric_limits<double>::min();
      scoreMoments.numberOfExaminees = scoreFrequencies.sum();
      scoreMoments.momentValues.setZero(4);

      Eigen::VectorXd deviations(numberOfScores);

      for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
        double score = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                            minimumScore,
                                                            scoreIncrement);

        int scoreFreq = scoreFrequencies(scoreLocation);

        deviations(scoreLocation) = score;

        if (scoreFreq >= 1) {
          scoreMoments.minimumObservedScore = std::min(score, scoreMoments.minimumObservedScore);
          scoreMoments.maximumObservedScore = std::max(score, scoreMoments.maximumObservedScore);
        }

        scoreMoments.momentValues(0) += score * static_cast<double>(scoreFreq);
      }

      scoreMoments.momentValues(0) /= static_cast<double>(scoreMoments.numberOfExaminees);

      deviations = deviations - Eigen::VectorXd::Constant(deviations.size(), scoreMoments.momentValues(0));

      for (size_t powCoeff = 2; powCoeff <= 4; powCoeff++) {
        scoreMoments.momentValues(powCoeff - 1) = (deviations.array().pow(powCoeff) * scoreFrequencies.array().cast<double>()).sum() /
                        static_cast<double>(scoreMoments.numberOfExaminees);

      }
      
      scoreMoments.momentValues(1) = std::sqrt(scoreMoments.momentValues(1));
      scoreMoments.momentValues(2) /= std::pow(scoreMoments.momentValues(1), 3);
      scoreMoments.momentValues(3) /= std::pow(scoreMoments.momentValues(1), 4);

      return scoreMoments;
    }
  } // namespace Structures
} // namespace EquatingRecipes