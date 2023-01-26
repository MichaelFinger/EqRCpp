#include <cmath>
#include <equating_recipes/structures/moments.hpp>

namespace EquatingRecipes {
  namespace Structures {
    Moments Moments::getScoreMoments(const Eigen::VectorXd& scores) {
      Moments scoreMoments;

      scoreMoments.numberOfPersons = scores.size();

      scoreMoments.minimumObservedScore = scores.minCoeff();
      scoreMoments.maximumObservedScore = scores.maxCoeff();

      scoreMoments.moments.setZero(4);
      scoreMoments.moments(0) = scores.mean();

      Eigen::VectorXd meanVector = Eigen::VectorXd::Constant(scores.size(), scoreMoments.moments(0));
      Eigen::VectorXd deviations = scores - meanVector;

      double variance = deviations.array().pow(2).mean();
      double skewness = deviations.array().pow(3).mean();
      double kurtosis = deviations.array().pow(4).mean();

      scoreMoments.moments(1) = std::sqrt(variance);
      scoreMoments.moments(2) = skewness / std::pow(scoreMoments.moments(1), 3);
      scoreMoments.moments(3) = kurtosis / std::pow(scoreMoments.moments(1), 3);

      return scoreMoments;
    }

    Moments Moments::getScoreMoments(const std::map<double, int>& scoreFreqDist) {
      Moments scoreMoments;

      scoreMoments.numberOfPersons = 0;
      scoreMoments.minimumObservedScore = (*(scoreFreqDist.begin())).first;
      scoreMoments.maximumObservedScore = (*(scoreFreqDist.end())).first;
      scoreMoments.moments.setZero(4);

      std::for_each(scoreFreqDist.begin(),
                    scoreFreqDist.end(),
                    [&](const std::pair<double, int>& entry) {
                      double scoreValue = entry.first;
                      int scoreFreq = entry.second; 

                      scoreMoments.numberOfPersons += scoreFreq;
                      scoreMoments.moments(0) += scoreValue * static_cast<double>(scoreFreq);
                    });

      scoreMoments.moments(0) /= static_cast<double>(scoreMoments.numberOfPersons);

      Eigen::VectorXd deviations(scoreFreqDist.size());
      size_t index = 0;
      std::for_each(scoreFreqDist.begin(),
                    scoreFreqDist.end(),
                    [&](const std::pair<double, int>& entry) {
                      double scoreValue = entry.first;
                      int scoreFreq = entry.second; 

                      deviations(index) = (scoreValue - scoreMoments.moments(0)) * static_cast<double>(scoreFreq);

                      ++index;
                    });

      double variance = deviations.array().pow(2).mean();
      double skewness = deviations.array().pow(3).mean();
      double kurtosis = deviations.array().pow(4).mean();

      scoreMoments.moments(1) = std::sqrt(variance);
      scoreMoments.moments(2) = skewness / std::pow(scoreMoments.moments(1), 3);
      scoreMoments.moments(3) = kurtosis / std::pow(scoreMoments.moments(1), 3);

      return scoreMoments;
    }
  } // namespace Structures
} // namespace EquatingRecipes