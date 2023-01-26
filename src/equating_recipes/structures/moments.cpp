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

      Eigen::VectorXd deviations = scores - Eigen::VectorXd::Constant(scoreMoments.moments(0));
      
      double variance = deviations.pow(2).mean();
      double skewness = deviations.pow(3).mean();
      double kurtosis = deviations.pow(4).mean();

      scoreMoments.moments(1) = std::sqrt(variance);
      scoreMoments.moments(2) = skewness / std::pow(scoreMoments.moments(1), 3);
      scoreMoments.moments(3) = kurtosis / std::pow(scoreMoments.moments(1), 3);

      return scoreMoments;
    }
  }
}