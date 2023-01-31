#ifndef LSAT6_RAW_SCORE_FREQ_DIST_HPP
#define LSAT6_RAW_SCORE_FREQ_DIST_HPP

#include <Eigen/Core>

namespace Tests {
  namespace Fixtures {
    struct LSAT6 {
      static Eigen::VectorXd rawScoreFrequencyDistribution() {
        Eigen::VectorXd freqDist(6);

        freqDist(0) = 3;
        freqDist(1) = 20;
        freqDist(2) = 85;
        freqDist(3) = 237;
        freqDist(4) = 357;
        freqDist(5) = 298;

        return freqDist;
      }

      static Eigen::VectorXd rawScoreRelativeFrequencyDistribution() {
        Eigen::VectorXd freqDist = LSAT6::rawScoreFrequencyDistribution();

        Eigen::VectorXd relativeFreqDist = freqDist.cast<double>() / static_cast<double>(freqDist.sum());

        return relativeFreqDist;
      }
    };
  } // namespace Fixtures
} // namespace Tests

#endif