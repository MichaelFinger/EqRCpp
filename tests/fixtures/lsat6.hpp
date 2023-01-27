#ifndef LSAT6_RAW_SCORE_FREQ_DIST_HPP
#define LSAT6_RAW_SCORE_FREQ_DIST_HPP

#include <map>

struct LSAT6 {
  static std::map<double, int> rawScoreFrequencyDistribution() {
    std::map<double, int> freqDist {
       { 0, 3 },
       { 1, 20 },
       { 2, 85 },
       { 3, 237 },
       { 4, 357 },
       { 5, 298 }
    };

    return freqDist;
  }
};

#endif