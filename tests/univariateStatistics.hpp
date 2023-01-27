#ifndef TESTS_UNIVARIATE_STATISTICS_HPP
#define TESTS_UNIVARIATE_STATISTICS_HPP

#include <iostream>
#include <map>
#include "fixtures/lsat6.hpp"
#include <equating_recipes/structures/univariate_statistics.hpp>

namespace Tests {
  struct UnivariateStatistics {
    void run() {
      std::map<double, int> lsat6FreqDist = LSAT6::rawScoreFrequencyDistribution();

      std::for_each(lsat6FreqDist.begin(),
                    lsat6FreqDist.end(),
                    [&](const std::pair<double, int>& entry) {
                      std::cout << entry.first << ", " << entry.second << "\n";
                    });

      EquatingRecipes::Structures::UnivariateStatistics univariateStatistics =
          EquatingRecipes::Structures::UnivariateStatistics::create(LSAT6::rawScoreFrequencyDistribution(),
                                                                    0,
                                                                    5,
                                                                    1,
                                                                    "X");

      std::cout << univariateStatistics.toString() << "\n";
    }
  };
} // namespace Tests
#endif