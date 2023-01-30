#ifndef TESTS_MOMENTS_HPP
#define TESTS_MOMENTS_HPP

#include <iostream>
#include <Eigen/Core>

#include "fixtures/lsat6.hpp"
#include <equating_recipes/structures/moments.hpp>

namespace Tests {
  struct Moments {
    void run() {
      Eigen::VectorXi lsat6FreqDist = Tests::Fixtures::LSAT6::rawScoreFrequencyDistribution();

      Eigen::VectorXd lsat6RelativeFreqDist = Tests::Fixtures::LSAT6::rawScoreRelativeFrequencyDistribution();

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Structures::Moments::fromScoreFrequencies(lsat6FreqDist,
                                                                                                                0.0,
                                                                                                                5.0,
                                                                                                                1.0);

      std::cout << moments.toString() << "\n";
    }
  };
} // namespace Tests

#endif