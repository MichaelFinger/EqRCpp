#ifndef TESTS_BETA_BIMONIAL_HPP
#define TESTS_BETA_BIMONIAL_HPP

#include <iostream>

#include <Eigen/Core>
#include <fmt/core.h>

#include "fixtures/actmathfreq.hpp"
#include <equating_recipes/structures/all_structures.hpp>
#include <equating_recipes/beta_binomial.hpp>
#include <equating_recipes/utilities.hpp>

namespace Tests {
  struct BetaBinomial {
    void run() {
      Tests::Fixtures::ACTMathFreqData actMathFreqData = Tests::Fixtures::ACTMathFreq::getFreqDist();

      Eigen::VectorXd betaParEstsX = Tests::Fixtures::ACTMathFreq::get4ParameterBetaEstimatesX();

      EquatingRecipes::BetaBinomial betaBinomial;

      EquatingRecipes::Utilities scoreStatistics;
      EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX = scoreStatistics.univariateFromScoreFrequencies(actMathFreqData.freqX,
                                                                                                                               0,
                                                                                                                               40,
                                                                                                                               1,
                                                                                                                               "X");

      EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomialSmoothing = betaBinomial.betaBinomialSmoothing(univariateStatisticsX,
                                                                                                                    4,
                                                                                                                    0.0);

      std::cout << betaBinomialSmoothing.toLongString() << "\n";
    }
  };
} // namespace Tests
#endif