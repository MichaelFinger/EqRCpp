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

      Eigen::VectorXd obsScoreDensityX;

      betaBinomial.observedDensity(actMathFreqData.rawScores.size(),
                                   actMathFreqData.freqX.sum(),
                                   betaParEstsX,
                                   obsScoreDensityX);

      std::cout << "Observed Score Density for X: " << EquatingRecipes::Utilities::vectorXdToString(obsScoreDensityX, false) << "\n";
    }
  };
} // namespace Tests
#endif