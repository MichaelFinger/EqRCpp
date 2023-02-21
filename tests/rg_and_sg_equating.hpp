#ifndef TESTS_RG_AND_SG_HPP
#define TESTS_RG_AND_SG_HPP

#include <iostream>

#include <Eigen/Core>
#include <fmt/core.h>

#include "fixtures/actmathfreq.hpp"
#include <equating_recipes/structures/all_structures.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/rg_and_sg_equating.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>

namespace Tests {
  struct RGAndSG {
    void run() {
      Tests::Fixtures::ACTMathFreqData actMathFreqData = Tests::Fixtures::ACTMathFreq::getFreqDist();
      double maximumScore = actMathFreqData.rawScores.maxCoeff();

      // std::cout << EquatingRecipes::Utilities::vectorXdToString(actMathFreqData.rawScores, false) << "\n"
      //           << EquatingRecipes::Utilities::vectorXdToString(actMathFreqData.freqX, false) << "\n"
      //           << EquatingRecipes::Utilities::vectorXdToString(actMathFreqData.freqY, false) << "\n";

      EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsX = EquatingRecipes::Utilities::univariateFromScoreFrequencies(actMathFreqData.freqX,
                                                                                                                                                 0,
                                                                                                                                                 maximumScore,
                                                                                                                                                 1,
                                                                                                                                                 "X");

      EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsY = EquatingRecipes::Utilities::univariateFromScoreFrequencies(actMathFreqData.freqY,
                                                                                                                                                 0,
                                                                                                                                                 maximumScore,
                                                                                                                                                 1,
                                                                                                                                                 "Y");

      // // actMathFreqData

      EquatingRecipes::Structures::PData pData;

      EquatingRecipes::RandomAndSingleGroupEquating rgSgEquating;
      EquatingRecipes::Structures::EquatedRawScoreResults results;
      rgSgEquating.randomGroupEquating(EquatingRecipes::Structures::Design::RandomGroups,
                                       EquatingRecipes::Structures::Method::LINEAR,
                                       EquatingRecipes::Structures::Smoothing::NO,
                                       univariateStatisticsX,
                                       univariateStatisticsY,
                                       0,
                                       pData,
                                       results);

      std::cout << results.toString() << "\n";
    }
  };
} // namespace Tests

#endif