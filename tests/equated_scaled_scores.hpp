#ifndef TESTS_EQUATED_SCALED_SCORES_HPP
#define TESTS_EQUATED_SCALED_SCORES_HPP

#include <iostream>
#include <Eigen/Core>

#include "fixtures/lsat6.hpp"

#include <equating_recipes/equated_scaled_scores.hpp>
#include <equating_recipes/structures/all_structures.hpp>
#include <equating_recipes/utilities.hpp>

namespace Tests {
  struct EquatedScaledScores {
    void run() {
      EquatingRecipes::EquatedScaledScores equatedScaledScoresWrapper;
      // EquatingRecipes::Structures::EquatedScaledScoresResults results = equatedScaledScoresWrapper.run();      
    }
  };
}

#endif