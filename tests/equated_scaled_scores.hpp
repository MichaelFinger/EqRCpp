#ifndef TESTS_EQUATED_SCALED_SCORES_HPP
#define TESTS_EQUATED_SCALED_SCORES_HPP

#include <iostream>
#include <Eigen/Core>

#include "fixtures/lsat6.hpp"

#include <equating_recipes/wrappers/equated_scaled_scores.hpp>
#include <equating_recipes/structures/equated_scaled_scores_input.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/utilities.hpp>

namespace Tests {
  struct EquatedScaledScores {
    void run() {
      EquatingRecipes::Structures::EquatedScaledScoresInput input;

      EquatingRecipes::Wrappers::EquatedScaledScores equatedScaledScoresWrapper;
      EquatingRecipes::Structures::EquatedScaledScoresResults results = equatedScaledScoresWrapper.run(input);
      
    }
  };
}

#endif