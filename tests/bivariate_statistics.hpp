#ifndef TESTS_BIVARIATE_STATISTICS_HPP
#define TESTS_BIVARIATE_STATISTICS_HPP

#include <iostream>

#include <Eigen/Core>
#include <fmt/core.h>

#include "fixtures/mondatx.hpp"
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>

namespace Tests {
  struct BivariateStatistics {
    void run() {
      Eigen::MatrixXd jointRawScores = Tests::Fixtures::MondatX::jointRawScores();

      EquatingRecipes::Structures::BivariateStatistics bivariateStatistics =
         EquatingRecipes::Structures::BivariateStatistics::create(jointRawScores,
         0,
         36,
         1,
         0,
         12,
         1,
         "X",
         "Y");

      std::cout << bivariateStatistics.toString() << "\n";
    }
  };
} // namespace Tests
#endif