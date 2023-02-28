// #include "univariate_statistics.hpp"
// #include "bivariate_statistics.hpp"
// #include "moments.hpp"
// #include "equated_scaled_scores.hpp"
// #include "matrix_decomp.hpp"
// #include "beta_binomial.hpp"
// #include "rg_and_sg_equating.hpp"

#include "solving.hpp"

#include <equating_recipes/analytic_standard_errors.hpp>
#include <equating_recipes/beta_binomial.hpp>
#include <equating_recipes/bootstrap.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/cg_no_smoothing.hpp>
#include <equating_recipes/cubic_spline.hpp>
#include <equating_recipes/equated_scaled_scores.hpp>
#include <equating_recipes/irt_equating.hpp>
#include <equating_recipes/irt_scale_transformation.hpp>
#include <equating_recipes/log_linear_equating.hpp>
#include <equating_recipes/results_document.hpp>
#include <equating_recipes/rg_and_sg_equating.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/haebara_function.hpp>
#include <equating_recipes/stocking_lord_function.hpp>
#include <equating_recipes/continuized_log_linear_equating.hpp>
#include <equating_recipes/kernel_equating.hpp>

#include <iostream>
#include <Eigen/Dense>

int main(int argc, char const *argv[]) {
  // Tests::UnivariateStatistics univariateStatisticsTest;  
  // univariateStatisticsTest.run();

  // Tests::BivariateStatistics bivariateStatisticsTest;
  // bivariateStatisticsTest.run();

  // Tests::Moments moments;
  // moments.run();

  // Tests::EquatedScaledScores equatedScaledScores;
  // equatedScaledScores.run();

  // Tests::MatrixDecomp matrixDecomp;
  // matrixDecomp.run();

  // Tests::BetaBinomial betaBinomial;
  // betaBinomial.run();

  // Tests::RGAndSG rgAndSG;
  // rgAndSG.run();

  Tests::Solving solving;
  solving.run();

  return 0;
}
