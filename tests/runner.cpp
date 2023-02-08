// #include "univariate_statistics.hpp"
// #include "bivariate_statistics.hpp"
// #include "moments.hpp"
// #include "equated_scaled_scores.hpp"
// #include "matrix_decomp.hpp"
// #include "beta_binomial.hpp"
#include "rg_and_sg.hpp"

int main(int argc, char const *argv[])
{
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

  Tests::RGAndSG rgAndSG;
  rgAndSG.run();

  return 0;
}
