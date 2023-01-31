#include "univariate_statistics.hpp"
#include "bivariate_statistics.hpp"
#include "moments.hpp"
#include "equated_scaled_scores.hpp"

int main(int argc, char const *argv[])
{
  // Tests::UnivariateStatistics univariateStatisticsTest;  
  // univariateStatisticsTest.run();

  // Tests::BivariateStatistics bivariateStatisticsTest;
  // bivariateStatisticsTest.run();

  // Tests::Moments moments;
  // moments.run();

  Tests::EquatedScaledScores equatedScaledScores;
  equatedScaledScores.run();

  return 0;
}
