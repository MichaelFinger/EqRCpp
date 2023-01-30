#include "univariate_statistics.hpp"
#include "bivariate_statistics.hpp"
#include "moments.hpp"

int main(int argc, char const *argv[])
{
  // Tests::UnivariateStatistics univariateStatisticsTest;  
  // univariateStatisticsTest.run();

  // Tests::BivariateStatistics bivariateStatisticsTest;
  // bivariateStatisticsTest.run();

  Tests::Moments moments;
  moments.run();

  return 0;
}
