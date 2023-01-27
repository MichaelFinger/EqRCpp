#include "univariate_statistics.hpp"
#include "bivariate_statistics.hpp"

int main(int argc, char const *argv[])
{
  // Tests::UnivariateStatistics univariateStatisticsTest;  
  // univariateStatisticsTest.run();

  Tests::BivariateStatistics bivariateStatisticsTest;
  bivariateStatisticsTest.run();

  return 0;
}
