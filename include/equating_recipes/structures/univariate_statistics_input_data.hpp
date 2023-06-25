#ifndef STRUCTURES_UNIVARIATE_STATISTICS_INPUT_DATA_HPP
#define STRUCTURES_UNIVARIATE_STATISTICS_INPUT_DATA_HPP

#include <string>
#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct UnivariateStatisticsInputData {
      std::string variableName;
      Eigen::VectorXd scoreFrequencies;
      double minimumScore;
      double maximumScore;
      double scoreIncrement;
      std::string id;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif