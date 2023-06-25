/* 
  From Source: ERutilities.h 
  Original Struct: BSTATS
  Description: raw-score statistics for a bivariate distribution 
*/

#ifndef STRUCTURES_BIVARIATE_STATISTICS_HPP
#define STRUCTURES_BIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <fmt/core.h>

#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct BivariateStatistics {
      std::string datasetName;
      std::string rowVariableName;
      std::string columnVariableName;
      EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsRow;
      EquatingRecipes::Structures::UnivariateStatistics univariateStatisticsColumn;
      int numberOfExaminees;
      Eigen::MatrixXi bivariateFreqDist;
      Eigen::MatrixXd bivariateFreqDistDouble;
      double covariance;
      double correlation;
      Eigen::MatrixXd bivariateProportions;

      // std::string toString() {
      //   std::string value = "Univariate Statistics: Row Scores\n";
      //   value.append(this->univariateStatisticsRow.toString());
      //   value.append("\n");

      //   value.append("Univariate Statistics: Column Scores\n");
      //   value.append(this->univariateStatisticsColumn.toString());
      //   value.append("\n");

      //   value.append(fmt::format("Number of Examinees: {}\n", this->numberOfExaminees));

      //   value.append(fmt::format("Bivariate Frequency Distribution:\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXiToString(this->bivariateFreqDist)));
      //   value.append(fmt::format("Bivariate Frequency Distribution (double):\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXdToString(this->bivariateFreqDistDouble)));
      //   value.append(fmt::format("Bivariate Proportions:\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXdToString(this->bivariateProportions)));
      //   value.append(fmt::format("Sum of Bivariate Proportions: {}\n", this->bivariateProportions.sum()));
      //   value.append(fmt::format("Covariance: {}\n", this->covariance));
      //   value.append(fmt::format("Correlation: {}\n", this->correlation));

      //   return value;
      // }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif