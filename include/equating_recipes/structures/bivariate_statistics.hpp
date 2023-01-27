/* 
  From Source: ERutilities.h 
  Original Struct: BSTATS
  Description: raw-score statistics for a bivariate distribution 

  From Source: ERutilities.h, ERutilities.c
  Original Method: ReadRawGet_BSTATS
  Description: Reads raw data and get or assign all elements for bivariate struct s
*/

#ifndef STRUCTURES_BIVARIATE_STATISTICS_HPP
#define STRUCTURES_BIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>

#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct BivariateStatistics {
      UnivariateStatistics rowScoreStatistics;
      UnivariateStatistics columnScoreStatistics;

      int numberOfExaminees;
      Eigen::MatrixXi bivariateFreqDist;
      Eigen::MatrixXd bivariateFreqDistDouble;
      double covariance;
      double correlation;
      Eigen::MatrixXd bivariateProportions;

      static BivariateStatistics create(const Eigen::MatrixXd& scores,
                                        const double& minimumRowScore,
                                        const double& maximumRowScore,
                                        const double& rowScoreIncrement,
                                        const double& minimumColumnScore,
                                        const double& maximumColumnScore,
                                        const double& columnScoreIncrement,
                                        const std::string& rowScoreId,
                                        const std::string& columnScoreId);
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif