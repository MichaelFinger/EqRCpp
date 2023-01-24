/* 
  From Source: ERutilities.h 
  Description: raw-score statistics for a bivariate distribution 
*/

#ifndef STRUCTURES_BIVARIATE_STATISTICS_HPP
#define STRUCTURES_BIVARIATE_STATISTICS_HPP

#include <Eigen/Core>

#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct BivariateStatistics {
      // char fname[100]; // name of input file 
      // int n;           // number of examinees 
      
      // var 1 -- rows of bfd[][] 
      UnivariateStatistics rowMarginalStatistics;
      
      // var 2 -- columns of bfd[][] 
      UnivariateStatistics columnMarginalStatistics;

      Eigen::MatrixXi bivariateFreqDist;        // bivariate freq dist bfd[0][0]...bfd[ns1-1][ns2-1]       
      Eigen::MatrixXd bivariateFreqDistDouble;  // double version of bfd[][] 
      double covariance;                        // covariance   
      double correlation;                       // correlation   
      double bivariateProportions;              // biv proportions[var1][var2]   
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif