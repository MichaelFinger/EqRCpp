/* 
  From Source: ERutilities.h 
  Original Struct: BOOT_ESS_RESULTS
  Description: equated scale-score results for bootstrap
*/

#ifndef STRUCTURES_BOOTSTRAP_EQUATED_SCALED_SCORES_RESULTS_HPP
#define STRUCTURES_BOOTSTRAP_EQUATED_SCALED_SCORES_RESULTS_HPP

#include <string>
#include <Eigen/Core>
#include <fmt/core.h>

namespace EquatingRecipes {
  namespace Structures {
    struct BootstrapEquatedScaledScoresResults {
      Eigen::MatrixXd unroundedScaledScoresSumsAndMeans;           /* for matrix of sum and mn for unrounded scale scores */
      Eigen::MatrixXd unroundedScaledScoresSumSquaresAndSDs;       /* for matrix of sum2 and sd for unrounded scale scores */
      Eigen::MatrixXd unroundedScaledScoresBoostrapStandardErrors; /* overall bootstrap se's for unrounded scale scales */
      Eigen::MatrixXd roundedScaledScoresSumsAndMeans;             /* for matrix of sum and mean for rounded scale scores */
      Eigen::MatrixXd roundedScaledScoresSumSquaresAndSDs;         /* for matrix of sum2 and sd for rounded scale scores */
      Eigen::MatrixXd roundedScaledScoresBoostrapStandardErrors;   /* overall bootstrap se's for rounded scale scales */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif