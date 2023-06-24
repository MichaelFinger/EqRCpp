/* 
  From Source: ERutilities.h 
  Original Struct: BOOT_ESS_RESULTS
  Description: equated scale-score results for bootstrap
*/

#ifndef IMPLEMENTATION_BOOTSTRAP_EQUATED_SCALED_SCORES_RESULTS_HPP
#define IMPLEMENTATION_BOOTSTRAP_EQUATED_SCALED_SCORES_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct BootstrapEquatedScaledScoresResults {
      Eigen::MatrixXd unroundedScaledScoresSumsAndMeans;           /* for matrix of sum and mn for unrounded scale scores */
      Eigen::MatrixXd unroundedScaledScoresSumSquaresAndSDs;       /* for matrix of sum2 and sd for unrounded scale scores */
      Eigen::VectorXd unroundedScaledScoresBoostrapStandardErrors; /* overall bootstrap se's for unrounded scale scales */
      Eigen::MatrixXd roundedScaledScoresSumsAndMeans;             /* for matrix of sum and mean for rounded scale scores */
      Eigen::MatrixXd roundedScaledScoresSumSquaresAndSDs;         /* for matrix of sum2 and sd for rounded scale scores */
      Eigen::VectorXd roundedScaledScoresBoostrapStandardErrors;   /* overall bootstrap se's for rounded scale scales */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif