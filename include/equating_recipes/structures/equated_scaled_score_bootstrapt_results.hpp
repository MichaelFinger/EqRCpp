/* 
  From Source: ERutilities.h 
  Original Struct: BOOT_ESS_RESULTS
  Description: equated scaled score results for bootstrap
*/

#ifndef IMPLEMENTATION_EQUATED_SCALED_SCORE_BOOTSTRAP_RESULTS_HPP
#define IMPLEMENTATION_EQUATED_SCALED_SCORE_BOOTSTRAP_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedScaledScoreBootstrapResults {
      Eigen::MatrixXd sumAndMeanUnroundedScaledScores;                    /* for matrix of sum and mn for unrounded scale scores */
      Eigen::MatrixXd sumSquaredAndSDUnroundedScaledScores;              /* for matrix of sum2 and sd for unrounded scale scores */
      Eigen::MatrixXd unroundedScaledScoresBootstrapStandardErrors;         /* overall bootstrap se's for unrounded scale scales */ 
      Eigen::MatrixXd sumAndMeanRoundedScaledScores;                      /* for matrix of sum and mean for rounded scale scores */
      Eigen::MatrixXd sumSquaredAndSDRoundedScaledScores;                  /* for matrix of sum2 and sd for rounded scale scores */
      Eigen::MatrixXd roundedScaledScoresBootstrapStandardErrors;             /* overall bootstrap se's for rounded scale scales */
    };
  }
}

#endif