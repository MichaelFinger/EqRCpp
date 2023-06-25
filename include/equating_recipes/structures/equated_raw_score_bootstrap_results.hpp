/* 
  From Source: ERutilities.h 
  Original Struct: BOOT_ERAW_RESULTS
  Description: equated raw score results for bootstrap
*/

#ifndef STRUCTURES_EQUATED_RAW_SCORE_BOOTSTRAP_RESULTS_HPP
#define STRUCTURES_EQUATED_RAW_SCORE_BOOTSTRAP_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedRawScoreBootstrapResults {
      Eigen::MatrixXd sumAndMeanRawScores;                    /* for matrix of sum and mean for scores */
      Eigen::MatrixXd sumSquaredAndSDRawScores;                /* for matrix of sum2 and sd for scores */
      Eigen::MatrixXd rawScoresBootstrapStandardErrors;                      /* overall bootstrap se's */
    };
  }
}

#endif