/* 
  From Source: ERutilities.h 
  Original Struct: BOOT_ERAW_RESULTS
  Description: equated raw-score results for bootstrap
*/

#ifndef IMPLEMENTATION_BOOTSTRAP_EQUATED_RAW_SCORE_RESULTS_HPP
#define IMPLEMENTATION_BOOTSTRAP_EQUATED_RAW_SCORE_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct BootstrapEquatedRawScoreResults {
      Eigen::MatrixXd sumAndMeanScores;                  /* for matrix of sum and mean for scores */
      Eigen::MatrixXd sumSquareAndSDScores;               /* for matrix of sum2 and sd for scores */
      Eigen::MatrixXd bootstrapStandardErrors;                          /* overall bootstrap se's */
    };
  }
}

#endif