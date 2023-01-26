/* 
  From Source: ERutilities.h 
  Original Struct: ESS_RESULTS
  Description: equated scaled sore results
*/

#ifndef STRUCTURES_EQUATED_SCALED_SCORE_RESULTS_HPP
#define STRUCTURES_EQUATED_SCALED_SCORE_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedScaledScoreResults {
      Eigen::MatrixXd unroundedEquatedScaledScores;                       /* unrounded equated scale scores */
      Eigen::MatrixXd roundedEquatedScaledScores;                           /* rounded equated scale scores */
      Eigen::MatrixXd unroundedEquatedScaledScoreMoments;     /* moments for equated unrounded scale scores */
      Eigen::MatrixXd roundedEquatedScaledScoreMoments;         /* moments for equated rounded scale scores */
    };
  }
}

#endif