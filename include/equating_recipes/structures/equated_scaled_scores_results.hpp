/* 
  From Source: ERutilities.h 
  Original Struct: ESS_RESULTS
  Description: equated scaled score results
*/

#ifndef IMPLEMENTATION_EQUATED_SCALED_SCORES_RESULTS_HPP
#define IMPLEMENTATION_EQUATED_SCALED_SCORES_RESULTS_HPP

#include <string>
#include <Eigen/Core>
#include <fmt/core.h>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedScaledScoresResults {
      Eigen::MatrixXd unroundedEquatedScaledScores;                       /* unrounded equated scale scores */
      Eigen::MatrixXd roundedEquatedScaledScores;                           /* rounded equated scale scores */
      Eigen::MatrixXd unroundedEquatedScaledScoreMoments;     /* moments for equated unrounded scale scores */
      Eigen::MatrixXd roundedEquatedScaledScoreMoments;         /* moments for equated rounded scale scores */

      void configure(size_t numberOfScores, size_t numberOfMethods) {
        unroundedEquatedScaledScoreMoments.setZero(numberOfMethods, numberOfScores);
        roundedEquatedScaledScoreMoments.setZero(numberOfMethods, numberOfScores);
        unroundedEquatedScaledScoreMoments.setZero(numberOfMethods, 4);
        roundedEquatedScaledScoreMoments.setZero(numberOfMethods, 4);
      }
    };
  }
}

#endif