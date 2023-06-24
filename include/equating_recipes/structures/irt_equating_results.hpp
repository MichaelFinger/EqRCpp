/* 
  From Source: IRTeq.h 
  Original Struct: RawTruObsEquiv
  Description: Structure that stores IRT true-score and observed-score equating results
*/

#ifndef IMPLEMENTATION_IRT_EQUATING_RESULTS_HPP
#define IMPLEMENTATION_IRT_EQUATING_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTEquatingResults {
      size_t numberOfRawScoreCategoriesNewForm;                                       /* Number of new form raw score categories */
      double minimumTrueScoreNewForm;                                          /* lower limit of true test score on the new form */
      double minimumTrueScoreOldForm;                                          /* lower limit of true test score on the old form */
      Eigen::VectorXd thetaEquivalentFormXScore;                                             /* Theta-equivalent of Form X score */
      Eigen::VectorXd unroundedEquatedTrueScore;                  /* Unrounded raw-to-raw conversion for IRT true score equating */
      Eigen::VectorXd roundedEquatedTrueScore;                      /* Rounded raw-to-raw conversion for IRT true score equating */
      Eigen::VectorXd unroundedEquatedObservedScore;              /*  Unrounded raw-to-raw conversion for IRT obs score equating */
      Eigen::VectorXd roundedEquatedObservedScore;                   /* Rounded raw-to-raw conversion for IRT obs score equating */
      Eigen::VectorXd momentsEquatedTrueScores;                       /* new form equated (true score) score (unrounded) moments */
      Eigen::VectorXd momentsEquatedObservedScores;                   /* new form equated (obs. score) score (unrounded) moments */
    };
  }
}

#endif