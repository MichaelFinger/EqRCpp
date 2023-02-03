/* 
  From Source: ERutilities.h 
  Original Struct: ERAW_RESULTS
  Description: equated raw score results
*/

#ifndef STRUCTURES_EQUATED_RAW_SCORE_RESULTS_HPP
#define STRUCTURES_EQUATED_RAW_SCORE_RESULTS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedRawScoreResults {
      Eigen::VectorXd xSyntheticPopulationMean; /* mean for x for synthetic pop */
      Eigen::VectorXd ySyntheticPopulationMean; /* mean for y for synthetic pop */
      Eigen::VectorXd xSyntheticPopulationSD;   /* sd for x for synthetic pop */
      Eigen::VectorXd ySyntheticPopulationSD;   /* sd for y for synthetic pop */
      Eigen::VectorXd gammaPopulation1;         /* gamma for pop 1 */
      Eigen::VectorXd gammaPopulation2;         /* gamma for pop 2 */
      Eigen::VectorXd slope;                    /* slope */
      Eigen::VectorXd intercept;                /* intercept */
      Eigen::MatrixXd equatedRawScores;         /* equated raw scores */
      Eigen::MatrixXd equatedRawScoreMoments;   /* moments for equated raw scores */
      Eigen::MatrixXd relativeFreqDistsX;       /* rel FD for X and syn pop: [0] for FE, [1] for MFE */
      Eigen::MatrixXd relativeFreqDistsY;       /* rel FD for Y and syn pop: [0] for FE, [1] for MFE */
    };
  }
}

#endif