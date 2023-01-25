/* 
  From Source: ERutilities.h 
  Original Struct: ERAW_RESULTS
  Description: raw-score statistics for a univariate distribution 
*/

/*
truct ERAW_RESULTS{
  /* 
    equated raw-score results 
  */
  // double msx[4];                         /* mean for x for synthetic pop */
  // double msy[4];                         /* mean for y for synthetic pop */
  // double ssx[4];                           /* sd for x for synthetic pop */
  // double ssy[4];                           /* sd for y for synthetic pop */
  // double gamma1[4];                                   /* gamma for pop 1 */
  // double gamma2[4];                                   /* gamma for pop 2 */
  // double a[4];                                                  /* slope */
  // double b[4];                                              /* intercept */
  // double **eraw;                                   /* equated raw scores */
  // double **mts;                        /* moments for equated raw scores */
  // double **fxs;     /* rel FD for X and syn pop: [0] for FE, [1] for MFE */
  // double **gys;     /* rel FD for Y and syn pop: [0] for FE, [1] for MFE */

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
      Eigen::MatrixXd equatedRawScoreMoments    /* moments for equated raw scores */
      Eigen::MatrixXd relativeFreqDistsX;       /* rel FD for X and syn pop: [0] for FE, [1] for MFE */
      Eigen::MatrixXd relativeFreqDistsY;       /* rel FD for Y and syn pop: [0] for FE, [1] for MFE */
    };
  }
}

#endif