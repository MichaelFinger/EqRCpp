// struct IRT_TRUE{
//   /*
//     IRT true score equating:
// 	item parameters and other statitistics for a test form.
//   */
//   int num_items;		                      /* number of items on test */
//   int num_scores;                     /* number of test score categories */
//   int *num_icat;          /* pointer to numbers of item score categories */
//   int *num_choice;         /* pointer to numbers of choices for MC items */
//   double *right;                          /* poiner to right item scores */
//   double *wrong;                          /* poiner to wrong item scores */
//   double *a_par;                              /* pointer to a parameters */
//   double *b_par;                              /* pointer to b parameters */
//   double *c_par;                              /* pointer to c parameters */
//   double **d_par;                 /* pointer to pointers to d parameters */
//   double *scores;                              /* pointer to test scores */
//   double *theta;                 /* theta values corresponding to scores */
//   double min;                             /* minimal possible test score */
//             /* min may not equal to the actually minimal score scores[0] */
//   double max;                             /* maximum possible test score */
//   double inc;                                         /* score increment */
//   double chance;                          /* the chance level test score */
// };

/* 
  From Source: ERutilities.h 
  Original Struct: IRT_TRUE
  Description: IRT true score equating: item parameters and other statitistics for a test form.
*/

#ifndef STRUCTURES_IRT_TRUE_EQUATING_HPP
#define STRUCTURES_IRT_TRUE_EQUATING_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct IrtTrueScoreEquating {
      int numberOfItems;                           /* number of items on test */
      int numberOfTestScoreCategories;             /* number of test score categories */
      Eigen::VectorXi numberOfItemScoreCategories; /* pointer to numbers of item score categories */
      Eigen::VectorXi numberOfMCItemChoices;       /* pointer to numbers of choices for MC items */
      Eigen::VectorXi correctItemScores;           /* poiner to right item scores */
      Eigen::VectorXi incorrectItemScores;         /* poiner to wrong item scores */
      Eigen::VectorXd a;                           /* pointer to a parameters */
      Eigen::VectorXd b;                           /* pointer to b parameters */
      Eigen::VectorXd c;                           /* pointer to c parameters */
      Eigen::MatrixXd d;                           /* pointer to pointers to d parameters */
      Eigen::VectorXd testScores;                  /* pointer to test scores */
      Eigen::VectorXd theta;                       /* theta values corresponding to scores */
      double minimumObservableTestScore;           /* minimal possible test score; min may not equal to the actually minimal score scores[0] */
      double maximumObservableTestScore;           /* maximum possible test score */
      double scoreIncrement;                       /* score increment */
      double chanceLevelTestScore;                 /* the chance level test score */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif