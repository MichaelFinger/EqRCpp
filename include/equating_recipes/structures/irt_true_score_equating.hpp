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
    struct IRTTrueScoreEquating {
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