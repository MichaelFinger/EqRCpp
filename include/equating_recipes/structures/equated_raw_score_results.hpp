/* 
  From Source: ERutilities.h 
  Original Struct: ERAW_RESULTS
  Description: equated raw score results
*/

#ifndef IMPLEMENTATION_EQUATED_RAW_SCORE_RESULTS_HPP
#define IMPLEMENTATION_EQUATED_RAW_SCORE_RESULTS_HPP

#include <string>
#include <Eigen/Core>
#include <fmt/core.h>
#include <fmt/format.h>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedRawScoreResults {
      Eigen::VectorXd xSyntheticPopulationMean = Eigen::VectorXd(4); /* mean for x for synthetic pop */
      Eigen::VectorXd ySyntheticPopulationMean = Eigen::VectorXd(4); /* mean for y for synthetic pop */
      Eigen::VectorXd xSyntheticPopulationSD = Eigen::VectorXd(4);   /* sd for x for synthetic pop */
      Eigen::VectorXd ySyntheticPopulationSD = Eigen::VectorXd(4);   /* sd for y for synthetic pop */
      Eigen::VectorXd gammaPopulation1 = Eigen::VectorXd(4);         /* gamma for pop 1 */
      Eigen::VectorXd gammaPopulation2 = Eigen::VectorXd(4);         /* gamma for pop 2 */
      Eigen::VectorXd slope = Eigen::VectorXd(4);                    /* slope */
      Eigen::VectorXd intercept = Eigen::VectorXd(4);                /* intercept */
      Eigen::MatrixXd equatedRawScores;            /* equated raw scores */
      Eigen::MatrixXd equatedRawScoreMoments;      /* moments for equated raw scores */
      Eigen::MatrixXd relativeFreqDistsX;          /* rel FD for X and syn pop: [0] for FE, [1] for MFE */
      Eigen::MatrixXd relativeFreqDistsY;          /* rel FD for Y and syn pop: [0] for FE, [1] for MFE */

      // std::string toString() {
      //   std::string msg = "Equated Raw Score Results\n";
        
      //   msg.append(fmt::format("xSyntheticPopulationMean:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(xSyntheticPopulationMean, false)));
      //   msg.append(fmt::format("ySyntheticPopulationMean:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(ySyntheticPopulationMean, false)));
      //   msg.append(fmt::format("xSyntheticPopulationSD:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(xSyntheticPopulationSD, false)));
      //   msg.append(fmt::format("ySyntheticPopulationSD:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(ySyntheticPopulationSD, false)));
      //   msg.append(fmt::format("gammaPopulation1:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(gammaPopulation1, false)));
      //   msg.append(fmt::format("gammaPopulation2:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(gammaPopulation2, false)));
      //   msg.append(fmt::format("slope:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(slope, false)));
      //   msg.append(fmt::format("intercept:\n{}\n", EquatingRecipes::Implementation::Utilities::vectorXdToString(intercept, false)));
      //   msg.append(fmt::format("equatedRawScores:\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXdToString(equatedRawScores.transpose())));
      //   msg.append(fmt::format("equatedRawScoreMoments:\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXdToString(equatedRawScoreMoments.transpose())));
      //   msg.append(fmt::format("relativeFreqDistsX:\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXdToString(relativeFreqDistsX)));
      //   msg.append(fmt::format("relativeFreqDistsY:\n{}\n", EquatingRecipes::Implementation::Utilities::matrixXdToString(relativeFreqDistsY)));

      //   return msg;
      // }
    };
  }
}

#endif