/* 
  From Source: ERutilities.h 
  Original Struct: ERAW_RESULTS
  Description: equated raw score results
*/

#ifndef STRUCTURES_EQUATED_RAW_SCORE_RESULTS_HPP
#define STRUCTURES_EQUATED_RAW_SCORE_RESULTS_HPP

#include <string>
#include <Eigen/Core>
#include <fmt/core.h>
#include <fmt/format.h>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct EquatedRawScoreResults {
      Eigen::VectorXd xSyntheticPopulationMean(4); /* mean for x for synthetic pop */
      Eigen::VectorXd ySyntheticPopulationMean(4); /* mean for y for synthetic pop */
      Eigen::VectorXd xSyntheticPopulationSD(4);   /* sd for x for synthetic pop */
      Eigen::VectorXd ySyntheticPopulationSD(4);   /* sd for y for synthetic pop */
      Eigen::VectorXd gammaPopulation1(4);         /* gamma for pop 1 */
      Eigen::VectorXd gammaPopulation2(4);         /* gamma for pop 2 */
      Eigen::VectorXd slope(4);                    /* slope */
      Eigen::VectorXd intercept(4);                /* intercept */
      Eigen::MatrixXd equatedRawScores;            /* equated raw scores */
      Eigen::MatrixXd equatedRawScoreMoments;      /* moments for equated raw scores */
      Eigen::MatrixXd relativeFreqDistsX;          /* rel FD for X and syn pop: [0] for FE, [1] for MFE */
      Eigen::MatrixXd relativeFreqDistsY;          /* rel FD for Y and syn pop: [0] for FE, [1] for MFE */

      std::string toString() {
        std::string msg = "Equated Raw Score Results\n";
        
        msg.append(fmt::format("xSyntheticPopulationMean:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(xSyntheticPopulationMean, false)));
        msg.append(fmt::format("ySyntheticPopulationMean:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(ySyntheticPopulationMean, false)));
        msg.append(fmt::format("xSyntheticPopulationSD:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(xSyntheticPopulationSD, false)));
        msg.append(fmt::format("ySyntheticPopulationSD:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(ySyntheticPopulationSD, false)));
        msg.append(fmt::format("gammaPopulation1:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(gammaPopulation1, false)));
        msg.append(fmt::format("gammaPopulation2:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(gammaPopulation2, false)));
        msg.append(fmt::format("slope:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(slope, false)));
        msg.append(fmt::format("intercept:\n{}\n", EquatingRecipes::Utilities::vectorXdToString(intercept, false)));
        msg.append(fmt::format("equatedRawScores:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(equatedRawScores.transpose())));
        msg.append(fmt::format("equatedRawScoreMoments:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(equatedRawScoreMoments.transpose())));
        msg.append(fmt::format("relativeFreqDistsX:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(relativeFreqDistsX)));
        msg.append(fmt::format("relativeFreqDistsY:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(relativeFreqDistsY)));

        return msg;
      }
    };
  }
}

#endif