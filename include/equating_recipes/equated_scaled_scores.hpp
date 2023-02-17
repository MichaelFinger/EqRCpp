/* 
  From Source: ERutilities.h, ERutilities.c
  Original Struct: Wrapper_ESS
  Description: Wrapper used for getting equated scale scores essu[][] and essr[][]
    Assigns or computes values for all variables in struct ESS_RESULTS s
    Can be used for any design or methodology.
*/

#ifndef WRAPPERS_EQUATED_SCALED_SCORES_HPP
#define WRAPPERS_EQUATED_SCALED_SCORES_HPP

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <Eigen/Core>

#include <equating_recipes/score_statistics.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  struct EquatedScaledScores {
    void run(EquatingRecipes::Structures::PData& pData,
             const EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults,
             const double& minimumRawScoreYct,
             const double& maximumRawScoreYct,
             const double& scoreIncrementYct,
             const size_t& roundToNumberOfDecimalPlaces,
             const int& lowestObservableRoundedScaledScore,
             const int& highestObservableRoundedScaledScore,
             const EquatingRecipes::Structures::RawToScaledScoreTable& rawToScaledScoreTable,
             EquatingRecipes::Structures::EquatedScaledScoresResults& results) {
      size_t numberOfScores = EquatingRecipes::Utilities::getNumberOfScores(pData.mininumScoreX,
                                                                         pData.maximumScoreX,
                                                                         pData.scoreIncrementX);

      pData.minimumRawScoreYct = minimumRawScoreYct;
      pData.maximumRawScoreYct = maximumRawScoreYct;
      pData.scoreIncrementYct = scoreIncrementYct;
      pData.roundToNumberOfDecimalPlaces = roundToNumberOfDecimalPlaces;
      pData.lowestObservableRoundedScaledScore = lowestObservableRoundedScaledScore;
      pData.highestObservableRoundedScaledScore = highestObservableRoundedScaledScore;
      pData.rawToScaledScoreTable = rawToScaledScoreTable;

      if (pData.bootstrapReplicationNumber <= 1) {
        // actual equating or 1st bootstrap replication
        results.configure(numberOfScores, pData.methods.size());
      }

      for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
        Eigen::VectorXd unroundedEquatedScaledScores = Eigen::VectorXd::Zero(results.unroundedEquatedScaledScores.rows());
        Eigen::VectorXd roundedEquatedScaledScores = Eigen::VectorXd::Zero(results.roundedEquatedScaledScores.rows());

        EquatingRecipes::Utilities::getEquatedScaledScores(pData.mininumScoreX,
                                                           pData.maximumScoreX,
                                                           pData.scoreIncrementX,
                                                           pData.minimumRawScoreYct,
                                                           pData.maximumRawScoreYct,
                                                           pData.scoreIncrementYct,
                                                           equatedRawScoreResults.equatedRawScores.col(methodIndex),
                                                           rawToScaledScoreTable,
                                                           pData.roundToNumberOfDecimalPlaces,
                                                           pData.lowestObservableRoundedScaledScore,
                                                           pData.highestObservableRoundedScaledScore,
                                                           unroundedEquatedScaledScores,
                                                           roundedEquatedScaledScores);

        results.unroundedEquatedScaledScores.col(methodIndex) = unroundedEquatedScaledScores;
        results.roundedEquatedScaledScores.col(methodIndex) = roundedEquatedScaledScores;
      }

      /* compute moments:  Note that when inc==1, essu[*][min-minp+1] is the
        unrounded scale score associated with fdx[0], where fdx[0]
        is associated with scores ranging from min to max.
        Recall that essu[*][0] is the unrounded scale score
        associated with minp-inc/2 = minp-.5 when inc = 1.  In this example,
        min-minp+1 = loc(min,minp,inc) + 1
      */

      if (pData.scoreFrequenciesX.size() >= 1) {
        EquatingRecipes::Structures::Moments moments;

        for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
          moments = EquatingRecipes::ScoreStatistics::momentsFromScoreFrequencies(results.unroundedEquatedScaledScores.col(methodIndex),
                                                                                  pData.scoreFrequenciesX);

          results.unroundedEquatedScaledScoreMoments.col(methodIndex) = moments.momentValues;

          moments = EquatingRecipes::ScoreStatistics::momentsFromScoreFrequencies(results.roundedEquatedScaledScores.col(methodIndex),
                                                                                  pData.scoreFrequenciesX);

          results.roundedEquatedScaledScoreMoments.col(methodIndex) = moments.momentValues;
        }
      }
    }
  };
} // namespace EquatingRecipes

#endif