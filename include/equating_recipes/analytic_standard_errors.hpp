/* 
  From Source: Analytic_SEs.h, Analytic_SEs.c
  Original Struct: SE_EPequate
  Description: Estimated standard errors for equipercentile equivalents 
  on scale of y (see Kolen & Brennan, 2004, p. 248).
*/

#ifndef STRUCTURES_ANALYTIC_STANDARD_ERRORS_HPP
#define STRUCTURES_ANALYTIC_STANDARD_ERRORS_HPP

#include <cmath>
#include <Eigen/Core>

namespace EquatingRecipes {
  struct AnalyticStandardErrors {
    Eigen::VectorXd calculate(const size_t& numberOfRawScoreCategoriesY,
                              const size_t& numberOfExamineesY,
                              const Eigen::VectorXd& cumulativeRelativeFreqDistY,
                              const size_t& numberOfRawScoreCategoriesX,
                              const double& scoreIncrementX,
                              const size_t& numberOfExamineesX,
                              const Eigen::VectorXd& percentileRankDistX) {
      Eigen::VectorXd standardErrors(numberOfRawScoreCategoriesX);

      for (size_t scoreLocationX = 0; scoreLocationX < numberOfRawScoreCategoriesX; scoreLocationX++) {
        double percentileRankProportion = percentileRankDistX(scoreLocationX) / 100.0;

        if (percentileRankProportion >= 1.0) {
          standardErrors(scoreLocationX) = 0.0;
        } else {
          size_t scoreLocationY;

          for (scoreLocationY = 1; scoreLocationY < numberOfRawScoreCategoriesY; scoreLocationY++) {
            if (cumulativeRelativeFreqDistY(scoreLocationY) > percentileRankProportion) {
              break;
            }
          }

          if (cumulativeRelativeFreqDistY(scoreLocationY) != cumulativeRelativeFreqDistY(scoreLocationY - 1)) {
            double relativeFrequencyY = cumulativeRelativeFreqDistY(scoreLocationY) - cumulativeRelativeFreqDistY(scoreLocationY - 1);

            standardErrors(scoreLocationX) = (1.0 / std::pow(relativeFrequencyY, 2)) *
                                             (percentileRankProportion * (1.0 - percentileRankProportion) *
                                                  static_cast<double>(numberOfExamineesX + numberOfExamineesY) /
                                                  static_cast<double>(numberOfExamineesX * numberOfExamineesY) -
                                              (cumulativeRelativeFreqDistY(scoreLocationY) - percentileRankProportion) * (percentileRankProportion - cumulativeRelativeFreqDistY(scoreLocationY - 1)) /
                                                  (static_cast<double>(numberOfExamineesY) * relativeFrequencyY));

            standardErrors(scoreLocationX) = standardErrors(scoreLocationX) > 0.0 ? std::sqrt(standardErrors(scoreLocationX)) : 0.0;

            /* The following line handles non-unit increments */
            standardErrors(scoreLocationX) = scoreIncrementX * standardErrors(scoreLocationX);
          } else {
            standardErrors(scoreLocationX) = 0.0;
          }
        }
      }

      return standardErrors;
    }
  };
} // namespace EquatingRecipes

#endif