/* 
  From Source: ERutilities.h 
  Original Struct: USTATS
  Description: raw-score statistics for a univariate distribution 
*/

#ifndef STRUCTURES_UNIVARIATE_STATISTICS_HPP
#define STRUCTURES_UNIVARIATE_STATISTICS_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct UnivariateStatistics {
      // std::string inputFilename;               // name of input file
      // std::string id;                          // single character id
      int numberOfExaminees;                      // number of examinees
      double minimumObservedScore;                // min score in data
      double maximumObservedScore;                // max score in data
      double freqDistMinScore;                    // min score for fd[]
      double freqDistMaxScore;                    // max score for fd[]
      double adjacentScoresIncrement;             // increment between adjacent scores
      int numberOfScores;                         // number of scores (or categories)
      Eigen::VectorXi freqDist;                   // freq dist fd[0]...fd[ns-1]
      Eigen::VectorXd freqDistDouble;             // double version of fd[]
      Eigen::VectorXi cumulativeFreqDist;         // cum freq dist
      Eigen::VectorXd relativeFreqDist;           // relative freq dist
      Eigen::VectorXd cumulativeRelativeFreqDist; // cum relative freq dist
      Eigen::VectorXd percentileRankDist;         // percentile rank dist
      Eigen::VectorXd moments;                    // moments: mean, sd, skew, kurt
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif