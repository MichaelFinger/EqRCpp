/* 
  From Source: ERutilities.h 
  Original Struct: USTATS
  Description: raw-score statistics for a univariate distribution 

  From Source: ERutilities.h, ERutilities.c
  Original Method: ReadFdGet_USTATS
  Description: Read frequency distribution from file fp

  From Source: ERutilities.h, ERutilities.c
  Original Method: ReadRawGet_USTATS
  Description: Read raw data and get or assign all elements for struct s
*/

#ifndef STRUCTURES_UNIVARIATE_STATISTICS_HPP
#define STRUCTURES_UNIVARIATE_STATISTICS_HPP

#include <string>
#include <Eigen/Core>
#include <equating_recipes/structures/moments.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct UnivariateStatistics {
      // std::string inputFilename;                 // name of input file
      std::string id;                               // single character id
      int numberOfExaminees;                        // number of examinees
      double minimumScore;                          // min score in data
      double maximumScore;                          // max score in data
      double freqDistMinimumScore;                  // min score for fd[]
      double freqDistMaximumScore;                  // max score for fd[]
      double scoreIncrement;                        // increment between adjacent scores
      int numberOfScores;                           // number of scores (or categories)
      Eigen::VectorXi freqDist;                     // freq dist fd[0]...fd[ns-1]
      Eigen::VectorXd freqDistDouble;               // double version of fd[]
      Eigen::VectorXi cumulativeFreqDist;           // cum freq dist
      Eigen::VectorXd relativeFreqDist;             // relative freq dist
      Eigen::VectorXd cumulativeRelativeFreqDist;   // cum relative freq dist
      Eigen::VectorXd percentileRankDist;           // percentile rank dist
      Eigen::VectorXd momentValues;                 // moments: mean, sd, skew, kurt

      void configure(const double& minimumScore,
                     const double& maximumScore,
                     const double& scoreIncrement);

      static UnivariateStatistics buildFromScoreFrequencies(const Eigen::VectorXi& scoreFrequencies,
                                         const double& minimumScore,
                                         const double& maximumScore,
                                         const double& scoreIncrement,
                                         const std::string& id);

      static UnivariateStatistics buildFromScores(const Eigen::VectorXd& scores,
                                         const double& minimumScore,
                                         const double& maximumScore,
                                         const double& scoreIncrement,
                                         const std::string& id);

      std::string toString();
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif