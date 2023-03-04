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

#include <algorithm>
#include <string>

#include <Eigen/Core>
#include <fmt/core.h>

namespace EquatingRecipes {
  namespace Structures {
    struct UnivariateStatistics {
      std::string id;                             // single character id
      int numberOfExaminees;                      // number of examinees
      double minimumScore;                        // min score in data
      double maximumScore;                        // max score in data
      double freqDistMinimumScore;                // min score for fd[]
      double freqDistMaximumScore;                // max score for fd[]
      double scoreIncrement;                      // increment between adjacent scores
      size_t numberOfScores;                         // number of scores (or categories)
      Eigen::VectorXd freqDist;                   // freq dist fd[0]...fd[ns-1]
      Eigen::VectorXd freqDistDouble;             // double version of fd[]
      Eigen::VectorXd cumulativeFreqDist;         // cum freq dist
      Eigen::VectorXd relativeFreqDist;           // relative freq dist
      Eigen::VectorXd cumulativeRelativeFreqDist; // cum relative freq dist
      Eigen::VectorXd percentileRankDist;         // percentile rank dist
      Eigen::VectorXd momentValues;               // moments: mean, sd, skew, kurt
      
      void configure(const double& minimumScore,
                     const double& maximumScore,
                     const double& scoreIncrement,
                     const size_t& numberOfScores) {
        this->minimumScore = minimumScore;
        this->maximumScore = maximumScore;
        this->scoreIncrement = scoreIncrement;
        this->numberOfScores = numberOfScores;

        this->freqDist.setZero(this->numberOfScores);
        this->freqDistDouble.setZero(this->numberOfScores);
        this->relativeFreqDist.setZero(this->numberOfScores);
        this->cumulativeFreqDist.setZero(this->numberOfScores);
        this->cumulativeRelativeFreqDist.setZero(this->numberOfScores);
        this->percentileRankDist.setZero(this->numberOfScores);
        this->momentValues.setZero(4);
      }

      // std::string toString() {
      //   std::string msg = "";

      //   msg.append(fmt::format("Score Variable ID: {}\n", this->id));
      //   msg.append(fmt::format("Number of examinees: {}\n", this->numberOfExaminees));
      //   msg.append(fmt::format("min score in data: {}\n", minimumScore));
      //   msg.append(fmt::format("max score in data: {}\n", maximumScore));
      //   msg.append(fmt::format("min score for fd[]: {}\n", freqDistMinimumScore));
      //   msg.append(fmt::format("max score for fd[]: {}\n", freqDistMaximumScore));
      //   msg.append(fmt::format("increment between adjacent scores: {}\n", scoreIncrement));
      //   msg.append(fmt::format("number of scores (or categories): {}\n", numberOfScores));
      //   msg.append(fmt::format("freq dist fd[0]...fd[ns-1]: {}\n", EquatingRecipes::Utilities::vectorXdToString(freqDist, false)));
      //   msg.append(fmt::format("double version of fd[]: {}\n", EquatingRecipes::Utilities::vectorXdToString(freqDistDouble, false)));
      //   msg.append(fmt::format("cum freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(cumulativeFreqDist, false)));
      //   msg.append(fmt::format("relative freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(relativeFreqDist, false)));
      //   msg.append(fmt::format("cum relative freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(cumulativeRelativeFreqDist, false)));
      //   msg.append(fmt::format("percentile rank dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(percentileRankDist, false)));
      //   msg.append(fmt::format("moments: mean, sd, skew, kurt: {}\n", EquatingRecipes::Utilities::vectorXdToString(momentValues, false)));

      //   return msg;
      // }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif