#include <algorithm>
#include <set>
#include <string>
#include <fmt/core.h>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/moments.hpp>

namespace EquatingRecipes {
  namespace Structures {
    UnivariateStatistics UnivariateStatistics::create(const std::map<double, int>& scoreFreqDist,
                                                      const double& minimumScore,
                                                      const double& maximumScore,
                                                      const double& scoreIncrement,
                                                      const std::string& id) {
      UnivariateStatistics univariateStatistics;

      univariateStatistics.id = id;
      univariateStatistics.minimumScore = minimumScore;
      univariateStatistics.maximumScore = maximumScore;
      univariateStatistics.adjacentScoresIncrement = scoreIncrement;
      univariateStatistics.numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore,
                                                                                       maximumScore,
                                                                                       scoreIncrement);

      univariateStatistics.freqDist.setZero(univariateStatistics.numberOfScores);
      univariateStatistics.freqDistDouble.setZero(univariateStatistics.numberOfScores);
      univariateStatistics.relativeFreqDist.setZero(univariateStatistics.numberOfScores);
      univariateStatistics.cumulativeFreqDist.setZero(univariateStatistics.numberOfScores);
      univariateStatistics.cumulativeRelativeFreqDist.setZero(univariateStatistics.numberOfScores);
      univariateStatistics.percentileRankDist.setZero(univariateStatistics.numberOfScores);

      univariateStatistics.moments.setZero(4);

      univariateStatistics.numberOfExaminees = 0;

      std::set<size_t> scoreIndicesWithNonzeroFreq;

      int cumulativeFreq = 0;

      std::for_each(scoreFreqDist.begin(),
                    scoreFreqDist.end(),
                    [&](const std::pair<double, int>& entry) {
                      double scoreValue = entry.first;
                      int scoreFreq = entry.second;

                      size_t scoreIndex = EquatingRecipes::Utilities::getScoreLocation(scoreValue,
                                                                                       minimumScore,
                                                                                       scoreIncrement);

                      if (scoreFreq > 0) {
                        scoreIndicesWithNonzeroFreq.insert(scoreIndex);
                      }

                      cumulativeFreq += scoreFreq;

                      univariateStatistics.freqDist(scoreIndex) = scoreFreq;
                      univariateStatistics.cumulativeFreqDist(scoreIndex) = cumulativeFreq;
                      univariateStatistics.numberOfExaminees += scoreFreq;
                    });

      univariateStatistics.freqDistDouble = univariateStatistics.freqDist.cast<double>();

      univariateStatistics.relativeFreqDist = univariateStatistics.freqDistDouble /
                                              static_cast<double>(univariateStatistics.numberOfExaminees);

      univariateStatistics.cumulativeRelativeFreqDist = univariateStatistics.cumulativeFreqDist.cast<double>() /
                                                        static_cast<double>(univariateStatistics.numberOfExaminees);

      size_t minimumFreqDistScoreIndex = *(scoreIndicesWithNonzeroFreq.begin());
      size_t maximumFreqDistScoreIndex = *(scoreIndicesWithNonzeroFreq.rbegin());

      univariateStatistics.freqDistMinimumScore = EquatingRecipes::Utilities::getScore(minimumFreqDistScoreIndex,
                                                                                       minimumScore,
                                                                                       scoreIncrement);

      univariateStatistics.freqDistMaximumScore = EquatingRecipes::Utilities::getScore(maximumFreqDistScoreIndex,
                                                                                       minimumScore,
                                                                                       scoreIncrement);

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Structures::Moments::getScoreMoments(scoreFreqDist);

      univariateStatistics.moments = moments.momentValues;

      std::for_each(scoreFreqDist.begin(),
                    scoreFreqDist.end(),
                    [&](const std::pair<double, int>& entry) {
                      double scoreValue = entry.first;
                      int scoreFreq = entry.second;
                      size_t scoreLocation = EquatingRecipes::Utilities::getScoreLocation(scoreValue,
                                                                                          minimumScore,
                                                                                          scoreIncrement);

                      univariateStatistics.percentileRankDist(scoreLocation) = EquatingRecipes::Utilities::percentileRank(minimumScore,
                                                                                                                          maximumScore,
                                                                                                                          scoreIncrement,
                                                                                                                          univariateStatistics.cumulativeRelativeFreqDist,
                                                                                                                          scoreValue);
                    });

      return univariateStatistics;
    }

    UnivariateStatistics UnivariateStatistics::create(const Eigen::VectorXd& scores,
                                                      const double& minimumScore,
                                                      const double& maximumScore,
                                                      const double& scoreIncrement,
                                                      const std::string& id) {
      std::map<double, int> freqDist;

      UnivariateStatistics univariateStatistics = UnivariateStatistics::create(freqDist,
                                                                               minimumScore,
                                                                               maximumScore,
                                                                               scoreIncrement,
                                                                               id);

      return univariateStatistics;
    }

    std::string UnivariateStatistics::toString() {
      std::string msg = "";

      msg.append(fmt::format("Score Variable ID: {}\n", this->id));
      msg.append(fmt::format("Number of examinees: {}\n", this->numberOfExaminees));
      msg.append(fmt::format("min score in data: {}\n", minimumScore));
      msg.append(fmt::format("max score in data: {}\n", maximumScore));
      msg.append(fmt::format("min score for fd[]: {}\n", freqDistMinimumScore));
      msg.append(fmt::format("max score for fd[]: {}\n", freqDistMaximumScore));
      msg.append(fmt::format("increment between adjacent scores: {}\n", adjacentScoresIncrement));
      msg.append(fmt::format("number of scores (or categories): {}\n", numberOfScores));
      msg.append(fmt::format("freq dist fd[0]...fd[ns-1]: {}\n", EquatingRecipes::Utilities::vectorXiToString(freqDist, false)));
      msg.append(fmt::format("double version of fd[]: {}\n", EquatingRecipes::Utilities::vectorXdToString(freqDistDouble, false)));
      msg.append(fmt::format("cum freq dist: {}\n", EquatingRecipes::Utilities::vectorXiToString(cumulativeFreqDist, false)));
      msg.append(fmt::format("relative freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(relativeFreqDist, false)));
      msg.append(fmt::format("cum relative freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(cumulativeRelativeFreqDist, false)));
      msg.append(fmt::format("percentile rank dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(percentileRankDist, false)));
      msg.append(fmt::format("moments: mean, sd, skew, kurt: {}\n", EquatingRecipes::Utilities::vectorXdToString(moments, false)));

      return msg;
    }
  } // namespace Structures
} // namespace EquatingRecipes