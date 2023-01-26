#include <algorithm>
#include <set>
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
      univariateStatistics.numberOfScores = EquatingRecipes::Utilities::numberOfScores(maximumScore,
                                                                                       minimumScore,
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
                      // univariateStatistics.freqDistDouble(scoreIndex) = static_cast<double>(scoreFreq);
                      univariateStatistics.cumulativeFreqDist(scoreIndex) = cumulativeFreq;

                      univariateStatistics.numberOfExaminees += scoreFreq;
                    });

      univariateStatistics.freqDistDouble = univariateStatistics.freqDist.cast<double>();

      univariateStatistics.relativeFreqDist = univariateStatistics.freqDistDouble /
                                              static_cast<double>(univariateStatistics.numberOfExaminees);

      univariateStatistics.cumulativeRelativeFreqDist = univariateStatistics.cumulativeFreqDist.cast<double>() /
                                                        static_cast<double>(univariateStatistics.numberOfExaminees);

      size_t minimumFreqDistScoreIndex = *(scoreIndicesWithNonzeroFreq.begin());
      size_t maximumFreqDistScoreIndex = *(scoreIndicesWithNonzeroFreq.end());

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
  } // namespace Structures
} // namespace EquatingRecipes