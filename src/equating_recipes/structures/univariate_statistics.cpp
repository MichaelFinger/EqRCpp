#include <algorithm>
#include <string>
#include <fmt/core.h>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/moments.hpp>

namespace EquatingRecipes {
  namespace Structures {
    void UnivariateStatistics::configure(const double& minimumScore,
                                         const double& maximumScore,
                                         const double& scoreIncrement) {
      this->minimumScore = minimumScore;
      this->maximumScore = maximumScore;
      this->scoreIncrement = scoreIncrement;
      this->numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore,
                                                                        maximumScore,
                                                                        scoreIncrement);

      this->freqDist.setZero(this->numberOfScores);
      this->freqDistDouble.setZero(this->numberOfScores);
      this->relativeFreqDist.setZero(this->numberOfScores);
      this->cumulativeFreqDist.setZero(this->numberOfScores);
      this->cumulativeRelativeFreqDist.setZero(this->numberOfScores);
      this->percentileRankDist.setZero(this->numberOfScores);
    }

    UnivariateStatistics UnivariateStatistics::create(const Eigen::VectorXi& scoreFrequencies,
                                                      const double& minimumScore,
                                                      const double& maximumScore,
                                                      const double& scoreIncrement,
                                                      const std::string& id) {
      UnivariateStatistics univariateStatistics;

      /*
      Eigen::VectorXi cumulativeFreqDist;           // cum freq dist
      Eigen::VectorXd cumulativeRelativeFreqDist;   // cum relative freq dist
      Eigen::VectorXd percentileRankDist;           // percentile rank dist
      */

      univariateStatistics.id = id;
      univariateStatistics.numberOfExaminees = 0;
      univariateStatistics.configure(minimumScore,
                                     maximumScore,
                                     scoreIncrement);

      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Structures::Moments::getScoreMoments(scoreFrequencies,
                                                                                                           minimumScore,
                                                                                                           maximumScore,
                                                                                                           scoreIncrement);

      univariateStatistics.numberOfExaminees = moments.numberOfExaminees;
      univariateStatistics.freqDistMinimumScore = moments.minimumObservedScore;
      univariateStatistics.freqDistMaximumScore = moments.maximumObservedScore;
      univariateStatistics.momentValues = moments.momentValues;

      univariateStatistics.freqDist = scoreFrequencies;
      univariateStatistics.freqDistDouble = univariateStatistics.freqDist.cast<double>();
      univariateStatistics.relativeFreqDist = univariateStatistics.freqDistDouble /
                                              static_cast<double>(univariateStatistics.numberOfExaminees);

      univariateStatistics.cumulativeFreqDist(0) = scoreFrequencies(0);
      for (size_t index = 1; index < scoreFrequencies.size(); index++) {
        univariateStatistics.cumulativeFreqDist(index) = univariateStatistics.cumulativeFreqDist(index - 1) + scoreFrequencies(index);
      }

      univariateStatistics.cumulativeRelativeFreqDist = univariateStatistics.cumulativeFreqDist.cast<double>() /
                                                        static_cast<double>(univariateStatistics.numberOfExaminees);

      univariateStatistics.percentileRankDist = EquatingRecipes::Utilities::percentileRanks(minimumScore,
                                                                                            maximumScore,
                                                                                            scoreIncrement,
                                                                                            univariateStatistics.cumulativeRelativeFreqDist);

      return univariateStatistics;
    }

    UnivariateStatistics UnivariateStatistics::create(const Eigen::VectorXd& scores,
                                                      const double& minimumScore,
                                                      const double& maximumScore,
                                                      const double& scoreIncrement,
                                                      const std::string& id) {
      Eigen::VectorXi freqDist = EquatingRecipes::Utilities::getRawScoreFrequencyDistribution(scores,
                                                                                                    minimumScore,
                                                                                                    maximumScore,
                                                                                                    scoreIncrement,
                                                                                                    true);

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
      msg.append(fmt::format("increment between adjacent scores: {}\n", scoreIncrement));
      msg.append(fmt::format("number of scores (or categories): {}\n", numberOfScores));
      msg.append(fmt::format("freq dist fd[0]...fd[ns-1]: {}\n", EquatingRecipes::Utilities::vectorXiToString(freqDist, false)));
      msg.append(fmt::format("double version of fd[]: {}\n", EquatingRecipes::Utilities::vectorXdToString(freqDistDouble, false)));
      msg.append(fmt::format("cum freq dist: {}\n", EquatingRecipes::Utilities::vectorXiToString(cumulativeFreqDist, false)));
      msg.append(fmt::format("relative freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(relativeFreqDist, false)));
      msg.append(fmt::format("cum relative freq dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(cumulativeRelativeFreqDist, false)));
      msg.append(fmt::format("percentile rank dist: {}\n", EquatingRecipes::Utilities::vectorXdToString(percentileRankDist, false)));
      msg.append(fmt::format("moments: mean, sd, skew, kurt: {}\n", EquatingRecipes::Utilities::vectorXdToString(momentValues, false)));

      return msg;
    }
  } // namespace Structures
} // namespace EquatingRecipes