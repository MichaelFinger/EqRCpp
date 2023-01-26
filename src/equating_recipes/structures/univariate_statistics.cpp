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
      univariateStatistics.minimumObservedScore = minimumScore;
      univariateStatistics.maximumObservedScore = maximumScore;
      univariateStatistics.adjacentScoresIncrement = scoreIncrement;
      univariateStatistics.numberOfScores = EquatingRecipes::Utilities::numberOfScores(maximumScore,
                                                                                       minimumScore,
                                                                                       scoreIncrement);

      univariateStatistics.freqDist.setZero(univariateStatistics.numberOfScores - 1);
      univariateStatistics.freqDistDouble.setZero(univariateStatistics.numberOfScores - 1);
      univariateStatistics.cumulativeFreqDist.setZero(univariateStatistics.numberOfScores - 1);
      univariateStatistics.cumulativeRelativeFreqDist.setZero(univariateStatistics.numberOfScores - 1);
      univariateStatistics.percentileRankDist.setZero(univariateStatistics.numberOfScores - 1);

      univariateStatistics.numberOfExaminees = 0;

      std::set<size_t> scoreIndicesWithNonzeroFreq;

      std::for_each(scoreFreqDist.begin(),
                    scoreFreqDist.end(),
                    [&](const std::pair<double, int>& entry) {
                      double scoreValue = entry.first;
                      int scoreFreq = entry.second;

                      size_t scoreIndex = EquatingRecipes::Utilities::scoreLocation(scoreValue,
                                                                                    minimumScore,
                                                                                    scoreIncrement);

                      if (scoreFreq > 0) {
                        scoreIndicesWithNonzeroFreq.insert(scoreIndex);
                      }

                      univariateStatistics.freqDist(scoreIndex) = scoreFreq;
                      univariateStatistics.freqDistDouble(scoreIndex) = static_cast<double>(scoreFreq);
                      univariateStatistics.numberOfExaminees += scoreFreq;
                    });

      univariateStatistics.freqDistMinimumScore = EquatingRecipes::Utilities::getScore(*(scoreIndicesWithNonzeroFreq.begin()),
                                                                                       minimumScore,
                                                                                       scoreIncrement);

      univariateStatistics.freqDistMaximumScore = EquatingRecipes::Utilities::getScore(*(scoreIndicesWithNonzeroFreq.end()),
                                                                                       minimumScore,
                                                                                       scoreIncrement);

      // for(i=0;i<=s->ns-1;i++) if(s->fd[i]) break; 
      //   s->mind = score(i,min,inc);
      //   for(i=s->ns-1;i>=0;i--) if(s->fd[i]) break; 
      //   s->maxd = score(i,min,inc);
      //   if((s->n = MomentsFromFD(min, max, inc, NULL, s->fd, s->mts)) != n)
      //     runerror("\nError somewhere in ReadFdGet_USTATS()");

      //   s->cfd[0] = s->fd[0];
      //   for(i=1;i<=loc(max,min,inc);i++) s->cfd[i] = s->cfd[i-1] + s->fd[i];
      //   for(i=0;i<=s->ns-1;i++) s->rfd[i] = s->fd[i]/((double) n);
      //   cum_rel_freqs(min,max,inc,s->rfd,s->crfd);

      //   for(i=0;i<=loc(max,min,inc);i++){
      //     s->prd[i] = perc_rank(min,max,inc,s->crfd,score(i,min,inc));
      //     s->dbl_fd[i] = s->fd[i];
      //   }

      

      return univariateStatistics;
    }

    UnivariateStatistics UnivariateStatistics::create(const Eigen::VectorXd& scores,
                                                      const double& minimumScore,
                                                      const double& maximumScore,
                                                      const double& scoreIncrement,
                                                      const std::string& id) {
    }
  } // namespace Structures
} // namespace EquatingRecipes