/* 
  From Source: ERutilities.h, ERutilities.c
  Original Method: ReadRawGet_BSTATS
  Description: Reads raw data and get or assign all elements for bivariate struct s

  Original Method: ReadRawGet_moments 
  Description: Compute moments from raw scores in file;
    scores need not be integers or positive;
    frequency distribution NOT computed or output
    assumes space allocated for moments[4];
    assumes data are whitespaced delimited
*/

#ifndef SCORE_STATISTICS_HPP
#define SCORE_STATISTICS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <Eigen/Core>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  struct ScoreStatistics {
    static EquatingRecipes::Structures::BivariateStatistics bivariate(const Eigen::MatrixXd& scores,
                                                                      const double& minimumRowScore,
                                                                      const double& maximumRowScore,
                                                                      const double& rowScoreIncrement,
                                                                      const double& minimumColumnScore,
                                                                      const double& maximumColumnScore,
                                                                      const double& columnScoreIncrement,
                                                                      const std::string& rowScoreId,
                                                                      const std::string& columnScoreId) {
      EquatingRecipes::Structures::BivariateStatistics bivariateStatistics;

      size_t rowScoreColumnIndex = 0;
      size_t columnScoreColumnIndex = 1;

      bivariateStatistics.univariateStatisticsRow = univariateFromScores(scores.col(rowScoreColumnIndex),
                                                                    minimumRowScore,
                                                                    maximumRowScore,
                                                                    rowScoreIncrement,
                                                                    rowScoreId);

      bivariateStatistics.univariateStatisticsColumn = univariateFromScores(scores.col(columnScoreColumnIndex),
                                                                       minimumColumnScore,
                                                                       maximumColumnScore,
                                                                       columnScoreIncrement,
                                                                       columnScoreId);

      bivariateStatistics.bivariateFreqDist.setZero(bivariateStatistics.univariateStatisticsRow.numberOfScores,
                                                    bivariateStatistics.univariateStatisticsColumn.numberOfScores);
      bivariateStatistics.bivariateFreqDistDouble.setZero(bivariateStatistics.univariateStatisticsRow.numberOfScores,
                                                    bivariateStatistics.univariateStatisticsColumn.numberOfScores);
      bivariateStatistics.bivariateProportions.setZero(bivariateStatistics.univariateStatisticsRow.numberOfScores,
                                                    bivariateStatistics.univariateStatisticsColumn.numberOfScores);

      bivariateStatistics.numberOfExaminees = scores.rows();

      bivariateStatistics.covariance = 0.0;

      for (size_t rowIndex = 0; rowIndex < scores.rows(); ++rowIndex) {
        double rowScore = scores(rowIndex, rowScoreColumnIndex);
        double columnScore = scores(rowIndex, columnScoreColumnIndex);

        size_t rowScoreLocation = EquatingRecipes::Utilities::getScoreLocation(rowScore,
                                                                               minimumRowScore,
                                                                               rowScoreIncrement);

        size_t columnScoreLocation = EquatingRecipes::Utilities::getScoreLocation(columnScore,
                                                                                  minimumColumnScore,
                                                                                  columnScoreIncrement);

        bivariateStatistics.univariateStatisticsRow.freqDist(rowScoreLocation)++;
        bivariateStatistics.univariateStatisticsColumn.freqDist(columnScoreLocation)++;
        bivariateStatistics.bivariateFreqDist(rowScoreLocation, columnScoreLocation)++;

        bivariateStatistics.covariance += rowScore * columnScore;
      }

      bivariateStatistics.bivariateFreqDistDouble = bivariateStatistics.bivariateFreqDist.cast<double>();
      bivariateStatistics.bivariateProportions = bivariateStatistics.bivariateFreqDistDouble / static_cast<double>(bivariateStatistics.numberOfExaminees);

      bivariateStatistics.covariance = (bivariateStatistics.covariance / static_cast<double>(bivariateStatistics.numberOfExaminees)) -
                                       (bivariateStatistics.univariateStatisticsRow.momentValues(0) * bivariateStatistics.univariateStatisticsColumn.momentValues(0));

      bivariateStatistics.correlation = bivariateStatistics.covariance /
                                        (bivariateStatistics.univariateStatisticsRow.momentValues(1) * bivariateStatistics.univariateStatisticsColumn.momentValues(1));

      return bivariateStatistics;
    }

    static EquatingRecipes::Structures::UnivariateStatistics univariateFromScoreFrequencies(const Eigen::VectorXd& scoreFrequencies,
                                                                                            const double& minimumScore,
                                                                                            const double& maximumScore,
                                                                                            const double& scoreIncrement,
                                                                                            const std::string& id) {
      EquatingRecipes::Structures::UnivariateStatistics univariateStatistics;

      univariateStatistics.id = id;
      univariateStatistics.numberOfExaminees = 0;
      univariateStatistics.configure(minimumScore,
                                     maximumScore,
                                     scoreIncrement);

      EquatingRecipes::Structures::Moments moments = momentsFromScoreFrequencies(scoreFrequencies,
                                                                                 minimumScore,
                                                                                 maximumScore,
                                                                                 scoreIncrement);

      univariateStatistics.numberOfExaminees = scoreFrequencies.sum();
      univariateStatistics.freqDistMinimumScore = EquatingRecipes::Utilities::getFirstObservedScore(scoreFrequencies,
                                                                                                    minimumScore,
                                                                                                    maximumScore,
                                                                                                    scoreIncrement,
                                                                                                    true);
      univariateStatistics.freqDistMaximumScore = EquatingRecipes::Utilities::getFirstObservedScore(scoreFrequencies,
                                                                                                    minimumScore,
                                                                                                    maximumScore,
                                                                                                    scoreIncrement,
                                                                                                    false);
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

    static EquatingRecipes::Structures::UnivariateStatistics univariateFromScores(const Eigen::VectorXd& scores,
                                                                                  const double& minimumScore,
                                                                                  const double& maximumScore,
                                                                                  const double& scoreIncrement,
                                                                                  const std::string& id) {
      Eigen::VectorXd freqDist = EquatingRecipes::Utilities::getRawScoreFrequencyDistribution(scores,
                                                                                              minimumScore,
                                                                                              maximumScore,
                                                                                              scoreIncrement,
                                                                                              true);

      EquatingRecipes::Structures::UnivariateStatistics univariateStatistics = univariateFromScoreFrequencies(freqDist,
                                                                                                              minimumScore,
                                                                                                              maximumScore,
                                                                                                              scoreIncrement,
                                                                                                              id);

      return univariateStatistics;
    }

    static EquatingRecipes::Structures::Moments momentsFromScores(const Eigen::VectorXd& scores) {
      EquatingRecipes::Structures::Moments scoreMoments;

      scoreMoments.momentValues.setZero(4);
      scoreMoments.momentValues(0) = scores.mean();

      Eigen::VectorXd meanVector = Eigen::VectorXd::Constant(scores.size(), scoreMoments.momentValues(0));
      Eigen::VectorXd deviations = scores - meanVector;

      double variance = deviations.array().pow(2).mean();
      double skewness = deviations.array().pow(3).mean();
      double kurtosis = deviations.array().pow(4).mean();

      scoreMoments.momentValues(1) = std::sqrt(variance);
      scoreMoments.momentValues(2) = skewness / std::pow(scoreMoments.momentValues(1), 3);
      scoreMoments.momentValues(3) = kurtosis / std::pow(scoreMoments.momentValues(1), 4);

      return scoreMoments;
    }

    static EquatingRecipes::Structures::Moments momentsFromScoreFrequencies(const Eigen::VectorXd& scoreFrequencies,
                                                                            const double& minimumScore,
                                                                            const double& maximumScore,
                                                                            const double& scoreIncrement) {
      size_t numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore,
                                                                         maximumScore,
                                                                         scoreIncrement);

      Eigen::VectorXd scores(numberOfScores);

      for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
        scores(scoreLocation) = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                     minimumScore,
                                                                     scoreIncrement);
      }

      EquatingRecipes::Structures::Moments scoreMoments = momentsFromScoreFrequencies(scores,
                                                                                      scoreFrequencies);

      return scoreMoments;
    }

    static EquatingRecipes::Structures::Moments momentsFromScoreFrequencies(const Eigen::VectorXd& scores,
                                                                            const Eigen::VectorXd& scoreFrequencies) {
      EquatingRecipes::Structures::Moments scoreMoments;
      scoreMoments.momentValues.setZero(4);

      size_t numberOfScores = scores.size();
      double numberOfExaminees = scoreFrequencies.sum();

      scoreMoments.momentValues(0) = scores.cwiseProduct(scoreFrequencies).sum() / numberOfExaminees;

      Eigen::VectorXd deviations = scores - Eigen::VectorXd::Constant(scores.size(), scoreMoments.momentValues(0));

      for (size_t powCoeff = 2; powCoeff <= 4; powCoeff++) {
        scoreMoments.momentValues(powCoeff - 1) = (deviations.array().pow(powCoeff)).cwiseProduct(scoreFrequencies.array()).sum() /
                                                  numberOfExaminees;
      }

      scoreMoments.momentValues(1) = std::sqrt(scoreMoments.momentValues(1));
      scoreMoments.momentValues(2) /= std::pow(scoreMoments.momentValues(1), 3);
      scoreMoments.momentValues(3) /= std::pow(scoreMoments.momentValues(1), 4);

      return scoreMoments;
    }
  };
} // namespace EquatingRecipes

#endif