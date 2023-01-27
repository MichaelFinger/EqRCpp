#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  namespace Structures {
    BivariateStatistics BivariateStatistics::create(const Eigen::MatrixXd& scores,
                                                    const double& minimumRowScore,
                                                    const double& maximumRowScore,
                                                    const double& rowScoreIncrement,
                                                    const double& minimumColumnScore,
                                                    const double& maximumColumnScore,
                                                    const double& columnScoreIncrement,
                                                    const std::string& rowScoreId,
                                                    const std::string& columnScoreId) {
      BivariateStatistics bivariateStatistics;

      size_t rowScoreColumnIndex = 0;
      size_t columnScoreColumnIndex = 1;

      bivariateStatistics.rowScoreStatistics.id = rowScoreId;
      bivariateStatistics.rowScoreStatistics.configure(minimumRowScore,
                                                       maximumRowScore,
                                                       rowScoreIncrement);
      
      bivariateStatistics.columnScoreStatistics.id = columnScoreId;
      bivariateStatistics.columnScoreStatistics.configure(minimumColumnScore,
                                                          maximumColumnScore,
                                                          columnScoreIncrement);
      
      bivariateStatistics.bivariateFreqDist.setZero(bivariateStatistics.rowScoreStatistics.numberOfScores,
                                                    bivariateStatistics.columnScoreStatistics.numberOfScores);
      bivariateStatistics.bivariateFreqDistDouble.setZero(bivariateStatistics.rowScoreStatistics.numberOfScores,
                                                          bivariateStatistics.columnScoreStatistics.numberOfScores);
      bivariateStatistics.bivariateProportions.setZero(bivariateStatistics.rowScoreStatistics.numberOfScores,
                                                       bivariateStatistics.columnScoreStatistics.numberOfScores);

      
      bivariateStatistics.rowScoreStatistics.freqDistMinimumScore = scores.col(rowScoreColumnIndex).minCoeff();
      bivariateStatistics.rowScoreStatistics.freqDistMaximumScore = scores.col(rowScoreColumnIndex).maxCoeff();

      bivariateStatistics.columnScoreStatistics.freqDistMinimumScore = scores.col(columnScoreColumnIndex).minCoeff();
      bivariateStatistics.columnScoreStatistics.freqDistMaximumScore = scores.col(columnScoreColumnIndex).maxCoeff();

      bivariateStatistics.rowScoreStatistics.numberOfExaminees = scores.rows();
      bivariateStatistics.columnScoreStatistics.numberOfExaminees = scores.rows();
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

        bivariateStatistics.rowScoreStatistics.freqDist(rowScoreLocation)++;
        bivariateStatistics.columnScoreStatistics.freqDist(columnScoreLocation)++;
        bivariateStatistics.bivariateFreqDist(rowScoreLocation, columnScoreLocation)++;

        bivariateStatistics.covariance += rowScore * columnScore;
      }

      EquatingRecipes::Structures::Moments rowScoreMoments = EquatingRecipes::Structures::Moments::getScoreMoments(scores.col(rowScoreColumnIndex));
      
      bivariateStatistics.rowScoreStatistics.freqDistMinimumScore = rowScoreMoments.minimumObservedScore;
      bivariateStatistics.rowScoreStatistics.freqDistMaximumScore = rowScoreMoments.maximumObservedScore;
      bivariateStatistics.rowScoreStatistics.momentValues = rowScoreMoments.momentValues;

      EquatingRecipes::Structures::Moments columnScoreMoments = EquatingRecipes::Structures::Moments::getScoreMoments(scores.col(columnScoreColumnIndex));

      bivariateStatistics.columnScoreStatistics.freqDistMinimumScore = columnScoreMoments.minimumObservedScore;
      bivariateStatistics.columnScoreStatistics.freqDistMaximumScore = columnScoreMoments.maximumObservedScore;
      bivariateStatistics.columnScoreStatistics.momentValues = columnScoreMoments.momentValues;

      return bivariateStatistics;
    }
  } // namespace Structures
} // namespace EquatingRecipes