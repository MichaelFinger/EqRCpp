#include <fmt/core.h>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
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

      bivariateStatistics.rowScoreStatistics = EquatingRecipes::Structures::UnivariateStatistics::buildFromScores(scores.col(rowScoreColumnIndex),
                                                                                                                  minimumRowScore,
                                                                                                                  maximumRowScore,
                                                                                                                  rowScoreIncrement,
                                                                                                                  rowScoreId);

      bivariateStatistics.columnScoreStatistics = EquatingRecipes::Structures::UnivariateStatistics::buildFromScores(scores.col(columnScoreColumnIndex),
                                                                                                                     minimumColumnScore,
                                                                                                                     maximumColumnScore,
                                                                                                                     columnScoreIncrement,
                                                                                                                     columnScoreId);

      bivariateStatistics.bivariateFreqDist.setZero(bivariateStatistics.rowScoreStatistics.numberOfScores,
                                                    bivariateStatistics.columnScoreStatistics.numberOfScores);
      bivariateStatistics.bivariateFreqDistDouble.setZero(bivariateStatistics.rowScoreStatistics.numberOfScores,
                                                          bivariateStatistics.columnScoreStatistics.numberOfScores);
      bivariateStatistics.bivariateProportions.setZero(bivariateStatistics.rowScoreStatistics.numberOfScores,
                                                       bivariateStatistics.columnScoreStatistics.numberOfScores);

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

      bivariateStatistics.bivariateFreqDistDouble = bivariateStatistics.bivariateFreqDist.cast<double>();
      bivariateStatistics.bivariateProportions = bivariateStatistics.bivariateFreqDistDouble / static_cast<double>(bivariateStatistics.numberOfExaminees);

      bivariateStatistics.covariance = (bivariateStatistics.covariance / static_cast<double>(bivariateStatistics.numberOfExaminees)) -
                                       (bivariateStatistics.rowScoreStatistics.momentValues(0) * bivariateStatistics.columnScoreStatistics.momentValues(0));

      bivariateStatistics.correlation = bivariateStatistics.covariance /
                                        (bivariateStatistics.rowScoreStatistics.momentValues(1) * bivariateStatistics.columnScoreStatistics.momentValues(1));

      return bivariateStatistics;
    }

    std::string BivariateStatistics::toString() {
      std::string value = "Univariate Statistics: Row Scores\n";
      value.append(this->rowScoreStatistics.toString());
      value.append("\n");

      value.append("Univariate Statistics: Column Scores\n");
      value.append(this->columnScoreStatistics.toString());
      value.append("\n");

      value.append(fmt::format("Number of Examinees: {}\n", this->numberOfExaminees));

      value.append(fmt::format("Bivariate Frequency Distribution:\n{}\n", EquatingRecipes::Utilities::matrixXiToString(this->bivariateFreqDist)));
      value.append(fmt::format("Bivariate Frequency Distribution (double):\n{}\n", EquatingRecipes::Utilities::matrixXdToString(this->bivariateFreqDistDouble)));
      value.append(fmt::format("Bivariate Proportions:\n{}\n", EquatingRecipes::Utilities::matrixXdToString(this->bivariateProportions)));
      value.append(fmt::format("Sum of Bivariate Proportions: {}\n", this->bivariateProportions.sum()));
      value.append(fmt::format("Covariance: {}\n", this->covariance));
      value.append(fmt::format("Correlation: {}\n", this->correlation));

      return value;
    }
  } // namespace Structures
} // namespace EquatingRecipes