#include <fmt/core.h>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  size_t Utilities::getScoreLocation(const double& score,
                                     const double& minimumScore,
                                     const double& scoreIncrement) {
    size_t location = static_cast<size_t>((score - minimumScore) / scoreIncrement + 0.5);

    return location;
  }

  size_t Utilities::numberOfScores(const double& minimumScore,
                                   const double& maximumScore,
                                   const double& scoreIncrement) {
    size_t nScores = getScoreLocation(maximumScore, minimumScore, scoreIncrement) + 1;

    return nScores;
  }

  double Utilities::getScore(const size_t& scoreLocation,
                             const double& minimumScore,
                             const double& scoreIncrement) {
    double scoreValue = minimumScore + static_cast<double>(scoreLocation) * scoreIncrement;

    return scoreValue;
  }

  Eigen::VectorXd Utilities::cumulativeRelativeFreqDist(const double& minimumScore,
                                                        const double& maximumScore,
                                                        const double& scoreIncrement,
                                                        const Eigen::VectorXd& relativeFreqDist) {
    size_t maximumScoreLocation = Utilities::getScoreLocation(maximumScore,
                                                              minimumScore,
                                                              scoreIncrement);

    Eigen::VectorXd cumRelFreqDist = Eigen::VectorXd::Zero(maximumScoreLocation + 1);

    cumRelFreqDist(0) = relativeFreqDist(0);

    for (size_t freqDistIndex = 1; freqDistIndex < maximumScoreLocation; ++freqDistIndex) {
      cumRelFreqDist(freqDistIndex) = cumRelFreqDist(freqDistIndex - 1) + relativeFreqDist(freqDistIndex);
    }

    return cumRelFreqDist;
  }

  double Utilities::percentileRank(const double& minimumScore,
                                   const double& maximumScore,
                                   const double& scoreIncrement,
                                   const Eigen::VectorXd& cumulativeRelativeFreqDist,
                                   const double& score) {
    double percRank;

    if (score < minimumScore - scoreIncrement / 2.0) {
      percRank = 0.0;
    } else if (score < minimumScore + scoreIncrement / 2.0) {
      percRank = 100 * ((score - (minimumScore - scoreIncrement / 2.0)) / scoreIncrement) * cumulativeRelativeFreqDist[0];
    } else if (score >= maximumScore + scoreIncrement / 2.0) {
      percRank = 100.0;
    } else {
      size_t maxScoreLocation = Utilities::getScoreLocation(maximumScore, minimumScore, scoreIncrement);

      size_t scoreLocation;
      double xStar;

      for (scoreLocation = 1; scoreLocation <= maxScoreLocation; ++scoreLocation) {
        xStar = Utilities::getScore(scoreLocation, minimumScore, scoreIncrement);
        if (score < xStar + scoreIncrement / 2.0) {
          break;
        }
      }

      percRank = 100 * (cumulativeRelativeFreqDist(scoreLocation - 1) +
                        ((score - (xStar - scoreIncrement / 2.0)) / scoreIncrement) *
                            (cumulativeRelativeFreqDist(scoreLocation) - cumulativeRelativeFreqDist(scoreLocation - 1)));
    }

    return percRank;
  }

  double Utilities::interpolate(const double& score,
                                const int& numberOfScoreCategories,
                                const Eigen::VectorXd& frequencies) {
    int truncatedScore = static_cast<int>(score);
    size_t maximumScore = numberOfScoreCategories - 1;
    double value;

    if (score <= 0.0) { /* extrapolate low */
      value = frequencies(0) + score * (frequencies(1) - frequencies(0));
    } else if (score >= static_cast<double>(maximumScore)) { /* extrapolate high */
      value = frequencies(maximumScore) + (score - maximumScore) * (frequencies(maximumScore) - frequencies(maximumScore - 1));
    } else { /* interpolate */
      value = frequencies(truncatedScore) + (score - truncatedScore) * (frequencies(truncatedScore + 1) - frequencies(truncatedScore));
    }

    /* an extrapolated value can be less than 0.  Hence ... */
    double result = (value > 0) ? value : 1.0e-10;

    return result;
  }

  //----------------------------------------------------------------------------------------------------
  // Custom Function Written for EqRCpp
  //----------------------------------------------------------------------------------------------------
  Eigen::VectorXi Utilities::getRawScoreFrequencyDistribution(const Eigen::VectorXd& rawScores,
                                                              const double& minimumScore,
                                                              const double& maximumScore,
                                                              const double& scoreIncrement,
                                                              const bool& includeRawScoresWithZeroFrequency) {
    size_t minimumScoreLocation = Utilities::getScoreLocation(minimumScore, minimumScore, scoreIncrement);
    size_t maximumScoreLocation = Utilities::getScoreLocation(maximumScore, minimumScore, scoreIncrement);
    size_t numberOfScores = Utilities::numberOfScores(minimumScore, maximumScore, scoreIncrement);

    Eigen::VectorXi freqDist = Eigen::VectorXi::Zero(numberOfScores);

    for (size_t scoreLocation = 0; scoreLocation < numberOfScores; ++scoreLocation) {
      double score = Utilities::getScore(scoreLocation,
                                         minimumScore,
                                         scoreIncrement);

      int scoreFreq = static_cast<int>(rawScores.cwiseEqual(score).count());

      if (scoreFreq >= 1 || (scoreFreq == 0 && includeRawScoresWithZeroFrequency)) {
        freqDist(scoreLocation) = scoreFreq;
      }
    }

    return freqDist;
  }

  std::string Utilities::vectorXiToString(const Eigen::VectorXi& vec,
                                          const bool& asRow) {
    std::string value = "";
    std::string sep = asRow ? ", " : "\n";

    value.append(fmt::format("{}", vec(0)));

    for (Eigen::Index index = 1; index < vec.size(); ++index) {
      value.append(fmt::format("{}{}", sep, vec(index)));
    }

    return value;
  }

  std::string Utilities::vectorXdToString(const Eigen::VectorXd& vec,
                                          const bool& asRow) {
    std::string value = "";
    std::string sep = asRow ? ", " : "\n";

    value.append(fmt::format("{}", vec(0)));

    for (Eigen::Index index = 1; index < vec.size(); ++index) {
      value.append(fmt::format("{}{}", sep, vec(index)));
    }

    return value;
  }

  std::string Utilities::matrixXiToString(const Eigen::MatrixXi& mat) {
    std::string value = "";

    for (Eigen::Index rowIndex = 0; rowIndex < mat.rows(); ++rowIndex) {
      std::string rowStr = fmt::format("{}", mat(rowIndex, 0));

      for (Eigen::Index columnIndex = 1; columnIndex < mat.cols(); ++columnIndex) {
        rowStr.append(fmt::format(", {}", mat(rowIndex, columnIndex)));
      }

      value.append(fmt::format("{}\n", rowStr));
    }

    return value;
  }

  std::string Utilities::matrixXdToString(const Eigen::MatrixXd& mat) {
    std::string value = "";

    for (Eigen::Index rowIndex = 0; rowIndex < mat.rows(); ++rowIndex) {
      std::string rowStr = fmt::format("{}", mat(rowIndex, 0));

      for (Eigen::Index columnIndex = 1; columnIndex < mat.cols(); ++columnIndex) {
        rowStr.append(fmt::format(", {}", mat(rowIndex, columnIndex)));
      }

      value.append(fmt::format("{}\n", rowStr));
    }

    return value;
  }

  Eigen::VectorXd Utilities::percentileRanks(const double& minimumScore,
                                             const double& maximumScore,
                                             const double& scoreIncrement,
                                             const Eigen::VectorXd& cumulativeRelativeFreqDist) {
    size_t numberOfScores = Utilities::numberOfScores(minimumScore,
                                                      maximumScore,
                                                      scoreIncrement);

    Eigen::VectorXd percRanks = Eigen::VectorXd::Zero(numberOfScores);

    for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
      double score = Utilities::getScore(scoreLocation,
                                         minimumScore,
                                         scoreIncrement);

      percRanks(scoreLocation) = Utilities::percentileRank(minimumScore,
                                                           maximumScore,
                                                           scoreIncrement,
                                                           cumulativeRelativeFreqDist,
                                                           score);
    }

    return percRanks;
  }
} // namespace EquatingRecipes