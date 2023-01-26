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

    if (score < minimumScore - scoreIncrement / 2) {
      percRank = 0.0;
    } else if (score < minimumScore + scoreIncrement / 2) {
      percRank = 100 * ((score - (minimumScore - scoreIncrement / 2)) / scoreIncrement) * cumulativeRelativeFreqDist[0];
    } else if (score >= maximumScore + scoreIncrement / 2) {
      percRank = 100.0;
    } else {
      size_t maxScoreLocation = getScoreLocation(maximumScore, minimumScore, scoreIncrement);

      size_t scoreIndex;
      double xStar;

      for (scoreIndex = 1; scoreIndex <= maxScoreLocation; ++scoreIndex) {
        double xStar = getScore(scoreIndex, minimumScore, scoreIncrement);
        if (score < xStar + scoreIncrement / 2) {
          break;
        }
      }

      percRank = 100 * (cumulativeRelativeFreqDist(scoreIndex) +
                        ((score - (xStar - scoreIncrement / 2)) / scoreIncrement) *
                            (cumulativeRelativeFreqDist(scoreIndex) - cumulativeRelativeFreqDist(scoreIndex - 1)));
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
  std::map<double, int> Utilities::getRawScoreFrequencyDistribution(const Eigen::VectorXd& rawScores,
                                                                    const double& minimumScore,
                                                                    const double& maximumScore,
                                                                    const double& scoreIncrement,
                                                                    const bool& includeRawScoresWithZeroFrequency) {
    std::map<double, int> freqDist;

    size_t minimumScoreLocation = Utilities::getScoreLocation(minimumScore, minimumScore, scoreIncrement);
    size_t maximumScoreLocation = Utilities::getScoreLocation(maximumScore, minimumScore, scoreIncrement);

    for (size_t location = minimumScoreLocation; location <= maximumScoreLocation; ++location) {
      double score = Utilities::getScore(location,
                                                 minimumScore,
                                                 scoreIncrement);

      int scoreFreq = static_cast<int>(rawScores.cwiseEqual(score).count());

      if (scoreFreq >= 1 || (scoreFreq == 0 && includeRawScoresWithZeroFrequency)) {
        freqDist[score] = scoreFreq;
      }
    }

    return freqDist;
  }
} // namespace EquatingRecipes