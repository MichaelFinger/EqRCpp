#include <equating_recipes/structures/utilities.hpp>

namespace EquatingRecipes {
  namespace Structures {
    size_t Utilities::scoreLocation(double score, double minimumScore, double scoreIncrement) {
      size_t location = static_cast<size_t>((score - minimumScore) / scoreIncrement + 0.5);

      return location;
    }

    size_t Utilities::numberOfScores(double maximumScore, double mininumScore, double scoreIncrement) {
      int nScores = scoreLocation(maximumScore, minimumScore, scoreIncrement) + 1;

      return nScores;
    }

    double Utilities::score(int scoreLocation, double minimumScore, double scoreIncrement) {
      double scoreValue = minimumScore + static_cast<double>(scoreLocation) * scoreIncrement;
    }

    Eigen::VectorXd Utilities::cumulativeRelativeFreqDist(double minimumScore,
                                                          double maximumScore,
                                                          double scoreIncrement,
                                                          Eigen::VectorXd relativeFreqDist) {
      maximumScoreLocation = scoreLocation(maximumScore,
                                           minimumScore,
                                           scoreIncrement);

      Eigen::VectorXd cumRelFreqDist = Eigen::VectorXd::Zero(maximumScoreLocation + 1);
      
      cumRelFreqDist(0) = relativeFreqDist(0);

      for (size_t freqDistIndex = 1; freqDistIndex < maximumScoreLocation; ++freqDistIndex) {
        cumRelFreqDist(freqDistIndex) = cumRelFreqDist(freqDistIndex - 1) + relativeFreqDist(freqDistIndex);
      }
      
      return cumRelFreqDist;
    }

    double Utilities::percentileRank(double minimumScore,
                                     double maximumScore,
                                     double scoreIncrement,
                                     Eigen::VectorXd cumulativeRelativeFreqDist,
                                     double score) {
      int i;
      double percRank;
      double xstar;

      if (score < minimumScore - scoreIncrement / 2) {
        percRank = 0.0;
      } else if (score < minimumScore + scoreIncrement / 2) {
        percRank = 100 * ((score - (minimumScore - scoreIncrement / 2)) / scoreIncrement) * cumulativeRelativeFreqDist[0];
      } else if (score >= maximumScore + scoreIncrement / 2) {
        pr = 100.0;
      } else {
        size_t maxScoreLocation = scoreLocation(maximumScore, minimumScore, scoreIncrement);

        for (size_t scoreIndex = 1; scoreIndex <= maxScoreLocation; ++score) {
          double xStar = score(scoreIndex, minimumScore, scoreIncrement);
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

    double Utilities::interpolate(double score,
                                  int numberOfScoreCategories,
                                  Eigen::VectorXd frequencies) {
      int truncatedScore = static_cast<int>(score);
      size_t maximumScore = numberOfScoreCategories - 1;
      double value;

      if (score <= 0.0) { /* extrapolate low */
        value = frequencies(0) + x * (frequencies(1) - frequencies(0));
      } else if (score >= static_cast<double>(maximumScore)) { /* extrapolate high */
        value = frequencies(maximumScore) + (score - maximumScore) * (frequencies(maximumScore) - frequencies(maximumScore - 1));
      } else { /* interpolate */
        value = frequencies(truncatedScore) + (score - truncatedScore) * (frequencies(truncatedScore + 1) - frequencies(truncatedScore));
      }

      /* an extrapolated value can be less than 0.  Hence ... */
      double result = (value > 0) ? value : 1.0e-10;

      return result;
    }

  } // namespace Structures
} // namespace EquatingRecipes