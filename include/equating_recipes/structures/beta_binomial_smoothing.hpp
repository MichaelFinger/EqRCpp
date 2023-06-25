/* 
  From Source: ERutilities.h 
  Original Struct: BB_SMOOTH
  Description: beta-binomial smoothing: 2 or 4 parameter beta with binomial or compound binomial error
*/

#ifndef STRUCTURES_BETA_BINOMIAL_SMOOTHING_HPP
#define STRUCTURES_BETA_BINOMIAL_SMOOTHING_HPP

#include <string>

#include <Eigen/Core>
#include <fmt/core.h>
#include <fmt/format.h>

#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct BetaBinomialSmoothing {
      size_t numberOfItems;                                       /* number of items on test */
      size_t numberOfExaminees;                                   /* sample size */
      size_t numberOfParameters;                                  /* number of parameters (2 or 4) */
      double reliablilty;                                         /* reliability -- almost always Kr20 */
      double lordK;                                               /* Lord's k for approximation of compound binomial */
      Eigen::VectorXd betaParameters = Eigen::VectorXd(4);        /* parameters of true score distribution */
      Eigen::VectorXd rawScoreMoments = Eigen::VectorXd(4);       /* Raw score moments */
      Eigen::VectorXd fittedRawScoreMoments = Eigen::VectorXd(4); /* Fitted raw score moments */
      Eigen::VectorXd trueScoreMoments = Eigen::VectorXd(4);      /* True score moments */
      double likelihoodRatioChiSq;                                /* likelihood ratio chi-square for fitted dist */
      double pearsonChiSq;                                        /* Pearson chi-square for fitted dist */
      size_t numberOfMomentsFit;                                  /* number of moments fit */
      Eigen::VectorXd fittedRawScoreDensity;                      /* fitted raw score dist (proportions, density) */
      Eigen::VectorXd fittedRawScoreCumulativeRelativeFreqDist;   /* cum rel freq dist for fitted dist */
      Eigen::VectorXd fittedRawScorePercentileRankDist;           /* percentile rank dist for fitted dist */
      EquatingRecipes::Structures::UnivariateStatistics rawScoreUnivariateStatistics;

      std::string toString() {
        std::string integerFormat = "6d";

        std::string msg = "Beta-Binomial Smoothing:  ";
        msg.append(fmt::format("{:1d}-parameter beta with ", numberOfParameters));
        msg.append(fmt::format("{}binomial error model", (lordK > 0) ? "compound " : ""));
        msg.append(fmt::format("\n\n       Number of examinees: {:6d}", numberOfExaminees));
        msg.append(fmt::format("\n           Number of items: {:6d}", numberOfItems));
        msg.append(fmt::format("\nReliability (usually KR20): {:12.5f}", reliablilty));
        msg.append(fmt::format("\n                  Lord's k: {:12.5f}", lordK));
        msg.append(fmt::format("\n     Number of momemts fit: {:6d}", numberOfParameters));

        msg.append("\n\n***Parameter Estimates for Beta Distribution***");
        msg.append(fmt::format("\n\n      alpha: {:12.5f}", betaParameters(0)));
        msg.append(fmt::format("\n       beta: {:12.5f}", betaParameters(1)));
        msg.append(fmt::format("\nlower limit: {:12.5f}", betaParameters(2)));
        msg.append(fmt::format("\nupper limit: {:12.5f}", betaParameters(3)));

        msg.append("\n\n***Moments***\n\n              Raw      Fitted-Raw              True\n");
        msg.append(fmt::format("\nMean {:12.5f}    {:12.5f}      {:12.5f}", rawScoreMoments(0), fittedRawScoreMoments(0), trueScoreMoments(0)));
        msg.append(fmt::format("\nS.D. {:12.5f}    {:12.5f}      {:12.5f}", rawScoreMoments(1), fittedRawScoreMoments(1), trueScoreMoments(1)));
        msg.append(fmt::format("\nSkew {:12.5f}    {:12.5f}      {:12.5f}", rawScoreMoments(2), fittedRawScoreMoments(2), trueScoreMoments(2)));
        msg.append(fmt::format("\nKurt {:12.5f}    {:12.5f}      {:12.5f}", rawScoreMoments(3), fittedRawScoreMoments(3), trueScoreMoments(3)));

        msg.append("\n\n***Chi-Square Statistics for Fitted Distribution***");
        msg.append(fmt::format("\n\nLikelihood ratio: {:12.5f}", likelihoodRatioChiSq));
        msg.append(fmt::format("\n         Pearson: {:12.5f}", pearsonChiSq));

        msg.append("\n\n***Frequencies and Proportions***");
        msg.append("\n\nScore   freq  fitted-freq          prop  fitted-prop  fitted-crfd   fitted-prd\n");

        for (size_t scoreLocation = 0; scoreLocation <= numberOfItems; scoreLocation++) {
          msg.append(fmt::format("\n{:5d} {:6.0f} {:12.5f}       {:12.5f} {:12.5f} {:12.5f} {:12.5f}",
                                 scoreLocation,
                                 rawScoreUnivariateStatistics.freqDistDouble(scoreLocation),
                                 static_cast<double>(numberOfExaminees) * fittedRawScoreDensity(scoreLocation),
                                 rawScoreUnivariateStatistics.relativeFreqDist(scoreLocation),
                                 fittedRawScoreDensity(scoreLocation),
                                 fittedRawScoreCumulativeRelativeFreqDist(scoreLocation),
                                 fittedRawScorePercentileRankDist(scoreLocation)));
        }

        msg.append("\n\n");

        for (size_t index = 1; index <= 78; index++) {
          msg.append("*");
        }

        return msg;
      }

      std::string toLongString() {
        std::string integerFormat = "6d";

        std::string msg = "Beta-Binomial Smoothing:  ";
        msg.append(fmt::format("{:1d}-parameter beta with ", numberOfParameters));
        msg.append(fmt::format("{}binomial error model", (lordK > 0) ? "compound " : ""));
        msg.append(fmt::format("\n\n       Number of examinees: {:16d}", numberOfExaminees));
        msg.append(fmt::format("\n           Number of items: {:6d}", numberOfItems));
        msg.append(fmt::format("\nReliability (usually KR20): {:30.15f}", reliablilty));
        msg.append(fmt::format("\n                  Lord's k: {:30.15f}", lordK));
        msg.append(fmt::format("\n     Number of momemts fit: {:6d}", numberOfParameters));

        msg.append("\n\n***Parameter Estimates for Beta Distribution***");
        msg.append(fmt::format("\n\n      alpha: {:30.15f}", betaParameters(0)));
        msg.append(fmt::format("\n       beta: {:30.15f}", betaParameters(1)));
        msg.append(fmt::format("\nlower limit: {:30.15f}", betaParameters(2)));
        msg.append(fmt::format("\nupper limit: {:30.15f}", betaParameters(3)));

        msg.append("\n\n***Moments***\n\n              Raw      Fitted-Raw              True\n");
        msg.append(fmt::format("\nMean {:30.15f}    {:30.15f}      {:30.15f}", rawScoreMoments(0), fittedRawScoreMoments(0), trueScoreMoments(0)));
        msg.append(fmt::format("\nS.D. {:30.15f}    {:30.15f}      {:30.15f}", rawScoreMoments(1), fittedRawScoreMoments(1), trueScoreMoments(1)));
        msg.append(fmt::format("\nSkew {:30.15f}    {:30.15f}      {:30.15f}", rawScoreMoments(2), fittedRawScoreMoments(2), trueScoreMoments(2)));
        msg.append(fmt::format("\nKurt {:30.15f}    {:30.15f}      {:30.15f}", rawScoreMoments(3), fittedRawScoreMoments(3), trueScoreMoments(3)));

        msg.append("\n\n***Chi-Square Statistics for Fitted Distribution***");
        msg.append(fmt::format("\n\nLikelihood ratio: {:30.15f}", likelihoodRatioChiSq));
        msg.append(fmt::format("\n         Pearson: {:30.15f}", pearsonChiSq));

        msg.append("\n\n***Frequencies and Proportions***");
        msg.append("\n\nScore             freq                    fitted-freq                                 prop                    fitted-prop                    fitted-crfd                    fitted-prd\n");

        for (size_t scoreLocation = 0; scoreLocation <= numberOfItems; scoreLocation++) {
          msg.append(fmt::format("\n{:5d} {:16.0f} {:30.15f}       {:30.15f} {:30.15f} {:30.15f} {:30.15f}",
                                 scoreLocation,
                                 rawScoreUnivariateStatistics.freqDistDouble(scoreLocation),
                                 static_cast<double>(numberOfExaminees) * fittedRawScoreDensity(scoreLocation),
                                 rawScoreUnivariateStatistics.relativeFreqDist(scoreLocation),
                                 fittedRawScoreDensity(scoreLocation),
                                 fittedRawScoreCumulativeRelativeFreqDist(scoreLocation),
                                 fittedRawScorePercentileRankDist(scoreLocation)));
        }

        msg.append("\n\n");

        for (size_t index = 1; index <= 78; index++) {
          msg.append("*");
        }

        return msg;
      }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif