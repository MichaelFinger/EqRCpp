/* 
  From Source: BetaBinomial.h, BetaBinomial.c
  Original Struct: ITEM_SPEC
  Description: Contains functions used for:
      (a) beta binomial smoothing--there are three types:
            2 parameter beta, binomial errors
            4 parameter beta, binomial errors
            4 parameter bete, compound binomial errorsin computing 
      (b) equating beta-binomial smoothed x to scale of 
          beta-binomial smoothed y
*/

#ifndef BETA_BINOMIAL_HPP
#define BETA_BINOMIAL_HPP

#include <Eigen/Dense>

#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/score_statistics.hpp>

namespace EquatingRecipes {
  class BetaBinomial {
  public:
    EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomialSmoothing(const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsX,
                                                                             const size_t& numberOfParameters,
                                                                             const double& reliability) {
      EquatingRecipes::Structures::BetaBinomialSmoothing results = smoothBetaBinomial(univariateStatisticsX.numberOfExaminees,
                                                                                      univariateStatisticsX.numberOfScores - 1,
                                                                                      univariateStatisticsX.freqDistDouble,
                                                                                      univariateStatisticsX.momentValues,
                                                                                      numberOfParameters,
                                                                                      reliability);

      return results;
    }

    EquatingRecipes::Structures::EquatedRawScoreResults randomGroupsEquipercentileEquating(const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsNewForm,
                                                                                           const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsOldForm,
                                                                                           const EquatingRecipes::Structures::BetaBinomialSmoothing& betaBinomialSmoothingNewForm,
                                                                                           const EquatingRecipes::Structures::BetaBinomialSmoothing& betaBinomialSmoothingOldForm,
                                                                                           size_t bootstrapReplicationNumber) {
      EquatingRecipes::Structures::EquatedRawScoreResults results;
      return results;
    }

  private:
    EquatingRecipes::Structures::BetaBinomialSmoothing smoothBetaBinomial(const size_t& numberOfExaminees,
                                                                          const size_t& numberOfItems,
                                                                          const Eigen::VectorXd& freqDist,
                                                                          const Eigen::VectorXd& rawScoreNoments,
                                                                          const size_t& numberOfParameters,
                                                                          const double& reliability) {
      EquatingRecipes::Structures::BetaBinomialSmoothing results;

      results.numberOfExaminees = numberOfExaminees;
      results.numberOfItems = numberOfItems;
      results.rawScoreMoments = rawScoreNoments;
      results.numberOfParameters = numberOfParameters;
      results.reliablilty = reliability;

      results.fittedRawScoreDist.setZero(numberOfItems);
      results.fittedRawScoreCumulativeRelativeFreqDist.setZero(numberOfItems);
      results.fittedRawScorePercentileRankDist.setZero(numberOfItems);

      if (numberOfParameters == 2) {
        twoParameterBetaSmooting(freqDist, results);
      } else {
        if (reliability > 0.0) {
          results.lordK = std::max(calculateLordK(reliability, numberOfItems, rawScoreNoments), 0.0);
        } else {
          results.lordK = 0.0;
        }

        fourParameterBetaSmooting(freqDist, results);
      }

      return results;
    }

    short fourParameterBetaSmooting(const Eigen::VectorXd& rawScoreFreqDist,
                                    EquatingRecipes::Structures::BetaBinomialSmoothing& betaFitResults) {
      short numberOfMomentsFit = 1;                  /* number of moments fit to obtain estimates */
      size_t numberOfExaminees;                      /* number of persons */
      Eigen::VectorXd smoothedProportions;           /* smoothed proportions */
      Eigen::VectorXd smoothedFrequencies;           /* smoothed frequencies */
      Eigen::VectorXd noncentralTrueScoreMoments(4); /* non-central true score moments */

      size_t numberOfItems;                /* number of items on test */
      size_t numberOfCategoriesLRChiSqare; /* number of categories used for chi-square */
      size_t numberOfCategoriesLRPValue;   /* number of categories used for chi-square p-value */

      numberOfExaminees = betaFitResults.numberOfExaminees;
      numberOfItems = betaFitResults.numberOfItems;
      smoothedProportions = betaFitResults.fittedRawScoreDist;
      smoothedFrequencies.setZero(numberOfItems);

      /* calculate true score moments */

      noncentralTrueScoreMoments = betaMoments(numberOfItems,
                                               betaFitResults.likelihoodRatioChiSq,
                                               betaFitResults.rawScoreMoments,
                                               betaFitResults.trueScoreMoments);

      /* calculate parameters of beta true score distribution  */

      numberOfMomentsFit = calculateBetaParameters(numberOfItems,
                                                   betaFitResults.trueScoreMoments,
                                                   noncentralTrueScoreMoments,
                                                   betaFitResults.betaParameters);

      /* calculate observed score density */

      if (observedDensity(numberOfItems + 1,
                          numberOfExaminees,
                          betaFitResults.betaParameters,
                          smoothedFrequencies)) {
        /* Assign uniform distribution and exit */
        smoothedFrequencies.setConstant(1.0 / static_cast<double>(numberOfItems + 1));
        numberOfMomentsFit = 0;
      } else {
        /* adjust for compound binomial if k != 0 */
        if (betaFitResults.lordK != 0.0) {
          observedDensityK(numberOfItems,
                           betaFitResults.lordK,
                           smoothedFrequencies);
        }

        /* compute smoothed proportions */
        smoothedProportions = smoothedFrequencies / static_cast<double>(numberOfExaminees);

        /* compute chi-square values */
        betaFitResults.likelihoodRatioChiSq = likelihoodRatioChiSquare(numberOfItems + 1,
                                                                       rawScoreFreqDist,
                                                                       smoothedFrequencies,
                                                                       numberOfCategoriesLRChiSqare);
        betaFitResults.pValueChiSq = chiSquarePValue(numberOfItems + 1,
                                                     0.0,
                                                     rawScoreFreqDist,
                                                     smoothedFrequencies,
                                                     numberOfCategoriesLRPValue);

        /* calculate fitted observed score moments */
        EquatingRecipes::ScoreStatistics scoreStatistics;
        EquatingRecipes::Structures::Moments moments = scoreStatistics.momentsFromScoreFrequencies(smoothedFrequencies,
                                                                                                   0,
                                                                                                   numberOfItems,
                                                                                                   1.0);
        betaFitResults.fittedRawScoreMoments = moments.momentValues;
      }

      betaFitResults.numberOfMomentsFit = numberOfMomentsFit;

      return numberOfMomentsFit >= 0 ? 0 : -1;
    }

    short twoParameterBetaSmooting(const Eigen::VectorXd& rawScoreFreqDist,
                                   EquatingRecipes::Structures::BetaBinomialSmoothing& betaFitResults) {
      short numberOfMomentsFit = 0;

      size_t numberOfItems;                /* number of items on test */
      size_t numberOfCategoriesLRChiSqare; /* number of categories used for chi-square */
      size_t numberOfCategoriesLRPValue;   /* number of categories used for chi-square p-value */

      size_t numberOfExaminees;                      /* number of persons */
      Eigen::VectorXd smoothedProportions;           /* pointer to smoothed proportions */
      Eigen::VectorXd smoothedFrequencies;           /* pointer to smoothed frequencies */
      Eigen::VectorXd parameterEstimates;            /* pointer to model parameter estimates */
      Eigen::VectorXd noncentralTrueScoreMoments(4); /* non-central true score moments */

      numberOfExaminees = betaFitResults.numberOfExaminees;
      numberOfItems = betaFitResults.numberOfItems;
      smoothedProportions = betaFitResults.fittedRawScoreDist;
      smoothedFrequencies.setZero(numberOfItems);
      betaFitResults.lordK = 0.0; /* binomial error distribution */

      /* calculate true score moments */
      noncentralTrueScoreMoments = betaMoments(numberOfItems,
                                               betaFitResults.lordK,
                                               betaFitResults.rawScoreMoments,
                                               betaFitResults.trueScoreMoments);

      /* calculate parameters of beta true score distribution  */

      parameterEstimates = betaFitResults.betaParameters;
      parameterEstimates(2) = 0.0;
      parameterEstimates(3) = 1.0;

      if (!estimateNegativeHyperGeometricParameters(numberOfItems,
                                                    betaFitResults.trueScoreMoments,
                                                    parameterEstimates)) {
        /* if two moments not fit, set alpha=1.0 and fit mean */
        parameterEstimates(2) = 0.0;
        parameterEstimates(3) = 1.0;
        parameterEstimates(0) = 1.0;
        parameterEstimates(1) = (static_cast<double>(numberOfItems) / betaFitResults.rawScoreMoments(0)) - parameterEstimates(0);

        if (parameterEstimates(1) > 0.0) {
          numberOfMomentsFit = 1;
        } else {
          /* if mean not fit return uniform distribution on 0 to 1.0 */
          parameterEstimates(1) = 1.0;
          numberOfMomentsFit = 0;
        }
      } else {
        numberOfMomentsFit = 2;
      }

      betaFitResults.betaParameters = parameterEstimates;

      /* calculate observed score density */

      if (observedDensity(numberOfItems + 1,
                          numberOfExaminees,
                          betaFitResults.betaParameters,
                          smoothedFrequencies)) {
        /* Assign uniform distribution and exit */
        smoothedProportions.setConstant(1.0 / static_cast<double>(numberOfItems + 1));
        numberOfMomentsFit = 0;
      } else {
        /* compute smoothed proportions */

        smoothedProportions = smoothedFrequencies / static_cast<double>(numberOfExaminees);

        /* compute chi-square values */
        betaFitResults.likelihoodRatioChiSq = likelihoodRatioChiSquare(numberOfItems + 1,
                                                                       rawScoreFreqDist,
                                                                       smoothedFrequencies,
                                                                       numberOfCategoriesLRChiSqare);

        betaFitResults.pValueChiSq = chiSquarePValue(numberOfItems + 1,
                                                     0.0,
                                                     rawScoreFreqDist,
                                                     smoothedFrequencies,
                                                     numberOfCategoriesLRPValue);

        /* calculate fitted observed score moments */

        ScoreStatistics scoreStatistics;
        EquatingRecipes::Structures::Moments moments = scoreStatistics.momentsFromScoreFrequencies(smoothedProportions,
                                                                                                   0.0,
                                                                                                   numberOfItems,
                                                                                                   1.0);
        betaFitResults.fittedRawScoreMoments = moments.momentValues;
      }

      betaFitResults.numberOfMomentsFit = numberOfMomentsFit;

      if (betaFitResults.numberOfMomentsFit < 0) {
        return -1;
      } else {
        return 0;
      }
    }

    // void Print_BB(FILE* fp, char tt[], struct USTATS* x, struct BB_SMOOTH* s);
    // void Print_RB(FILE* fp, char tt[], struct PDATA* inall, struct ERAW_RESULTS* r);

    short calculateBetaParameters(const size_t& numberOfItems, const Eigen::VectorXd& moments, const Eigen::VectorXd& nctmoment, Eigen::VectorXd& parameterEstimates) {
      short result = 0;
      return result;
    }

    bool estimateNegativeHyperGeometricParameters(const size_t& numberOfItems, const Eigen::VectorXd& moments, Eigen::VectorXd& parameterEstimates) {
      bool result = false;
      return result;
    }

    short calculateBetaParameterMM(const size_t& numberOfItems, const Eigen::VectorXd& moments, Eigen::VectorXd& parameterEstimates) {
      short result = 0;
      return result;
    }

    short calculateBetaParameterLS(const size_t& numberOfItems, const Eigen::VectorXd& moments, const Eigen::VectorXd& nctmoment, Eigen::VectorXd& parameterEstimates) {
      short result = 0;
      return result;
    }

    short findUpper(const Eigen::VectorXd& parameterEstimates, const Eigen::VectorXd& tmoment) {
      short result = 0;
      return result;
    }

    short findLower(const Eigen::VectorXd& parameterEstimates, const Eigen::VectorXd& tmoment) {
      short result = 0;
      return result;
    }

    double calcKurtosis(const Eigen::VectorXd& parameterEstimates) {
      double result = 0;
      return result;
    }

    short kurtosisFunctionD(const double& x, const size_t& numberOfItems, const Eigen::VectorXd& parameterEstimates, const Eigen::VectorXd& tmoment, const Eigen::VectorXd& nctmoment, Eigen::VectorXd& kurtosis) {
      short result = 0;
      return result;
    }

    short kurtosisFunctionUpper(const double& x, const size_t& numberOfItems, const Eigen::VectorXd& parameterEstimates, const Eigen::VectorXd& tmoment, const Eigen::VectorXd& nctmoment, Eigen::VectorXd& kurtosis) {
      short result = 0;
      return result;
    }

    Eigen::VectorXd betaMoments(const size_t& n, const double& k, const Eigen::VectorXd& rmoment, const Eigen::VectorXd& tmoment) {
      Eigen::VectorXd nctmoment;

      return nctmoment;
    }

    bool observedDensity(const size_t& n, const size_t& nsamp, const Eigen::VectorXd& beta, Eigen::VectorXd& scounts) {
      short result = 0;
      return result;
    }

    void observedDensityK(const size_t& n, const double& k, Eigen::VectorXd& scounts) {}

    void calculateM3(const size_t& n, const Eigen::VectorXd& beta, Eigen::MatrixXd& m3) {}

    short calculateM24(const size_t& n, const double& para, Eigen::VectorXd& m) {
      short result = 0;
      return result;
    }

    void calculateM15(const size_t& n, const double& para, Eigen::VectorXd& m) {}

    double likelihoodRatioChiSquare(const size_t& n, const Eigen::VectorXd& rawc, const Eigen::VectorXd& fitc, size_t& ncat) {
      double result = 0;
      return result;
    }

    double chiSquarePValue(const size_t& n, const double& minexpcount, const Eigen::VectorXd& rawc, const Eigen::VectorXd& fitc, size_t& ncat) {
      double result = 0;
      return result;
    }

    double calculateLordK(const double& kr20, const size_t numberOfItems, const Eigen::VectorXd& rmoment) {
      double result = 0;
      return result;
    }
  };
} // namespace EquatingRecipes

#endif