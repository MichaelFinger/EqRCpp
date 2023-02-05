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

#include <iostream>
#include <string>
#include <fmt/core.h>

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

    // private:
    size_t maximumNumberOfIterations = 20; /* maximum number of iterations for computing upper limit */
    double KACC = 0.001;                   /* accuracy for computing lower limit in CalcBetaParaLS */

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

      betaMoments(numberOfItems,
                  betaFitResults.likelihoodRatioChiSq,
                  betaFitResults.rawScoreMoments,
                  betaFitResults.trueScoreMoments,
                  noncentralTrueScoreMoments);

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
      betaMoments(numberOfItems,
                  betaFitResults.lordK,
                  betaFitResults.rawScoreMoments,
                  betaFitResults.trueScoreMoments,
                  noncentralTrueScoreMoments);

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

    short calculateBetaParameters(const size_t& numberOfItems,
                                  const Eigen::VectorXd& trueScoreMoments,
                                  const Eigen::VectorXd& noncentralTrueScoreMoments,
                                  Eigen::VectorXd& parameterEstimates) {
      short numberOfMomentsFit = 0;

      /* 	Only try to fit more than one moment if the true score
		  s.d. is greater than zero */
      if (trueScoreMoments(1) > 0.0) {
        /* try to fit all four moments */
        if (calculateBetaParameterMM(numberOfItems,
                                     trueScoreMoments,
                                     parameterEstimates)) {
          numberOfMomentsFit = 4;
        }
        /* if four moments not fit try to fit three moments 
        such that the squared difference in the estimated 
        and observed kurtosis is as small as possible */
        else if (calculateBetaParameterLS(numberOfItems,
                                          trueScoreMoments,
                                          noncentralTrueScoreMoments,
                                          parameterEstimates)) {
          numberOfMomentsFit = 3;
        }
        /* if three moments not fit set lower = 0.0 and upper = 1.0 
        and fit first two moments */
        else {
          parameterEstimates(2) = 0.0;
          parameterEstimates(3) = 1.0;

          if (estimateNegativeHyperGeometricParameters(numberOfItems,
                                                       trueScoreMoments,
                                                       parameterEstimates)) {
            numberOfMomentsFit = 2;
          }
        }
      } else {
        /* if two moments not fit, set alpha=1.0 and fit mean */
        parameterEstimates(2) = 0.0;
        parameterEstimates(3) = 1.0;
        parameterEstimates(0) = 1.0;
        parameterEstimates(1) = (static_cast<double>(numberOfItems) / trueScoreMoments(0)) - parameterEstimates(0);

        if (parameterEstimates(1) > 0.0) {
          numberOfMomentsFit = 1;
        } else {
          /* if mean not fit return uniform distribution on 0 to 1.0 */
          parameterEstimates(1) = 1.0;
          numberOfMomentsFit = 0;
        }
      }

      return numberOfMomentsFit;
    }

    /*
        Estimates parameters of 4-parameter beta binomial model
        by method of moments. Formula used is from pages 40 and 41 of
        Johnson, N.L., & Kotz, S.  Continuous univariate distributions - volume 2.
        NOTE: formula for alpha and beta on page 41 is wrong but the correct
              formula can be obtained from expression for "r" on page 41 (which
              is correct) and expression for kurtosis on page 40.

        Input
          nitems = number of items on test
          moment = true score moments

        Output
          para = method of moments estimates of parameters

        Function calls other than C or NR utilities: None
                                                      
        B. A. Hanson with updates by R. L. Brennan

        Date of last revision: 6/30/08 

      */
    bool calculateBetaParameterMM(const size_t& numberOfItems, const Eigen::VectorXd& moments, Eigen::VectorXd& parameterEstimates) {
      bool result = false;

      double b2 = moments(3);
      double b1 = std::pow(moments(2), 2);

      double r = 6.0 * (b2 - b1 - 1.0);
      r /= 6.0 + 3.0 * b1 - 2.0 * b2;

      double rad = (r + 2.0) * (r + 3.0) * b2;
      rad -= 3.0 * (r + 1.0) * (r - 6.0);
      rad = (24.0 * (r + 1.0)) / rad;

      if (rad > 1.0) {
        /* cannot compute alpha and beta */
        return false;
      }

      rad = std::sqrt(1.0 - rad);

      r /= 2.0;
      /* compute values of alpha and beta */
      double a;
      double b;

      if (moments(2) < 0.0) {
        a = r * (1.0 + rad);
        b = r * (1.0 - rad);
      } else {
        a = r * (1.0 - rad);
        b = r * (1.0 + rad);
      }

      if (a <= 0.0 || b <= 0.0) {
        return false;
      }

      /* compute values of lower and upper limit parameters */

      double bma = (a + b) * std::sqrt(a + b + 1);
      bma *= moments(1);
      bma /= std::sqrt(a * b);

      /* compute lower limit */
      double l = -bma * (a / (a + b));
      l += moments(0);

      /* compute upper limit */
      double u = bma + l;

      // /* assign values for output */
      parameterEstimates(0) = a;
      parameterEstimates(1) = b;
      parameterEstimates(2) = l / static_cast<double>(numberOfItems);
      parameterEstimates(3) = u / static_cast<double>(numberOfItems);

      if (parameterEstimates(2) < 0.0 || parameterEstimates(3) > 1.0) {
        return false;
      } else {
        return true;
      }
    }

    /*
      Find value of lower limit of 4-parameter beta compound binomial model
      that minimizes squared difference between predicted and observed kurtosis, 
      such that observed and predicted mean, variance and skewness are equal.
      Used when method of moments fails to produce acceptable values
      of all four parameters.

      When the solution that minimizes the squared difference in kurtosis
      is not the solution with the lower limit = 0 nor the solution with
      the upper limit = 1 it may be missed since an initial solution in this
      case is found by a rough grid search.

      Input
        nitems    = number of items 
        trueScoreMoments   = true score moments (mean, s.d., skewness, kurtosis)
        noncentralTrueScoreMoments = non-central true score moments

      Output
        para = method of moments estimates of parameters

      Returns 1 if solution found, otherwise returns 0.
        Accuracy of solution is +-KACC.

      Function calls other than C or NR utilities: 
        Kurtfuncd()
        KurtfuncUpper()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool calculateBetaParameterLS(const size_t& numberOfItems,
                                  const Eigen::VectorXd& trueScoreMoments,
                                  const Eigen::VectorXd& noncentralTrueScoreMoments,
                                  Eigen::VectorXd& parameterEstimates) {
      bool result = false;

      // double kurt,     /* squared difference of fitted and observed kurtosis */
      //     ll, tll;     /* holds lower limits */
      // double tpara[4]; /* temporary space for beta parameters */
      // double delta, rkurt, lkurt;
      // bool existr, existl;

      /* Initialize squared kurtosis difference to be infinity.
        If no valid solution is found kurt will never be less than
        this. HUGE_VAL is a macro that should be declared in <math.h>
        for any compiler following the ANSI C standard. It represents
        the largest representable floating point number or the floating
        point representation of infinity. */

      double kurt = std::numeric_limits<double>::infinity();
      double lkurt;

      /* find solution for lower limit of 0 */
      bool existl = kurtosisFunctionD(0.0,
                                      numberOfItems,
                                      parameterEstimates,
                                      trueScoreMoments,
                                      noncentralTrueScoreMoments,
                                      lkurt);
      if (existl) {
        kurt = lkurt;
      }

      /* find solution for upper limit of 1 */
      Eigen::VectorXd tempParameterEstimates;
      double rkurt;
      bool existr = kurtosisFunctionUpper(1.0,
                                          numberOfItems,
                                          tempParameterEstimates,
                                          trueScoreMoments,
                                          noncentralTrueScoreMoments,
                                          rkurt);

      if (existr && rkurt < kurt) {
        kurt = rkurt;

        parameterEstimates = tempParameterEstimates;
      }

      /* attempt to find a solution with a lower squared difference
		  in kurtosis than two solutions found above */

      double delta = .05;
      double ll = delta;
      double tll;
      existl = false; /* existl will flag whether a better solution is found */

      while (ll < .55) {
        if (kurtosisFunctionD(ll,
                              numberOfItems,
                              tempParameterEstimates,
                              trueScoreMoments,
                              noncentralTrueScoreMoments,
                              lkurt)) {
          if (lkurt < kurt) { /* a better solution is found if kurtl < kurt */
            kurt = lkurt;
            parameterEstimates = tempParameterEstimates;
            tll = ll;
            existl = true;
          }
        }
        ll += delta;
      }

      if (kurt == std::numeric_limits<double>::infinity()) {
        /* no valid solution has been found */

        result = false;
      } else if (!existl) {
        /* if no better solution than lower limit = 0 or
        upper limit = 1 has been found then return */

        return true;
      } else {
        /* loop to find solution to accuracy of KACC */

        delta = .01;
        ll = tll;
        while (delta >= this->KACC) {
          /* evaluate function at points to left and right 
          of current solution "ll" */

          existr = kurtosisFunctionD(ll + delta,
                                     numberOfItems,
                                     parameterEstimates,
                                     trueScoreMoments,
                                     noncentralTrueScoreMoments,
                                     rkurt);

          existl = kurtosisFunctionD(ll - delta,
                                     numberOfItems,
                                     parameterEstimates,
                                     trueScoreMoments,
                                     noncentralTrueScoreMoments,
                                     lkurt);

          if (existr && (rkurt < kurt)) { /* step to the right */
            ll += delta;
            while (((ll + delta) < 1.0) &&
                   kurtosisFunctionD(ll + delta,
                                     numberOfItems,
                                     parameterEstimates,
                                     trueScoreMoments,
                                     noncentralTrueScoreMoments,
                                     rkurt) &&
                   (rkurt < kurt)) {
              ll += delta;
              kurt = rkurt;
            }
          } else if (existl && (lkurt < kurt)) { /* step to the left */
            ll -= delta;
            while (((ll - delta) > 0.0) &&
                   kurtosisFunctionD(ll - delta,
                                     numberOfItems,
                                     parameterEstimates,
                                     trueScoreMoments,
                                     noncentralTrueScoreMoments,
                                     lkurt) &&
                   (lkurt < kurt)) {
              ll -= delta;
              kurt = lkurt;
            }
          }

          delta /= 10.0;
        }

        /* calculate final values to return */
        result = kurtosisFunctionD(ll,
                                   numberOfItems,
                                   parameterEstimates,
                                   trueScoreMoments,
                                   noncentralTrueScoreMoments,
                                   kurt);
      }

      return result;
    }

    /*
      Find upper limit of 4-parameter beta-binomial model
      that for input value of lower limit produces moments
      that match up to third moment.  

      Input
        para = parameters of beta binomial dist. (para(2) is input)
        tmoment = non-central true score moments

      Ooutput
        para = upper limit of beta binomial dist. (para(3) is output)  

      Returns 1 if successful; otherwise returns 0.

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool findUpper(Eigen::VectorXd& parameterEstimates, const Eigen::VectorXd& trueScoreMoments) {
      bool result = false;

      double lowerLimit = parameterEstimates(2);

      double m1 = trueScoreMoments(0); // first central moment;
      double m2 = trueScoreMoments(1); // second central moment;
      double m3 = trueScoreMoments(2); // third central moment;

      /* calculate upper limit */

      /* numerator */
      double fr1 = m1 * m1 * (m2 * lowerLimit - 2.0 * m3);
      fr1 += m2 * m2 * (m1 - 2.0 * lowerLimit);
      fr1 += m3 * (m2 + m1 * lowerLimit);
      double numerator = fr1;

      /* denominator */
      fr1 = m1 * m1 * (2.0 * m1 * lowerLimit - m2);
      fr1 += m2 * (2.0 * m2 - 3.0 * m1 * lowerLimit);
      fr1 += m3 * (lowerLimit - m1);

      parameterEstimates(3) = numerator / fr1;

      if (parameterEstimates(3) <= parameterEstimates(2) || parameterEstimates(3) > 1.0) {
        result = false;
      } else {
        result = true;
      }

      return result;
    }

    /*
      Find lower limit of 4-parameter beta-binomial model
      that for input value of upper limit produces moments
      that match up to third moment. 

      Input
        para = parameters of beta binomial dist (para(3) is input)
        tmoment = non-central true score moments

      Output
        para = lower limit  (para(2) is output)

      Returns 1 if successful; otherwise returns 0.

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool findLower(Eigen::VectorXd& parameterEstimates, const Eigen::VectorXd& trueScoreMoments) {
      bool result = false;

      double upperLimit = parameterEstimates(3);

      double m1 = trueScoreMoments(0); // first central moment;
      double m2 = trueScoreMoments(1); // second central moment;
      double m3 = trueScoreMoments(2); // third central moment;

      /* calculate lower limit */

      /* numerator */
      double fr1 = m1 * m1 * (m2 * upperLimit - 2.0 * m3);
      fr1 += m2 * m2 * (m1 - 2.0 * upperLimit);
      fr1 += m3 * (m2 + m1 * upperLimit);
      double numerator = fr1;

      /* denominator */
      fr1 = m1 * m1 * (2.0 * m1 * upperLimit - m2);
      fr1 += m2 * (2.0 * m2 - 3.0 * m1 * upperLimit);
      fr1 += m3 * (upperLimit - m1);

      parameterEstimates(2) = numerator / fr1;

      if (parameterEstimates(3) <= parameterEstimates(2) || parameterEstimates(2) < 0.0) {
        result = false;
      } else {
        result = true;
      }

      return result;
    }

    /*
      Estimate parameters alpha and beta of negative 
      hypergeometric distribution, given mean, sd, min, and max.

      Input
        nitems = number of items 
        moment = true score moments (only first and second used)
        para = last two elements contain minimum and maximum

      Output
        para = first two elements contain alpha and beta

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool estimateNegativeHyperGeometricParameters(const size_t& numberOfItems,
                                                  const Eigen::VectorXd& trueScoreMoments,
                                                  Eigen::VectorXd& parameterEstimates) {
      bool result = false;

      /* difference between minimum and maximum */
      double maximumMinusMinimum = static_cast<double>(numberOfItems) * (parameterEstimates(3) - parameterEstimates(2));

      double standardizedMean = (trueScoreMoments(0) - static_cast<double>(numberOfItems) * parameterEstimates(2)) / maximumMinusMinimum;
      double standardizedVariance = std::pow(trueScoreMoments(1), 2) / std::pow(maximumMinusMinimum, 2);

      parameterEstimates(0) = std::pow(standardizedMean, 2) * (1.0 - standardizedMean);
      parameterEstimates(0) /= standardizedVariance;
      parameterEstimates(0) -= standardizedMean;

      parameterEstimates(1) = standardizedMean * (1.0 - standardizedMean);
      parameterEstimates(1) /= standardizedVariance;
      parameterEstimates(1) -= 1.0;
      parameterEstimates(1) -= parameterEstimates(0);

      if (parameterEstimates(0) > 0.0 && parameterEstimates(1) > 0.0) {
        result = true;
      } else {
        result = false;
      }

      return result;
    }

    /*
      Returns value of kurtosis give parameters alpha and beta 
      of 4 parameter beta distribution.

    Input
      para - parameters of beta distribution

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    double calculateKurtosis(const Eigen::VectorXd& parameterEstimates) {
      double kurtosis;

      /* compute components used in calculations */
      double a = parameterEstimates(0);
      double b = parameterEstimates(1);
      double k1 = 3.0 * (a + b + 1.0);
      double k2 = 2.0 * (a + b) * (a + b);
      double k3 = a * b;
      double k4 = a + b - 6.0;
      double k5 = a + b + 2.0;
      double k6 = a + b + 3.0;

      /* compute kurtosis */
      kurtosis = (k1 * (k2 + k3 * k4)) / (k3 * k5 * k6);

      return kurtosis;
    }

    /*
      Using input lower limit find upper limit that matches skewness,
      and using these parameters compute squared difference of
      predicted and observed kurtosis.

      Input
        x = lower limit (0 - 1 scale)
        nitems = number of items 
        tmoment = true score moments (mean, s.d., skewness, kurtosis)
        noncentralTrueScoreMoments = non-central true score moments

      Output
        kurt = squared difference in observed and predicted kurtosis
        para = parameters of four-parameter beta binomial

      Returns 1 if successful; otherwise returns 0.

      Function calls other than C or NR utilities: 
        FindUpper()
        EstNegHypGeo()
        CalcKurt()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool kurtosisFunctionD(const double& lowerLimit,
                           const size_t& numberOfItems,
                           Eigen::VectorXd& parameterEstimates,
                           const Eigen::VectorXd& trueScoreMoments,
                           const Eigen::VectorXd& noncentralTrueScoreMoments,
                           double& kurtosis) {
      bool result = false;

      kurtosis = -1.0;
      parameterEstimates.setZero(4);

      parameterEstimates(2) = lowerLimit;

      if (!findUpper(parameterEstimates, noncentralTrueScoreMoments)) {
        result = false;
      } else if (!estimateNegativeHyperGeometricParameters(numberOfItems,
                                                           trueScoreMoments,
                                                           parameterEstimates)) {
        result = false;
      } else {
        kurtosis = calculateKurtosis(parameterEstimates);
        kurtosis = std::pow(trueScoreMoments(3) - kurtosis, 2);

        result = true;
      }

      return result;
    }

    /*
      Using input upper limit find lower limit that matches skewness,
      and using these parameters compute squared difference of
      predicted and observed kurtosis.

      Input
        x = upper limit (0 - 1 scale)
        nitems = number of items 
        para = beta-binomial parameters (lower limit used for input)
        tmoment = true score moments (mean, s.d., skewness, kurtosis)
        noncentralTrueScoreMoments = non-central true score moments

      Output
        kurt = squared difference in observed and predicted kurtosis

      Returns 1 if successful; otherwise returns 0.

      Function calls other than C or NR utilities: 
        FindLower()
        EstNegHypGeo()
        CalcKurt()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool kurtosisFunctionUpper(const double& upperLimit,
                               const size_t& numberOfItems,
                               Eigen::VectorXd& parameterEstimates,
                               const Eigen::VectorXd& trueScoreMoments,
                               const Eigen::VectorXd& noncentralTrueScoreMoments,
                               double& kurtosis) {
      bool result = false;

      kurtosis = -1.0;
      parameterEstimates(3) = upperLimit;

      if (!findLower(parameterEstimates,
                     noncentralTrueScoreMoments)) {
        result = false;
      } else if (!estimateNegativeHyperGeometricParameters(numberOfItems,
                                                           trueScoreMoments,
                                                           parameterEstimates)) {
        result = false;
      } else {
        kurtosis = calculateKurtosis(parameterEstimates);
        kurtosis = std::pow(trueScoreMoments(3) - kurtosis, 2);

        result = true;
      }

      return result;
    }

    /*
      Given raw score mean, s.d., skewness and kurtosis,
      compute true score mean, s.d., skewness and kurtosis

      Input
        n = number of items
        k = Lord's k
        rawScoreMoments = raw score mean, s.d., skewness and kurtosis

      Output
        trueScoreMoments = True score mean, s.d., skewness and kurtosis. 
                  If true score variance is negative,
                  then s.d., skew and kurt are set to zero.
        noncentralTrueScoreMoments = Non-central true score moments

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void betaMoments(const size_t& numberOfItems,
                     const double& lordK,
                     const Eigen::VectorXd& rawScoreMoments,
                     Eigen::VectorXd& trueScoreMoments,
                     Eigen::VectorXd& noncentralTrueScoreMoments) {
      Eigen::VectorXd centralMoments(4);    /* central moments */
      Eigen::VectorXd noncentralMoments(4); /* non-central moments */
      Eigen::VectorXd factorialMoments(4);  /* factoral moments */

      trueScoreMoments.setZero(4);
      noncentralTrueScoreMoments.setZero(4);

      double prod2, prod3, n2, dn, dnm2, kr; /* hold temporary results */

      /* compute raw score central moments */

      centralMoments(2) = rawScoreMoments(2) * rawScoreMoments(1) * rawScoreMoments(1) * rawScoreMoments(1);
      centralMoments(1) = rawScoreMoments(1) * rawScoreMoments(1); /* variance */
      centralMoments(3) = rawScoreMoments(3) * centralMoments(1) * centralMoments(1);

      /* compute raw score non-central moments */

      prod2 = rawScoreMoments(0) * rawScoreMoments(0);
      noncentralMoments(1) = centralMoments(1) + prod2;
      prod3 = prod2 * rawScoreMoments(0);
      noncentralMoments(2) = centralMoments(2) + 3.0 * centralMoments(1) * rawScoreMoments(0) + prod3;
      noncentralMoments(3) = centralMoments(3) + 4.0 * rawScoreMoments(0) * centralMoments(2) +
                             6.0 * prod2 * centralMoments(1) + prod3 * rawScoreMoments(0);

      /* compute raw factoral moments */

      factorialMoments(1) = noncentralMoments(1) - rawScoreMoments(0);
      factorialMoments(2) = noncentralMoments(2) - 3.0 * noncentralMoments(1) + 2.0 * rawScoreMoments(0);
      factorialMoments(3) = noncentralMoments(3) - 6.0 * noncentralMoments(2) + 11.0 * noncentralMoments(1) -
                            6.0 * rawScoreMoments(0);

      /* compute true score non-central moments */

      /* first moment */
      dn = static_cast<double>(numberOfItems);
      n2 = static_cast<double>(numberOfItems) * static_cast<double>(numberOfItems - 1);

      /* second moment */
      dnm2 = 1.0;
      kr = lordK * 2.0;
      noncentralMoments(0) = rawScoreMoments(0) / dn;
      noncentralTrueScoreMoments(0) = noncentralMoments(0);
      noncentralMoments(1) = (factorialMoments(1) / dnm2 + kr * noncentralMoments(0)) / (n2 + kr);
      noncentralTrueScoreMoments(1) = noncentralMoments(1);

      /* third moment */
      dnm2 = dn - 2.0;
      kr = lordK * 6.0;
      noncentralMoments(2) = (factorialMoments(2) / dnm2 + kr * noncentralMoments(1)) / (n2 + kr);
      noncentralTrueScoreMoments(2) = noncentralMoments(2);
      /* fourth moment */
      dnm2 *= dn - 3.0;
      kr = lordK * 12.0;
      noncentralMoments(3) = (factorialMoments(3) / dnm2 + kr * noncentralMoments(2)) / (n2 + kr);
      noncentralTrueScoreMoments(3) = noncentralMoments(3);

      /* put true score moments on observed score scale */

      dn *= dn;
      noncentralMoments(1) *= dn;
      dn *= static_cast<double>(numberOfItems);
      noncentralMoments(2) *= dn;
      dn *= static_cast<double>(numberOfItems);
      noncentralMoments(3) *= dn;

      /* compute true score central moments */

      prod2 = rawScoreMoments(0) * rawScoreMoments(0);
      centralMoments(1) = noncentralMoments(1) - prod2;
      prod3 = prod2 * rawScoreMoments(0);
      centralMoments(2) = noncentralMoments(2) - 3.0 * noncentralMoments(1) * rawScoreMoments(0) + 2.0 * prod3;
      centralMoments(3) = noncentralMoments(3) - 4.0 * rawScoreMoments(0) * noncentralMoments(2) +
                          6.0 * prod2 * noncentralMoments(1) - 3.0 * prod3 * rawScoreMoments(0);

      /* compute true score mean, s.d., skewness and kurtosis */

      trueScoreMoments(0) = rawScoreMoments(0);

      /* The central true score moment may be negative, if so assign
	    the true score s.d., skewness and kurtosis to zero. */
      if (centralMoments(1) > 0.0) {
        trueScoreMoments(1) = std::sqrt(centralMoments(1));
        trueScoreMoments(2) = centralMoments(2) / (trueScoreMoments(1) * trueScoreMoments(1) * trueScoreMoments(1));
        trueScoreMoments(3) = centralMoments(3) / (centralMoments(1) * centralMoments(1));
      } else {
        trueScoreMoments(1) = 0.0;
        trueScoreMoments(2) = 0.0;
        trueScoreMoments(3) = 0.0;
      }
    }

    /*
      Calculate observed density given parameters of 4-parameter
      beta-binomial distribution. 

      Input
        n = number of score points
        nsamp = sample size
        beta = four parameters of beta distribution 
              (alpha, beta, lower limit, upper limit)

      Output
        scounts = counts of observed score distribution

      Returns 0 if no error, returns -1 if
        memory could not be allocated, -2 if numerical error 
        (overflow or underflow).

      Function calls other than C or NR utilities: 
        CalcM15()
        CalcM24()
        CalcM3()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool observedDensity(const size_t& numberOfScores,
                         const size_t& numberOfExaminees,
                         const Eigen::VectorXd& betaParameterEstimates,
                         Eigen::VectorXd& observedScoreDensity) {
      /* difference in high and low cutoffs for beta */
      double differenceInBetaCutoffs = betaParameterEstimates(3) - betaParameterEstimates(2);
      double densityAdjustment = std::pow(differenceInBetaCutoffs, static_cast<double>(numberOfScores - 1));
      if (densityAdjustment <= 0.0) {
        return false; // underflow
      }

      Eigen::MatrixXd m1;
      Eigen::MatrixXd m2;
      Eigen::MatrixXd m3;
      Eigen::MatrixXd m4;
      Eigen::MatrixXd m5;
      Eigen::VectorXd prod1(numberOfScores);
      Eigen::VectorXd prod2(numberOfScores);

      /* calculate m1 */
      calculateM1(numberOfScores,
                   betaParameterEstimates(0),
                   m1);

      std::cout << "M1: " << EquatingRecipes::Utilities::matrixXdToString(m1) << "\n";

      /* calculate m5 */
      calculateM5(numberOfScores,
                   betaParameterEstimates(1),
                   m5);

      std::cout << "M5: " << EquatingRecipes::Utilities::matrixXdToString(m5) << "\n";

      /* calculate m2 */
      if (!calculateM24(numberOfScores,
                        betaParameterEstimates(2) / differenceInBetaCutoffs,
                        m2)) {
        return false;
      }

      std::cout << "M2: " << EquatingRecipes::Utilities::matrixXdToString(m2) << "\n";


      /* calculate m4 */
      if (!calculateM24(numberOfScores,
                        (1.0 - betaParameterEstimates(3)) / differenceInBetaCutoffs,
                        m4)) {
        return false;
      }

      std::cout << "M4: " << EquatingRecipes::Utilities::matrixXdToString(m4) << "\n";


      /* calculate M3 */
      calculateM3(numberOfScores,
                  betaParameterEstimates,
                  m3);

      std::cout << "M3: " << EquatingRecipes::Utilities::matrixXdToString(m3) << "\n";
      

      /* loop to calculate beta-binomial frequencies */
      observedScoreDensity.setZero(numberOfScores);

      // for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
      //   /* calculate first matrix product - M1 * M2 */
      //   for (size_t innerScoreLocation = scoreLocation + 1; innerScoreLocation >= 0; innerScoreLocation--) {
      //     prod1(innerScoreLocation) = m1(numberOfScores - 1 - scoreLocation) *
      //                                 m2(innerScoreLocation);
      //   }

      //   /* calculate second matrix product - result of first product * M3 */
      //   for (size_t innerScoreLocation = 0; innerScoreLocation < numberOfScores - scoreLocation; innerScoreLocation++) {
      //     prod2(innerScoreLocation) = prod1.segment(0, scoreLocation + 1).dot(m3.row(innerScoreLocation).segment(0, scoreLocation + 1));
      //   }

      //   /* calculate third matrix product - result of second product * M4 */
      //   for (size_t innerScoreLocation = numberOfScores - scoreLocation; innerScoreLocation >= 0; innerScoreLocation--) {
      //     prod1(innerScoreLocation) = prod2(innerScoreLocation) * m4(innerScoreLocation);
      //   }

      //   /* calculate fourth matrix product - result of third product * M5 */
      //   observedScoreDensity(scoreLocation) = densityAdjustment * prod1.segment(0, numberOfScores - scoreLocation).dot(m5.segment(scoreLocation, numberOfScores));
      // }

      Eigen::MatrixXd matProd = m1 * m2 * m3 * m4 * m5;
      Eigen::Diagonal<Eigen::MatrixXd> matProdDiag = matProd.diagonal();

      double uMinusLFactorial = 1;

      for (size_t scoreLocation = 1; scoreLocation < numberOfScores; scoreLocation++) {
        uMinusLFactorial *= (betaParameterEstimates(3) - betaParameterEstimates(2));
        matProdDiag(scoreLocation) *= uMinusLFactorial;
      }

      std::cout << "MProduct: " << EquatingRecipes::Utilities::vectorXdToString(matProdDiag, false) << "\n";

      return true;
    }

    /*
      Adjusts observed density calculated by "ObsDensity" for
      compound binomial error using Lord's k.  Uses equations (53)-(55)
      in Lord (1965,Psychometrika,239-270) page 268.

      Input
        n = number of items
        k = Lord's k
        scounts = fitted counts of observed score distribution
                  computed by ObsDensity()

      Output
        scounts = adjusted fitted counts of observed score distribution

      Function calls other than C or NR utilities: None
                                      
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    Eigen::VectorXd observedDensityK(const size_t& numberOfItems,
                                     const double& lordK,
                                     const Eigen::VectorXd& fittedObservedScoreDist) {
      double d2, pxm1, px, pxp1, dsp1, kn2, dn;

      Eigen::VectorXd adjustedFittedObservedScoreDist = fittedObservedScoreDist;

      dn = static_cast<double>(numberOfItems);
      kn2 = lordK / (dn * static_cast<double>(numberOfItems - 1));

      /* raw score of 0 */
      pxp1 = (dn - 1.0) * adjustedFittedObservedScoreDist(1);
      adjustedFittedObservedScoreDist(0) -= kn2 * pxp1;

      px = pxp1;
      pxm1 = 0.0;

      /* raw scores 1 through n-1 */
      for (size_t scoreLocation = 1; scoreLocation < numberOfItems; scoreLocation++) {
        dsp1 = static_cast<double>(scoreLocation + 1);

        pxp1 = (dsp1) * (dn - dsp1) * adjustedFittedObservedScoreDist(scoreLocation + 1);
        d2 = pxm1 - 2.0 * px + pxp1;

        adjustedFittedObservedScoreDist(scoreLocation) -= kn2 * d2;
        pxm1 = px;
        px = pxp1;
      }

      /* raw score n */
      adjustedFittedObservedScoreDist(numberOfItems) -= kn2 * pxm1;

      return adjustedFittedObservedScoreDist;
    }

    /*
      Calculate matrix M3

      Input
        n = number of score points
        beta = parameters of 4-parameter beta true score distribution

      Output
        m3 = matrix used for calculation of raw score distribution

      Function calls other than C or NR utilities: None
                                      
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void calculateM3(const size_t& numberOfScores, const Eigen::VectorXd& betaParameters, Eigen::MatrixXd& m3) {
      // register double *pc;
      // register double fr;	   /* floating point register for temp storage */
      // register int i;
      // double apb;	       /* sum of shape parameters of beta distribution */
      // int j;
      // int nitems;

      double apb = betaParameters(0) + betaParameters(1);
      size_t numberOfItems = numberOfScores - 1;

      m3.setZero(numberOfScores, numberOfScores);

      /* compute first elements of all m3[j] */
      m3(0, numberOfScores - 1) = 1.0;

      for (size_t scoreLocation = 1; scoreLocation < numberOfScores; scoreLocation++) {
        double fr = static_cast<double>(scoreLocation);

        m3(0, numberOfScores - scoreLocation - 1) = (fr / (apb + fr - 1.0)) * m3(0, numberOfScores - scoreLocation);        
      }

      /* compute other elements */
      for (size_t columnIndex = 1; columnIndex < numberOfScores; columnIndex++) {        
        for (size_t rowIndex = 1; rowIndex <= columnIndex; rowIndex++) {        
          double fr = apb - static_cast<double>(rowIndex) + static_cast<double>(columnIndex);
          fr /= static_cast<double>(numberOfItems - rowIndex + 1);
          fr *= m3(rowIndex - 1, numberOfScores - columnIndex - 1);
          
          m3(rowIndex, numberOfScores - columnIndex - 1) = fr;
        }
      }
    }

    /*
      Calculates (diagonal) matrices M2 and M4.

      Input
        n = number of score points
        para = parameter from beta distribution to use in calculation

      Output
        m = matrix to calculate

      Return 0 if no error, -2 if underflow

      Function calls other than C or NR utilities: None
                                      
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    bool calculateM24(const size_t& numberOfScores, const double& betaParameter, Eigen::MatrixXd& m) {
      size_t numberOfItems = numberOfScores - 1;

      m.setZero(numberOfScores, numberOfScores);

      m(0, 0) = 1.0;

      for (size_t scoreLocation = 1; scoreLocation < numberOfScores; scoreLocation++) {
        double fr = (static_cast<double>(numberOfItems - (scoreLocation - 1))) / (static_cast<double>(scoreLocation));
        fr *= betaParameter;
        fr *= m(scoreLocation - 1, scoreLocation - 1);

        if (fr == HUGE_VAL) {
          return false;
        }

        m(scoreLocation, scoreLocation) = fr;
      }

      return true;
    }

    /*
      Calculate matrices M1 and M5

      Input
        n = number of score points
        para = parameter from beta distribution to use in calculation

      Output
        m = matrix to calculate

      Function calls other than C or NR utilities: None
                                      
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void calculateM15(const size_t& numberOfScores, const double& betaParameter, Eigen::MatrixXd& m) {
      m.setZero(numberOfScores, numberOfScores);

      m(numberOfScores - 1) = 1.0;

      for (size_t scoreLocation = 1; scoreLocation < numberOfScores; scoreLocation++) {
        double di = static_cast<double>(scoreLocation);

        m(numberOfScores - 1 - scoreLocation) = ((betaParameter + (di - 1.0)) / di) * m(numberOfScores - scoreLocation);
      }
    }

    void calculateM1(const size_t& numberOfScores, const double& betaParameter, Eigen::MatrixXd& m) {
      m.setZero(numberOfScores, numberOfScores);

      for (size_t rowIndex = 0; rowIndex < numberOfScores; rowIndex++) {        
        for (size_t columnIndex = 0; columnIndex <= rowIndex; columnIndex++) {        
          if (rowIndex == columnIndex) {
            m(rowIndex, columnIndex) = 1.0;
          } else {
            m(rowIndex, columnIndex) = m(rowIndex - 1, columnIndex) * 
              (betaParameter + static_cast<double>(rowIndex - 1)) / static_cast<double>(rowIndex - columnIndex);
          }
        }
      }
    }

    void calculateM5(const size_t& numberOfScores, const double& betaParameter, Eigen::MatrixXd& m) {
      Eigen::MatrixXd tempMat;
      calculateM1(numberOfScores, betaParameter, tempMat);

      m.setZero(numberOfScores, numberOfScores);

      for (size_t rowIndex = 0; rowIndex < numberOfScores; rowIndex++) {        
        for (size_t columnIndex = 0; columnIndex <= rowIndex; columnIndex++) {        
          m(numberOfScores - 1 - rowIndex, columnIndex) = tempMat(rowIndex, columnIndex);
        }
      }
    }

    double likelihoodRatioChiSquare(const size_t& n, const Eigen::VectorXd& rawc, const Eigen::VectorXd& fitc, size_t& ncat) {
      double result = 0;
      return result;
    }

    double chiSquarePValue(const size_t& n, const double& minexpcount, const Eigen::VectorXd& rawc, const Eigen::VectorXd& fitc, size_t& ncat) {
      double result = 0;
      return result;
    }

    double calculateLordK(const double& kr20, const size_t numberOfItems, const Eigen::VectorXd& rawScoreMoments) {
      double result = 0;
      return result;
    }
  };
} // namespace EquatingRecipes

#endif