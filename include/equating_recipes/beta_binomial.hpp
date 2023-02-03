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
      short result = 0;
      return result;
    }

    short twoParameterBetaSmooting(const Eigen::VectorXd& rawScoreFreqDist,
                                   EquatingRecipes::Structures::BetaBinomialSmoothing& betaFitResults) {
      short result = 0;
      return result;
    }

    // void Print_BB(FILE* fp, char tt[], struct USTATS* x, struct BB_SMOOTH* s);
    // void Print_RB(FILE* fp, char tt[], struct PDATA* inall, struct ERAW_RESULTS* r);

    short calculateBetaParameters(const size_t& numberOfItems, const Eigen::VectorXd& moments, const Eigen::VectorXd& nctmoment, Eigen::VectorXd& parameterEstimates) {
      short result = 0;
      return result;
    }

    short estimateNegativeHyperGeometricParameters(const size_t& numberOfItems, const Eigen::VectorXd& moments, Eigen::VectorXd& parameterEstimates) {
      short result = 0;
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

    void betaMoments(const size_t& n, const double& k, const Eigen::VectorXd& rmoment, const Eigen::VectorXd& tmoment, Eigen::VectorXd& nctmoment) {
    }

    short observedDensity(const size_t& n, const size_t& nsamp, const Eigen::VectorXd& beta, const Eigen::VectorXd& scounts) {
      short result = 0;
      return result;
    }

    void observedDensityK(const size_t& n, const double& k, const Eigen::VectorXd& scounts) {}

    void calculateM3(const size_t& n, const Eigen::VectorXd& beta, Eigen::MatrixXd& m3) {}

    short calculateM24(const size_t& n, const double& para, Eigen::VectorXd& m) {
      short result = 0;
      return result;
    }

    void calculateM15(const size_t& n, const double& para, Eigen::VectorXd& m) {}

    double likelihoodRatioChiSquare(const size_t& n, const Eigen::VectorXd& rawc, const Eigen::VectorXd& fitc, const Eigen::VectorXi& ncat) {
      double result = 0;
      return result;
    }

    double chiSquarePValue(const size_t& n, const double& minexpcount, const Eigen::VectorXd& rawc, const Eigen::VectorXd& fitc, const Eigen::VectorXi& ncat) {
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