/* 
  From Source: ERutilities.h 
  Original Struct: BB_SMOOTH
  Description: beta-binomial smoothing: 2 or 4 parameter beta with binomial or compound binomial error
*/

#ifndef STRUCTURES_BETA_BINOMIAL_SMOOTHING_HPP
#define STRUCTURES_BETA_BINOMIAL_SMOOTHING_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct BetaBinomialSmoothing {
      int numberOfItems;		                                                           /* number of items on test */
      int numberOfExaminees;                                                                       /* sample size */
      int numberOfParameters;                                                    /* number of parameters (2 or 4) */
      double reliablilty;                                                    /* reliability -- almost always Kr20 */
      double lordK;		                                         /* Lord's k for approximation of compound binomial */
      Eigen::VectorXd betaParameters;	                 /* parameters of true score distribution */
      Eigen::VectorXd rawScoreMoments;	                                                     /* Raw score moments */
      Eigen::VectorXd fittedRawScoreMoments;	                                        /* Fitted raw score moments */
      Eigen::VectorXd trueScoreMoments;	                                                    /* True score moments */
      double likelihoodRatioChiSq;		                             /* likelihood ratio chi-square for fitted dist */
      double pValueChiSq;		                                                /* Pearson chi-square for fitted dist */
      short numberOfMomentsFit;	                                                         /* number of moments fit */
      Eigen::VectorXd fittedRawScoreDist;	                      /* pointer to fitted raw score dist (proportions) */
      Eigen::VectorXd fittedRawScoreCumulativeRelativeFreqDist;   /* pointer to cum rel freq dist for fitted dist */
      Eigen::VectorXd fittedRawScorePercentileRankDist;         /*pointer to percentile rank dist for fitted dist */
    };
  }
}

#endif