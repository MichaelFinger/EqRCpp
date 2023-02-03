/* 
  From Source: ERutilities.h 
  Original Struct: ULL_SMOOTH
  Description: univariate log-linear smoothing
*/

#ifndef STRUCTURES_UNIVARIATE_LOG_LINEAR_SMOOTHING_HPP
#define STRUCTURES_UNIVARIATE_LOG_LINEAR_SMOOTHING_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct UnivariateLogLinearSmoothing {
      int numberOfPersons;   /* number of persons */
      int numberOfScoreCategories;            /* number of score categories */
      double mininumRawScore;        /* minimum raw score */
      double rawScoreIncrement;        /* increment in raw scores */
      int degreesOfSmoothing;             /* number of degrees of smoothing */
      Eigen::MatrixXd rawScoreDesignMatrix;    /* design matrix for raw scores */
      Eigen::MatrixXd solutionDesignMatrix;        /* design matrix used for solution */
      Eigen::VectorXd actualFrequencies;       /* actual frqeuncies */
      Eigen::VectorXd fittedFrequencies;       /* fitted frequencies */
      Eigen::VectorXd betaCoefficients;      /* coefficients */
      double cllNormalizingConstant;         /* normalizing constant used in CLL */
      Eigen::VectorXd bActualMoments;     /* actual moments for B */
      Eigen::VectorXd bFittedMoments;     /* fitted moments for B */
      Eigen::VectorXd actualCentralMoments; /* central moments for actual and */
      Eigen::VectorXd fittedCentralMoments; /* central moments for fitted and */
      int numberOfIterations;           /* number of iterations to convergence */
      int maximumNumberOfIterators;       /* maximum number of iterations allowed */
      bool useRelativeCriterionComparison;         /* comparison for criterion (0-->absolute; 1-->relative) */
      bool useBRawAndCentralMoments;         /* 0--> use B mts for crit; 1--> use B_raw and cent mts */
      bool useScalingForBDesignMatrix;         /* 0 --> no scaling for B design matrix; 1--> scaling  */
      double convergenceCriterion;       /* convergence criterion */
      double likelihoodRatioChiSq;    /* likelihood-ratio chi-square */
      int numberOfZeros;         /* # 0's (see iteration()) */
      Eigen::VectorXd fittedRawScoreDist;   /* pointer to fitted raw score dist (proportions) */
      Eigen::VectorXd fittedRawScoreCumulativeRelativeDist;      /* pointer to cum rel freq dist for fitted dist */
      Eigen::VectorXd fittedRawScorePercentileRankDist;       /*pointer to percentile rank dist for fitted dist */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif