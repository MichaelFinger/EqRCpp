/* 
  From Source: ERutilities.h 
  Original Struct: BLL_SMOOTH
  Description: Bivariate log-linear smoothing, where bivariate distribution that is
    smoothed is (rows = non-common items) x (cols = common items);
	  i.e. u by v.  nct[] is the row-major version of this matrix for 
	  actual frequencies; mct[] is the corresponding matrix for fitted
	  frequencies.  bfd[][]is for the x by v matrix of 
	  fitted frequencies; subsequent elements are marginals for this 
	  matrix.  For an external anchor, u by v is identical to x by v.

	  Note that notation is in terms of x/u/v (implicitly for pop 1); 
	  same logic applies to y/u/v (implicitly for pop 2). 
  
	  For a single-group design, the "anchor" is effectively external 
	  in the sense that x=u, v plays the role of y, and Wrapper_SL 
	  will put x (rows) on scale of y (cols).   
*/

#ifndef STRUCTURES_BIVARIATE_LOG_LINEAR_SMOOTHING_HPP
#define STRUCTURES_BIVARIATE_LOG_LINEAR_SMOOTHING_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct BivariateLogLinearSmoothing {
      bool useInternalAnchor;                                             /* = 0 (external); = 1 (internal) */
      int numberOfPersons;                                                /* number of persons */
      int nonCommonItemsNumberOfScoreCategories;                          /* number of score categories for u (non-common items) */
      double nonCommonItemsMinimumRawScore;                               /* minimum raw score for u */
      double nonCommonItemsRawScoreIncrement;                             /* increment in raw scores for u */
      int commonItemsNumberOfScoreCategories;                             /* number of score categories for v (common items) */
      double commonItemsMinimumRawScore;                                  /* minimum raw score for v */
      double commonItemsRawScoreIncrement;                                /* increment in raw scores for v */
      int totalNumberOfScoreCategories;                                   /* total number of score categories */
      int nonCommonItemsDegreesOfSmoothing;                               /* number of degrees of smoothing for u */
      int commonItemsDegreesOfSmoothing;                                  /* number of degrees of smoothing for v */
      int numberOfCrossProductMomemts;                                    /* number of (u,v) cross product moments for smoothing */
      Eigen::MatrixXi crossProductMoments;                                /* cross-product moments for smoothing */
      int numberOfColumnsInDesignMatrix;                                  /* number of columns in design matrix */
      Eigen::MatrixXd rawScoreDesignMatrix;                               /* design matrix for raw scores */
      Eigen::MatrixXd solutionDesignMatrix;                               /* design matrix used for solution */
      Eigen::VectorXd actualFrequencies;                                  /* actual frqeuncies */
      Eigen::VectorXd fittedFrequencies;                                  /* fitted frequencies */
      Eigen::VectorXd betaCoefficients;                                   /* coefficients */
      double cllNormalizingConstant;                                      /* normalizing constant used in CLL */
      Eigen::VectorXd bActualMoments;                                     /* actual moments for B */
      Eigen::VectorXd bFittedMoments;                                     /* fitted moments for B */
      Eigen::VectorXd actualCentralMoments;                               /* central moments for actual using B_raw */
      Eigen::VectorXd fittedCentralMoments;                               /* central moments for fitted using B_raw */
      int numberOfIterations;                                             /* number of iterations to convergence */
      int maximumNumberOfIterations;                                      /* maximum number of iterations allowed */
      bool useRelativeCriterionComparison;                                /* comparison for criterion (0-->absolute; 1-->relative) */
      bool useRawAndCentralMoments;                                       /* 0--> use B mts for crit; 1--> use B_raw and cent mts */
      bool useScalingForBDeisngMatrix;                                    /* 0 --> no scaling for B design matrix; 1--> scaling  */
      double convergenceCriterion;                                        /* convergence criterion */
      double likelihoodRatioChiSq;                                        /* likelihood-ratio chi-square */
      int numberOfZeros;                                                  /* # 0's (see iteration()) */
      int nsx;                                                            /* = nsu + nsv -1 if internal anchor; = nsu if ext anchor  */
      double xMinimumRawScore;                                            /* minimum raw score for x */
      double xRawScoreIncrement;                                          /* increment in raw scores for x */
      Eigen::MatrixXd fittedBivariateFreqDist;                            /* fitted biv freq dist for x by v (see note below) */
      Eigen::VectorXd xFittedFrequences;                                  /* pointer to fitted frequencies for x */
      Eigen::VectorXd xFittedRawScoreDensity;                             /* pointer to fitted raw score props for x */
      Eigen::VectorXd xFittedRawScoreCumulativeRelativeFreqDist;          /* pointer to cum rel freq dist for fit. dist for x */
      Eigen::VectorXd xFittedRawScorePercentileRankDist;                  /* pointer to PR dist for fitted dist for x */
      Eigen::VectorXd commonItemFittedFrequences;                         /* pointer to fitted frequencies for v */
      Eigen::VectorXd commonItemFittedRawScoreDensity;                    /* pointer to fitted raw score proportions for v */
      Eigen::VectorXd commonItemFittedRawScoreCumulativeRelativeFreqDist; /* pointer to cum rel freq dist for fitted dist for v */
      Eigen::VectorXd commonItemFittedRawScorePercentileRankDist;         /* pointer to percentile rank dist for fitted dist for v */
      Eigen::MatrixXd xVFittedBivariateRelativeFreqDist;                  /* fitted biv rel freq dist for x by v
                                                                             Note: if external anchor, bfd[][] has dimensions nsu x nsv; 
                                                                               if internal anchor, bfd[][] has dimensions nsx x nsv, where
                                                                               nsx = nsu + nsv - 1 (bfd[][] includes structural zeros);
                                                                               density[], crfd[], and prd[] are associated with rows of 
                                                                               bfd[][]; density_v[], crfd_v[], and prd_v[] are associated 
                                                                               with columns of bfd[][]; */
      Eigen::VectorXd cumulativeRelativeFreqDistVector;                   /* cum rel fd as a row-major vector from bfd[][]; used for bootstrap */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif