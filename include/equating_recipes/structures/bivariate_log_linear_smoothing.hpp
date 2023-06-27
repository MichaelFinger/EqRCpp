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
#include <equating_recipes/structures/criterion_comparison_type.hpp>
#include <equating_recipes/structures/design_matrix_type.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct BivariateLogLinearSmoothing {
      bool isInternalAnchor;                                                          /* = 0 (external); = 1 (internal) */
      size_t numberOfExaminees;                                                       /* number of persons */
      size_t numberOfScoresU;                                                         /* number of score categories for u (non-common items) */
      double minimumRawScoreU;                                                        /* minimum raw score for u */
      double scoreIncrementU;                                                         /* increment in raw scores for u */
      size_t numberOfScoresV;                                                         /* number of score categories for v (common items) */
      double minimumRawScoreV;                                                        /* minimum raw score for v */
      double scoreIncrementV;                                                         /* increment in raw scores for v */
      size_t totalNumberOfScores;                                                     /* total number of score categories */
      size_t numberOfDegreesOfSmoothingU;                                             /* number of degrees of smoothing for u */
      size_t numberOfDegreesOfSmoothingV;                                             /* number of degrees of smoothing for v */
      size_t numberOfCrossProductMoments;                                             /* number of (u,v) cross product moments for smoothing */
      Eigen::MatrixXi crossProductMoments;                                            /* cross-product moments for smoothing, zero-offset matrix designating cross-
                                                                                         product moments.  Example: let cuv = 3,
                                                                                         and the desired cross-product moments be 
                                                                                         (u^1)*(v^1), (u^1)*(v^2), and (u^2)*(v^1). 
                                                                                         Then cpm[0] = {1,1}; cpm[1] = {1,2}; and
                                                                                         cpm[2] = {2,1}. */
      size_t numberOfColumnsInDesignMatrix;                                           /* number of columns in design matrix */
      Eigen::MatrixXd rawScoreDesignMatrix;                                           /* design matrix for raw scores */
      Eigen::MatrixXd solutionDesignMatrix;                                           /* design matrix used for solution */
      Eigen::VectorXd observedFrequencies;                                            /* actual frqeuncies */
      Eigen::VectorXd fittedFrequencies;                                              /* fitted frequencies */
      Eigen::VectorXd betaCoefficients;                                               /* coefficients */
      double cllNormalizingConstant;                                                  /* normalizing constant used in CLL */
      Eigen::VectorXd solutionDesignMatrixObservedMoments;                            /* actual moments for B */
      Eigen::VectorXd solutionDesignMatrixFittedMoments;                              /* fitted moments for B */
      Eigen::VectorXd rawScoreDesignMatrixObservedCentralMoments;                     /* central moments for actual using B_raw */
      Eigen::VectorXd rawScoreDesignMatrixFittedCentralMoments;                       /* central moments for fitted using B_raw */
      size_t numberOfIterations;                                                      /* number of iterations to convergence */
      size_t maximumNumberOfIterations;                                               /* maximum number of iterations allowed */
      EquatingRecipes::Structures::CriterionComparisonType criterionComparisonType;   /* comparison for criterion (0-->absolute; 1-->relative) */
      EquatingRecipes::Structures::DesignMatrixType designMatrixType;                 /* 0--> use B mts for crit; 1--> use B_raw and cent mts */
      bool useScalingForBDeisngMatrix;                                                /* 0 --> no scaling for B design matrix; 1--> scaling  */
      double convergenceCriterion;                                                    /* convergence criterion */
      double likelihoodRatioChiSquare;                                                /* likelihood-ratio chi-square */
      size_t numberOfZeros;                                                           /* # 0's (see iteration()) */
      size_t numberOfScoresX;                                                         /* = nsu + nsv -1 if internal anchor; = nsu if ext anchor  */
      double minimumRawScoreX;                                                        /* minimum raw score for x */
      double scoreIncrementX;                                                         /* increment in raw scores for x */
      Eigen::MatrixXd fittedBivariateFreqDist;                                        /* fitted biv freq dist for x by v (see note below) */
      Eigen::VectorXd fittedFrequencesX;                                              /* pointer to fitted frequencies for x */
      Eigen::VectorXd fittedRawScoreDensityX;                                         /* pointer to fitted raw score props for x */
      Eigen::VectorXd fittedRawScoreCumulativeRelativeFreqDistX;                      /* pointer to cum rel freq dist for fit. dist for x */
      Eigen::VectorXd fittedRawScorePercentileRankDistX;                              /* pointer to PR dist for fitted dist for x */
      Eigen::VectorXd fittedFrequencesV;                                              /* pointer to fitted frequencies for v */
      Eigen::VectorXd fittedRawScoreDensityV;                                         /* pointer to fitted raw score proportions for v */
      Eigen::VectorXd fittedRawScoreCumulativeRelativeFreqDistV;                      /* pointer to cum rel freq dist for fitted dist for v */
      Eigen::VectorXd fittedRawScorePercentileRankDistV;                              /* pointer to percentile rank dist for fitted dist for v */
      Eigen::MatrixXd fittedBivariateRelativeFreqDistXV;                              /* fitted biv rel freq dist for x by v
                                                                                                 Note: if external anchor, bfd[][] has dimensions nsu x nsv; 
                                                                                                   if internal anchor, bfd[][] has dimensions nsx x nsv, where
                                                                                                   nsx = nsu + nsv - 1 (bfd[][] includes structural zeros);
                                                                                                   density[], crfd[], and prd[] are associated with rows of 
                                                                                                   bfd[][]; density_v[], crfd_v[], and prd_v[] are associated 
                                                                                                   with columns of bfd[][]; */
      Eigen::VectorXd cumulativeRelativeFreqDistRowMajorVector;                       /* cum rel fd as a row-major vector from bfd[][]; used for bootstrap */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif