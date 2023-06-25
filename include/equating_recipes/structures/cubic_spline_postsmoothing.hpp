/* 
  From Source: ERutilities.h 
  Original Struct: CS_SMOOTH
  Description: cubic-spline postsmoothing; need one of these structures in
 	  InputParameters structure for putting x on scale of y, and one for
 	  InputParameters structure for putting y on scale of x.
*/

#ifndef STRUCTURES_CUBIC_SPLINE_POSTSMOOTHING_HPP
#define STRUCTURES_CUBIC_SPLINE_POSTSMOOTHING_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct CubicSplinePostsmoothing {
      size_t numberOfScores;                                 /* number of score categories */
      double smoothingParameter;                             /* smoothing or "flexibility" parameter */
      double percentileRankLowestScore;                      /* percentile rank for lowest score that is smoothed */
      double percentileRankHighestScore;                     /* percentile rank for highest score that is smoothed */
      double lowestSmoothedPseudoRawScorePercentileRank;     /* lowest (pseudo raw) score that is smoothed */
      double higestSmoothedPseudoRawScorePercentileRank;     /* highest (pseudo raw) score that is smoothed */
      size_t boundedNumberOfScores;                          /* high-low+1 = number of score categories in [low,high] */
      Eigen::MatrixXd equipercentileEquivalents;             /* equipercentile equivalents: eeq[ns] */
      Eigen::MatrixXd standardErrors;                        /* standard errors: se[ns] */
      Eigen::MatrixXd coefficients;                          /* vector containing a, b, c, d coeffs: cmat[4*ns]; Note--dimensioned for maximum */
      Eigen::MatrixXd cubicSplineSmoothedEquivalents;        /* cubic-spline smoothed equivalents including interpolated values: eeqs[ns] */
      Eigen::MatrixXd inverseCubicSplineSmoothedEquivalents; /* inverse of eeqs[]; computed only for Y to X; number of score categories is ns for X, NOT ns for Y: inv[ns for X] */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif