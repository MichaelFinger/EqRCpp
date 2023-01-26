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
      int numberOfScoreCategories;                           /* number of score categories */
      double smoothingParameter;                             /* smoothing or "flexibility" parameter */
      double lowestSmoothedScorePercentileRank;              /* percentile rank for lowest score that is smoothed */
      double highestSmoothedScorePercentileRank;             /* percentile rank for highest score that is smoothed */
      int lowestSmoothedPseudoRawScorePercentileRank;        /* lowest (pseudo raw) score that is smoothed */
      int higestSmoothedPseudoRawScorePercentileRank;        /* highest (pseudo raw) score that is smoothed */
      int numberOfSmoothedPseudoRawScoreCategories;          /* high-low+1 = number of score categories in [low,high] */
      Eigen::VectorXd equipercentileEquivalents;             /* equipercentile equivalents: eeq[ns] */
      Eigen::VectorXd standardErrors;                        /* standard errors: se[ns] */
      Eigen::VectorXd coefficients;                          /* vector containing a, b, c, d coeffs: cmat[4*ns]; Note--dimensioned for maximum */
      Eigen::VectorXd cubicSplineSmoothedEquivalents;        /* cubic-spline smoothed equivalents including interpolated values: eeqs[ns] */
      Eigen::VectorXd inverseCubicSplineSmoothedEquivalents; /* inverse of eeqs[]; computed only for Y to X; number of score categories is ns for X, NOT ns for Y: inv[ns for X] */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif