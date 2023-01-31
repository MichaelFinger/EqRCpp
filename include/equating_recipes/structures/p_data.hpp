/* 
  From Source: ERutilities.h
  Original Method: PDATA 
  Description: structure that contains all input for a particular 
    design/method/smoothing.  Used extensively in the argument lists for
    for Wrapper and Print functions 
*/

#ifndef STRUCTURES_PDATA_HPP
#define STRUCTURES_PDATA_HPP

#include <vector>
#include <Eigen/Core>

#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/cubic_spline_postsmoothing.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/irt_input.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct PData {
      UnivariateStatistics summaryRawDataX;     /* structure for summary raw data for x and v */
      UnivariateStatistics summaryRawDataY;     /* structure for summary raw data for y and v */
      BivariateStatistics summaryRawDataXV;     /* structure for summary raw data for x and v */
      BivariateStatistics summaryRawDataYV;     /* structure for summary raw data for y and v */
      BivariateStatistics summaryRawDataXY;     /* structure for summary raw data for x and y */
      Design design;
      Method method;
      Smoothing smoothing;

      double weightSyntheticPopulation1;                               /* weight for synthetic pop 1 */
      bool anchorIsInternal;                          /* = false (external); = true (internal) */
      double rreliabilityCommonItemsPopulation1;                    /* reliability of common items in pop 1 */
      double rreliabilityCommonItemsPopulation2;                    /* reliability of common items in pop 2 */
      std::vector<Method> methods;
      double mininumScoreX;                                         /* min score for x */
      double maximumScoreX;                                         /* max score for x */
      double scoreIncrementX;                 /* increment between adjacent scores for x */
      Eigen::VectorXd scoreFrequenciesX;                /* fd for new form x */
      RawToScaledScoreTable rawToScaledScoreTable /* conversion table for Y */
      size_t numberOfBootstrapReplications;                      /* number of replications for bootstrap */
      size_t bootstrapReplicationNumber = 0;     /* rep number for bootstrap; set to 0 for actual equating */
      size_t numberOfDecimalPlacesToRound = 0;
      BetaBinomalSmoothing betaBinomalSmoothingX;  /* structure for beta binomial smoothing for x */
      BetaBinomalSmoothing betaBinomalSmoothingY;  /* structure for beta binomial smoothing for y */ 
      UnivariateLogLinearSmoothing univariateLogLinearSmoothingX; /* structure for univ log-lin smoothing for x */ 
      UnivariateLogLinearSmoothing univariateLogLinearSmoothingY; /* structure for univ log-lin smoothing for y */
      BivariateLogLinearSmoothing bivariateLogLinearSmoothingXV; /* struc for biv log-lin smoothing for x & v */ 
      BivariateLogLinearSmoothing bivariateLogLinearSmoothingYV; /* struc for biv log-lin smoothing for y & v */
      BivariateLogLinearSmoothing bivariateLogLinearSmoothingXY; /* struc for biv log-lin smoothing for x & y */
      CubicSplinePostsmoothing cubicSplinePostsmoothing;      /* structure for cubic-spline postsmoothing */
      struct IRT_INPUT *IRT_Input;                /* structure for IRT input */
    }
  }
}

#endif