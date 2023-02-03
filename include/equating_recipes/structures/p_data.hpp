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
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
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
      EquatingRecipes::Structures::UnivariateStatistics summaryRawDataX;     /* structure for summary raw data for x and v */
      EquatingRecipes::Structures::UnivariateStatistics summaryRawDataY;     /* structure for summary raw data for y and v */
      EquatingRecipes::Structures::BivariateStatistics summaryRawDataXV;     /* structure for summary raw data for x and v */
      EquatingRecipes::Structures::BivariateStatistics summaryRawDataYV;     /* structure for summary raw data for y and v */
      EquatingRecipes::Structures::BivariateStatistics summaryRawDataXY;     /* structure for summary raw data for x and y */
      EquatingRecipes::Structures::Design design;
      EquatingRecipes::Structures::Method method;
      EquatingRecipes::Structures::Smoothing smoothing;

      double minimumRawScoreYct;
      double maximumRawScoreYct;
      double scoreIncrementYct;
      int lowestObservableRoundedScaledScore;                                       /* lowest possible rounded scale score */
      int highestObservableRoundedScaledScore;                                     /* highest possible rounded scale score */

      double weightSyntheticPopulation1;                                                     /* weight for synthetic pop 1 */
      bool anchorIsInternal;                                                      /* = false (external); = true (internal) */
      double rreliabilityCommonItemsPopulation1;                                   /* reliability of common items in pop 1 */
      double rreliabilityCommonItemsPopulation2;                                   /* reliability of common items in pop 2 */
      std::vector<Method> methods;
      double mininumScoreX;                                                                             /* min score for x */
      double maximumScoreX;                                                                             /* max score for x */
      double scoreIncrementX;                                                   /* increment between adjacent scores for x */
      Eigen::VectorXd scoreFrequenciesX;                                                              /* fd for new form x */
      EquatingRecipes::Structures::RawToScaledScoreTable rawToScaledScoreTable;                                               /* conversion table for Y */
      size_t numberOfBootstrapReplications;                                        /* number of replications for bootstrap */
      size_t bootstrapReplicationNumber = 0;                     /* rep number for bootstrap; set to 0 for actual equating */
      size_t roundToNumberOfDecimalPlaces = 0;
      EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomalSmoothingX;                           /* structure for beta binomial smoothing for x */
      EquatingRecipes::Structures::BetaBinomialSmoothing betaBinomalSmoothingY;                           /* structure for beta binomial smoothing for y */ 
      EquatingRecipes::Structures::UnivariateLogLinearSmoothing univariateLogLinearSmoothingX;           /* structure for univ log-lin smoothing for x */ 
      EquatingRecipes::Structures::UnivariateLogLinearSmoothing univariateLogLinearSmoothingY;           /* structure for univ log-lin smoothing for y */
      EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingXV;            /* struc for biv log-lin smoothing for x & v */ 
      EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingYV;            /* struc for biv log-lin smoothing for y & v */
      EquatingRecipes::Structures::BivariateLogLinearSmoothing bivariateLogLinearSmoothingXY;            /* struc for biv log-lin smoothing for x & y */
      EquatingRecipes::Structures::CubicSplinePostsmoothing cubicSplinePostsmoothing;                    /* structure for cubic-spline postsmoothing */
      EquatingRecipes::Structures::IRTInput *IRT_Input;                                          /* structure for IRT input */
    };
  }
}

#endif