/* 
  From Source: ERutilities.h
  Original Method: PDATA 
  Description: structure that contains all input for a particular 
    design/method/smoothing.  Used extensively in the argument lists for
    for Wrapper and Print functions 
*/

#ifndef IMPLEMENTATION_PDATA_HPP
#define IMPLEMENTATION_PDATA_HPP

#include <limits>
#include <string>
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
      EquatingRecipes::Structures::Design design = EquatingRecipes::Structures::Design::NOT_SPECIFIED;
      EquatingRecipes::Structures::Method method = EquatingRecipes::Structures::Method::NOT_SPECIFIED;
      EquatingRecipes::Structures::Smoothing smoothing = EquatingRecipes::Structures::Smoothing::NOT_SPECIFIED;
      std::string methodCode = "";
      double minimumRawScoreYct = std::numeric_limits<double>::quiet_NaN();
      double maximumRawScoreYct = std::numeric_limits<double>::quiet_NaN();
      double scoreIncrementYct = std::numeric_limits<double>::quiet_NaN();
      int lowestObservableRoundedScaledScore = std::numeric_limits<int>::quiet_NaN();   /* lowest possible rounded scale score */
      int highestObservableRoundedScaledScore = std::numeric_limits<int>::quiet_NaN();  /* highest possible rounded scale score */
      double weightSyntheticPopulation1 = std::numeric_limits<double>::quiet_NaN();        /* weight for synthetic pop 1 */
      bool isInternalAnchor = false;                    /* = false (external); = true (internal) */
      double reliabilityCommonItemsPopulation1 = std::numeric_limits<double>::quiet_NaN(); /* reliability of common items in pop 1 */
      double reliabilityCommonItemsPopulation2 = std::numeric_limits<double>::quiet_NaN(); /* reliability of common items in pop 2 */
      std::vector<std::string> methods;
      double minimumScoreX = std::numeric_limits<double>::quiet_NaN();              /* min score for x */
      double maximumScoreX = std::numeric_limits<double>::quiet_NaN();              /* max score for x */
      double scoreIncrementX = std::numeric_limits<double>::quiet_NaN();            /* increment between adjacent scores for x */
      Eigen::VectorXd scoreFrequenciesX; /* fd for new form x */
      size_t numberOfExaminees = 0;      
      size_t numberOfBootstrapReplications = 0;                                     /* number of replications for bootstrap */
      size_t bootstrapReplicationNumber = 0;                                    /* rep number for bootstrap; set to 0 for actual equating */
      size_t roundToNumberOfDecimalPlaces = 0;

      std::optional<EquatingRecipes::Structures::UnivariateStatistics> summaryRawDataX; /* structure for summary raw data for x and v */
      std::optional<EquatingRecipes::Structures::UnivariateStatistics> summaryRawDataY; /* structure for summary raw data for y and v */
      std::optional<EquatingRecipes::Structures::BivariateStatistics> summaryRawDataXV; /* structure for summary raw data for x and v */
      std::optional<EquatingRecipes::Structures::BivariateStatistics> summaryRawDataYV; /* structure for summary raw data for y and v */
      std::optional<EquatingRecipes::Structures::BivariateStatistics> summaryRawDataXY; /* structure for summary raw data for x and y */

      std::optional<EquatingRecipes::Structures::RawToScaledScoreTable> rawToScaledScoreTable; /* conversion table for Y */
      std::optional<EquatingRecipes::Structures::BetaBinomialSmoothing> betaBinomalSmoothingX;                /* structure for beta binomial smoothing for x */
      std::optional<EquatingRecipes::Structures::BetaBinomialSmoothing> betaBinomalSmoothingY;                /* structure for beta binomial smoothing for y */
      std::optional<EquatingRecipes::Structures::UnivariateLogLinearSmoothing> univariateLogLinearSmoothingX; /* structure for univ log-lin smoothing for x */
      std::optional<EquatingRecipes::Structures::UnivariateLogLinearSmoothing> univariateLogLinearSmoothingY; /* structure for univ log-lin smoothing for y */
      std::optional<EquatingRecipes::Structures::BivariateLogLinearSmoothing> bivariateLogLinearSmoothingXV;  /* struc for biv log-lin smoothing for x & v */
      std::optional<EquatingRecipes::Structures::BivariateLogLinearSmoothing> bivariateLogLinearSmoothingYV;  /* struc for biv log-lin smoothing for y & v */
      std::optional<EquatingRecipes::Structures::BivariateLogLinearSmoothing> bivariateLogLinearSmoothingXY;  /* struc for biv log-lin smoothing for x & y */
      std::optional<EquatingRecipes::Structures::CubicSplinePostsmoothing> cubicSplinePostsmoothing;          /* structure for cubic-spline postsmoothing */
      std::optional<EquatingRecipes::Structures::IRTInput> irtInput;                                          /* structure for IRT input */
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif