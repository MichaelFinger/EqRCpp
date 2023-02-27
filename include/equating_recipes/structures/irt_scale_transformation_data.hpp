/* 
  From Source: IRTst.h
  Original Struct: IRTstControl
  Description: 
*/

#ifndef STRUCTURES_IRT_SCALE_TRANSFORMATION_DATA_HPP
#define STRUCTURES_IRT_SCALE_TRANSFORMATION_DATA_HPP

#include <optional>
#include <Eigen/Core>

#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/form_type.hpp>
#include <equating_recipes/structures/irt_equating_results.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/irt_scale_transformation_item_results.hpp>
#include <equating_recipes/structures/irt_scale_tranformation_method.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/symmetry.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTScaleTransformationData {
      double minimumRawScoreNewForm;
      double maximumRawScoreNewForm;
      double rawScoreIncrementNewForm;
      
      double minimumRawScoreOldForm;
      double maximumRawScoreOldForm;
      double rawScoreIncrementOldForm;
      
      EquatingRecipes::Structures::Quadrature quadratureNewForm;
      EquatingRecipes::Structures::Quadrature quadratureOldForm;
      
      std::vector<EquatingRecipes::Structures::ItemSpecification> newItems;
      std::vector<EquatingRecipes::Structures::ItemSpecification> oldItems;
      std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems;

      bool runHaebara;
      EquatingRecipes::Structures::Symmetry haebaraSymmetryOption;
      bool haebaraFunctionStandardization;
      std::optional<double> haebaraSlopeStartingValue;
      std::optional<double> haebaraInterceptStartingValue;

      bool runStockingLord;
      EquatingRecipes::Structures::Symmetry stockingLordSymmetryOption;
      bool stockingLordFunctionStandardization;
      std::optional<double> stockingLordSlopeStartingValue;
      std::optional<double> stockingLordInterceptStartingValue;

      std::optional<double> haebaraSlope;
      std::optional<double> haebaraIntercept;
      
      std::optional<double> stockingLordSlope;
      std::optional<double> stockingLordIntercept;
      
      double meanMeanSlope;
      double meanMeanIntercept;
      
      double meanSigmaSlope;
      double meanSigmaIntercept;

      EquatingRecipes::Structures::Quadrature transformedQuadratureNewForm;
      std::vector<EquatingRecipes::Structures::IRTScaleTransformationItemResults> itemResultsNewForm;

      EquatingRecipes::Structures::IRTScaleTranformationMethod irtScaleTranformationMethod;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif