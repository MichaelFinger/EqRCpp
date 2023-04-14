/* 
  From Source: IRTst.h
  Original Struct: IRTstControl
  Description: 
*/

#ifndef STRUCTURES_IRT_SCALE_TRANSFORMATION_DATA_HPP
#define STRUCTURES_IRT_SCALE_TRANSFORMATION_DATA_HPP

#include <map>
#include <optional>
#include <Eigen/Core>

#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/structures/irt_scale_transformation_item_results.hpp>
#include <equating_recipes/structures/irt_scale_transformation_method.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTScaleTransformationData {
      // double minimumRawScoreNewForm;
      // double maximumRawScoreNewForm;
      // double rawScoreIncrementNewForm;

      // double minimumRawScoreOldForm;
      // double maximumRawScoreOldForm;
      // double rawScoreIncrementOldForm;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, EquatingRecipes::Structures::Quadrature> quadratureNewForm;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, EquatingRecipes::Structures::Quadrature> quadratureOldForm;

      std::vector<EquatingRecipes::Structures::ItemSpecification> newItems;
      std::vector<EquatingRecipes::Structures::ItemSpecification> oldItems;
      std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems;

      EquatingRecipes::Structures::Symmetry haebaraSymmetryOption;
      bool haebaraFunctionStandardization;

      EquatingRecipes::Structures::Symmetry stockingLordSymmetryOption;
      bool stockingLordFunctionStandardization;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> slopeStartingValue;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> interceptStartingValue;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> slopeEstimate;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> interceptEstimate;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, EquatingRecipes::Structures::Quadrature> transformedQuadratureNewForm;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, std::vector<EquatingRecipes::Structures::IRTScaleTransformationItemResults>> itemResultsNewForm;

      std::vector<EquatingRecipes::Structures::IRTScaleTransformationMethod> irtScaleTranformationMethods;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif