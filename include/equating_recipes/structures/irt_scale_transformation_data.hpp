/* 
  From Source: IRTst.h
  Original Struct: IRTstControl
  Description: 
*/

#ifndef STRUCTURES_IRT_SCALE_TRANSFORMATION_DATA_HPP
#define STRUCTURES_IRT_SCALE_TRANSFORMATION_DATA_HPP

#include <map>
#include <optional>
#include <set>
#include <Eigen/Core>

#include <equating_recipes/structures/common_item_specification.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/quadrature.hpp>
#include <equating_recipes/structures/symmetry.hpp>
#include <equating_recipes/structures/irt_scale_transformation_item_results.hpp>
#include <equating_recipes/structures/irt_scale_transformation_method.hpp>
#include <equating_recipes/structures/optimization_results.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTScaleTransformationData {
      // double minimumRawScoreNewForm;
      // double maximumRawScoreNewForm;
      // double rawScoreIncrementNewForm;

      // double minimumRawScoreOldForm;
      // double maximumRawScoreOldForm;
      // double rawScoreIncrementOldForm;

      // Input
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, int> maximumNumberOfIterations;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> maximumAbsoluteChangeInFunctionValue;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> maximumRelativeChangeInFunctionValue;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> maximumAbsoluteChangeInParameterValues;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> maximumRelativeChangeInParameterValues;

      EquatingRecipes::Structures::Quadrature quadratureNewForm;
      EquatingRecipes::Structures::Quadrature quadratureOldForm;

      std::vector<EquatingRecipes::Structures::ItemSpecification> newItems;
      std::vector<EquatingRecipes::Structures::ItemSpecification> oldItems;
      std::vector<EquatingRecipes::Structures::CommonItemSpecification> commonItems;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, EquatingRecipes::Structures::Symmetry> symmetryOptions;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, bool> standardizations;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> slopeStartingValue;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> interceptStartingValue;

      std::set<EquatingRecipes::Structures::IRTScaleTransformationMethod> irtScaleTranformationMethods;

      // Results
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> slopeEstimate;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, double> interceptEstimate;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, EquatingRecipes::Structures::Quadrature> transformedQuadratureNewForm;
      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, std::vector<EquatingRecipes::Structures::IRTScaleTransformationItemResults>> itemResultsNewForm;

      std::map<EquatingRecipes::Structures::IRTScaleTransformationMethod, EquatingRecipes::Structures::OptimizationResults> optimizationResults;
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif