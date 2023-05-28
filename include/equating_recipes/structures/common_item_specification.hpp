/* 
  From Source: IRTst.h
  Original Struct: ITEM_SPEC
  Description: 
*/

#ifndef STRUCTURES_COMMON_ITEM_SPECIFICATION_HPP
#define STRUCTURES_COMMON_ITEM_SPECIFICATION_HPP

#include <Eigen/Core>

#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/item_specification.hpp>

namespace EquatingRecipes {
  namespace Structures {
    class CommonItemSpecification {
    public:
      int newID;
      int oldID;
      size_t numberOfCategories;
      double scalingConstant;
      Eigen::VectorXd scoringFunctionValues;
      Eigen::VectorXd newA;
      Eigen::VectorXd newB;
      Eigen::VectorXd newC;
      Eigen::VectorXd newD;
      Eigen::VectorXd oldA;
      Eigen::VectorXd oldB;
      Eigen::VectorXd oldC;
      Eigen::VectorXd oldD;
      EquatingRecipes::Structures::IRTModel irtModel;

      static CommonItemSpecification build(const EquatingRecipes::Structures::ItemSpecification& oldItemSpec,
                                           const EquatingRecipes::Structures::ItemSpecification& newItemSpec) {
        CommonItemSpecification commonItemSpec;

        commonItemSpec.oldID = oldItemSpec.itemID;
        commonItemSpec.newID = newItemSpec.itemID;
        commonItemSpec.numberOfCategories = newItemSpec.numberOfCategories;
        commonItemSpec.scalingConstant = newItemSpec.scalingConstant;
        commonItemSpec.scoringFunctionValues = newItemSpec.scoringFunctionValues;
        commonItemSpec.irtModel = newItemSpec.irtModel;
        commonItemSpec.newA = newItemSpec.a;
        commonItemSpec.newB = newItemSpec.b;
        commonItemSpec.newC = newItemSpec.c;
        commonItemSpec.newD = newItemSpec.d;
        commonItemSpec.oldA = oldItemSpec.a;
        commonItemSpec.oldB = oldItemSpec.b;
        commonItemSpec.oldC = oldItemSpec.c;
        commonItemSpec.oldD = oldItemSpec.d;

        return commonItemSpec;
      }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif