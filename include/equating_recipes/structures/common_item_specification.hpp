/* 
  From Source: IRTst.h
  Original Struct: ITEM_SPEC
  Description: 
*/

#ifndef IMPLEMENTATION_COMMON_ITEM_SPECIFICATION_HPP
#define IMPLEMENTATION_COMMON_ITEM_SPECIFICATION_HPP

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

        commonItemSpec.newA.setZero(commonItemSpec.numberOfCategories + 1);
        commonItemSpec.newB.setZero(commonItemSpec.numberOfCategories + 1);
        commonItemSpec.newC.setZero(commonItemSpec.numberOfCategories + 1);
        commonItemSpec.newD.setZero(commonItemSpec.numberOfCategories + 1);

        commonItemSpec.oldA.setZero(commonItemSpec.numberOfCategories + 1);
        commonItemSpec.oldB.setZero(commonItemSpec.numberOfCategories + 1);
        commonItemSpec.oldC.setZero(commonItemSpec.numberOfCategories + 1);
        commonItemSpec.oldD.setZero(commonItemSpec.numberOfCategories + 1);
        
        switch (commonItemSpec.irtModel) {
          case EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC:
            commonItemSpec.newA(2) = newItemSpec.a(2);
            commonItemSpec.newB(2) = newItemSpec.b(2);
            commonItemSpec.newC(2) = newItemSpec.c(2);
            commonItemSpec.oldA(2) = oldItemSpec.a(2);
            commonItemSpec.oldB(2) = oldItemSpec.b(2);
            commonItemSpec.oldC(2) = oldItemSpec.c(2);

            break;
          case EquatingRecipes::Structures::IRTModel::GRADED_RESPONSE:
            commonItemSpec.newA(2) = newItemSpec.a(2);
            commonItemSpec.oldA(2) = oldItemSpec.b(2);

            for (size_t k = 2; k <= commonItemSpec.numberOfCategories; k++) {
              commonItemSpec.newB(k) = newItemSpec.b(k);
              commonItemSpec.oldB(k) = oldItemSpec.b(k);
            }

            break;
          case EquatingRecipes::Structures::IRTModel::PARTIAL_CREDIT:
            commonItemSpec.newA(2) = newItemSpec.a(2);
            commonItemSpec.newB(0) = newItemSpec.b(0);
            commonItemSpec.oldA(2) = oldItemSpec.a(2);
            commonItemSpec.oldB(0) = oldItemSpec.b(0);

            for (size_t k = 1; k <= commonItemSpec.numberOfCategories; k++) {
              commonItemSpec.newB(k) = newItemSpec.b(k);
              commonItemSpec.newD(k) = newItemSpec.d(k);
              commonItemSpec.oldB(k) = oldItemSpec.b(k);
              commonItemSpec.oldD(k) = oldItemSpec.d(k);
            }
            
            break;
          case EquatingRecipes::Structures::IRTModel::NOMINAL_RESPONSE:
            for(size_t k = 1; k <= commonItemSpec.numberOfCategories; k++) {
              commonItemSpec.newA(k) = newItemSpec.a(k);
              commonItemSpec.newC(k) = newItemSpec.c(k);
              commonItemSpec.oldA(k) = oldItemSpec.a(k);
              commonItemSpec.oldC(k) = oldItemSpec.c(k);
            }
            
            break;
          default:
            break;
        }

        return commonItemSpec;
      }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif