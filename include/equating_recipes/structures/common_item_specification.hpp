/* 
  From Source: IRTst.h
  Original Struct: ITEM_SPEC
  Description: 
*/

#ifndef STRUCTURES_COMMON_ITEM_SPECIFICATION_HPP
#define STRUCTURES_COMMON_ITEM_SPECIFICATION_HPP

#include <Eigen/Core>

#include <equating_recipes/structures/irt_model.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct CommonItemSpecification {
      int NewID;
      int OldID;
      size_t numberOfCategories;
      double scaleConstant;
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
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif