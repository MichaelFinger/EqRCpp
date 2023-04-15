/* 
  From Source: IRTst.h
  Original Struct: ITEM_SPEC
  Description: 
*/

#ifndef STRUCTURES_ITEM_SPECIFICATION_HPP
#define STRUCTURES_ITEM_SPECIFICATION_HPP

#include <Eigen/Core>

#include <equating_recipes/structures/irt_model.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct ItemSpecification {
      int itemID;
      unsigned long numberOfCategories;
      double scalingConstant;
      Eigen::VectorXd scoringFunctionValues;
      Eigen::VectorXd a;
      Eigen::VectorXd b;
      Eigen::VectorXd c;
      Eigen::VectorXd d;
      EquatingRecipes::Structures::IRTModel irtModel;
    };
  }
}

#endif