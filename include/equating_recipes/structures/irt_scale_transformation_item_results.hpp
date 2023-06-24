/* 
  From Source: IRTst.h
  Original Struct: IRTstControl
  Description: 
*/

#ifndef IMPLEMENTATION_IRT_SCALE_TRANSFORMATION_ITEM_RESULTS_HPP
#define IMPLEMENTATION_IRT_SCALE_TRANSFORMATION_ITEM_RESULTS_HPP

#include <Eigen/Core>

#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/item_specification.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTScaleTransformationItemResults {
      int itemID;
      EquatingRecipes::Structures::IRTModel irtModel;
      size_t numberOfCategories;
      double scalingConstant;
      Eigen::VectorXd scoringFunctionValues;
      Eigen::VectorXd transformedA;
      Eigen::VectorXd transformedB;
      Eigen::VectorXd transformedC;

      void configure(const EquatingRecipes::Structures::ItemSpecification& item) {
        this->irtModel = item.irtModel;
        this->itemID = item.itemID;
        this->numberOfCategories = item.numberOfCategories;
        this->scalingConstant = item.scalingConstant;
        this->scoringFunctionValues = item.scoringFunctionValues;

        if (item.a.size() >= 1) {
          this->transformedA.resize(item.a.size());
        }

        if (item.b.size() >= 1) {
          this->transformedB.resize(item.b.size());
        }

        if (item.c.size() >= 1) {
          this->transformedC.resize(item.c.size());
        }
      }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif