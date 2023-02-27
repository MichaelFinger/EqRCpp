/* 
  From Source: IRTeq.h
  Original Struct: IRT_INPUT
  Description: Structure that contains all input for IRT equating conducted using Wrapper_IRTeq(). 
               This struct is a member of the PDATA structure
*/

#ifndef STRUCTURES_IRT_INPUT_HPP
#define STRUCTURES_IRT_INPUT_HPP

#include <vector>

#include <Eigen/Core>

#include <equating_recipes/structures/irt_method.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/irt_scale_transformation_control.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_equating_results.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTInput {
      EquatingRecipes::Structures::IRTMethod method;                                                     /* 'T' for true score; 'O' for observed score; 'A' for both */
      Eigen::VectorXd newFormFrequencyDistribution;                                       /* Actual frequency distribution for new form */
      std::vector<EquatingRecipes::Structures::ItemSpecification> newItems;
      std::vector<EquatingRecipes::Structures::ItemSpecification> oldItems;
      EquatingRecipes::Structures::IRTScaleTransformationData irtScaleTransformationData;
  
      IRTFittedDistribution newFormIRTFittedDistribution;
      IRTFittedDistribution oldFormIRTFittedDistribution;
      IRTEquatingResults irtEquatingResults;
    };
  }
}

#endif