/* 
  From Source: IRTeq.h
  Original Struct: IRT_INPUT
  Description: Structure that contains all input for IRT equating conducted using Wrapper_IRTeq(). 
               This struct is a member of the PDATA structure
*/

#ifndef STRUCTURES_IRT_INPUT_HPP
#define STRUCTURES_IRT_INPUT_HPP

#include <Eigen/Core>

#include <equating_recipes/structures/irt_method.hpp>
#include <equating_recipes/structures/item_specification.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct IRTInput {
      IRTMethod method;      /* 'T' for true score; 'O' for observed score; 'A' for both */
      Eigen::VectorXd newFormFrequencyDistribution;                               /* Actual frequency distribution for new form */
  
  struct ItemSpec *NewItems;
  struct ItemSpec *OldItems;
  struct IRTstControl *stControl;
  struct RawFitDist *NewForm;
  struct RawFitDist *OldForm;
  struct RawTruObsEquiv *RawEq;
    };
  }
}

#endif