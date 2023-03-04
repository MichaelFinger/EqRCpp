/* 
  From Source: IRTst.h
  Original Struct: 
  Description: 
*/

#ifndef STRUCTURES_IRT_METHOD_HPP
#define STRUCTURES_IRT_METHOD_HPP

#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    enum class IRTMethod {
      NOT_SPECIFIED,
      TRUE_SCORE,
      OBSERVED_SCORE,
      TRUE_AND_OBSERVED_SCORE
    };
  }
}

#endif