/* 
  From Source: IRTst.h
  Original Struct: 
  Description: 
*/

#ifndef IMPLEMENTATION_IRT_METHOD_HPP
#define IMPLEMENTATION_IRT_METHOD_HPP

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