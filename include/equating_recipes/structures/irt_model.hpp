/* 
  From Source: IRTst.h
  Original Enum: ModelSpec
  Description: 
*/

#ifndef STRUCTURES_IRT_MODEL_HPP
#define STRUCTURES_IRT_MODEL_HPP

namespace EquatingRecipes {
  namespace Structures {
    enum class IRTModel {
      NOT_SPECIFIED,
      THREE_PARAMETER_LOGISTIC,
      GRADED_RESPONSE,
      PARTIAL_CREDIT,
      NOMINAL_RESPONSE
    };
  }
}

#endif