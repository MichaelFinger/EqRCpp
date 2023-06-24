/* 
  Based on Source: IRTst.c
*/

#ifndef IMPLEMENTATION_IRT_SCALE_TRANFORMATION_METHOD_HPP
#define IMPLEMENTATION_IRT_SCALE_TRANFORMATION_METHOD_HPP

namespace EquatingRecipes {
  namespace Structures {
    enum class IRTScaleTransformationMethod {
      NOT_SPECIFIED,
      HAEBARA,
      MEAN_MEAN,
      MEAN_SIGMA,
      STOCKING_LORD     
    };
  }
}

#endif