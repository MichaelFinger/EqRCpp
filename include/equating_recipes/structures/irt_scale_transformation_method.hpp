/* 
  Based on Source: IRTst.c
*/

#ifndef STRUCTURES_IRT_SCALE_TRANFORMATION_METHOD_HPP
#define STRUCTURES_IRT_SCALE_TRANFORMATION_METHOD_HPP

namespace EquatingRecipes {
  namespace Structures {
    enum class IRTScaleTransformationMethod {
      NONE,
      HAEBARA,
      MEAN_MEAN,
      MEAN_SIGMA,
      STOCKING_LORD     
    };
  }
}

#endif