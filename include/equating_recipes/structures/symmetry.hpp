#ifndef IMPLEMENTATION_SYMMETRY_HPP
#define IMPLEMENTATION_SYMMETRY_HPP

// Original coding:
//   enum symmetry {old_scale, new_scale, symmetric};

namespace EquatingRecipes {
  namespace Structures {
    enum class Symmetry {
      NOT_SPECIFIED,
      OLD_SCALE, 
      NEW_SCALE, 
      SYMMETRIC
    };
  }
}

#endif