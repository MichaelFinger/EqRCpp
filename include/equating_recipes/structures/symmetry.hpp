#ifndef STRUCTURES_SYMMETRY_HPP
#define STRUCTURES_SYMMETRY_HPP

// Original coding:
//   enum symmetry {old_scale, new_scale, symmetric};

namespace EquatingRecipes {
  namespace Structures {
    enum class Symmetry {
      OLD_SCALE, 
      NEW_SCALE, 
      SYMMETRIC
    };
  }
}

#endif