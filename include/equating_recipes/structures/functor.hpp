#ifndef STRUCTURES_FUNCTOR_HPP
#define STRUCTURES_FUNCTOR_HPP

#include <initializer_list>

namespace EquatingRecipes {
  namespace Structures {
    template<typename T = double>
    class Functor {
    public:
      void operator()(const T& x, T& functionValue, T& functionDerivative) = 0;
    };
  }
}

#endif