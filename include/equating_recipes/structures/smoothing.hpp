#ifndef STRUCTURES_SMOOTHING_HPP
#define STRUCTURES_SMOOTHING_HPP

// Original coding: 'N'= no; 'L' = loglin; 'B' = betabin; 'S' = cub spl; 'K' = kernel; 'Z' = CLL

namespace EquatingRecipes {
  namespace Structures {
    enum class Smoothing {
       NOT_SPECIFIED,
       LOG_LINEAR,
       BETA_BINOMIAL,
       CUBIC_SPLINE,
       KERNEL,
       CONTINUIZED_LOG_LINEAR_EQUATING
    };
  }
}

#endif