#ifndef IMPLEMENTATION_DESIGN_HPP
#define IMPLEMENTATION_DESIGN_HPP

// Original coding: 'R' = RG; 'S' = SG; 'C' = CINEG, 'B' (single group w/ counter balance)

namespace EquatingRecipes {
  namespace Structures {
    enum class Design {
      NOT_SPECIFIED,
      RANDOM_GROUPS,
      SINGLE_GROUP,
      COMMON_ITEN_NON_EQUIVALENT_GROUPS,
      SINGLE_GROUP_COUNTER_BALANCE      
    };
  }
}

#endif