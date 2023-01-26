#ifndef STRUCTURES_DESIGN_HPP
#define STRUCTURES_DESIGN_HPP

// Original coding: 'R' = RG; 'S' = SG; 'C' = CINEG

namespace EquatingRecipes {
  namespace Structures {
    enum class Design {
       RandomGroups,
       SingleGroup,
       CommonItenNonEquivalentGroups      
    };
  }
}

#endif