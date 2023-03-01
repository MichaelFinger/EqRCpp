#ifndef STRUCTURES_DESIGN_HPP
#define STRUCTURES_DESIGN_HPP

// Original coding: 'R' = RG; 'S' = SG; 'C' = CINEG, 'B' (single group w/ counter balance)

namespace EquatingRecipes {
  namespace Structures {
    enum class Design {
       RandomGroups,
       SingleGroup,
       CommonItenNonEquivalentGroups,
       SingleGroupCounterBalance      
    };
  }
}

#endif