#ifndef STRUCTURES_LOSS_SPEC_HPP
#define STRUCTURES_LOSS_SPEC_HPP

// Original coding:
//   enum LossSpec {mix_ha, mix_sl};

namespace EquatingRecipes {
  namespace Structures {
    enum class LossSpec {
      NOT_SPECIFIED,
      MIX_HA,
      MIX_SL
    };
  }
}

#endif