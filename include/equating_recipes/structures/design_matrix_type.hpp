#ifndef EQUATING_RECIPES_STRUCTURES_DESIGN_MATRIX_TYPE_HPP
#define EQUATING_RECIPES_STRUCTURES_DESIGN_MATRIX_TYPE_HPP

namespace EquatingRecipes {
  namespace Structures {
    // DesignMatrixType (Btype): type of design matrix and, hence, type of moments
    //   for criterion:
    //     SOLUTION --> use B (could be scaled or unscaled as indicated
    //                  in design_matrix()) --- see note below
    //     RAW_SCORE --> use B_raw and central moments based on it
    enum class DesignMatrixType {
      SOLUTION,
      RAW_SCORE
    };
  }
} // namespace EquatingRecipes

#endif