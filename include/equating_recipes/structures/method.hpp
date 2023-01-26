#ifndef STRUCTURES_METHOD_HPP
#define STRUCTURES_METHOD_HPP

// Original coding:
//   'M' = mean; 'L' = lin; 'E' equi (for RG and CG designs); 
// 	 'For CG design,
// 		  'E' = Freq esti FE) with Braun-Holland (BH-FE) 
//        'F' = Modified freq est (MFE) with Braun-Holland (BH-MFE) 
//        'G' = FE + BH-FE + MFE + BH-MFE
//        'C' = Chained
// 		  'H' = FE + BH-FE + Chained
//      'A' = FE + BH-FE + MFE + BH-MFE + Chained

namespace EquatingRecipes {
  namespace Structures {
    enum class Method {
      MEAN,
      LINEAR,
      EQUIPERCENTILE,
      FREQUENCY_ESTIMATION_BRUAN_HOLLAND,
      MODIFIED_FREQUENCY_ESTIMATION_BRUAN_HOLLAND,
      CHAINED
    };
  }
}

#endif