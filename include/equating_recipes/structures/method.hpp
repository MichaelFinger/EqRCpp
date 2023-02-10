#ifndef STRUCTURES_METHOD_HPP
#define STRUCTURES_METHOD_HPP

#include <map>
#include <string>

// Original coding:
//   'M' = mean; 'L' = lin; 'E' equi (for RG and CG designs); 
// 	 'For CG design,
// 		  'E' = Freq esti FE) with Braun-Holland (BH-FE) 
//      'F' = Modified freq est (MFE) with Braun-Holland (BH-MFE) 
//      'G' = FE + BH-FE + MFE + BH-MFE
//      'C' = Chained
// 		  'H' = FE + BH-FE + Chained
//      'A' = FE + BH-FE + MFE + BH-MFE + Chained

namespace EquatingRecipes {
  namespace Structures {
    enum class Method {
      MEAN,
      LINEAR,
      EQUIPERCENTILE,
      FE_BH,
      MFE_BH,
      FE_BH_MFE_BH,
      CHAINED,
      FE_BH_CHAINED,
      FE_BH_MFE_BH_CHAINED
    };
  }
}

#endif