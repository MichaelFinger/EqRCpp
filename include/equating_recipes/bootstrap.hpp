/* Bootstrap.c

   For different equating design/method/smoothing the only statements 
   that need to be changed are those between the 'start' and 'end' comments
   in Wrapper_Bootstrap
   
   NOTES:
   
     Boot_BSTATS and Boot_USTATS() use the following functions from NR:
       ran2() and sort()

     Parametric_boot_univ_BB(), Parametric_boot_univ_ULL(), 
	 and Parametric_boot_biv(),uses ran2()
       
     See comments for Equated_ss() for conventions used to find 
       (a) raw scores associated with locations in vectors
       (b) locations in vectors associated with raw scores
 
This file, which is part of Equating Recipes, is free softrware.
You can distribute it and/or modify it under the terms of the
GNU Lesser General Public License, version 3, as published by 
the Free Software Foundation.

This file is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY, WITHOUT EVEN THE IMPLIED WARRANTY OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License, version 3, for more details.

You should have received a copy of the GNU Lesser General Public 
License, version 3, in a ReadMe file distributed along with this
file.  If not, see <http://www.gnu.org/licenses/>   

Copyright 2009 
Center for Advanced Studies in Measurement and Assessment (CASMA)
University of Iowa
                      
*/   

#ifndef BOOTSTRAP_HPP
#define BOOTSTRAP_HPP

#include <equating_recipes/rg_and_sg_equating.hpp>
#include <equating_recipes/cg_no_smoothing.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/beta_binomial.hpp>
#include <equating_recipes/log_linear_equating.hpp>

#include <equating_recipes/score_statistics.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class Bootstrap {
  public:

  private:
  };
}

// void Wrapper_Bootstrap(struct PDATA *inall, int nrep, long *idum, 
//                       struct BOOT_ERAW_RESULTS *t, struct BOOT_ESS_RESULTS *u);
// void Boot_initialize_eraw(struct PDATA *inall, struct BOOT_ERAW_RESULTS *t);
// void Boot_USTATS(struct USTATS *x, long *idum, int rep, struct USTATS *xb);
// void Boot_BSTATS(struct BSTATS *xv, long *idum, int rep, struct BSTATS *s);
// void Boot_accumulate_eraw(struct PDATA *inall, struct ERAW_RESULTS *b, 
//                           struct BOOT_ERAW_RESULTS *t);
// void Boot_se_eraw(struct PDATA *inall, struct BOOT_ERAW_RESULTS *t);
// void Print_Boot_se_eraw(FILE *fp, char tt[], struct PDATA *inall,
//                         struct ERAW_RESULTS *r, 
//                         struct BOOT_ERAW_RESULTS *t, int mdiff);
// void Boot_initialize_ess(struct PDATA *inall, struct BOOT_ESS_RESULTS *u);
// void Boot_accumulate_ess(struct PDATA *inall, struct ESS_RESULTS *s,
//                          struct BOOT_ESS_RESULTS *u);
// void Boot_se_ess(struct PDATA *inall, struct BOOT_ESS_RESULTS *u);
// void Print_Boot_se_ess(FILE *fp, char tt[], struct PDATA *inall,
//                        struct ESS_RESULTS *s, 
//                        struct BOOT_ESS_RESULTS *u, int mdiff);
// void Parametric_boot_univ_BB(struct BB_SMOOTH *x, long *idum, int rep, 
//                              struct BB_SMOOTH *btx);
// void Parametric_boot_univ_ULL(struct ULL_SMOOTH *x, long *idum, int rep, 
//                               struct ULL_SMOOTH *btx);
// void Parametric_boot_biv(struct BLL_SMOOTH *xv, long *idum, int rep, 
//                          struct BLL_SMOOTH *btxv);

#endif