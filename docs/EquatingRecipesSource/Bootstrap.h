/* 
  Bootstrap.h --- header for Bootstrap.c
*/

#include "ERutilities.h"

#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

void Wrapper_Bootstrap(struct PDATA *inall, int nrep, long *idum, 
                      struct BOOT_ERAW_RESULTS *t, struct BOOT_ESS_RESULTS *u);
void Boot_initialize_eraw(struct PDATA *inall, struct BOOT_ERAW_RESULTS *t);
void Boot_USTATS(struct USTATS *x, long *idum, int rep, struct USTATS *xb);
void Boot_BSTATS(struct BSTATS *xv, long *idum, int rep, struct BSTATS *s);
void Boot_accumulate_eraw(struct PDATA *inall, struct ERAW_RESULTS *b, 
                          struct BOOT_ERAW_RESULTS *t);
void Boot_se_eraw(struct PDATA *inall, struct BOOT_ERAW_RESULTS *t);
void Print_Boot_se_eraw(FILE *fp, char tt[], struct PDATA *inall,
                        struct ERAW_RESULTS *r, 
                        struct BOOT_ERAW_RESULTS *t, int mdiff);
void Boot_initialize_ess(struct PDATA *inall, struct BOOT_ESS_RESULTS *u);
void Boot_accumulate_ess(struct PDATA *inall, struct ESS_RESULTS *s,
                         struct BOOT_ESS_RESULTS *u);
void Boot_se_ess(struct PDATA *inall, struct BOOT_ESS_RESULTS *u);
void Print_Boot_se_ess(FILE *fp, char tt[], struct PDATA *inall,
                       struct ESS_RESULTS *s, 
                       struct BOOT_ESS_RESULTS *u, int mdiff);
void Parametric_boot_univ_BB(struct BB_SMOOTH *x, long *idum, int rep, 
                             struct BB_SMOOTH *btx);
void Parametric_boot_univ_ULL(struct ULL_SMOOTH *x, long *idum, int rep, 
                              struct ULL_SMOOTH *btx);
void Parametric_boot_biv(struct BLL_SMOOTH *xv, long *idum, int rep, 
                         struct BLL_SMOOTH *btxv);

#endif 