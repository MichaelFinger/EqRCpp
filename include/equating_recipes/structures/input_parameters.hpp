// struct PDATA{
//   /* 
//     structure that contains all input for a particular 
//     design/method/smoothing.  Used extensively in the argument lists for
//     for Wrapper and Print functions 
//   */
//   char xfname[100];                        /* name of x or xv input file */
//   char yfname[100];                        /* name of y or yv input file */
//   char xyfname[100];                            /* name of xy input file */
//   struct USTATS *x;        /* structure for summary raw data for x and v */
//   struct USTATS *y;        /* structure for summary raw data for y and v */
//   struct BSTATS *xv;       /* structure for summary raw data for x and v */
//   struct BSTATS *yv;       /* structure for summary raw data for y and v */
//   struct BSTATS *xy;       /* structure for summary raw data for x and y */
//   char design;                       /* 'R' = RG; 'S' = SG; 'C' = CINEG  */
//   char method; /* 'M' = mean; 'L' = lin; 'E' equi (for RG and CG designs); 
// 			      'For CG design,
// 			      'E' = Freq esti FE) with Braun-Holland (BH-FE) 
//                   'F' = Modified freq est (MFE) with Braun-Holland (BH-MFE) 
//                   'G' = FE + BH-FE + MFE + BH-MFE
//                   'C' = Chained
// 			      'H' = FE + BH-FE + Chained
//                   'A' = FE + BH-FE + MFE + BH-MFE + Chained */
//   char smoothing; /* 'N'= no; 'L' = loglin; 'B' = betabin; 'S' = cub spl;
// 				                                 'K' = kernel; 'Z' = CLL */
//   double w1;                               /* weight for synthetic pop 1 */
//   int anchor;                          /* = 0 (external); = 1 (internal) */
//   double rv1;                    /* reliability of common items in pop 1 */
//   double rv2;                    /* reliability of common items in pop 2 */
//   int nm;                                           /* number of methods */
//   char **names;                                      /* names of methods */
//   double min;                                         /* min score for x */
//   double max;                                         /* max score for x */
//   double inc;                 /* increment between adjacent scores for x */
//   int *fdx;                                         /* fd for new form x */
//   int n;                           /* number of examinees for new form x */
//   double minp;      /* min raw score for yct-- see ReadSSConvTableForY() */
//   double maxp;      /* max raw score for yct-- see ReadSSConvTableForY() */
//   double incp;      /* raw score inc for yct-- see ReadSSConvTableForY() */
//   char nameyct[100];                  /* name of file containing yct[][] */
//   double **yct;                                /* conversion table for Y */
//   int round;                   /* if round = i, then round to i-th place */
//   int lprss;                      /* lowest possible rounded scale score */
//   int hprss;                     /* highest possible rounded scale score */
//   int nrep;                      /* number of replications for bootstrap */
//   int rep;     /* rep number for bootstrap; set to 0 for actual equating */
//   struct BB_SMOOTH *bbx;  /* structure for beta binomial smoothing for x */
//   struct BB_SMOOTH *bby;  /* structure for beta binomial smoothing for y */ 
//   struct ULL_SMOOTH *ullx; /* structure for univ log-lin smoothing for x */ 
//   struct ULL_SMOOTH *ully; /* structure for univ log-lin smoothing for y */
//   struct BLL_SMOOTH *bllxv; /* struc for biv log-lin smoothing for x & v */ 
//   struct BLL_SMOOTH *bllyv; /* struc for biv log-lin smoothing for y & v */
//   struct BLL_SMOOTH *bllxy; /* struc for biv log-lin smoothing for x & y */
//   struct CS_SMOOTH *cs;      /* structure for cubic-spline postsmoothing */
//   struct IRT_INPUT *IRT_Input;                /* structure for IRT input */
// };

/* 
  From Source: ERutilities.h 
  Original Struct: PDATA
  Description: structure that contains all input for a particular 
               design/method/smoothing.  Used extensively in the argument lists for
               for Wrapper and Print functions 
*/

#ifndef STRUCTURES_INPUT_PARAMETERS_HPP
#define STRUCTURES_INPUT_PARAMETERS_HPP

#include <string>
#include <vector>
#include <Eigen/Core>

#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/cubic_spline_post_smoothing.hpp>
#include <equating_recipes/structures/irt_input.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct InputParameters {
      std::string xInputFilename;                                                         /* name of x or xv input file */ 
      std::string yInputFilename;                                                         /* name of y or yv input file */ 
      std::string xyInputFilename;                                                             /* name of xy input file */ 
      EquatingRecipes::Structures::UnivariateStatistics xRawDataSummary;  /* structure for summary raw data for x and v */                                       
      EquatingRecipes::Structures::UnivariateStatistics yRawDataSummary;  /* structure for summary raw data for y and v */                                       
      EquatingRecipes::Structures::BivariateStatistics xvRawDataSummary;  /* structure for summary raw data for x and v */                                       
      EquatingRecipes::Structures::BivariateStatistics yvRawDataSummary;  /* structure for summary raw data for y and v */                                       
      EquatingRecipes::Structures::BivariateStatistics xyRawDataSummary;  /* structure for summary raw data for x and y */                                       
      EquatingRecipes::Structures::Design design;                                   /* 'R' = RG; 'S' = SG; 'C' = CINEG  */                
      std::vector<EquatingRecipes::Structures::Method> method;            /* 'M' = mean; 'L' = lin; 'E' equi (for RG and CG designs); 
                                                                  			      'For CG design,
                                                                  			      'E' = Freq esti FE) with Braun-Holland (BH-FE) 
                                                                                    'F' = Modified freq est (MFE) with Braun-Holland (BH-MFE) 
                                                                                    'G' = FE + BH-FE + MFE + BH-MFE
                                                                                    'C' = Chained
                                                                  			      'H' = FE + BH-FE + Chained
                                                                                    'A' = FE + BH-FE + MFE + BH-MFE + Chained */
      EquatingRecipes::Structures::Smoothing smoothing; /* 'N'= no; 'L' = loglin; 'B' = betabin; 'S' = cub spl; 'K' = kernel; 'Z' = CLL */

    };
  }
}

#endif