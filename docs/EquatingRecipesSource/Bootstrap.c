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
#include "NRutilities.h" 
#include "RGandSG_NoSmooth.h"
#include "CG_NoSmooth.h"
#include "CG_EquiEquate.h" 
#include "BetaBinomial.h"
#include "Bootstrap.h"
#include "LogLinear.h"
#include "ERutilities.h"
                    
/******************************************************************************/

void Wrapper_Bootstrap(struct PDATA *inall, int nrep, long *idum, 
                       struct BOOT_ERAW_RESULTS *t, struct BOOT_ESS_RESULTS *u)
/*
  Input
  
    struct PDATA inall = input for actual equating
    nrep = number of bootstrap replications
    idum = seed                         
      
  Output
  
    struct BOOT_ERAW_RESULTS t = boot results for raw scores
    struct BOOT_ESS_RESULTS u = boot results for scale scores 
                           (NULL indicates no bootstrap for scale scores)
                           
  NOTES:
    (a) struct PDATA inall is the input associated with the actual equating
        (starts out with inall->rep = 0); 
    (b) for any bootstrap replication, everything in struct PDATA inall
        remains the same except the USTATS and BSTATS structures 
        and the replication number
          (1) SG design requires one bootstrap struct BSTATS
          (2) RG design requires two bootstrap struct USTATS
          (3) CI design requires two bootstrap struct BSTATS 
    (c) For different equating design/method/smoothing the only statements 
        that need to be changed are those between 
        the '*****start...*****' and '*****end...******' comments   
       
  Function calls other than C or NR utilities:
    Boot_USTATS()
    Boot_BSTATS()  
    Wrapper_RN()
    Wrapper_SN()
    Wrapper_CN()
    Wrapper_RB()
	Wrapper_RL()    
	Wrapper_SL()
    Parametric_boot_univ_BB()
    Parametric_boot_univ_ULL()
	Parametric_boot_biv()
	Boot_initialize_eraw()
    Boot_accumulate_eraw()
    Boot_accumulate_ess()
    Boot_se_eraw()
    Boot_se_ess()
	Wrapper_ESS()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08           
*/
{
  int rep;
  struct BSTATS xvb,yvb;                /* for bootstrap reps and CI design */
  struct USTATS xb,yb;                  /* for bootstrap reps and RG design */
  struct BSTATS xyb;                    /* for bootstrap reps and SG design */
  struct ERAW_RESULTS br;    /* for bootstrap reps; not dependent on design */
  struct ESS_RESULTS bs;     /* for bootstrap reps; not dependent on design */
  struct BB_SMOOTH bt_bbx, bt_bby;                /* for boot reps with REB */
  struct ULL_SMOOTH bt_ullx, bt_ully;             /* for boot reps with REL */
  struct BLL_SMOOTH bt_bllxy;                     /* for boot reps with SEL */
  struct BLL_SMOOTH bt_bllxv, bt_bllyv;           /* for boot reps with CEL */
  
  inall->nrep = nrep; 
  Boot_initialize_eraw(inall,t);
  if(u != NULL) Boot_initialize_ess(inall,u);
  for(rep=1;rep<=nrep;rep++){
    inall->rep = rep;            /* rep always stored in struct PDATA inall */
       
    /***** start for different designs/methods/smoothing *****/

    /* Notes: (a) Since struct PDATA has the name "in" here,
                  the name of struct PDATA must be "in" everywhere; 
              (b) when smoothing=='N', storage allocation done in
                  boot_USTATS() and boot_BSTATS(); otherwise,
                  storage allocation done in Parametric_boot_univ_BB(),
				  Parametric_boot_univ_ULL(), or Parametric_boot_biv().
              (c) Drivers here are the same as those for actual
                  equating; here the last three parameters for all of them 
                  must be: rep,inall,&br */

     /* nonparametric bootstrap procedures (no smoothing involved) */
    
    if(inall->design == 'R' && inall->smoothing == 'N'){
      Boot_USTATS(inall->x,idum,rep,&xb);
      Boot_USTATS(inall->y,idum,rep,&yb);         
      Wrapper_RN('R',inall->method,'N',&xb,&yb,rep,inall,&br);
    }

    else if(inall->design == 'S' && inall->smoothing == 'N'){
      Boot_BSTATS(inall->xy,idum,rep,&xyb);      
      Wrapper_SN('S',inall->method,'N',&xyb,rep,inall,&br);
    }

    else if(inall->design == 'C' && inall->smoothing == 'N'){
      Boot_BSTATS(inall->xv,idum,rep,&xvb);
      Boot_BSTATS(inall->yv,idum,rep,&yvb);       
      Wrapper_CN('C',inall->method,'N',inall->w1,inall->anchor,
                inall->rv1,inall->rv2,&xvb,&yvb,rep,inall,&br);
    }

    /* parametric bootstrap procedures (involve smoothing) */

    else if(inall->design == 'R' && inall->smoothing == 'B'){ 
    
      Parametric_boot_univ_BB(inall->bbx, idum, rep, &bt_bbx);
      Parametric_boot_univ_BB(inall->bby, idum, rep, &bt_bby);     
      Wrapper_RB('R','E','B',inall->x, inall->y, &bt_bbx, &bt_bby, 
		         rep, inall, &br);
    }

	else if(inall->design == 'R' && inall->smoothing == 'L'){ 

      Parametric_boot_univ_ULL(inall->ullx, idum, rep, &bt_ullx);
      Parametric_boot_univ_ULL(inall->ully, idum, rep, &bt_ully);      
      Wrapper_RL('R','E','L',inall->x, inall->y, &bt_ullx, &bt_ully, 
		          rep, inall, &br);
    }

	else if(inall->design == 'S' && inall->smoothing == 'L'){
      Parametric_boot_biv(inall->bllxy,idum,rep,&bt_bllxy);
	  Wrapper_SL('S','E','L',inall->xy, &bt_bllxy, rep, inall, &br);
	}
       
    else if(inall->design == 'C' && inall->smoothing == 'L'){
      Parametric_boot_biv(inall->bllxv,idum,rep,&bt_bllxv);
	  Parametric_boot_biv(inall->bllyv,idum,rep,&bt_bllyv);
	  Wrapper_CL('C',inall->method,'L',inall->w1, inall->anchor, 
		         inall->rv1, inall->rv2,
		         inall->xv, inall->yv, &bt_bllxv, &bt_bllyv, 
				 rep, inall, &br);
	}
    else;
    
    /***** end for different designs/methods/smoothing *****/
                                        
    Boot_accumulate_eraw(inall,&br,t);
    if(u != NULL){
      Wrapper_ESS(inall,&br,inall->minp,inall->maxp,inall->incp,inall->nameyct,
                 inall->round,inall->lprss,inall->hprss,&bs);
      Boot_accumulate_ess(inall,&bs,u);
    }                                                  
  }                                                       /* end of rep loop */

  Boot_se_eraw(inall,t);
  if(u != NULL) Boot_se_ess(inall,u);  
  inall->rep = 0;                                      /* bootstrap concluded */

  return;
} 
                     
/******************************************************************************/

void Boot_initialize_eraw(struct PDATA *inall, struct BOOT_ERAW_RESULTS *t)
/*
  
  Input
    struct PDATA inall
      
  Output  
    struct BOOT_ERAW_RESULTS t
      double mn[][] = matrix of sum and mean for scores
      double sd[][] = matrix of sum2 and sd for scores
      double bse[] = overall bootstrap se's 

  Function calls other than C or NR utilities:
                                                
  R. L. Brennan

  Date of last revision: 6/30/08          
*/
{
  int i,j;
                                                          /* allocate storage */
  t->mn = dmatrix(0,inall->nm-1,0,loc(inall->max,inall->min,inall->inc));         
  t->sd = dmatrix(0,inall->nm-1,0,loc(inall->max,inall->min,inall->inc));    
  t->bse = dvector(0,inall->nm-1);      
                                                                /* initialize */
  for(i=0;i<inall->nm;i++){                                          
    t->bse[i] = 0.;
    for(j=0;j<=loc(inall->max,inall->min,inall->inc);j++) 
      t->mn[i][j] = t->sd[i][j] = 0.;
  }                           
}


/******************************************************************************/

void Boot_USTATS(struct USTATS *x, long *idum, int rep, struct USTATS *xb)
/*
  Based on USTATS x for actual equating, get a bootstrap sample 
  and store results in USTATS xb
  
  NOTE:  some initialization done here rather than in Boot_initialize-eraw()
         because otherwise there would have to be several versions of 
         Boot_initialize_eraw() for different designs. 
  
  Input: 
    struct USTATS x
    idum = seed
    rep = replication number
    
  Output
    struct USTATS xb:  xb must be different from x

  Function calls other than C or NR utilities:
    ran2()
    sort()
    cum_rel_freqs()
    perc_rank()
    score()
    MomentsFromFD()
    runerror()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08      
*/
{
  int i,
      lochigh = loc(x->max,x->min,x->inc),             /* loc of high score */
      sum = 0,                    /* cumulative value over cells in x->fd[] */
      k =0;                          /* counts number of elements in pvec[] */
  float *pvec,                          /* vector containing random persons */
      rnumber;                                             /* random number */
  
  if(rep == 1){            /* assign and allocate only at time of first rep */
    strcpy(xb->fname,"NONE");
    xb->n = x->n;  
    xb->min = x->min; 
    xb->max = x->max;
    xb->inc = x->inc; 
    xb->ns = x->ns;     
    xb->fd = ivector(0,lochigh); 
    xb->dbl_fd = dvector(0,lochigh);  
    xb->cfd = ivector(0,lochigh);  
    xb->rfd = dvector(0,lochigh); 
    xb->crfd = dvector(0,lochigh); 
    xb->prd = dvector(0,lochigh);       
  }
  
  pvec = vector(1,xb->n);              /* NR -- NOTE pvec starts at pvec[1] */    
  for(i=0;i<=lochigh;i++) xb->fd[i] = 0;                      /* initialize */
  /* ===== Before version 1.0 update ===== 
  for(i=1;i<=xb->n;i++)                   
    pvec[i] = (int) (xb->n*ran2(idum) + 1.);  
  sort(xb->n,pvec);
  */  
  for(i=1;i<=xb->n;i++){                       /* get random person numbers */  
    rnumber = er_random(idum); 
	pvec[i] = (int) (xb->n*rnumber + 1.); 
  } 
  er_sort(pvec,1,xb->n);   /* NR-- sort random person numbers from low to high */
  /* ===== End of version 1.0 update ===== */ 
  /* get bootstrap fd */
  for(i=0;i<=lochigh;i++){
    sum += x->fd[i];
    while ((int)(pvec[++k] + .00001) <= sum){
      xb->fd[i]++;
      if(k==xb->n) goto exit;
    }
    k--;                                                 /* gone 1 too far */
  }

  exit:  
  
  /* get actual minimum and maximum scores in bootstrap data */
                                
  for(i=0;i<=lochigh;i++)
    if(xb->fd[i]!=0){xb->mind = score(i,xb->min,xb->inc); break;}
  for(i=lochigh;i>=0;i--)
    if(xb->fd[i]!=0){xb->maxd = score(i,xb->min,xb->inc); break;}

  /* get cfd, rfd, crfd, prd in bootstrap data */

  for(i=1;i<=lochigh;i++) xb->cfd[i] = xb->cfd[i-1] + xb->fd[i];
  for(i=0;i<=lochigh;i++) xb->rfd[i] = (double) xb->fd[i]/xb->n;
  cum_rel_freqs(xb->min,xb->max,xb->inc,xb->rfd,xb->crfd);

  for(i=0;i<=lochigh;i++){
    xb->prd[i] = perc_rank(xb->min,xb->max,xb->inc,
                           xb->crfd,score(i,xb->min,xb->inc));
    xb->dbl_fd[i] = xb->fd[i];
  }
    
  /*  get moments */
                                           
  if(MomentsFromFD(xb->min, xb->max, xb->inc, NULL, 
      xb->fd, xb->mts) != xb->n)
      runerror("\nError somewhere in Boot_USTATS()");

  free_vector(pvec,1,xb->n);

}

/******************************************************************************/

void Boot_BSTATS(struct BSTATS *xv, long *idum, int rep, struct BSTATS *xvb)
/*
  Based on BSTATS xv for actual equating, get a bootstrap sample 
  and store results in BSTATS xvb
  
  NOTE:  some initialization done here rather than in Boot_initialize-eraw()
         because otherwise there would have to be several versions of 
         Boot_initialize_eraw() for different designs. 
  
  Input: 
    sort()
    struct BSTATS xv
    idum = seed
    rep = replication number
    
  Output
    struct BSTATS xvb:  xvb must be different from xv

  Function calls other than C or NR utilities:
    ran2()
    sort()
    cum_rel_freqs()
    perc_rank()
    score()
    MomentsFromFD()
    runerror()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08      
*/
{
  int i,j,v,
      loc1 = loc(xv->max1,xv->min1,xv->inc1), /* loc of high score in var 1 */
      loc2 = loc(xv->max2,xv->min2,xv->inc2), /* loc of high score in var 2 */
      sum = 0,                /* cumulative value over cells in xv->bfd[][] */
      k =0;                          /* counts number of elements in pvec[] */
  float *pvec,                          /* vector containing random persons */
	  rnumber;                                             /* random number */
  
  if(rep == 1){            /* assign and allocate only at time of first rep */
    strcpy(xvb->fname,"NONE");
    xvb->n = xv->n;  
    xvb->min1 = xv->min1; 
    xvb->max1 = xv->max1;
    xvb->inc1 = xv->inc1;
    xvb->ns1 = xv->ns1;  
    xvb->min2 = xv->min2;
    xvb->max2 = xv->max2;
    xvb->inc2 = xv->inc2;
    xvb->ns2 = xv->ns2;    
    xvb->fd1 = ivector(0,loc1);   
    xvb->fd2 = ivector(0,loc2); 
    xvb->bfd = imatrix(0,loc1,0,loc2);  

    xvb->dbl_fd1 = dvector(0,loc1);   
    xvb->dbl_fd2 = dvector(0,loc2); 
    xvb->dbl_bfd = dmatrix(0,loc1,0,loc2);  

    xvb->cfd1 = ivector(0,loc1);
    xvb->rfd1 = dvector(0,loc1); 
    xvb->crfd1 = dvector(0,loc1); 
    xvb->prd1 = dvector(0,loc1);   

    xvb->cfd2 = ivector(0,loc2);
    xvb->rfd2 = dvector(0,loc2); 
    xvb->crfd2 = dvector(0,loc2); 
    xvb->prd2 = dvector(0,loc2);
    xvb->bp12 = dmatrix(0,loc1,0,loc2);
  }
  
  pvec = vector(1,xvb->n);              /* NR -- NOTE pvec starts at pvec[1] */    
  xvb->cov = 0.;                                           /* initialize cov */
  for(i=0;i<=loc1;i++){ 
    xvb->fd1[i] = 0;                             /* initialize row marginals */
    for(j=0;j<=loc2;j++)
      xvb->bfd[i][j] = 0;                        /* initialize bfd[][] cells */
  }
  for(j=0;j<=loc2;j++) 
    xvb->fd2[j] = 0;                             /* initialize col marginals */
  /* ===== Before version 1.0 update =====  
  for(i=1;i<=xvb->n;i++)                         
    pvec[i] = (int) (xvb->n*ran2(idum) + 1.); 
  sort(xvb->n,pvec); 
  */
  for(i=1;i<=xvb->n;i++){                       /* get random person numbers */  
    rnumber = er_random(idum); 
    pvec[i] = (int) (xvb->n*rnumber + 1.); 
  } 
  er_sort(pvec,1,xvb->n);/* NR-- sort random person numbers from low to high */
  /* ===== End of version 1.0 update ===== */
  
  /* get bootstrap bivariate and marginal distributions */
  
  for(i=0;i<=loc1;i++)
    for(j=0;j<=loc2;j++){
      sum += xv->bfd[i][j];
      while ((int)(pvec[++k] + .00001) <= sum){
        xvb->bfd[i][j]++;
        xvb->fd1[i]++;
        xvb->fd2[j]++;
        xvb->cov += (score(i,xvb->min1,xvb->inc1))*
                    (score(j,xvb->min2,xvb->inc2));
        if(k==xvb->n) goto exit;
      }
      k--;                                                 /* gone 1 too far */
    }
  exit:  
  
 /*  get actual minimum and maximum scores in bootstrap data */
                                
  for(i=0;i<=loc1;i++)
    if(xvb->fd1[i]!=0){xvb->mind1 = score(i,xvb->min1,xvb->inc1); break;}
  for(i=loc1;i>=0;i--)
    if(xvb->fd1[i]!=0){xvb->maxd1 = score(i,xvb->min1,xvb->inc1); break;}
    
  for(j=0;j<=loc2;j++)
    if(xvb->fd2[j]!=0){xvb->mind2 = score(j,xvb->min2,xvb->inc2); break;}
  for(j=loc2;j>=0;j--)
    if(xvb->fd2[j]!=0){xvb->maxd2 = score(j,xvb->min2,xvb->inc2); break;}

  /* get cfd1,rfd1,crfd1,prd1, and cfd2,rfd2,crfd2,prd2 in boot data */

  xvb->cfd1[0] = xvb->fd1[0];
  for(i=1;i<=loc1;i++) xvb->cfd1[i] = xvb->cfd1[i-1] + xvb->fd1[i];
  for(i=0;i<=loc1;i++) xvb->rfd1[i] = (double) xvb->fd1[i]/xvb->n;
  cum_rel_freqs(xvb->min1,xvb->max1,xvb->inc1,xvb->rfd1,xvb->crfd1);

  for(i=0;i<=loc1;i++){
    xvb->prd1[i] = perc_rank(xvb->min1,xvb->max1,xvb->inc1,
                           xvb->crfd1,score(i,xvb->min1,xvb->inc1));
    xvb->dbl_fd1[i] = xvb->fd1[i];
  }

  xvb->cfd2[0] = xvb->fd2[0];
  for(i=1;i<=loc2;i++) xvb->cfd2[i] = xvb->cfd2[i-1] + xvb->fd2[i];
  for(i=0;i<=loc2;i++) xvb->rfd2[i] = (double) xvb->fd2[i]/xvb->n;
  cum_rel_freqs(xvb->min2,xvb->max2,xvb->inc2,xvb->rfd2,xvb->crfd2);
  
  for(i=0;i<=loc2;i++){
    xvb->prd2[i] = perc_rank(xvb->min2,xvb->max2,xvb->inc2,
                           xvb->crfd2,score(i,xvb->min2,xvb->inc2));
    xvb->dbl_fd2[i] = xvb->fd2[i];
  }

  /*  get moments, cov, and corr in bootstrap data */
                                           
  if(MomentsFromFD(xvb->min1, xvb->max1, xvb->inc1, NULL, 
      xvb->fd1, xvb->mts1) != xvb->n)
      runerror("\nError somewhere in Boot_BSTATS()");
  if(MomentsFromFD(xvb->min2, xvb->max2, xvb->inc2, NULL, 
      xvb->fd2, xvb->mts2) != xvb->n)
      runerror("\nError somewhere in Boot_BSTATS()");
  xvb->cov = xvb->cov/xvb->n -  (xvb->mts1[0])*(xvb->mts2[0]);
  xvb->corr = xvb->cov/(xvb->mts1[1] * xvb->mts2[1]);

  /* bivariate proportions for frequency estimation, as well as
     double version of xvb->bfd[][] */             

  for(v=0;v<=xvb->ns2-1;v++)        /* v indexes second variable here */                                         
    for(i=0;i<=xvb->ns1-1;i++){
      xvb->bp12[i][v] = (double) xvb->bfd[i][v]/xvb->n; 
      xvb->dbl_bfd[i][v] = xvb->bfd[i][v];
    }

  /* deallocate space */
  
  free_vector(pvec,1,xvb->n);

}

/***********************************************************************/

void Boot_accumulate_eraw(struct PDATA *inall, struct ERAW_RESULTS *b, 
                          struct BOOT_ERAW_RESULTS *t)
/*

  accumulate b->eraw[][] in t->mn[][]
  accumulate b->eraw[][]*b->eraw[][] in t->sd[][]
  
  
  Input
     struct PDATA inall: following variables are used
       nm = number of methods 
       min = min raw score 
       max = max raw score 
       inc = raw score increment
       
     struct ERAW_RESULTS b 
       double eraw[][] --- for a boot rep, not actual equating
  
  Output   
     struct BOOT_ERAW_RESULTS s
       double **mn: for matrix of sum and mean for scores 
       double **sd; for matrix of sum2 and sd for scores

  Function calls other than C or NR utilities: 
    loc()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08       
*/
{
  int i,j;
  
  for(i=0;i<inall->nm;i++)
    for(j=0;j<=loc(inall->max,inall->min,inall->inc);j++){ 
      t->mn[i][j] += b->eraw[i][j];      
      t->sd[i][j] += b->eraw[i][j]*b->eraw[i][j];
    }  

}

/***********************************************************************/

void Boot_se_eraw(struct PDATA *inall, struct BOOT_ERAW_RESULTS *t)
/*
  final computations for s->mn[][] and s->sd[][] == bootstrap se's
  
  Input
     struct PDATA inall: following variables are used
       nm = number of methods 
       min = min raw score
       max = max raw score
       inc = raw score increment
       fdx[] = fd for x
       nrep = number of replications
       n = sample size associated with fdx[] 
       
     struct BOOT_ERAW_RESULTS t
       mn[][] = for matrix of sum's for scores 
       sd[][] = for matrix of sum2's for scores 
  
  Output   
     struct BOOT_ERAW_RESULTS t
       mn[][] = for matrix of mean for scores 
       sd[][] = for matrix of sd for scores 
       bse[] = overall bootstrap se's (one for each method)

  NOTE: bootstrap se's require a divisor of nrep-1, not nrep;
        see "divisor" statement in code

  Function calls other than C or NR utilities: 
    loc()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08  
*/
{
  int i,j;
  
  for(i=0;i<inall->nm;i++){
    for(j=0;j<=loc(inall->max,inall->min,inall->inc);j++){ 
      t->mn[i][j] /= inall->nrep;      
      t->sd[i][j] = t->sd[i][j]/inall->nrep - t->mn[i][j]*t->mn[i][j];
      t->sd[i][j] *= (double) (inall->nrep)/(inall->nrep - 1); /* divisor */
      t->bse[i] += inall->fdx[j]*t->sd[i][j];
      t->sd[i][j] = (t->sd[i][j]>0.00001) ? sqrt(t->sd[i][j]) : 0.;
    }    
    t->bse[i] = sqrt(t->bse[i]/inall->n);
  } 

}

/***********************************************************************/

void Print_Boot_se_eraw(FILE *fp, char tt[], struct PDATA *inall,
                        struct ERAW_RESULTS *r, 
                        struct BOOT_ERAW_RESULTS *t, int mdiff)
/*
  Input
  
    fp = file pointer for output
    tt[] = user-identified title
 
    struct PDATA inall: only variables used:
      nrep = number of replications
      nm = number of methods 
      names[][] = names of methods 
      min = min raw score 
      max = max raw score
      inc = raw score increment
      fdx[] = frequency distribution for new form x
     
    struct ERAW_RESULTS r;  only variable:

      eraw[][] = actual equated raw scores  
    
    struct BOOT_ERAW_RESULTS t;  variables used:
  
      mn[][] = matrix of means for equated scores 
      sd[][] = matrix of sd's = bse's for equated scores
      bse[] = overall bootstrap se 
      
    mdiff: = 1 --> print eraw[][]-mn[][], which indicates the extent to which
    equated raw scores (using original data) are above or below the mean 
    of the equated raw scores for the nrep bootstrap replications.  Ideally
    these values should be 0. Large values may indicate need for more
    bootstrap replications.

  Function calls other than C or NR utilities: 
    score()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j; 
  
  fprintf(fp,"\n\n%s\n\n",tt);
  
  fprintf(fp,"Bootstrap standard errors for raw scores\n\n");
    
  fprintf(fp,"Number of bootstrap replications = %d\n\n",inall->nrep);
      
  fprintf(fp,"                   ");  
  for(j=0;j<=inall->nm-1;j++) 
    fprintf(fp,"     Method %d:  %s",j,inall->names[j]);
     
  fprintf(fp,"\n                   ");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"  ------------------------"); 
  fprintf(fp,"\n   Score (x)  fd(x)");
  for(j=0;j<=inall->nm-1;j++) fprintf(fp,"          eraw         bse");
  fprintf(fp,"\n");

  for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
    fprintf(fp,"\n%12.5f%7d",score(i,inall->min,inall->inc),inall->fdx[i]);
    for(j=0;j<=inall->nm-1;j++){     
      fprintf(fp,"  %12.5f",r->eraw[j][i]);
      if(t->sd[j][i] >0.)
        fprintf(fp,"%12.5f",t->sd[j][i]);
      else
        fprintf(fp,"         ***");        /* means bse==0 or no freq */
    }
  }
  
  fprintf(fp,"\n\n    Ave. bse       ");
    for(j=0;j<=inall->nm-1;j++){
      fprintf(fp,"              ");
      fprintf(fp,"%12.5f",t->bse[j]);
    }
    
  fprintf(fp,"\n");
    for(j=1;j<=19+26*inall->nm;j++) fprintf(fp,"-");
  fprintf(fp,"\n *** means bse = 0");
  fprintf(fp,"\n");
    for(j=1;j<=19+26*inall->nm;j++) fprintf(fp,"-");

  
  /*if mdiff==1, print eraw[][] - mn[][] */
  
  if(mdiff==1){
    fprintf(fp,"\n\n   Score (x)  fd(x)");
    for(j=0;j<=inall->nm-1;j++) fprintf(fp,"     eraw-mean            ");
    fprintf(fp,"\n");
    
    for(i=0;i<=loc(inall->max,inall->min,inall->inc);i++){
      fprintf(fp,"\n%12.5f%7d",score(i,inall->min,inall->inc),inall->fdx[i]);
      for(j=0;j<=inall->nm-1;j++){     
        fprintf(fp,"  %12.5f",r->eraw[j][i]-t->mn[j][i]);
        fprintf(fp,"            ");
      }
    }  
  
    fprintf(fp,"\n"); 
    for(j=1;j<=19+26*inall->nm;j++) fprintf(fp,"-");

  }

}

/*************************************************************************/

void Boot_initialize_ess(struct PDATA *inall, struct BOOT_ESS_RESULTS *u)
/*
  
  Input
    
    struct PDATA inall: variables used:
      nm = number of methods
      min = minimum raw score for x 
      max = maximum raw score for x
	  inc = increment in raw scores for x
   
   Output             
    
     struct BOOT_ESS_RESULTS u: variables used:
       int rep = replication number
       mnu[][] = matrix of sum and mean for unrounded scale scores
       sdu[][] = matrix of sum2 and sd for unrounded scale scores
       bseu[] = overall bootstrap se's for unrounded scale scales 
       mnr[][] = matrix of sum and mean for rounded scale scores
       sdr[[]] = matrix of sum2 and sd for rounded scale scores
       bser[] = overall bootstrap se's for rounded scale scales

  Function calls other than C or NR utilities: 
    loc()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j,
	  ns = nscores(inall->max,inall->min,inall->inc);
                                                 /* allocate storage */
  u->mnu =  dmatrix(0,inall->nm-1,0,ns-1);      
  u->sdu =  dmatrix(0,inall->nm-1,0,ns-1); 
  u->bseu = dvector(0,inall->nm-1);                      
  u->mnr =  dmatrix(0,inall->nm-1,0,ns-1);       
  u->sdr =  dmatrix(0,inall->nm-1,0,ns-1);
  u->bser = dvector(0,inall->nm-1);                              
                                                       /* initialize */
  for(i=0;i<inall->nm;i++){                                
    u->bseu[i] = 0.;
    u->bser[i] = 0.;
    for(j=0;j<ns;j++) 
      u->mnu[i][j] = u->sdu[i][j] = u->mnr[i][j] = u->sdr[i][j] = 0.;
  }                           
}

/***********************************************************************/

void Boot_accumulate_ess(struct PDATA *inall, struct ESS_RESULTS *s,
                         struct BOOT_ESS_RESULTS *u)
/*

  accumulate s->essu[][] and store in u->mnu[][]
  accumulate s->essu[][]*s->essu[][] and store in u->sdu[][]
  accumulate s->essr[][] and store in s->mnr[][]
  accumulate s->essr[][]*s->essr[][] and store in u->sdr[][]
  
  Input
    
    struct PDATA inall: variables used:
      nm = number of methods
      min = min raw score for x 
      max = max raw score for x 
      inc = raw score increment for x          
      round: if round = i, then round to i-th place
  
     struct ESS_RESULTS s (must be same struct as in Wrapper_ESS()
                           within the bootstrap loop)
 
       essu[][] = unrounded equated scale scores for a replication
       essr[][] = rounded equated scale scores for a replication
  
  Output
     
    struct Boot_ess_results u (must be same struct as in 
                               Boot_initialize_ess())
                             
    rep = replication number for bootstrap
    mnu[][] = sum and mean for unrounded scale scores
    sdu[][] = sum2 and sd for unrounded scale scores 
    mnr[][] = sum and mean for rounded scale scores
    sdr[][] = sum2 and sd for rounded scale scores

  Function calls other than C or NR utilities: None
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j,
	   ns = nscores(inall->max,inall->min,inall->inc);
  
  for(i=0;i<inall->nm;i++)
    for(j=0;j<ns;j++){ 
      u->mnu[i][j] += s->essu[i][j];      
      u->sdu[i][j] += s->essu[i][j]*s->essu[i][j];
      if(inall->round>0){
        u->mnr[i][j] += s->essr[i][j];      
        u->sdr[i][j] += s->essr[i][j]*s->essr[i][j];
      }
    }  

}

/***********************************************************************/

void Boot_se_ess(struct PDATA *inall, struct BOOT_ESS_RESULTS *u)
/*
  final computations for
    u->muu[][] and u->sdu[][] == u->bseu[][] = bootstrap unrounded se's 
    u->mur[][] and u->sdr[][] == u->bser[][] = bootstrap rounded se's
      
  Input
    
    struct PDATA inall: variables used:
      nm = number of methods
      min = min raw score for x 
      max = max raw score for x
      inc = raw score increment for x
      nrep = number of replications
      fdx[] = FD for old form x
      n = sample size associated with fdx[]
      round: if round = i, then round to i-th place
    
     struct BOOT_ESS_RESULTS u (must be same struct as in 
                                Boot_initialize_ess() and
                                Boot_accumulate_ess())
  
       mnu[][] = sums for unrounded scale scores
       sdu[][] = sum2's for unrounded scale scores 
       mnr[][] = sums for rounded scale scores
       sdr[][] = sum2's for rounded scale scores
       bser[] = overall bootstrap se's for rounded scale scales
       
   Output
     struct BOOT_ESS_RESULTS u (must be same struct as in 
                                Boot_initialize_ess() and
                                Boot_accumulate_ess())
  
       mnu[][] = means for unrounded scale scores
       sdu[][] = sd's for unrounded scale scores
       bseu[] = overall bootstrap se's for unrounded scale scales 
       mnr[][] = means for rounded scale scores
       sdr[][] = sd's for rounded scale scores
       bser[] = overall bootstrap se's for rounded scale scales

  NOTE: bootstrap se's require a divisor of nrep-1, not nrep;
        see "divisor" statements in code

  Function calls other than C or NR utilities:
    loc()
    score()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j,
	  ns = nscores(inall->max,inall->min,inall->inc);
  
  for(i=0;i<inall->nm;i++){

    /* for unrounded scale scores */

    for(j=0;j<ns;j++){
      u->mnu[i][j] /= inall->nrep;      
      u->sdu[i][j] = u->sdu[i][j]/inall->nrep - u->mnu[i][j]*u->mnu[i][j];
      u->sdu[i][j] *= (double) (inall->nrep)/(inall->nrep - 1); /* divisor */

      if (u->sdu[i][j]>0.00001) u->bseu[i] += inall->fdx[j]*u->sdu[i][j];        
                 
      u->sdu[i][j] = (u->sdu[i][j]>0.00001) ? sqrt(u->sdu[i][j]) : 0.;
	}    
    u->bseu[i] = sqrt(u->bseu[i]/inall->n);

	/* for rounded scale scores */
    
    if(inall->round>0){
      for(j=0;j<ns;j++){
        u->mnr[i][j] /= inall->nrep;      
        u->sdr[i][j] = u->sdr[i][j]/inall->nrep - u->mnr[i][j]*u->mnr[i][j];
        u->sdr[i][j] *= (double) (inall->nrep)/(inall->nrep - 1); /* divisor */
        
        if (u->sdr[i][j]>0.00001)u->bser[i] += inall->fdx[j]*u->sdr[i][j]; 
             
        u->sdr[i][j] = (u->sdr[i][j]>0.00001) ? sqrt(u->sdr[i][j]) : 0.;
      }    
      u->bser[i] = sqrt(u->bser[i]/inall->n);
	}
  } 

}  

/***********************************************************************/

void Print_Boot_se_ess(FILE *fp, char tt[], struct PDATA *inall,
                       struct ESS_RESULTS *s, 
                       struct BOOT_ESS_RESULTS *u, int mdiff)                       
/*
  prints bootstrap se's for unrounded equated scale scores.
  If inall->round>0 then boot se's also printed for 
  rounded equated scale scores
  
  Input
  
    fp = file pointer for output
    tt[] = user-identified title
 
    struct PDATA inall: only variables used:
 
      nm = number of methods 
      names[][] = names of methods 
      min = min raw score
      max = max raw score 
      inc = raw score increment
      fdx[] = frequency distribution for new form x
      n = sample size associated with fdx[] 
      nameyct[100] = name of file containing yct[][]
      yct[][] = conversion table for Y
      round =: if round = i, then round to i-th place
      nrep = number of replications
     
    struct ESS_RESULTS s;  (must be same as in Wrapper_ESS()
    for actual equating): variables:

      essu[][] = unrounded equated scale scores for actual equating
      essr[][] = rounded equated scale scores for actual equating  
    
    struct BOOT_ESS_RESULTS u;  variables used:
  
      mnu[][] = means for unrounded scale scores
      sdu[][] = sd's = bse's for unrounded scale scores
      bseu[] = overall bootstrap se's for unrounded scale scales 
      mnr[][] = means for rounded scale scores
      sdr[][] = sd's = bse's for rounded scale scores
      bser[] = overall bootstrap se's for rounded scale scales
      
    mdiff: = 1 --> print essu[][]-mnu[][]. Also print essr[][]-mnr[][] 
    if round>0.  essr[][]-mnr[][] indicates the extent to which
    equated scale scores (using original data) are above or below the mean 
    of the equated scale scores for the nrep bootstrap replications.  Ideally
    these values should be 0. Large values may indicate need for more
    bootstrap replications.

  Function calls other than C or NR utilities:
    score()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/   
{
  int j,k,m,
	  ns = nscores(inall->max,inall->min,inall->inc); 
  
  for(k=1;k<=(inall->round>0 ? 2 : 1);k++){   /* beginning of loop for k */
                                       /* k=1 (unrounded); k=2 (rounded) */
    fprintf(fp,"\n\n%s\n\n",tt);
  
    fprintf(fp,"Bootstrap standard errors for %s scale scores\n\n",
      (k==1) ? "unrounded" : "rounded");
    
    fprintf(fp,"Number of bootstrap replications = %d\n\n",inall->nrep);
  
    fprintf(fp,"Name of file containing conversion table for Y: %s\n\n",
      inall->nameyct);
      
    fprintf(fp,"                   ");   
    for(m=0;m<=inall->nm-1;m++) 
      fprintf(fp,"     Method %d:  %s",m,inall->names[m]);         
    fprintf(fp,"\n                   ");
    for(m=0;m<=inall->nm-1;m++) fprintf(fp,"  ------------------------"); 
    fprintf(fp,"\n   Score (x)  fd(x)");
    for(m=0;m<=inall->nm-1;m++) fprintf(fp,"          %s         bse",
        (k==1) ? "essu" : "essr");
    fprintf(fp,"\n");

    for(j=0;j<ns;j++){       
      fprintf(fp,"\n%12.5f",score(j,inall->min,inall->inc));
      fprintf(fp,"%7d",inall->fdx[j]); 
                             
      for(m=0;m<=inall->nm-1;m++){
        if(k==1){                                            /* unrounded */
          fprintf(fp,"  %12.5f",s->essu[m][j]);       
          if(u->sdu[m][j] >0.)
            fprintf(fp,"%12.5f",u->sdu[m][j]);
          else
            fprintf(fp,"         ***");        /* means bse==0 or no freq */
        }
        else{                                                  /* rounded */
          fprintf(fp,"  %12.0f",s->essr[m][j]);       
          if(u->sdr[m][j] >0.)
            fprintf(fp,"%12.5f",u->sdr[m][j]);
          else
            fprintf(fp,"         ***");        /* means bse==0 or no freq */        
        }
	  }	  
    }
  
    fprintf(fp,"\n\n    Ave. bse       ");
      for(m=0;m<=inall->nm-1;m++){
        fprintf(fp,"              ");
        if(k==1) fprintf(fp,"%12.5f",u->bseu[m]);
        else fprintf(fp,"%12.5f",u->bser[m]);
      }
    
    fprintf(fp,"\n");
      for(m=1;m<=121;m++) fprintf(fp,"-");
    fprintf(fp,"\n *** means bse = 0"); 
    fprintf(fp,"\n");
      for(m=1;m<=121;m++) fprintf(fp,"-");
    
    /* print difference between essu[][] and mnu[][]
       or essr[][] and mnr[][] if mdiff==1 */
  
    if(mdiff==1){
      fprintf(fp,"\n\n   Score (x)  fd(x)");
      for(m=0;m<=inall->nm-1;m++) fprintf(fp,"     %s            ",
        (k==1) ? "essu-mean" : "essr-mean");
      fprintf(fp,"\n");
      
      for(j=0;j<ns;j++){                             
        fprintf(fp,"\n%12.5f",score(j,inall->min,inall->inc));                    
        fprintf(fp,"%7d",inall->fdx[j]);      
        for(m=0;m<=inall->nm-1;m++){     
          if(k==1) fprintf(fp,"  %12.5f",s->essu[m][j] - u->mnu[m][j]);
          else fprintf(fp,"  %12.5f",s->essr[m][j] - u->mnr[m][j]);
          fprintf(fp,"            ");
        }		
	  }
	}
  
    fprintf(fp,"\n");
      for(m=1;m<=121;m++) fprintf(fp,"-");

  }                                                  /* end of loop for k */

}

/******************************************************************************/

void Parametric_boot_univ_BB(struct BB_SMOOTH *x, long *idum, int rep, 
                             struct BB_SMOOTH *btx)
/*
   Parametric bootstrap sample for a univariate distribution

   Calling sequence from Wrapper_Bootstrap():
     Parametric_boot_univ_BB(inall->bbx, idum, rep, bt_bbx)

   Input
     x    = struct BB_SMOOTH for actual fitted distribution;
	        taken as population distribution
     idum = seed
	 rep  = replication number

   Output
     btx = struct BB_SMOOTH for bootstrap sample from x; only elements
	       needed are btx->crfd and btx->prd

  Function calls other than C or NR utilities:
    ran2()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j,
	  np = x->num_persons,                           /* number of persons */
	  ns = x->num_items+1;                  /* number of score categories */
  float r;

  if(rep==1){
    btx->crfd = dvector(0,ns-1);            
    btx->prd = dvector(0,ns-1); 
  }

  for(j=0;j<ns;j++) btx->crfd[j] = 0.;                  /* initialization */

  /* ===== Before version 1.0 update =====
  for(i=1;i<=np;i++){
    r = ran2(idum);
    for(j=0;j<ns;j++)
      if(r < x->crfd[j]){
        btx->crfd[j]++;                            
        break;
      }
  }
  */
  for(i=1;i<=np;i++){  
    r = er_random(idum); 
    for(j=0;j<ns;j++)
      if(r < x->crfd[j]){
        btx->crfd[j]++;                           /* increment frequencies */
        break;
      }
  }
  /* ===== End of version 1.0 update ===== */

  /* here btx->crfd[] is frequencies; 
     next two statemenst makes btx->crfd[] cum rel freqs */
 
  for(j=0;j<ns;j++) btx->crfd[j] /= np;                  /* rel freq dist */
  for(j=1;j<ns;j++) btx->crfd[j] += btx->crfd[j-1];  /* cum rel freq dist */

  for(j=0;j<ns;j++)                                   /* percentile ranks */
    btx->prd[j] = perc_rank(0, ns-1, 1, btx->crfd, j); 

  return;
}

/*************************************************************************/

void Parametric_boot_univ_ULL(struct ULL_SMOOTH *x, long *idum, int rep, 
                              struct ULL_SMOOTH *btx)
/*
   Parametric bootstrap sample for a univariate distribution

   Calling sequence from Wrapper_Bootstrap():
     Parametric_boot_univ_ULL(inall->ullx, idum, rep, bt_ullx)

   Input
     x    = struct ULL_SMOOTH for actual fitted distribution;
	        taken as population distribution
     idum = seed
	 rep  = replication number

   Output
     btx = struct ULL_SMOOTH for bootstrap sample from x; only elements
	       needed are btx->crfd and btx->prd

  Function calls other than C or NR utilities:
    ran2()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j,
	  np = x->num_persons,
	  ns = x->ns;
  float r;

  if(rep==1){
    btx->crfd = dvector(0,ns-1);            
    btx->prd = dvector(0,ns-1); 
  }

  for(j=0;j<ns;j++) btx->crfd[j] = 0.;                  /* initialization */
  /* ===== Before version 1.0 update =====
  for(i=1;i<=np;i++){
    r = ran2(idum);
    for(j=0;j<ns;j++)
      if(r < x->crfd[j]){
        btx->crfd[j]++;                          
        break;
      }
  }
  */ 
  for(i=1;i<=np;i++){ 
    r = er_random(idum); 
    for(j=0;j<ns;j++)
      if(r < x->crfd[j]){
        btx->crfd[j]++;                          /* increment frequencies */
        break;
      }
  }
  /* ===== End of version 1.0 update ===== */

  /* here btx->crfd[] is frequencies; 
     next two statemenst makes btx->crfd[] cum rel freqs */
 
  for(j=0;j<ns;j++) btx->crfd[j] /= np;                  /* rel freq dist */
  for(j=1;j<ns;j++) btx->crfd[j] += btx->crfd[j-1];  /* cum rel freq dist */

  for(j=0;j<ns;j++)                                   /* percentile ranks */
    btx->prd[j] = perc_rank(0, ns-1, 1, btx->crfd, j); 

  return;
}
/**************************************************************************/

void Parametric_boot_biv(struct BLL_SMOOTH *xv, long *idum, int rep, 
                         struct BLL_SMOOTH *btxv)
/*
   Parametric bootstrap sample for a bivariate distribution

   Input
     xv      = BLL_SMOOTH for original data (smoothed);  
	           xv->brd[][] is taken as population distribution
     idum    = seed
	 rep     = replication number

   Output
     btxv   = struct BLL_SMOOTH for bootstrap sample; not
	          all data elements are needed

  Function calls other than C or NR utilities:
    ran2()
                                                
  R. L. Brennan

  Date of last revision: 6/30/08    
*/
{
  int i,j,
	  nsx = xv->nsx,
	  nsv = xv->nsv,
	  ns = nsx*nsv,                    /* total number of score categories */
	  np = xv->num_persons;
  float r;
  double *ptr;                                         /* pointer variable */
  
  if(rep==1){
    btxv->anchor = xv->anchor;           /* = 0 (external); = 1 (internal) */ 
    btxv->num_persons = xv->num_persons;              /* number of persons */

	btxv->nsx = xv->nsx;               /* number of score categories for x */
    btxv->minx = xv->minx;                      /* minimum raw score for x */
    btxv->incx = xv->incx;                /* increment in raw scores for x */
    btxv->nsv = xv->nsv;               /* number of score categories for v */
    btxv->minv = xv->minv;                      /* minimum raw score for v */
    btxv->incv = xv->incv;                /* increment in raw scores for v */
    btxv->ns = xv->ns;                 /* total number of score categories */

    btxv->bfd = dmatrix(0,nsx-1,0,nsv-1);                 /* bootstrap bfd */

    btxv->fd_x = dvector(0,nsx-1);        /* row marg fd for bootstrap bfd */ 
    btxv->density_x = dvector(0,nsx-1);  /* row marg rfd for bootstrap bfd */ 
    btxv->crfd_x = dvector(0,nsx-1);    /* row marginal crfd bootstrap bfd */ 
    btxv->prd_x = dvector(0,nsx-1);       /* row marginal PR bootstrap bfd */ 

    btxv->fd_v = dvector(0,nsv-1);        /* row marg fd for bootstrap bfd */ 
    btxv->density_v = dvector(0,nsv-1);  /* col marg rfd for bootstrap bfd */ 
    btxv->crfd_v = dvector(0,nsv-1);    /* col marginal crfd bootstrap bfd */ 
    btxv->prd_v = dvector(0,nsv-1);       /* col marginal PR bootstrap bfd */

	btxv->crfd_vector_bfd = dvector(0,ns-1);/* vector version of btxv->bfd */
	btxv->brfd = dmatrix(0,nsx-1,0,nsv-1); /* rel freq verion of btxv->bfd */
  }


  for(j=0;j<ns;j++) btxv->crfd_vector_bfd[j] = 0.;      /* initialization */
  
  /* In the following code, xv->crfd_vector_bfd[] is the crfd  for the 
     fitted actual data, while btxv->crfd_vector_bfd[] is the bootstrap fd 
	 at end of for loop */
  /* ===== Before version 1.0 update ===== 
  for(i=1;i<=np;i++){
    r = ran2(idum);
    for(j=0;j<ns;j++)
      if(r < xv->crfd_vector_bfd[j]){
        btxv->crfd_vector_bfd[j]++;              
        break;
      }
  }
  */
  for(i=1;i<=np;i++){ 
    r = er_random(idum); 
    for(j=0;j<ns;j++)
      if(r < xv->crfd_vector_bfd[j]){
        btxv->crfd_vector_bfd[j]++;              /* increment frequencies */
        break;
      }
  } 
  /* ===== End of version 1.0 update ===== */

  /* Here btxv->crfd_vector_bfd[] is bootstrap frequencies in vector format; 
     next statements put it in matrix format and store it in btxv->bfd[][] */
  
  ptr = btxv->crfd_vector_bfd;
  for(i=0;i<nsx;i++)
	for(j=0;j<nsv;j++)
	  btxv->bfd[i][j] = *ptr++;

  /* get row and col marginal fd, rfd, crfd, and PR */
                      
  for(j=0;j<nsv;j++) btxv->fd_v[j] = 0.;               /* initialization */

  for(i=0;i<nsx;i++){   
	btxv->fd_x[i] = 0.;                                /* initialization */
	for(j=0;j<nsv;j++){
      btxv->fd_x[i] += btxv->bfd[i][j];                        /* row fd */
	  btxv->fd_v[j] += btxv->bfd[i][j];                        /* col fd */
	}
	btxv->density_x[i] = btxv->fd_x[i]/np;                /* row density */
  }
  for(j=0;j<nsv;j++) 
	btxv->density_v[j] = btxv->fd_v[j]/np;                /* col density */

  cum_rel_freqs(0, nsx-1, 1, btxv->density_x, btxv->crfd_x); /* row crfd */
  for (i=0;i<nsx;i++)
    btxv->prd_x[i] = perc_rank(0, nsx-1, 1, btxv->crfd_x, i); /* row prd */

  cum_rel_freqs(0, nsv-1, 1, btxv->density_v, btxv->crfd_v); /* col crfd */
  for (j=0;j<nsv;j++)
    btxv->prd_v[j] = perc_rank(0, nsv-1, 1, btxv->crfd_v, j); /* col prd */

  /* following code gets rel fd version of btxv->bfd[]; i.e., 
     btxv->brfd[][] is the bootstrap rel freq biv dist for x by v; 
	 needed for FEorMFE_EE() */ 

  for(i=0;i<nsx;i++)
    for(j=0;j<nsv;j++) 
	  btxv->brfd[i][j] = btxv->bfd[i][j]/np;

  return;
}

/*************************************************************************/