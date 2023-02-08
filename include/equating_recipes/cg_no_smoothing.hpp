/*
  CG_NoSmooth.c 

  Wrapper, print, and linear functions for common-item non-equivalent groups 
  design with no smoothing.

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

#ifndef CG_NO_SMOOTHING_HPP
#define CG_NO_SMOOTHING_HPP

#include <equating_recipes/structures/all_structures.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>

namespace EquatingRecipes {
  struct CINEGDeisgnEquatingNoSmoothing {
    /*
      Wrapper for getting mean, linear, or equipercentile equating for CINEG design
      with no smoothing. Equipercentile equating includes frequency estimation with 
      Braun-Holland (linear) results, modified frequency estimation with 
      Braun-Holland (linear) results, and chained equipercentile equating

      Purposes of Wrapper_CN() are to:
          (i) allocate space for r->eraw[][] and r->mts[0][] 
        (ii) populate elements of inall, including the 
              determination of proportional weights, if requested
        (iii) get equating results and store them in r
        (iv) get moments
        
      Assumes that in xv, score 1 is for x and score 2 is for v
      Assumes that in yv, score 1 is for y and score 2 is for v
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep
      
      Input:
      
        design = 'C' (CINEG)
        method:  'M' = mean, 
                'L' = linear (Tucker, Levine-obs, Levine-true, chained)
                'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
                'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
                'G' = FE + BH-FE + MFE + BH-MFE
                'C' = Chained
          'H' = FE + BH-FE + Chained
                'A' = FE + BH-FE + MFE + BH-MFE + Chained
                  
        smoothing = 'N' (none)  
        w1 = weight for pop. 1 (associated with xv)
            [0,1] except that for any number outside this 
            range, proportional weights are used -- i.e.,
            w1 = xv->n/(xv->n + yv->n)
        anchor = 0 --> external; otherwise internal
        rv1 = reliability of common items for population 1 
              (set to 0 for all methods except 'F', 'G, and 'A')
        rv2 = reliability of common items for population 2
              (set to 0 for all methods except 'F', 'G, and 'A')
        xv = struct BSTATS
        yv = struct BSTATS 
        rep = replication number for bootstrap; should be set to 0
              for actual equating; 
      
        NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
        
      Output:
        
        struct PDATA inall:   populates selected values of inall 
        
        struct ERAW_RESULTS r          

          msx[] = means for x for synthetic pop 
          msy[] = means for y for synthetic pop 
          dssx[] = sd's for x for synthetic pop 
          ssy[] = sd's for y for synthetic pop 
          gamma1[] = gamma's for pop 1 
          gamma2[] = gamma's for pop 2 
          a[] = slopes
          b[] = intercepts
          eraw[][]:  equated raw scores
          mts[][]:  moments for equated raw scores  
          
          NOTE: vectors in ERAW_RESULTS for nm = 4 linear methods 
                give results for following methods: 0 - Tucker,
                1 - Levine Observed, 2 - Levine True, 3 - Chained  
              
          NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
                of equated scores; memory allocated here;
                eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
                because we are getting equated raw scores for x.
                eraw[][] stored in this "row"  manner so that 
                Equated_ss() can be used directly on 
                eraw[k] where k is the method number  
              
      NOTE: Whenever method differs, there must be different structures
            passed as struct PDATA and struct ERAW_RESULTS 
        
      NOTE: If Wrapper_CN() is called in a bootstrap loop, then in
            the calling function struct ERAW_RESULTS must be different
            from struct ERAW_RESULTS for the actual equating. 
                                                
      Function calls other than C or NR utilities:
        CI_LinEq() 
        FEorMFE_EE()
        Chained_EE()
        runerror()
                                                  
      R. L. Brennan

      Date of last revision: 6/30/08   
    */
    void run(const EquatingRecipes::Structures::Design& design,
             const EquatingRecipes::Structures::Method& method,
             const EquatingRecipes::Structures::Smoothing& smoothing,
             const double& population1Weight,
             const bool& isInternalAnchor,
             const double& reliabilityCommonItemsPopulation1,
             const double& reliabilityCommonItemsPopulation2,
             const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXV,
             const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsYV,
             const size_t& bootstrapReplicationNumber,
             EquatingRecipes::Structures::PData& pData,
             EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {

    }

    /*
      computes results for common-item linear equating
      of x to scale of y for methods:
        0. Tucker
        1. Levine Observed
        2. Levine True
        3. Chained
                  
      Kolen and Brennan (2004) notation used here
      
      Input:
          
          mnx1 = mean for x for pop 1    
          sdx1 = sd for x for pop 1    
          mnv1 = mean for v for pop 1    
          sdv1 = sd for v for pop 1     
          covxv1 = cov for x and v in pop 1  

          mny2 = mean for y in pop 2     
          sdy2 = sd for y in pop 2     
          mnv2 = mean for v in pop 2     
          sdv2 = sd for v in pop 2    
          covyv2 = cov for y and v in pop 2
          
          w1 = weight for pop 1
          anchor = 0 --> external; otherwise internal
          method = 'M' --> mean equating (i.e., slope = 1)
          min = min score for x 
          max = max score for x
          inc = increment for x
          fdx[] = fd for x 
      
      Output:
      
          msx[] = vector of means for x for synthetic pop (for the methods)
          msy[] = vector of means for y for synthetic pop (for the methods)
          ssx[] = vector of sd's for x for synthetic pop (for the methods) 
          ssy[] = vector of sd's for y for synthetic pop (for the methods) 
          gamma1[] = vector of gamma's for pop 1 (for the methods)
          gamma2[] = vector of gamma's for pop 2 (for the methods)      
          a[] = vector of slopes (for the methods) 
          b[] = vector of intercepts (for the methods)
          eraw[][] =  method (rows) by raw score (columns) matrix
                      of equated scores; see Wrapper_CLN()for more details
          mts[][] = method (rows) by moments (columns) matrix
                    
          NOTE:  it is assumed that storage already allocated for
                all vectors and eraw[][]

        Function calls other than C or NR utilities:
          CI_LinObsEq()
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08                                   
      */
    void commonItemLinearEquating(const double& meanXPop1,
                                  const double& sdXPop1,
                                  const double& meanVPop1,
                                  const double& sdVPop1,
                                  const double& covYVPop2,
                                  const double& population1Weight,
                                  const bool& isInternalAnchor,
                                  const EquatingRecipes::Structures::Method& method,
                                  const double& mininumScoreX,
                                  const double& maximumScoreX,
                                  const double& scoreIncrementX, 
                                  const Eigen::VectorXd freqDistX,
                                  Eigen::VectorXd& meanVectorXSynPop,
                                  Eigen::VectorXd& meanVectorYSynPop,
                                  Eigen::VectorXd& sdVectorXSynPop,
                                  Eigen::VectorXd& sdVectorYSynPop,
                                  Eigen::VectorXd& gammaVectorPop1,
                                  Eigen::VectorXd& gammaVectorPop2,
                                  Eigen::VectorXd& slopes,
                                  Eigen::VectorXd& intercepts,
                                  Eigen::MatrixXd& methodByEquatedRawScores,
                                  Eigen::MatrixXd& methodByMomemnts) {}

    /* 
      For a linear observed score equating method and a CINEG design,
      get synthetic pop means and sd's, and get slope (*a) and intercept (*b)  

      Kolen and Brennan (2004) notation used here
    
      Input:
        
        mnx1 = mean for x for pop 1    
        sdx1 = sd for x for pop 1    
        mnv1 = mean for v for pop 1    
        sdv1 = sd for v for pop 1     

        mny2 = mean for y in pop 2     
        sdy2 = sd for y in pop 2     
        mnv2 = mean for v in pop 2     
        sdv2 = sd for v in pop 2    
        
        w1 = weight for pop 1
        method = 'M' --> mean equating (i.e., slope = 1)
        gamma1 = gamma for pop 1 
        gamma2 = gamma for pop 2
    
    Output:

        msx = mean for x for synthetic pop 
        msy = mean for y for synthetic pop 
        ssx = sd for x for synthetic pop  
        ssy = sd for y for synthetic pop       
        a = slope  
        b = intercept 
      
      Function calls other than C or NR utilities: None
                                                  
      R. L. Brennan

      Date of last revision: 6/30/08         
    */
    void commonItemLinearObservedEquating(const double& meanXPop1,
                                          const double& sdXPop1,
                                          const double& meanVPop1,
                                          const double& sdVPop1,
                                          const double& meanYPop2,
                                          const double& sdYPop2,
                                          const double& meanVPop2,
                                          const double& sdVPop2,
                                          const double& population1Weight,
                                          const EquatingRecipes::Structures::Method& method,
                                          const double& gammaPop1,
                                          const double& gammaPop2,
                                          double& meanXSynPop,
                                          double& meanYSynPop,
                                          double& sdXSynPop,
                                          double& sdYSynPop,
                                          double& slope,
                                          double& intercept) {}
  };
} // namespace EquatingRecipes

#endif