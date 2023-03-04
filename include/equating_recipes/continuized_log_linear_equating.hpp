/*

CLL_Equate.c  File for continuized log-linear equating

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

#ifndef CONTINUIZED_LOG_LINEAR_EQUATING_HPP
#define CONTINUIZED_LOG_LINEAR_EQUATING_HPP

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/log_linear_equating.hpp>


namespace EquatingRecipes {
  class ContinuizedLogLinearEquating {
  public:
    /*
    Wrapper for doing CLL equating with RG design
      and log-linear smoothing smoothing
      
    Assumes that equating puts raw scores for x on scale of y
    
    NOTE: This function is used (unaltered) for both actual equating and 
          equating done in Wrapper_Bootstrap().  Distinguishing between the
          two is the purpose of the variable rep

    Input
    
      design = 'R'(random groups)
      method = 'E'(equipercentile)
      smoothing = 'Z' (CLL 'smoothing')  
      *x = pointer to struct USTATS (new form)
      *y = pointer to struct USTATS (old form)
      *ullx = pointer to struct BB_SMOOTH (new form)
      *ully = pointer to struct BB_SMOOTH (old form)
      rep = replication number for bootstrap; should be set to 0
            for actual equating;  
      
    Output
      
      struct PDATA *inall:   populates selected values of inall 
      
      struct ERAW_RESULTS *r: populates            

        **eraw: equated raw scores;          
                method (rows) by raw score (columns) matrix
                of equated scores. Here there is only one method.
                So, memory allocated for eraw[][] is: 
                eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                              (loc(x->max,x->min,x>-inc)]
                because we are getting equated raw scores for x 
        **mts:  moments for equated raw scores           
        
    NOTE: If Wrapper_RC() is called in a bootstrap loop,
          then in the calling function struct ERAW_RESULTS must
          be different from struct ERAW_RESULTS for the actual
          equating. 

    Function calls other than C or NR utilities:                   
      CLLEquateEG()
      MomentsFromFD()  
                                                  
    Tianyou Wang and Robert L. Brennan

    Date of last revision: 6/30/08       
  */
    void Wrapper_RC(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const EquatingRecipes::Structures::UnivariateStatistics& x,
                    const EquatingRecipes::Structures::UnivariateStatistics& y,
                    const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& ullx,
                    const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& ully,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {
      /* method name --- 10 characters; right justified */
      std::vector<std::string> names {"   y-equiv"};
      // double maxx, maxy, *scoresx, *scoresy;
      // int i;
      // double *parax, *paray;

      std::vector<double> parax;
      std::vector<double> paray;

      parax.resize(ullx.degreesOfSmoothing + 1);
      paray.resize(ully.degreesOfSmoothing + 1);

      parax[0] = ullx.cllNormalizingConstant;
      for (size_t i = 0; i < ullx.degreesOfSmoothing; i++) {
        parax[i + 1] = ullx.betaCoefficients(i);
      }

      paray[0] = ully.cllNormalizingConstant;
      for (size_t i = 0; i < ully.degreesOfSmoothing; i++) {
        paray[i + 1] = ully.betaCoefficients(i);
      }

      std::vector<double> scoresx; // = dvector(0, ullx->ns);
      std::vector<double> scoresy; // = dvector(0, ully->ns);

      scoresx.resize(ullx.numberOfScores + 1);
      scoresy.resize(ully.numberOfScores + 1);

      double maxx = ullx.mininumRawScore + static_cast<double>(ullx.numberOfScores - 1) * ullx.rawScoreIncrement;
      double maxy = ully.mininumRawScore + static_cast<double>(ully.numberOfScores - 1) * ully.rawScoreIncrement;
      for (size_t i = 0; i < ullx.numberOfScores; i++) {
        scoresx[i] = static_cast<double>(i);
      }

      for (size_t i = 0; i < ully.numberOfScores; i++) {
        scoresy[i] = static_cast<double>(i);
      }

      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating */
                                                            /* counting of replications done in Wrapper_Bootstrap() */

      /* allocation and assignments for struct PDATA inall
        Note that for every assignment of the form inall->(var) = x->(var)
        or inall->(var) = y->(var), values vary depending on whether x or y 
        is for actual equating or a bootstrap sample; all other values are 
        the same for the actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) { /* no assignment or stor alloc for bootstrap reps */
        pData.summaryRawDataX = x;
        pData.summaryRawDataY = y;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;

        pData.methods.push_back(names[0]); /* only one row/method, 0 */

        pData.mininumScoreX = x.minimumScore;
        pData.maximumScoreX = x.maximumScore;
        pData.scoreIncrementX = x.scoreIncrement;
        pData.scoreFrequenciesX = x.freqDistDouble;
        pData.numberOfExaminees = x.numberOfExaminees;

        pData.univariateLogLinearSmoothingX = ullx;
        pData.univariateLogLinearSmoothingY = ullx;
      }

      /* allocation and assignments for results */

      if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.mininumScoreX,
                                                                                   pData.scoreIncrementX);

        results.equatedRawScores.resize(pData.methods.size(), maximumScoreLocation + 1);
        results.equatedRawScoreMoments.resize(pData.methods.size(), 4);
      }

      /* Compute equating results */
      Eigen::VectorXd equatedRawScores = cllEquateEG(ullx.mininumRawScore,
                                                     maxx,
                                                     parax,
                                                     ully.mininumRawScore,
                                                     maxy,
                                                     paray,
                                                     scoresx);

      results.equatedRawScores.row(0) = equatedRawScores;

      /* get moments */
      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(results.equatedRawScores.row(0),
                                                                                                             pData.scoreFrequenciesX);

      results.equatedRawScoreMoments.row(0) = moments.momentValues;
    }

    /*
      Wrapper for doing CLL equating with SG design
      and log-linear smoothing.  
      
      NOTE: This is for the SG design in which x and y do not share any items in 
      common, which means that functionally this is the external anchor case.  
      The bivariate log-linear smoothing procedure needs to know this. So, when
      Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
      anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
      (with anchor set to 0) and Wrapper_SC() can still be used, but convergence 
      of the smoothing algorithm may be compromised because of dependencies between
      x and y.
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep

      Input
      
        design = 'S' (single group)
        method = 'E' (equipercentile)
        smoothing = 'Z' (CLL 'smoothing')  
        xy = struct BSTATS 
      bllxy = struct BLL_SMOOTH
        rep = replication number for bootstrap; should be set to 0
              for actual equating;  
        
        NOTE: it is assumed that the first variable
              in xy is indeed x and the second variable is y.
          For bllxy, data memebers with '_x' are for x,
          and data members with '_v' are for y.  This somewhat
          inconsistent notation arises because BLL_SMOOTH is
          usually used with the CG design in which the variables
          are more naturally designated x (or u) and v.

      Output
        
        struct PDATA inall   populates selected values of inall  

        struct ERAW_RESULTS *r: populates            

          **eraw: equated raw scores;          
                  method (rows) by raw score (columns) matrix
                  of equated scores. Here there is only one method.
                  So, memory allocated for eraw[][] is: 
                  eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                                (loc(x->max,x->min,x>-inc)]
                  because we are getting equated raw scores for x 
          **mts:  moments for equated raw scores            
          
      NOTE: If Wrapper_SC() is called in a bootstrap loop,
            then in the calling function struct ERAW_RESULTS must
            be different from struct ERAW_RESULTS for the actual
            equating. 
                                                
      Function calls other than C or NR utilities:
        CLLEquateSG()
        MomentsFromFD()  
                                                    
      Tianyou Wang

      Date of last revision: 6/30/08   
    */
    void Wrapper_SC(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const EquatingRecipes::Structures::BivariateStatistics& xy,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllxy,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {
      /* method names --- 10 characters; right justified */
      std::vector<std::string> names {"   y-equiv"};

      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating */
                                                            /* counting of replications done in Wrapper_Bootstrap() */

      /* Allocation and assignments for struct PDATA inall>
      Note that for every assignment of the form inall->(var) = x->(var)
      or inall->(var) = y->(var), values vary depending on whether x or y 
      is for actual equating or a bootstrap sample; all other values are 
      the same for the actual equating and a bootstrap sample */
      if (pData.bootstrapReplicationNumber == 0) { /* no assignment or stor alloc for bootstrap reps */
        pData.summaryRawDataXY = xy;
        pData.bivariateLogLinearSmoothingXY = bllxy;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        pData.isInternalAnchor = false; /* implicitly, anchor is external for biv log-linear
						                        smoothing with the SG design */
        pData.methods.push_back(names[0]);
        pData.mininumScoreX = xy.univariateStatisticsRow.minimumScore;
        pData.maximumScoreX = xy.univariateStatisticsRow.maximumScore;
        pData.scoreIncrementX = xy.univariateStatisticsRow.scoreIncrement;
        pData.scoreFrequenciesX = xy.univariateStatisticsRow.freqDistDouble;
        pData.numberOfExaminees = xy.numberOfExaminees;
      }

      /* allocation and assignments for r */

      if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.mininumScoreX,
                                                                                   pData.scoreIncrementX);

        results.equatedRawScores.resize(pData.methods.size(), maximumScoreLocation + 1);
        results.equatedRawScoreMoments.resize(pData.methods.size(), 4);
      }

      /* Compute equating results. Put x on scale of y.
        Note that in struct xy, '1' designates x and '2' designates y; 
        in struct bllxy, '_x' designates x and '_v' designates y. So: 
          xy->ns2 = number of score categories for y
        xy->min2 = minimum score for y
        xy->inc2 = increment for y
        bllxy->crfd_v = log-linear smoothed cum rel fd for y
        xy->ns1 = number of score categories for x
          bllxy->prd_x = log-linear smoothed PR dist for x
        r->eraw[0] = y equivalents for x (output) */

      Eigen::VectorXd equatedRawScores = cllEquateSG(bllxy);

      results.equatedRawScores.row(0) = equatedRawScores;

      /* get moments */
      EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(results.equatedRawScores.row(0),
                                                                                                             pData.scoreFrequenciesX);

      results.equatedRawScoreMoments.row(0) = moments.momentValues;

      return;
    }

    /*
      Wrapper for doing CLL equating for SG with counter-balance design 
      and log-linear smoothing.  
      
      NOTE: This is for the SG with counter balance design in which x and y 
      in both group 1 and group 2 do not share any items in 
      common, which means that functionally this is the external anchor case.  
      The bivariate log-linear smoothing procedure needs to know this. So, when
      Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
      anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
      (with anchor set to 0) and Wrapper_SL() can still be used, but convergence 
      of the smoothing algorithm may be compromised because of dependencies between
      x and y. 
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep

      Input
      
        design = 'B' (single group)
        method = 'E' (equipercentile)
        smoothing = 'Z' (CLL 'smoothing')  
      wtsx = weight of x for population 1
      wtsy = weight of y for population 2
        xy1 = struct BSTATS for group 1
        xy2 = struct BSTATS for group 2
      bllxy1 = struct BLL_SMOOTH for group 1
      bllxy2 = struct BLL_SMOOTH for group 2
        rep = replication number for bootstrap; should be set to 0
              for actual equating;  
        
        NOTE: it is assumed that the first variable
              in xy is indeed x and the second variable is y.
          For bllxy, data memebers with '_x' are for x,
          and data members with '_v' are for y.  This somewhat
          inconsistent notation arises because BLL_SMOOTH is
          usually used with the CG design in which the variables
          are more naturally designated x (or u) and v.

      Output
        
        struct PDATA inall   populates selected values of inall  

        struct ERAW_RESULTS *r: populates            

          **eraw: equated raw scores;          
                  method (rows) by raw score (columns) matrix
                  of equated scores. Here there is only one method.
                  So, memory allocated for eraw[][] is: 
                  eraw[0][[0 ... (nscores(x->max,x->min,x>-inc)-1) =
                                (loc(x->max,x->min,x>-inc)]
                  because we are getting equated raw scores for x 
          **mts:  moments for equated raw scores            
          
      NOTE: If Wrapper_BC() is called in a bootstrap loop,
            then in the calling function struct ERAW_RESULTS must
            be different from struct ERAW_RESULTS for the actual
            equating. 
                                                
      Function calls other than C or NR utilities:
        CLLEquateCB()
        MomentsFromFD()  
                                                    
      Tianyou Wang

      Date of last revision: 6/30/08   
    */
    void Wrapper_BC(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const double& wtsx,
                    const double& wtsy,
                    const EquatingRecipes::Structures::BivariateStatistics& xy1,
                    const EquatingRecipes::Structures::BivariateStatistics& xy2,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllxy1,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllxy2,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {
      /* method names --- 10 characters; right justified */
      std::vector<std::string> names {"   y-equiv"};

      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating */
                                                            /* counting of replications done in Wrapper_Bootstrap() */

      /* Allocation and assignments for struct PDATA inall>
        Note that for every assignment of the form inall->(var) = x->(var)
        or inall->(var) = y->(var), values vary depending on whether x or y 
        is for actual equating or a bootstrap sample; all other values are 
        the same for the actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) { /* no assignment or stor alloc for bootstrap reps */

        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        pData.isInternalAnchor = false; /* implicitly, anchor is external for biv log-linear
						                        smoothing with the CB design */

        pData.methods.push_back(names[0]); /* only one row/method, 0 */
      }

      /* allocation and assignments for r */

      if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                   pData.mininumScoreX,
                                                                                   pData.scoreIncrementX);
        results.equatedRawScores.resize(pData.methods.size(), maximumScoreLocation + 1);
        results.equatedRawScoreMoments.resize(pData.methods.size(), 4);
      }

      /* Compute equating results. Put x on scale of y.
        Note that in struct xy, '1' designates x and '2' designates y; 
        in struct bllxy, '_x' designates x and '_v' designates y. So: 
          xy->ns2 = number of score categories for y
        xy->min2 = minimum score for y
        xy->inc2 = increment for y
        bllxy->crfd_v = log-linear smoothed cum rel fd for y
        xy->ns1 = number of score categories for x
          bllxy->prd_x = log-linear smoothed PR dist for x
        r->eraw[0] = y equivalents for x (output) */

      Eigen::VectorXd equatedRawScores = cllEquateCB(bllxy1, bllxy2, wtsx, wtsy);

      results.equatedRawScores.row(0) = equatedRawScores;

      EquatingRecipes::Structures::Moments moments =
          EquatingRecipes::Utilities::momentsFromScoreFrequencies(results.equatedRawScores.row(0),
                                                                  pData.scoreFrequenciesX);

      results.equatedRawScoreMoments.row(0) = moments.momentValues;
    }

    /*
      Wrapper for CLL equating for CG design with log-linear smoothing. 
      Equipercentile equating includes frequency estimation with 
      Braun-Holland (linear) results, modified frequency estimation with 
      Braun-Holland (linear) results, and chained equipercentile equating
        
      Assumes that in xv, score 1 is for x and score 2 is for v
      Assumes that in yv, score 1 is for y and score 2 is for v
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep
      
      Input:
      
        design = 'C' (CINEG)

        method:  'E' = Frequency estimation (FE)
                'F' = Modified freq est (MFE) 
                'G' = FE +  MFE
                'C' = Chained
          'H' = FE +  Chained
                'A' = FE + MFE + Chained
                  
        smoothing = 'Z' (CLL 'smoothing') 

        w1 = weight for pop. 1 (associated with xv)
            [0,1] except that for any number outside this 
            range, proportional weights are used -- i.e.,
            w1 = xv->n/(xv->n + yv->n)
        anchor = 0 --> external; 1 --> internal
        rv1 = reliability of common items for population 1 
              (set to 0 for all methods except 'F', 'G, and 'A')
        rv2 = reliability of common items for population 2
              (set to 0 for all methods except 'F', 'G, and 'A')
        xv = struct BSTATS
        yv = struct BSTATS 
        bllxv = struct BLL_SMOOTH; uses brfd, prd_x,  and crfd_v 
      bllyv = struct BLL_SMOOTH; uses brfd, crfd_x, and crfd_v
              (Note that bllyv->crfd_x is really crfd for y in pop 2)
        rep = replication number for bootstrap; should be set to 0
              for actual equating; 
      
        NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be conducted 
        
      Output:
        
        struct PDATA inall:   populates selected values of inall 
        
        struct ERAW_RESULTS r          
    
          a[] = slopes for Braun-Holland
          b[] = intercepts for Braun-Holland
          eraw[][]:  equated raw scores
          mts[][]:  moments for equated raw scores   
              
          NOTE: eraw[][] is a method (rows) by raw score (columns) matrix
                of equated scores; memory allocated here;
                eraw[0...(nm-1)][[0...(loc(xv->max1,xv->min1,xv>-inc1)]                                      ]
                because we are getting equated raw scores for x.
                eraw[][] stored in this "row"  manner so that 
                Equated_ss() can be used directly on 
                eraw[k] where k is the method number  
              
      NOTE: Whenever method differs, there must be different structures
            passed as struct PDATA and struct ERAW_RESULTS 
        
      NOTE: If Wrapper_CC() is called in a bootstrap loop, then in
            the calling function struct ERAW_RESULTS must be different
            from struct ERAW_RESULTS for the actual equating. 
                                                
      Function calls other than C or NR utilities:
        CLLEquateNEATPS()
        CLLEquateNEATChn()
        runerror()
                                                  
      T. D. Wang

      Date of last revision: 6/30/08   
    */
    void Wrapper_CC(const EquatingRecipes::Structures::Design& design,
                    const EquatingRecipes::Structures::Method& method,
                    const EquatingRecipes::Structures::Smoothing& smoothing,
                    const double& w1,
                    const bool& isInternalAnchor,
                    const double& reliabilityCommonItemsPopulation1,
                    const double& reliabilityCommonItemsPopulation2,
                    const EquatingRecipes::Structures::BivariateStatistics& xv,
                    const EquatingRecipes::Structures::BivariateStatistics& yv,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllxv,
                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bllyv,
                    const size_t& replicationNumber,
                    EquatingRecipes::Structures::PData& pData,
                    EquatingRecipes::Structures::EquatedRawScoreResults& results) {
      std::vector<std::string> names {"        FE", "       MFE", "  ChainedE"};
      
      std::string methodCode = EquatingRecipes::Utilities::getMethodCode(method);

      pData.bootstrapReplicationNumber = replicationNumber; /* should be set to 0 for actual equating. */
                                                            /* Counting of replications done in Wrapper_Bootstrap(), 
                                                               which is why this statement cannot be in the if statement below */

      /* allocation and assignments for inall
        Note that for every assignment of the form inall->(var) = xv->(var)
        or inall->(var) = yv->(var) values vary depending on whether xv or yv
        is for actual equating or a bootstrap sample; all other values are 
        the same for the actual equating and a bootstrap sample */

      if (pData.bootstrapReplicationNumber == 0) { /* no assignment or stor alloc for bootstrap reps */
        pData.summaryRawDataXV = xv;
        pData.summaryRawDataYV = yv;
        pData.bivariateLogLinearSmoothingXV = bllxv;
        pData.bivariateLogLinearSmoothingYV = bllyv;
        pData.design = design;
        pData.method = method;
        pData.smoothing = smoothing;
        if (w1 < 0 || w1 > 1) {
          /* proportional wts if w1 outside [0,1] */
          pData.weightSyntheticPopulation1 = static_cast<double>(xv.numberOfExaminees) /
                                             static_cast<double>(xv.numberOfExaminees + yv.numberOfExaminees);
        } else {
          pData.weightSyntheticPopulation1 = w1;
        }

        pData.isInternalAnchor = isInternalAnchor;
        pData.reliabilityCommonItemsPopulation1 = reliabilityCommonItemsPopulation1;
        pData.reliabilityCommonItemsPopulation2 = reliabilityCommonItemsPopulation2;

        if ((methodCode == "F" || methodCode == "G" || methodCode == "A") &&
            (reliabilityCommonItemsPopulation1 == 0 || reliabilityCommonItemsPopulation2 == 0)) {
          throw std::runtime_error("\nMFE cannot be conducted since rv1 == 0 or rv2 == 0");
        }

        pData.methods.clear();

        if (methodCode == "E") { /* method == "E" */
          pData.methods.push_back(names[0]);

        } else if (methodCode == "F") { /* method == "F" */
          pData.methods.push_back(names[1]);

        } else if (methodCode == "G") { /* method == "G" */
          pData.methods.push_back(names[0]);
          pData.methods.push_back(names[1]);

        } else if (methodCode == "C") { /* method == "C" */
          pData.methods.push_back(names[2]);

        } else if (methodCode == "H") { /* method == "H" */
          pData.methods.push_back(names[0]);
          pData.methods.push_back(names[2]);

        } else { /* method == "A" */
          pData.methods = names;
        }

        pData.mininumScoreX = xv.univariateStatisticsRow.minimumScore;
        pData.maximumScoreX = xv.univariateStatisticsRow.maximumScore;
        pData.scoreIncrementX = xv.univariateStatisticsRow.scoreIncrement;
        pData.scoreFrequenciesX = xv.univariateStatisticsRow.freqDistDouble;
        pData.numberOfExaminees = xv.numberOfExaminees;
      }

      /* allocation and assignments for r */

      if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(xv.univariateStatisticsRow.maximumScore,
                                                                                   xv.univariateStatisticsRow.minimumScore,
                                                                                   xv.univariateStatisticsRow.scoreIncrement);
        results.equatedRawScores.resize(pData.methods.size(), maximumScoreLocation + 1);
        results.equatedRawScoreMoments.resize(pData.methods.size(), maximumScoreLocation + 1);
        results.relativeFreqDistsX.resize(1, maximumScoreLocation + 1);

        maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(yv.univariateStatisticsRow.maximumScore,
                                                                            yv.univariateStatisticsRow.minimumScore,
                                                                            yv.univariateStatisticsRow.scoreIncrement);

        results.relativeFreqDistsY.resize(1, maximumScoreLocation + 1);
      }

      /* Equipercentile results, including Braun-Holland (BH) linear. 
         Note: For FE syn densities are in fxs[0] and gys[0]
         For MFE syn densities are in fxs[1] and gys[1] */

      /* FE + BH-FE in positions 0 and 1*/

      if (methodCode == "E" || methodCode == "G" || methodCode == "A" || methodCode == "H") {
        Eigen::VectorXd equatedRawScores = cllEquateNEATPS(bllxv, bllyv, pData.weightSyntheticPopulation1);
        results.equatedRawScores.row(0) = equatedRawScores;
      }

      if (methodCode == "C" || methodCode == "A" || methodCode == "H") {
        Eigen::VectorXd equatedRawScores = cllEquateNEATChn(bllxv, bllyv);

        size_t methodIndex;
        if (methodCode == "C") {
          methodIndex = 0;
        } else if (methodCode == "A") {
          methodIndex = 2;
        } else {
          methodIndex = 1;
        }
        
        results.equatedRawScores.row(methodIndex) = equatedRawScores;
      }        

      /* get moments */

      for (size_t i = 0; i < pData.methods.size(); i++) {
        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(results.equatedRawScores.row(i),
          xv.univariateStatisticsRow.freqDistDouble);

        results.equatedRawScoreMoments.row(i) = moments.momentValues;
      }
    }

  private:
    /*--------------------------------------------------------------------------
      CLLEGPdf
      
      functionality

      Computes the continuous pdf based on the fitted loglinear model
      parameters for a discrete distribution
      
      Author: Tianyou Wang 10/29/2004.
      
      input
        min         lower limit of the distribution
        max         upper limit of the distribution
            npara       number of parameters
            para		a vector of parameters for the loglinear model
                        with a design matrix from polynominals of natural
                        basis.
        x           a particular score for which the pdf and cdf is 
                    generated   
        nc          normalizing constant 

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllEGPdf(const double& min,
                    const double& max,
                    const std::vector<double>& para,
                    const double& x,
                    const double& nc) {
      double pdf;

      pdf = expPolynomial(para, x);

      pdf /= nc;

      return pdf;
    }

    /*--------------------------------------------------------------------------
      CLLEGCdf
      
      functionality

      Computes the continuous pdf and cdf based on the fitted loglinear model
      parameters for a discrete distribution
      
      Author: Tianyou Wang 10/29/2004.
      
      input
        min         lower limit of the distribution
        max         upper limit of the distribution
            npara       number of parameters
            para		a vector of parameters for the loglinear model
                        with a design matrix from polynominals of natural
                        basis.
        x           a particular score for which the pdf and cdf is 
                    generated            
        nc          normalizing constant 

    
      output
        the function returns the smoothed cdf
    --------------------------------------------------------------------------*/
    double cllEGCdf(const double& min,
                    const double& max,
                    const std::vector<double>& para,
                    const double& x,
                    const double& nc) {
      double cdf;

      cdf = gaussianQuadrature64(min, x, para);

      cdf /= nc;

      return cdf;
    }

    double gaussianQuadrature(const double& a,
                              const double& b,
                              const std::vector<double>& para,
                              const std::vector<double>& x,
                              const std::vector<double>& w) {
      double xm = 0.5 * (b + a);
      double xr = 0.5 * (b - a);
      double s = 0.0;

      for (size_t j = 0; j < x.size(); j++) {
        double dx = xr * x[j];
        s += w[j] * (expPolynomial(para, (xm + dx)) + expPolynomial(para, (xm - dx)));
      }

      s *= xr;

      return s;
    }

    /*--------------------------------------------------------------------------
    GaussianQuadrature16
    
    functionality

    Computes the numerical integration using Gaussian quadrature with 
    16 points
    
    Author: Tianyou Wang 10/29/2004.
    
    input
          (*func)()   pointer to a function which is the integrand
          a   		lower limit of integral
          b   		upper limit of integral
      npara       number of parameters for the integrand
      para        vector of parameters for the integrand
  
    output
          The function returns the integrated value.
    --------------------------------------------------------------------------*/
    double gaussianQuadrature16(const double& a,
                                const double& b,
                                const std::vector<double>& para) {
      std::vector<double> x {0.09501250983764, 0.28160355077926, 0.45801677765723,
                             0.61787624440264, 0.755404408355, 0.86563120238783,
                             0.94457502307323, 0.98940093499165};
      std::vector<double> w {0.18945061045507, 0.18260341504492, 0.169156519395,
                             0.14959598881658, 0.12462897125553, 0.09515851168249,
                             0.06225352393865, 0.02715245941175};

      double s = gaussianQuadrature(a,
                                    b,
                                    para,
                                    x,
                                    w);

      return s;
    }

    /*--------------------------------------------------------------------------
    GaussianQuadrature32
    
    functionality

    Computes the numerical integration using Gaussian quadrature with 
    32 points
    
    Author: Tianyou Wang 10/29/2004.
    
    input
          (*func)()   pointer to a function which is the integrand
          a   		lower limit of integral
          b   		upper limit of integral
      npara       number of parameters for the integrand
      para        vector of parameters for the integrand
  
    output
          The function returns the integrated value.
    --------------------------------------------------------------------------*/
    double gaussianQuadrature32(const double& a,
                                const double& b,
                                const std::vector<double>& para) {
      std::vector<double> x {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double s = gaussianQuadrature(a,
                                    b,
                                    para,
                                    x,
                                    w);

      return s;
    }

    /*--------------------------------------------------------------------------
    GaussianQuadrature64
    
    functionality

    Computes the numerical integration using Gaussian quadrature with 
    64 points
    
    Author: Tianyou Wang 10/29/2004.
    
    input
          (*func)()   pointer to a function which is the integrand
          a   		lower limit of integral
          b   		upper limit of integral
      npara       number of parameters for the integrand
      para        vector of parameters for the integrand
  
    output
          The function returns the integrated value.
    --------------------------------------------------------------------------*/
    double gaussianQuadrature64(const double& a,
                                const double& b,
                                const std::vector<double>& para) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};
      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double s = gaussianQuadrature(a,
                                    b,
                                    para,
                                    x,
                                    w);

      return s;
    }

    double gaussianQuadrature64i(const double& a,
                                 const double& b,
                                 const std::vector<double>& para,
                                 const size_t& i) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};
      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double xm = 0.5 * (b + a);
      double xr = 0.5 * (b - a);
      double s = 0;

      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s += w[j] * (expPolynomialxi(para, (xm + dx), i) +
                     expPolynomialxi(para, (xm - dx), i));
      }

      s *= xr;

      return s;
    }

    /*--------------------------------------------------------------------------
    ExpPolynomial
    
    functionality

    Computes the exponential function of a polynomial (the fitted mean of 
    the loglinear model)
    
    Author: Tianyou Wang 10/29/2004.
    
    input
          npara       the number of polynomial coefficients
          *para  		the vector that contains the parameters
          x   		the plugged in value
  
    output
          The function returns exponential function of a polynomial.
    --------------------------------------------------------------------------*/
    double expPolynomial(const std::vector<double>& para,
                         const double& x) {
      double s = 0;

      for (size_t j = 0; j < para.size(); j++) {
        s += para[j] * std::pow(x, static_cast<double>(j));
      }

      s = std::exp(s);

      return s;
    }

    /*--------------------------------------------------------------------------
    ExpPolynomialxi
    
    functionality

    Computes the exponential function of a polynomial (the fitted mean of 
    the loglinear model)
    
    Author: Tianyou Wang 10/29/2004.
    
    input
          npara       the number of polynomial coefficients
          *para  		the vector that contains the parameters
          x   		the plugged in value
  
    output
          The function returns exponential function of a polynomial.
    --------------------------------------------------------------------------*/
    double expPolynomialxi(const std::vector<double>& para,
                           const double& x,
                           const size_t& i) {
      double s = 0;

      for (size_t j = 0; j < para.size(); j++)
        s += para[j] * std::pow(x, j);

      s = exp(s);
      s *= pow(x, static_cast<double>(i));

      return s;
    }

    /*------------------------------------------------------------------------------
      CalcLLContinuMoments	
      
      functionality

      calculates mean, sd, skewness and kurtosis for a continuous distribution.

      input -
            (*pdf)()    pointer to a function which is the pdf of the continuous 
                    distribution
            a   		lower limit of distribution
            b   		upper limit of distribution
        npara       number of parameters for the distribution
        para        vector of parameters for the distribution

      output -
        moments - mean, sd, skewness and kurtosis of distribution

    ------------------------------------------------------------------------------*/
    void calcLLContinuMoments(std::function<double(const double&, const double&, const std::vector<double>&, const double&)> pdf,
                              const double& a,
                              const double& b,
                              const std::vector<double>& para,
                              std::vector<double>& moments) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};
      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double xm = 0.5 * (b + a);
      double xr = 0.5 * (b - a);
      double s1 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s1 += w[j] * (pdf(a, b, para, (xm + dx)) * (xm + dx) +
                      pdf(a, b, para, (xm - dx)) * (xm - dx));
      }
      s1 *= xr;

      double s2 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s2 += w[j] * (pdf(a, b, para, (xm + dx)) * std::pow(xm + dx - s1, 2) +
                      pdf(a, b, para, (xm - dx)) * std::pow(xm - dx - s1, 2));
      }
      s2 *= xr;

      double s3 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s3 += w[j] * (pdf(a, b, para, (xm + dx)) * std::pow(xm + dx - s1, 3) +
                      pdf(a, b, para, (xm - dx)) * std::pow(xm - dx - s1, 3));
      }
      s3 *= xr / std::pow(s2, 1.5);

      double s4 = 0;
      for (size_t j = 0; j < 32; j++) {
        double dx = xr * x[j];
        s4 += w[j] * (pdf(a, b, para, (xm + dx)) * std::pow(xm + dx - s1, 4) +
                      pdf(a, b, para, (xm - dx)) * std::pow(xm - dx - s1, 4));
      }
      s4 *= xr / std::pow(s2, 2);

      moments = {s1,
                 std::sqrt(s2),
                 s3,
                 s4};
    }

    /*--------------------------------------------------------------------------
      CLLEquateEG
      
      functionality:

      Computes equating function based on continuized Log-linear cdf in Wang (2005). 
          
      author: Tianyou Wang 1/5/2005.
      
      input:
        minx        lower limit of the distribution for the new form
        maxx        upper limit of the distribution for the new form
        nparax      number of parameters for the new form
        paraxx		  a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the new form.
        miny        lower limit of the distribution for the old form
        maxy        upper limit of the distribution for the old form
        nparay      number of parameters for the old form
        paraxy		  a vector of parameters for the loglinear model
                    with a design matrix from polynominals of natural
                    basis for the old form.
        nDistCatx   Number of discrete score categories for the new form
        scoresx     vector containing the discrete scores for the new form

      output:
        Equatedx   a vector containing the equated score 
    --------------------------------------------------------------------------*/
    Eigen::VectorXd cllEquateEG(const double& minx,
                                const double& maxx,
                                const std::vector<double>& parax,
                                const double& miny,
                                const double& maxy,
                                const std::vector<double>& paray,
                                const std::vector<double>& scoresx) {
      double ncx = gaussianQuadrature64(minx, maxx, parax);
      double ncy = gaussianQuadrature64(miny, maxy, paray);
      Eigen::VectorXd equatedx;

      equatedx.resize(scoresx.size());

      std::for_each(scoresx.begin(),
                    scoresx.end(),
                    [&](const double& x) {
                      double cdfx = cllEGCdf(minx, maxx, parax, x, ncx);
                      double eqautedScore = cllInverseCdf(miny, maxy, paray, cdfx, ncy);

                      equatedx(eqautedScore);
                    });

      return equatedx;
    }

    /*--------------------------------------------------------------------------
      CLLInverseCdf
      
      functionality:

      Computes the inverse of the cdf for the continuized log-linear cdf in Wang 
      (2005). 
      
      author: Tianyou Wang 1/5/2005.
      
      input:
        min         lower limit of the distribution
        max         upper limit of the distribution
            npara       number of parameters
            para		a vector of parameters for the loglinear model
                        with a design matrix from polynominals of natural basis.
        cdf         a particular cdf for which the score is found
    
      output:
        The function returns the inverse of cdf
    --------------------------------------------------------------------------*/
    double cllInverseCdf(const double& min,
                         const double& max,
                         const std::vector<double>& para,
                         const double& cdf,
                         const double& nc) {
      double eps = .000001;

      double lb = min;
      double ub = max;
      double cdfl = cllEGCdf(min, max, para, lb, nc);
      double cdfu = cllEGCdf(min, max, para, ub, nc);
      double half;

      if (cdf < cdfl) {
        half = min;
      } else if (cdf > cdfu) {
        half = max;
      } else {
        for (size_t iter = 1; iter <= 200; iter++) {
          half = 0.5 * (lb + ub);
          double cdfhalf = cllEGCdf(min, max, para, half, nc);
          double absdif = std::abs(cdf - cdfhalf);

          if (absdif < eps) {
            break;
          } else if (cdfhalf < cdf) {
            lb = half;
          } else {
            ub = half;
          }
        }
      }

      return half;
    }

    /*--------------------------------------------------------------------------
      CLLBivPdf
      
      functionality

      Computes the continuous bivariate pdf based on the fitted bivariate 
      loglinear model	parameters for a discrete distribution
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
        x           a particular score for which the pdf and cdf is 
                    generated            
        y           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllBivPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                     const double& x,
                     const double& y) {
      // static double integ;
      // double minx, maxx, miny, maxy, pdf;

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;

      double integ = bivGaussianQuadrature64(bivar, minx, maxx, miny, maxy);
      double pdf = bivExpPolynomial(bivar, x, y) / integ;

      return pdf;
    }

    /*--------------------------------------------------------------------------
      CLLBivCdf
      
      functionality

      Computes the continuous bivariate cdf based on the fitted bivariate 
      loglinear model	parameters for a discrete distribution
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
        x           a particular score for which the pdf and cdf is 
                    generated            
        y           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllBivCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                     const double& x,
                     const double& y,
                     const double& nc) {
      // double minx, maxx, miny, maxy, cdf;

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;

      double cdf = bivGaussianQuadrature32(bivar, minx, x, miny, y);

      cdf /= nc; /* bivar->num_persons; */

      return cdf;
    }

    /*--------------------------------------------------------------------------
      CLLMargYPdf
      
      functionality

      Computes the continuous marginal pdf for Y based on the fitted bivariate 
      loglinear model	parameters for a discrete distribution
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
        y           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllMargYPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing bivar,
                       const double& y,
                       const double& nc) {
      // int j;
      // double minx, maxx, miny, maxy;
      // double xr, xm, dx, s;
      std::vector<double> x {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;
      double xm = 0.5 * (maxx + minx);
      double xr = 0.5 * (maxx - minx);
      double s = 0;
      for (size_t j = 0; j < 16; j++) {
        double dx = xr * x[j];
        s += w[j] * (bivExpPolynomial(bivar, xm + dx, y) + bivExpPolynomial(bivar, xm - dx, y));
      }

      s *= xr;

      s /= nc;

      return s;
    }

    /*--------------------------------------------------------------------------
      CLLMargXPdf
      
      functionality

      Computes the continuous marginal pdf for X based on the fitted bivariate 
      loglinear model	parameters for a discrete distribution
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
        x           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllMargXPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                       const double& x,
                       const double& nc) {
      // int j;
      // double minx, maxx, miny, maxy;
      // double yr, ym, dy, s;
      std::vector<double> y {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double minx = bivar.minimumRawScoreX - 0.5;
      double maxx = bivar.minimumRawScoreX + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double miny = bivar.minimumRawScoreV - 0.5;
      double maxy = bivar.minimumRawScoreV + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;
      double ym = 0.5 * (maxy + miny);
      double yr = 0.5 * (maxy - miny);
      double s = 0;
      for (size_t j = 0; j < 16; j++) {
        double dy = yr * y[j];
        s += w[j] * (bivExpPolynomial(bivar, x, ym + dy) + bivExpPolynomial(bivar, x, ym - dy));
      }

      s *= yr;

      s /= nc;

      return s;
    }

    /*--------------------------------------------------------------------------
      BivExpPolynomial
      
      functionality

      Computes the bivariate exponential function of a polynomial (the fitted mean of 
      the loglinear model)
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
            x   		the plugged in value for X
            y   		the plugged in value for Y
    
      output
            The function returns exponential function of a polynomial.
      --------------------------------------------------------------------------*/
    double bivExpPolynomial(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                            const double& x,
                            const double& y) {
      // int j, nx, ny, nxbyy;
      // double s = 0;

      double xValue = x;
      double yValue = y;

      /* when internal anchor, f(x,v) = f(u,v) for u = x - v and u is within valid range */
      if (bivar.isInternalAnchor) {
        xValue -= yValue;

        if (xValue < bivar.minimumRawScoreU) {
          return 0;
        } else if (xValue > (bivar.minimumRawScoreU + static_cast<double>(bivar.numberOfScoresU - 1) * bivar.scoreIncrementU)) {
          return 0;
        }
      }

      size_t nx = bivar.numberOfDegreesOfSmoothingU;
      size_t ny = bivar.numberOfDegreesOfSmoothingV;
      size_t nxbyy = bivar.numberOfCrossProductMoments;

      double s = bivar.cllNormalizingConstant;

      for (size_t j = 1; j <= nx; j++) {
        s += bivar.betaCoefficients(j - 1) * std::pow(xValue, static_cast<double>(j));
      }

      for (size_t j = 1; j <= ny; j++) {
        s += bivar.betaCoefficients(nx + j - 1) * std::pow(yValue, static_cast<double>(j));
      }

      for (size_t j = 1; j <= nxbyy; j++) {
        s += bivar.betaCoefficients(nx + ny + j - 1) *
             std::pow(x, bivar.crossProductMoments(j - 1, 0)) *
             std::pow(y, bivar.crossProductMoments(j - 1, 1));
      }

      s = std::exp(s);

      return s;
    }

    /*--------------------------------------------------------------------------
      BivGaussianQuadrature64
      
      functionality

      Computes the tow-dimensional numerical integration using Gaussian quadrature with 
      64 points
      
      Author: Tianyou Wang 4/29/2007.
      
      input
            (*func)()   pointer to a function which is the integrand
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
            ax   		lower limit of x variable
            bx   		upper limit of x variable
            ay   		lower limit of y variable
            by   		upper limit of y variable
    
      output
            The function returns the integrated value.
      --------------------------------------------------------------------------*/
    double bivGaussianQuadrature64(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                                   const double& ax,
                                   const double& bx,
                                   const double& ay,
                                   const double& by) {
      std::vector<double> x {0.02435029266342, 0.07299312178780, 0.12146281929612,
                             0.16964442042399, 0.21742364374001, 0.26468716220877, 0.31132287199021,
                             0.35722015833767, 0.40227015796399, 0.44636601725346, 0.48940314570705,
                             0.53127946401989, 0.57189564620263, 0.61115535517239, 0.64896547125466,
                             0.68523631305423, 0.71988185017161, 0.75281990726053, 0.78397235894334,
                             0.81326531512280, 0.84062929625258, 0.86599939815409, 0.88931544599511,
                             0.91052213707850, 0.92956917213194, 0.94641137485840, 0.96100879965205,
                             0.97332682778991, 0.98333625388463, 0.99101337147674, 0.99634011677196,
                             0.99930504173577};

      std::vector<double> w {0.04869095700914, 0.04857546744150, 0.04834476223480,
                             0.04799938859646, 0.04754016571483, 0.04696818281621, 0.04628479658131,
                             0.04549162792742, 0.04459055816376, 0.04358372452932, 0.04247351512365,
                             0.04126256324262, 0.03995374113272, 0.03855015317862, 0.03705512854024,
                             0.03547221325688, 0.03380516183714, 0.03205792835485, 0.03023465707240,
                             0.02833967261426, 0.02637746971505, 0.02435270256871, 0.02227017380838,
                             0.02013482315353, 0.01795171577570, 0.01572603047602, 0.01346304789672,
                             0.01116813946013, 0.00884675982636, 0.00650445796898, 0.00414703326056,
                             0.00178328072170};

      double xm = 0.5 * (bx + ax);
      double xr = 0.5 * (bx - ax);
      double ym = 0.5 * (by + ay);
      double yr = 0.5 * (by - ay);
      double si = 0;

      for (size_t i = 0; i < 32; i++) {
        double dx = xr * x[i];
        double sj1 = 0;
        double sj2 = 0;
        for (size_t j = 0; j < 32; j++) {
          double dy = yr * x[j];
          sj1 += yr * w[j] * (bivExpPolynomial(bivar, (xm + dx), (ym + dy)) + bivExpPolynomial(bivar, (xm + dx), (ym - dy)));
        // }
        // for (size_t j = 0; j < 32; j++) {
        //   double dy = yr * x[j];
          sj2 += yr * w[j] * (bivExpPolynomial(bivar, (xm - dx), (ym + dy)) + bivExpPolynomial(bivar, (xm - dx), (ym - dy)));
        }
        si += w[i] * (sj1 + sj2);
      }

      si *= xr;

      return si;
    }

    /*--------------------------------------------------------------------------
      BivGaussianQuadrature32
      
      functionality

      Computes the tow-dimensional numerical integration using Gaussian quadrature with 
      64 points
      
      Author: Tianyou Wang 4/29/2007.
      
      input
            (*func)()   pointer to a function which is the integrand
        bivar       the structure that contain bivariate distribution parameters
                    defined in "MLBivLogLin.h"
            ax   		lower limit of x variable
            bx   		upper limit of x variable
            ay   		lower limit of y variable
            by   		upper limit of y variable
    
      output
            The function returns the integrated value.
    --------------------------------------------------------------------------*/
    double bivGaussianQuadrature32(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar,
                                   const double& ax,
                                   const double& bx,
                                   const double& ay,
                                   const double& by) {
      // int i, j;
      // double xr, xm, dx, yr, ym, dy, si, sj1, sj2;
      std::vector<double> x {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};

      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double xm = 0.5 * (bx + ax);
      double xr = 0.5 * (bx - ax);
      double ym = 0.5 * (by + ay);
      double yr = 0.5 * (by - ay);
      double si = 0;

      for (size_t i = 0; i < 16; i++) {
        double dx = xr * x[i];
        double sj1 = 0;
        double sj2 = 0;
        for (size_t j = 0; j < 16; j++) {
          double dy = yr * x[j];
          sj1 += yr * w[j] * (bivExpPolynomial(bivar, (xm + dx), (ym + dy)) + bivExpPolynomial(bivar, (xm + dx), (ym - dy)));
          // }
          // for (size_t j = 0; j < 16; j++) {
          // double dy = yr * x[j];
          sj2 += yr * w[j] * (bivExpPolynomial(bivar, (xm - dx), (ym + dy)) + bivExpPolynomial(bivar, (xm - dx), (ym - dy)));
        }
        si += w[i] * (sj1 + sj2);
      }

      si *= xr;

      return si;
    }

    /*--------------------------------------------------------------------------
      CLLNEATPSMargPdf
      
      functionality

      Computes the continuous marginal pdf for X for the synthetic population
      using the basic assumption of the frequency estimation based on the fitted 
      bivariate loglinear model parameters for a discrete distribution
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar1      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 1
        bivar2      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 2
        wts         weights for population 1 in the synthetic population
        x           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllNEATPSMargPdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                            const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2,
                            const std::vector<double>& fv,
                            const double& wts, const double& x) {
      // int j;
      // static double nc;
      // double minv, maxv, minx, maxx;
      // double vr, vm, dv, s, fxv, fx1, fx2;
      std::vector<double> v {0.04830766568774, 0.14447196158280, 0.23928736225214,
                             0.33186860228213, 0.42135127613064, 0.50689990893223,
                             0.58771575724076, 0.66304426693022, 0.73218211874029,
                             0.79448379596794, 0.84936761373257, 0.89632115576605,
                             0.93490607593774, 0.96476225558751, 0.98561151154527,
                             0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double minv = bivar1.minimumRawScoreV - 0.5;
      double maxv = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      double minx = bivar1.minimumRawScoreX - 0.5;
      double maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double vm = 0.5 * (maxv + minv);
      double vr = 0.5 * (maxv - minv);
      double nc = bivGaussianQuadrature64(bivar1, minx, maxx, minv, maxv);
      double s = 0;
      for (size_t j = 0; j < 16; j++) {
        double dv = vr * v[j];

        double fxv = cllBivPdf(bivar1, x, vm + dv) / fv[j] * fv[j + 32] +
                     cllBivPdf(bivar1, x, vm - dv) / fv[j + 16] * fv[j + 48];

        s += w[j] * fxv;
      }

      double fx2 = s * vr;

      double fx1 = cllMargXPdf(bivar1, x, nc);

      double pdf = fx1 * wts + fx2 * (1 - wts);

      return pdf;
    }

    /*--------------------------------------------------------------------------
      CLLNEATPSMargCdf
      
      functionality

      Computes the continuous marginal cdf for X for the synthetic population
      using the basic assumption of the frequency estimation based on the fitted 
      bivariate loglinear model parameters for a discrete distribution
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar1      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 1
        bivar2      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 2
        wts         weights for population 1 in the synthetic population
        x           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        the function returns the smoothed pdf
    --------------------------------------------------------------------------*/
    double cllNEATPSMargCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                            const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2,
                            const double& wts,
                            const double& x) {
      // int j;
      // static double nc1, nc2;
      // double minv, maxv, minx, maxx, miny, maxy;
      // double vr, vm, dv, xr, xm, dx, fx, s;
      std::vector<double> fv; /* an array containing 4 arrays of V distribution at quadrature points,
					     1: f1v plus, 2: f1v minxs, 3: f2v plus, 4: f2v minxs, */
      fv.resize(4 * 16);

      std::vector<double> xx {0.04830766568774, 0.14447196158280, 0.23928736225214,
                              0.33186860228213, 0.42135127613064, 0.50689990893223,
                              0.58771575724076, 0.66304426693022, 0.73218211874029,
                              0.79448379596794, 0.84936761373257, 0.89632115576605,
                              0.93490607593774, 0.96476225558751, 0.98561151154527,
                              0.99726386184948};
      std::vector<double> w {0.09654008851473, 0.09563872007927, 0.09384439908080,
                             0.09117387869576, 0.08765209300440, 0.08331192422695,
                             0.07819389578707, 0.07234579410885, 0.06582222277636,
                             0.05868409347854, 0.05099805926238, 0.04283589802223,
                             0.03427386291302, 0.02539206530926, 0.01627439473091,
                             0.00701861000947};

      double minx = bivar1.minimumRawScoreX - 0.5;
      double maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double miny = bivar2.minimumRawScoreX - 0.5;
      double maxy = bivar2.minimumRawScoreX + bivar2.scoreIncrementX * static_cast<double>(bivar2.numberOfScoresX - 1) + 0.5;
      double xm = 0.5 * (x + minx);
      double xr = 0.5 * (x - minx);
      double minv = bivar1.minimumRawScoreV - 0.5;
      double maxv = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      double vm = 0.5 * (maxv + minv);
      double vr = 0.5 * (maxv - minv);
      double s = 0;
      double nc1 = bivGaussianQuadrature64(bivar1, minx, maxx, minv, maxv);
      double nc2 = bivGaussianQuadrature64(bivar2, miny, maxy, minv, maxv);

      for (size_t j = 0; j < 16; j++) {
        double dv = vr * xx[j];
        fv[j] = cllMargYPdf(bivar1, vm + dv, nc1);
        fv[j + 16] = cllMargYPdf(bivar1, vm - dv, nc1);
        fv[j + 32] = cllMargYPdf(bivar2, vm + dv, nc2);
        fv[j + 48] = cllMargYPdf(bivar2, vm - dv, nc2);
      }

      for (size_t j = 0; j < 16; j++) {
        double dx = xr * xx[j];
        double fx = cllNEATPSMargPdf(bivar1, bivar2, fv, wts, xm + dx) +
                    cllNEATPSMargPdf(bivar1, bivar2, fv, wts, xm - dx);
        s += w[j] * fx;
      }

      s *= xr;

      return s;
    }

    /*--------------------------------------------------------------------------
	CLLNEATPSInverseCdf
	
	functionality

	Computes the continuous marginal cdf for X for the synthetic population
	using the basic assumption of the frequency estimation based on the fitted 
	bivariate loglinear model parameters for a discrete distribution
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wts         weights for population 1 in the synthetic population
		x           a particular score for which the pdf and cdf is 
		            generated            

 
	output
 		the function returns the smoothed pdf
--------------------------------------------------------------------------*/

    double cllNEATPSInverseCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                               const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2,
                               const double& wts,
                               const double& cdf) {
      // int niter = 0;
      // double ub, lb, half, cdfu, cdfl, cdfhalf;
      // double absdif
      double eps = .00001;

      double lb = bivar1.minimumRawScoreX - 0.5;
      double ub = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double cdfl = cllNEATPSMargCdf(bivar1, bivar2, wts, lb);
      double cdfu = cllNEATPSMargCdf(bivar1, bivar2, wts, ub);

      double half;

      if (cdf < cdfl) {
        half = lb;
      } else if (cdf > cdfu) {
        half = ub;
      } else {
        for (size_t iter = 1; iter <= 200; iter++) {
          half = 0.5 * (lb + ub);
          double cdfhalf = cllNEATPSMargCdf(bivar1, bivar2, wts, half);
          double absdif = std::abs(cdf - cdfhalf);
          if (absdif < eps) {
            break;
          } else if (cdfhalf < cdf) {
            lb = half;
          } else {
            ub = half;
          }
        }
      }

      return half;
    }

    /*--------------------------------------------------------------------------
      CLLEquateNEATPS
      
      functionality

      Computes the equating fucntion for the NEAT design using the poststratification
      (frequency estimation) + CLL method
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar1      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 1
        bivar2      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 2
        wts         weight for the X group in the synthetic population            

    
      output
        Equatedx    The equating function 
    --------------------------------------------------------------------------*/
    // int
    Eigen::VectorXd cllEquateNEATPS(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                                    const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2,
                                    const double& wts) {
      Eigen::VectorXd equatedx(bivar1.numberOfScoresX);

      for (size_t i = 0; i < bivar1.numberOfScoresX; i++) {
        double cdfx = cllNEATPSMargCdf(bivar1, bivar2, wts, bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(i));
        equatedx(i) = cllNEATPSInverseCdf(bivar2, bivar1, 1.0 - wts, cdfx);
      }

      return equatedx;
    }

    /*--------------------------------------------------------------------------
      CLLMargInverseYCdf
      
      functionality

      Computes the inverse of continuous marginal cdf for Y for a bivariate distribution.
      
      Author: Tianyou Wang 4/29/2007.
      
      input
        bivar1      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 1
        bivar2      the structure that contain bivariate distribution parameters
                    (defined in "MLBivLogLin.h") for population 2
        ycdf        the cdf for y that is used to compute the inverse            

    
      output
        the function returns y score correspdonding the the ycdf
    --------------------------------------------------------------------------*/
    double cllMargInverseYCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                              const double& ycdf,
                              const double& nc) {
      // int niter = 0;
      // double ub, lb, half, cdfu, cdfl, cdfhalf;
      // double absdif;
      double eps = .00001;
      // double maxx;

      double lb = bivar1.minimumRawScoreV - 0.5;
      double ub = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      double maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double cdfl = cllBivCdf(bivar1, maxx, lb, nc);
      double cdfu = cllBivCdf(bivar1, maxx, ub, nc);

      double half;

      if (ycdf < cdfl) {
        half = lb;
      } else if (ycdf > cdfu) {
        half = ub;
      } else {
        for (size_t iter = 1; iter <= 200; iter++) {
          half = .5 * (lb + ub);

          double cdfhalf = cllBivCdf(bivar1, maxx, half, nc);

          double absdif = std::abs(ycdf - cdfhalf);
          if (absdif < eps) {
            break;
          } else if (cdfhalf < ycdf) {
            lb = half;
          } else {
            ub = half;
          }
        }
      }

      return half;
    }

    /*--------------------------------------------------------------------------
	CLLMargInverseXCdf
	
	functionality

	Computes the inverse of continuous marginal cdf for X for a bivariate distribution.
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		xcdf        the cdf for x that is used to compute the inverse            

 
	output
 		the function returns y score correspdonding the the ycdf
--------------------------------------------------------------------------*/
    double cllMargInverseXCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                              const double& xcdf,
                              const double& nc) {
      // int niter = 0;
      // double ub, lb, half, cdfu, cdfl, cdfhalf;
      // double absdif, eps = .00001;
      // double maxy;
      double eps = .00001;
      double half;

      double lb = bivar1.minimumRawScoreX - 0.5;
      double ub = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double maxy = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      double cdfl = cllBivCdf(bivar1, lb, maxy, nc);
      double cdfu = cllBivCdf(bivar1, ub, maxy, nc);

      if (xcdf < cdfl) {
        half = lb;
      } else if (xcdf > cdfu) {
        half = ub;
      } else {
        for (size_t iter = 1; iter <= 200; iter++) {
          half = 0.5 * (lb + ub);

          double cdfhalf = cllBivCdf(bivar1, half, maxy, nc);
          double absdif = std::abs(xcdf - cdfhalf);

          if (absdif < eps) {
            break;
          } else if (cdfhalf < xcdf) {
            lb = half;
          } else {
            ub = half;
          }
        }
      }

      return half;
    }

    /*--------------------------------------------------------------------------
	CLLMargInverseCBYCdf
	
	functionality

	Computes the inverse of continuous marginal cdf for Y for the weighted
	bivariate distribution of the CB design.
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
		wtsy        the weight for population 1 and Y
		ycdf        the cdf for y that is used to compute the inverse            

 
	output
 		the function returns y score correspdonding the the ycdf
--------------------------------------------------------------------------*/

    double cllMargInverseCBYCdf(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                                const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2,
                                const double& wtsy,
                                const double& ycdf,
                                const double& nc1,
                                const double& nc2) {
      // int niter = 0;
      // double ub, lb, half, cdfu, cdfl, cdfhalf;
      // double absdif, eps = .00001;
      // double maxx;

      double eps = .00001;
      double half;

      double lb = bivar1.minimumRawScoreV - 0.5;
      double ub = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      double maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double cdfl = wtsy * cllBivCdf(bivar1, maxx, lb, nc1) + (1 - wtsy) * cllBivCdf(bivar2, maxx, lb, nc2);
      double cdfu = wtsy * cllBivCdf(bivar1, maxx, ub, nc1) + (1 - wtsy) * cllBivCdf(bivar2, maxx, ub, nc2);

      if (ycdf < cdfl) {
        half = lb;
      } else if (ycdf > cdfu) {
        half = ub;
      } else {
        for (size_t iter = 1; iter <= 200; iter++) {
          half = 0.5 * (lb + ub);
          double cdfhalf = wtsy * cllBivCdf(bivar1, maxx, half, nc1) + (1 - wtsy) * cllBivCdf(bivar2, maxx, half, nc2);
          double absdif = std::abs(ycdf - cdfhalf);

          if (absdif < eps) {
            break;
          } else if (cdfhalf < ycdf) {
            lb = half;
          } else {
            ub = half;
          }
        }
      }

      return half;
    }

    /*--------------------------------------------------------------------------
	CLLEquateNEATChn
	
	functionality

	Computes the equating fucntion for the NEAT design using the chained 
	equipercentile + CLL method
	
	Author: Tianyou Wang 4/29/2007.
	
	input
		bivar1      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 1
		bivar2      the structure that contain bivariate distribution parameters
		            (defined in "MLBivLogLin.h") for population 2
 
	output
 		Equatedx    The equating function 
--------------------------------------------------------------------------*/
    // int
    Eigen::VectorXd cllEquateNEATChn(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                          const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2) {
      double minx = bivar1.minimumRawScoreX - 0.5;
      double maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double miny = bivar2.minimumRawScoreX - 0.5;
      double maxy = bivar2.minimumRawScoreX + bivar2.scoreIncrementX * static_cast<double>(bivar2.numberOfScoresX - 1) + 0.5;
      double minv = bivar1.minimumRawScoreV - 0.5;
      double maxv = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      double nc1 = bivGaussianQuadrature64(bivar1, minx, maxx, minv, maxv);
      double nc2 = bivGaussianQuadrature64(bivar2, miny, maxy, minv, maxv);
      
      maxv = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      maxy = bivar2.minimumRawScoreX + bivar2.scoreIncrementX * static_cast<double>(bivar2.numberOfScoresX - 1) + 0.5;

      Eigen::VectorXd equatedx(bivar1.numberOfScoresX);

      for (size_t i = 0; i < bivar1.numberOfScoresX; i++) {
        double cdfx = cllBivCdf(bivar1, bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(i), maxv, nc1);
        double eqtempv = cllMargInverseYCdf(bivar1, cdfx, nc1);
        double cdfv = cllBivCdf(bivar2, maxx, eqtempv, nc2);
        equatedx(i) = cllMargInverseXCdf(bivar2, cdfv, nc2);
      }

      return equatedx;
    }

    /*--------------------------------------------------------------------------
      CLLEquateSG
      
      functionality

      Computes the equating fucntion for the single group design using the  
      equipercentile + CLL method
      
      Author: Tianyou Wang 7/10/2007.
      
      input
        bivar      the structure that contain bivariate distribution of X and 
                  Y (defined in "MLBivLogLin.h") 
    
      output
        Equatedx    The equating function 
    --------------------------------------------------------------------------*/
    // int
    Eigen::VectorXd cllEquateSG(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar) {
      double minx = bivar.minimumRawScoreX;
      double miny = bivar.minimumRawScoreV;
      double maxx = minx + static_cast<double>(bivar.numberOfScoresX - 1) * bivar.scoreIncrementX + 0.5;
      double maxy = miny + static_cast<double>(bivar.numberOfScoresV - 1) * bivar.scoreIncrementV + 0.5;
      double nc = bivGaussianQuadrature64(bivar, minx, maxx, miny, maxy);

      Eigen::VectorXd equatedx(bivar.numberOfScoresX);

      for (size_t i = 0; i < bivar.numberOfScoresX; i++) {
        double cdfx = cllBivCdf(bivar, bivar.minimumRawScoreX + bivar.scoreIncrementX * static_cast<double>(i), maxy, nc);
        equatedx(i) = cllMargInverseYCdf(bivar, cdfx, nc);
      }

      return equatedx;
    }

    /*--------------------------------------------------------------------------
      CLLEquateCB
      
      functionality

      Computes the equating fucntion for the counter balance design using the  
      equipercentile + CLL method
      
      Author: Tianyou Wang 7/10/2007.
      
      input
        bivar1      the structure that contain bivariate distribution (X and Y) 
                    (defined in "MLBivLogLin.h") for population 1
        bivar2      the structure that contain bivariate distribution (X and Y)
                    (defined in "MLBivLogLin.h") for population 2
          wtsx        weight for population 1 and X
        wtsy        weight for population 1 and Y
    
      output
        Equatedx    The equating function 
    --------------------------------------------------------------------------*/
    // int
    Eigen::VectorXd cllEquateCB(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar1,
                                const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivar2,
                                const double& wtsx,
                                const double& wtsy) {
      // int i, nDistCatx;
      // double maxx, maxy;
      // double cdfx, nc1, nc2;

      double maxx = bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(bivar1.numberOfScoresX - 1) + 0.5;
      double maxy = bivar1.minimumRawScoreV + bivar1.scoreIncrementV * static_cast<double>(bivar1.numberOfScoresV - 1) + 0.5;
      // nDistCatx = bivar1->nsx;
      double nc1 = bivGaussianQuadrature32(bivar1, bivar1.minimumRawScoreU, maxx, bivar1.minimumRawScoreV, maxy);
      double nc2 = bivGaussianQuadrature32(bivar2, bivar1.minimumRawScoreX, maxx, bivar1.minimumRawScoreV, maxy);

      Eigen::VectorXd equatedx(bivar1.numberOfScoresX);

      for (size_t i = 0; i < bivar1.numberOfScoresX; i++) {
        double cdfx = wtsx * cllBivCdf(bivar1, bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(i), maxy, nc1) +
                      (1 - wtsx) * cllBivCdf(bivar2, bivar1.minimumRawScoreX + bivar1.scoreIncrementX * static_cast<double>(i), maxy, nc2);
        equatedx(i) = cllMargInverseCBYCdf(bivar1, bivar2, wtsy, cdfx, nc1, nc2);
      }

      return equatedx;
    }

    /*--------------------------------------------------------------------------
      CLLSEEEG
      
      functionality:

      Computes equating function based on continuized Log-linear cdf in Wang 	(2005). 
          
      author: Tianyou Wang 2/22/2006.
      
      input:
            nparax      number of parameters for the new form
            paraxx		a vector of parameters for the loglinear model
                        with a design matrix from polynominals of natural
                        basis for the new form.
            nparay      number of parameters for the old form
            paraxy		a vector of parameters for the loglinear model
                        with a design matrix from polynominals of natural
                        basis for the old form.
        nCatx		Number of discrete score categories for the new form
        Bx          design matrix from the new form loglinear model

      output:
        SEE        a vector containing the standard error of equating
    --------------------------------------------------------------------------*/
    Eigen::VectorXd cllSEEEG(const long& npx,
                             const std::vector<double>& parax,
                             const std::vector<double>& scoresx,
                             const long& npy,
                             const std::vector<double>& paray,
                             const std::vector<double>& scoresy) {
      // int i, j;
      // double cdfx, pdfy;
      // double **Sigmamx, **Sigmamy, *px, *py, **Sigmabetax, **Sigmabetay, **Bx, **By;
      // double *eYbetax, *eYbetay;
      // double integx, integx2, integy, integy2, minx, maxx, miny, maxy, Equatedx, integ1, integ2, integ3;

      Eigen::MatrixXd Sigmamx(scoresx.size(), scoresx.size()); // Sigmamx = dmatrix(0, nCatx - 1, 0, nCatx - 1);
      Eigen::MatrixXd Sigmamy(scoresy.size(), scoresy.size()); // Sigmamy = dmatrix(0, nCaty - 1, 0, nCaty - 1);
      Eigen::MatrixXd Sigmabetax(parax.size(), parax.size());  // Sigmabetax = dmatrix(0, nparax - 1, 0, nparax - 1);
      Eigen::MatrixXd Sigmabetay(paray.size(), paray.size());  // Sigmabetay = dmatrix(0, nparay - 1, 0, nparay - 1);
      Eigen::MatrixXd Bx(scoresx.size(), parax.size());        // Bx = dmatrix(0, nCatx - 1, 0, nparax - 1);
      Eigen::MatrixXd By(scoresy.size(), paray.size());        // By = dmatrix(0, nCaty - 1, 0, nparay - 1);
      Eigen::MatrixXd px(scoresx.size(), 1);                   // px = dvector(0, nCatx - 1);
      Eigen::MatrixXd py(scoresy.size(), 1);                   // py = dvector(0, nCaty - 1);
      Eigen::MatrixXd eYbetax(parax.size(), 1);                // eYbetax = dvector(0, nparax - 1);
      Eigen::MatrixXd eYbetay(paray.size(), 1);                // eYbetay = dvector(0, nparay - 1);
      Eigen::VectorXd SEE(scoresx.size());

      for (size_t i = 0; i < parax.size(); i++) {
        for (size_t j = 0; j < scoresx.size(); j++) {
          if (j == 0) {
            Bx(j, i) = 0.0;
          } else {
            Bx(j, i) = std::pow(static_cast<double>(j), static_cast<double>(i));
          }
        }
      }
      Bx(0, 0) = 1.0;

      for (size_t i = 0; i < paray.size(); i++) {
        for (size_t j = 0; j < scoresy.size(); j++) {
          if (j == 0) {
            By(j, i) = 0;
          } else {
            By(j, i) = pow(static_cast<double>(j), static_cast<double>(i));
          }
        }
      }
      By(0, 0) = 1;

      for (size_t i = 0; i < scoresx.size(); i++) {
        for (size_t j = 0; j < scoresx.size(); j++) {
          Sigmamx(i, j) = -1.0 * px(i) * px(j) / npx;
          if (i == j) {
            Sigmamx(i, j) += px(i);
          }
        }
      }

      for (size_t i = 0; i < scoresy.size(); i++) {
        for (size_t j = 0; j < scoresy.size(); j++) {
          Sigmamy(i, j) = -1.0 * py(i) * py(j) / npy;
          if (i == j) {
            Sigmamy(i, j) += py(i);
          }
        }
      }

      Sigmabetax = Bx.transpose() * Sigmamx * Bx;
      Sigmabetay = By.transpose() * Sigmamy * By;

      Eigen::MatrixXd InverseSigmabetax = Sigmabetax.inverse();
      Eigen::MatrixXd InverseSigmabetay = Sigmabetay.inverse();

      double minx = scoresx.front() - 0.5;
      double maxx = scoresx.back() + 0.5;

      double miny = scoresy.front() - 0.5;
      double maxy = scoresy.back() + 0.5;
      double integx = gaussianQuadrature64(minx, maxx, parax);
      double integx2 = integx * integx;
      double integy = gaussianQuadrature64(miny, maxy, paray);
      double integy2 = integy * integy;

      for (size_t i = 0; i < scoresx.size(); i++) {
        double cdfx = cllEGCdf(minx, maxx, parax, scoresx[i], integx);
        double equatedx = cllInverseCdf(miny, maxy, paray, cdfx, integy);
        double pdfy = cllEGPdf(miny, maxy, paray, equatedx, integy);
        double integ2 = gaussianQuadrature64(minx, scoresx[i], parax);

        for (size_t j = 0; j < parax.size(); j++) {
          double integ1 = gaussianQuadrature64i(minx, scoresx[i], parax, j);
          double integ3 = gaussianQuadrature64i(minx, maxx, parax, j);
          eYbetax(j) = integ1 * integx - integ2 * integ3;
          eYbetax(j) /= integx2;
          eYbetax(j) /= pdfy;
        }

        integ2 = gaussianQuadrature64(miny, equatedx, paray);
        for (size_t j = 0; j < paray.size(); j++) {
          double integ1 = gaussianQuadrature64i(miny, equatedx, paray, j);
          double integ3 = gaussianQuadrature64i(miny, maxy, paray, j);
          eYbetay(j) = integ2 * integ3 - integ1 * integy;
          eYbetay(j) /= integy2;
          eYbetay(j) /= pdfy;
        }

        SEE(i) = (eYbetax.transpose() * Sigmabetax * eYbetax).coeff(0, 0);  // vtAv0(nparax, eYbetax, Sigmabetax);
        SEE(i) += (eYbetay.transpose() * Sigmabetay * eYbetay).coeff(0, 0); // vtAv0(nparay, eYbetay, Sigmabetay);
        SEE(i) = sqrt(SEE(i));
      }

      return SEE;
    }
  };
} // namespace EquatingRecipes

#endif