/*
kernel_Equate.c   code for kernel equating

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

#ifndef KERNEL_EQUATING_HPP
#define KERNEL_EQUATING_HPP

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/QR>
#include <boost/math/distributions/normal.hpp>

#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/smoothing.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/structures/univariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>

#include <equating_recipes/implementation/cg_equipercentile_equating.hpp>
#include <equating_recipes/implementation/log_linear_equating.hpp>
#include <equating_recipes/implementation/utilities.hpp>

namespace EquatingRecipes {
  namespace Implementation {
    class KernelEquating {
    public:
      /*
			Wrapper for doing kernel equating with RG design
				and log-linear smoothing smoothing
				
			Assumes that equating puts raw scores for x on scale of y
			
			NOTE: This function is used (unaltered) for both actual equating and 
						equating done in Wrapper_Bootstrap().  Distinguishing between the
						two is the purpose of the variable rep

			Input
			
				design = 'R'(random groups)
				method = 'E'(equipercentile equating)
				smoothing = 'K' (kernel 'smoothing')  
				*x = pointer to struct USTATS (new form)
				*y = pointer to struct USTATS (old form)
				*ullx = pointer to struct ULL_SMOOTH (new form)
				*ully = pointer to struct ULL_SMOOTH (old form)
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
					
			NOTE: If Wrapper_RK() is called in a bootstrap loop,
						then in the calling function struct ERAW_RESULTS must
						be different from struct ERAW_RESULTS for the actual
						equating. 

			Function calls other than C or NR utilities:                   
				KernelEquateSEERG()
				MomentsFromFD()  
																										
			Tianyou Wang

			Date of last revision: 6/30/08       
		*/
      void runWithRGDesign(const EquatingRecipes::Structures::Design& design,
                           const EquatingRecipes::Structures::Method& method,
                           const EquatingRecipes::Structures::Smoothing& smoothing,
                           const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsX,
                           const EquatingRecipes::Structures::UnivariateStatistics& univariateStatisticsY,
                           const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& univariateLogLinearSmoothingX,
                           const EquatingRecipes::Structures::UnivariateLogLinearSmoothing& univariateLogLinearSmoothingY,
                           const size_t& replicationNumber,
                           EquatingRecipes::Structures::PData& pData,
                           EquatingRecipes::Structures::EquatedRawScoreResults& results) {
        /* method name --- 10 characters; right justified */
        std::vector<std::string> methodNames {"   y-equiv"};

        Eigen::VectorXd scoresX(0, univariateLogLinearSmoothingX.numberOfScores);
        Eigen::VectorXd scoresY(0, univariateLogLinearSmoothingY.numberOfScores);

        for (size_t scoreIndex = 0; scoreIndex < univariateLogLinearSmoothingX.numberOfScores; scoreIndex++) {
          scoresX(scoreIndex) = univariateLogLinearSmoothingX.minimumRawScore +
                                static_cast<double>(scoreIndex) * univariateLogLinearSmoothingX.rawScoreIncrement;
        }

        Eigen::VectorXd scoresY(univariateLogLinearSmoothingY.numberOfScores);
        for (size_t scoreIndex = 0; scoreIndex < univariateLogLinearSmoothingY.numberOfScores; scoreIndex++) {
          scoresY(scoreIndex) = univariateLogLinearSmoothingY.minimumRawScore +
                                static_cast<double>(scoreIndex) * univariateLogLinearSmoothingY.rawScoreIncrement;
        }

        size_t maximumScoreLocation = EquatingRecipes::Implementation::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                                   pData.minimumScoreX,
                                                                                                   pData.scoreIncrementX);

        pData.bootstrapReplicationNumber = replicationNumber;

        /* allocation and assignments for struct PDATA inall
          Note that for every assignment of the form inall->(var) = x->(var)
          or inall->(var) = y->(var), values vary depending on whether x or y
          is for actual equating or a bootstrap sample; all other values are
          the same for the actual equating and a bootstrap sample */

        if (pData.bootstrapReplicationNumber == 0) {
          pData.summaryRawDataX = univariateStatisticsX;
          pData.summaryRawDataY = univariateStatisticsY;
          pData.design = design;
          pData.method = method;
          pData.smoothing = smoothing;
          pData.methods = methodNames;
          pData.minimumScoreX = univariateStatisticsX.minimumScore;
          pData.maximumScoreX = univariateStatisticsX.maximumScore;
          pData.scoreIncrementX = univariateStatisticsX.scoreIncrement;
          pData.scoreFrequenciesX = univariateStatisticsX.freqDistDouble;
          pData.numberOfExaminees = univariateStatisticsX.numberOfExaminees;

          pData.univariateLogLinearSmoothingX = univariateLogLinearSmoothingX;
          pData.univariateLogLinearSmoothingY = univariateLogLinearSmoothingY;
        }

        /* allocation and assignments for r */

        if (pData.bootstrapReplicationNumber <= 1) { /* no storage allocation for bootstrap reps >1 */
          results.equatedRawScores.resize(1, maximumScoreLocation);
          results.equatedRawScoreMoments.resize(1, 4);
          results.equatingStandardErrors.resize(1, maximumScoreLocation + 1);
        }

        /* Compute equating results */
        Eigen::VectorXd equatedRawScores(maximumScoreLocation + 1);
        Eigen::VectorXd equatingStandardErrors(maximumScoreLocation + 1);

        kernelEquateSEERG(univariateLogLinearSmoothingX.numberOfScores,
                          univariateLogLinearSmoothingX.degreesOfSmoothing,
                          univariateLogLinearSmoothingX.numberOfExaminees,
                          scoresX,
                          univariateLogLinearSmoothingX.fittedRawScoreDist,
                          univariateLogLinearSmoothingY.numberOfScores,
                          univariateLogLinearSmoothingY.degreesOfSmoothing,
                          univariateLogLinearSmoothingY.numberOfExaminees,
                          scoresY,
                          univariateLogLinearSmoothingY.fittedRawScoreDist,
                          equatedRawScores,
                          equatingStandardErrors);

        results.equatedRawScores.row(0) = equatedRawScores;
        results.equatingStandardErrors.row(0) = equatingStandardErrors;

        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Implementation::Utilities::momentsFromScoreFrequencies(results.equatedRawScores.row(0),
                                                                                                                               pData.scoreFrequenciesX);

        results.equatedRawScoreMoments.row(0) = moments.momentValues;
      }

      /*
			Wrapper for doing kernel equating with SG design
			and log-linear smoothing.  
			
			NOTE: This is for the SG design in which x and y do not share any items in 
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
			
				design = 'S' (single group)
				method = 'E' (equipercentile equating)
				smoothing = 'K' (kernel 'smoothing')  
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
					
			NOTE: If Wrapper_SK() is called in a bootstrap loop,
						then in the calling function struct ERAW_RESULTS must
						be different from struct ERAW_RESULTS for the actual
						equating. 
																								
			Function calls other than C or NR utilities:
				KernelEquateSG()
				MomentsFromFD()  
																										
			Tianyou Wang

			Date of last revision: 6/30/08   
		*/
      void runWithSGDesign(const EquatingRecipes::Structures::Design& design,
                           const EquatingRecipes::Structures::Method& method,
                           const EquatingRecipes::Structures::Smoothing& smoothing,
                           const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXY,
                           const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingXY,
                           const size_t& replicationNumber,
                           EquatingRecipes::Structures::PData& pData,
                           EquatingRecipes::Structures::EquatedRawScoreResults& results) {
        /* method names --- 10 characters; right justified */
        std::vector<std::string> methodNames {"   y-equiv"};

        pData.bootstrapReplicationNumber = replicationNumber;

        /* Allocation and assignments for struct PDATA inall>
          Note that for every assignment of the form inall->(var) = x->(var)
          or inall->(var) = y->(var), values vary depending on whether x or y
          is for actual equating or a bootstrap sample; all other values are
          the same for the actual equating and a bootstrap sample */

        if (pData.bootstrapReplicationNumber == 0) {
          pData.summaryRawDataXY = bivariateStatisticsXY;
          pData.bivariateLogLinearSmoothingXY = bivariateLogLinearSmoothingXY;
          pData.design = design;
          pData.method = method;
          pData.smoothing = smoothing;
          pData.isInternalAnchor = false; /* implicitly, anchor is external for biv log-linear smoothing with the SG design */
          pData.methods = methodNames;
          pData.minimumScoreX = bivariateStatisticsXY.univariateStatisticsRow.minimumScore;
          pData.maximumScoreX = bivariateStatisticsXY.univariateStatisticsRow.maximumScore;
          pData.scoreIncrementX = bivariateStatisticsXY.univariateStatisticsRow.scoreIncrement;
          pData.scoreFrequenciesX = bivariateStatisticsXY.univariateStatisticsRow.freqDistDouble;
          pData.numberOfExaminees = bivariateStatisticsXY.numberOfExaminees;
        }

        /* allocation and assignments for r */

        size_t maximumScoreLocation = EquatingRecipes::Implementation::Utilities::getScoreLocation(pData.maximumScoreX,
                                                                                                   pData.minimumScoreX,
                                                                                                   pData.scoreIncrementX);

        if (pData.bootstrapReplicationNumber <= 1) {
          results.equatedRawScores.resize(1, maximumScoreLocation + 1);
          results.equatedRawScoreMoments.resize(1, 4);
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
        Eigen::VectorXd equatedScores(maximumScoreLocation + 1);
        kernelEquateSG(bivariateLogLinearSmoothingXY, equatedScores);
        results.equatedRawScoreMoments.row(0) = equatedScores;

        /* get moments */
        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Implementation::Utilities::momentsFromScoreFrequencies(equatedScores,
                                                                                                                               pData.scoreFrequenciesX);
      }

      /*
      Wrapper for conitnuized log-linear equating for CG design with log-linear smoothing. 
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

        method:  'E' = Frequency estimation (FE) with Braun-Holland (BH) under FE
                'F' = Modified freq est (MFE) with Braun-Holland (BH) under MFE
                'G' = FE + BH-FE + MFE + BH-MFE
                'C' = Chained
          'H' = FE + BH-FE + Chained
                'A' = FE + BH-FE + MFE + BH-MFE + Chained
                  
        smoothing = 'K' (kernel 'smoothing') 

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
      
        NOTE: if rv1 == 0 or rv2 == 0, then MFE cannot be obtained 
        
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
        
      NOTE: If Wrapper_CK() is called in a bootstrap loop, then in
            the calling function struct ERAW_RESULTS must be different
            from struct ERAW_RESULTS for the actual equating. 
                                                
      Function calls other than C or NR utilities:
        KernelEquateNEATPS()
        KernelEquateNEATChn()
        runerror()
                                                  
      Tianyou Wang

      Date of last revision: 6/30/08   
    */
      void runWithCGDesign(const EquatingRecipes::Structures::Design& design,
                           const EquatingRecipes::Structures::Method& method,
                           const EquatingRecipes::Structures::Smoothing& smoothing,
                           const double& w1,
                           const bool& isInternalAnchor,
                           const double& reliabilityCommonItemsPopulation1,
                           const double& reliabilityCommonItemsPopulation2,
                           const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXV,
                           const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsYV,
                           const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingXV,
                           const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingYV,
                           const size_t& replicationNumber,
                           EquatingRecipes::Structures::PData& pData,
                           EquatingRecipes::Structures::EquatedRawScoreResults& results) {
        std::vector<std::string> methods {"        FE",
                                          "       MFE",
                                          "  ChainedE"};

        pData.bootstrapReplicationNumber = replicationNumber;

        std::string methodCode = EquatingRecipes::Implementation::Utilities::getMethodCode(method);

        /* allocation and assignments for inall
           Note that for every assignment of the form inall->(var) = xv->(var)
           or inall->(var) = yv->(var) values vary depending on whether xv or yv
           is for actual equating or a bootstrap sample; all other values are
           the same for the actual equating and a bootstrap sample */

        if (pData.bootstrapReplicationNumber == 0) {
          pData.summaryRawDataXV = bivariateStatisticsXV;
          pData.summaryRawDataYV = bivariateStatisticsYV;
          pData.bivariateLogLinearSmoothingXV = bivariateLogLinearSmoothingXV;
          pData.bivariateLogLinearSmoothingYV = bivariateLogLinearSmoothingYV;
          pData.design = design;
          pData.method = method;
          pData.smoothing = smoothing;

          if (w1 < 0.0 || w1 > 1.0) {
            pData.weightSyntheticPopulation1 = static_cast<double>(bivariateStatisticsXV.numberOfExaminees) /
                                               static_cast<double>(bivariateStatisticsXV.numberOfExaminees + bivariateStatisticsYV.numberOfExaminees);
          } else {
            pData.weightSyntheticPopulation1 = w1;
          }

          pData.isInternalAnchor = isInternalAnchor;
          pData.reliabilityCommonItemsPopulation1 = reliabilityCommonItemsPopulation1;
          pData.reliabilityCommonItemsPopulation2 = reliabilityCommonItemsPopulation2;

          if ((methodCode == "F" || methodCode == "G" || methodCode == "A") &&
              (reliabilityCommonItemsPopulation1 == 0 || reliabilityCommonItemsPopulation2 == 0)) {
            throw std::runtime_error("MFE cannot be conducted since rv1 == 0 or rv2 == 0");
          }

          if (methodCode == "E") {
            pData.methods.push_back(methods[0]);
          } else if (methodCode == "F") {
            pData.methods.push_back(methods[1]);
          } else if (methodCode == "G") {
            pData.methods.push_back(methods[0]);
            pData.methods.push_back(methods[1]);
          } else if (methodCode == "C") {
            pData.methods.push_back(methods[2]);
          } else if (methodCode == "H") {
            pData.methods.push_back(methods[0]);
            pData.methods.push_back(methods[2]);
          } else {
            // method == 'A'
            pData.methods = methods;
          }

          pData.minimumScoreX = bivariateStatisticsXV.univariateStatisticsRow.minimumScore;
          pData.maximumScoreX = bivariateStatisticsXV.univariateStatisticsRow.maximumScore;
          pData.scoreIncrementX = bivariateStatisticsXV.univariateStatisticsRow.scoreIncrement;
          pData.scoreFrequenciesX = bivariateStatisticsXV.univariateStatisticsRow.freqDistDouble;
          pData.numberOfExaminees = bivariateStatisticsXV.numberOfExaminees;
        }

        size_t maximumScoreLocationX = EquatingRecipes::Implementation::Utilities::getScoreLocation(bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                                                                                    bivariateStatisticsXV.univariateStatisticsRow.maximumScore,
                                                                                                    bivariateStatisticsXV.univariateStatisticsRow.scoreIncrement);

        size_t maximumScoreLocationY = EquatingRecipes::Implementation::Utilities::getScoreLocation(bivariateStatisticsYV.univariateStatisticsRow.maximumScore,
                                                                                                    bivariateStatisticsYV.univariateStatisticsRow.minimumScore,
                                                                                                    bivariateStatisticsYV.univariateStatisticsRow.scoreIncrement);

        /* allocation and assignments for r */
        if (pData.bootstrapReplicationNumber <= 1) {
          results.equatedRawScores.resize(pData.methods.size(), maximumScoreLocationX + 1);
          results.equatedRawScoreMoments.resize(pData.methods.size(), 4);
          results.relativeFreqDistsX.resize(2, maximumScoreLocationX + 1);
          results.relativeFreqDistsY.resize(2, maximumScoreLocationY + 1);
        }

        /* Equipercentile results, including Braun-Holland (BH) linear.
          Note: For FE syn densities are in fxs[0] and gys[0]
          For MFE syn densities are in fxs[1] and gys[1]
          For BH under FE, slope in a[0] and intercept in b[0]
          For BH under MFE, slope in a[1] and intercept in b[1] */

        /* FE + BH-FE in positions 0 and 1*/
        Eigen::VectorXd equatedRawScores(maximumScoreLocationX + 1);
        if (methodCode == "E" || methodCode == "G" || methodCode == "A" || methodCode == "H") {
          kernelEquateNEATPS(bivariateLogLinearSmoothingXV,
                             bivariateLogLinearSmoothingYV,
                             pData.weightSyntheticPopulation1,
                             equatedRawScores);

          results.equatedRawScores.row(0) = equatedRawScores;
        }

        size_t methodIndex;
        bool isNEATChainedMethod = false;
        if (methodCode == "C") {
          methodIndex = 0;
          isNEATChainedMethod = true;
        } else if (methodCode == "A") {
          methodIndex = 2;
          isNEATChainedMethod = true;
        } else if (methodCode == "H") {
          methodIndex = 1;
          isNEATChainedMethod = true;
        }

        if (isNEATChainedMethod) {
          equatedRawScores = results.equatedRawScores.row(methodIndex);

          kernelEquateNEATChn(bivariateLogLinearSmoothingXV,
                              bivariateLogLinearSmoothingYV,
                              equatedRawScores);

          results.equatedRawScores.row(methodIndex) = equatedRawScores;
        }

        /* get moments */
        for (methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
          equatedRawScores = results.equatedRawScores.row(methodIndex);
          EquatingRecipes::Structures::Moments moments = EquatingRecipes::Implementation::Utilities::momentsFromScoreFrequencies(equatedRawScores,
                                                                                                                                 pData.scoreFrequenciesX);
          results.equatedRawScoreMoments.row(methodIndex) = moments.momentValues;
        }
      }

      /*
      Wrapper for doing kernel equating with SG with counter-balance design 
      and log-linear smoothing.  
      
      NOTE: This is for the SG with counter balance design in which x and y 
      in both group 1 and group 2 do not share any items in 
      common, which means that functionally this is the external anchor case.  
      The bivariate log-linear smoothing procedure needs to know this. So, when
      Wrapper_Smooth_BLL() is called (as it must be prior to calling Wrapper_SL()),
      anchor must be set to 0. If x and y share common items, Wrapper_Smooth_BLL()
      (with anchor set to 0) and Wrapper_SL() can still be used, but convergence 
      of the smoothing algorithm may be compromised because of dependencies between
      x and y. (To date my experience does not suggest this is a problem.)
        
      Assumes that equating puts raw scores for x on scale of y
      
      NOTE: This function is used (unaltered) for both actual equating and 
            equating done in Wrapper_Bootstrap().  Distinguishing between the
            two is the purpose of the variable rep

      Input
      
        design = 'B' (single group w/ counter balance)
        method = 'E' (equipercentile equating)
        smoothing = 'K' (kernel 'smoothing')  
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
          
      NOTE: If Wrapper_SK() is called in a bootstrap loop,
            then in the calling function struct ERAW_RESULTS must
            be different from struct ERAW_RESULTS for the actual
            equating. 
                                                
      Function calls other than C or NR utilities:
        KernelEquateSEECB()
        MomentsFromFD()  
                                                    
      Tianyou Wang

      Date of last revision: 6/30/08   
    */
      void runWithSGCounterBalanceDesign(const EquatingRecipes::Structures::Design& design,
                                         const EquatingRecipes::Structures::Method& method,
                                         const EquatingRecipes::Structures::Smoothing& smoothing,
                                         const double& wtsX,
                                         const double& wtsY,
                                         const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXY1,
                                         const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXY2,
                                         const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingXY1,
                                         const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingXY2,
                                         const size_t& replicationNumber,
                                         EquatingRecipes::Structures::PData& pData,
                                         EquatingRecipes::Structures::EquatedRawScoreResults& results) {
        std::vector<std::string> methods {"   y-equiv"};

        Eigen::VectorXd wts(2);
        wts(0) = wtsX;
        wts(1) = wtsY;

        pData.bootstrapReplicationNumber = replicationNumber;

        if (pData.bootstrapReplicationNumber == 0) {
          pData.design = design;
          pData.method = method;
          pData.smoothing = smoothing;
          pData.isInternalAnchor = false; /* implicitly, anchor is external for biv log-linear smoothing with the SG design */
          pData.methods = methods;
        }

        /* allocation and assignments for r */

        size_t maximumScoreLocation = EquatingRecipes::Implementation::Utilities::getNumberOfScores(pData.maximumScoreX,
                                                                                                    pData.minimumScoreX,
                                                                                                    pData.scoreIncrementX);

        if (pData.bootstrapReplicationNumber <= 1) {
          results.equatedRawScores.resize(1, maximumScoreLocation + 1);
          results.equatedRawScoreMoments.resize(1, 4);
          results.equatingStandardErrors.resize(1, maximumScoreLocation + 1);
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
        Eigen::VectorXd equatedRawScores(maximumScoreLocation + 1);
        Eigen::VectorXd equatingStandardErrors(maximumScoreLocation + 1);

        kernelEquateSEECB(bivariateLogLinearSmoothingXY1,
                          bivariateLogLinearSmoothingXY2,
                          wts,
                          equatedRawScores,
                          equatingStandardErrors);

        results.equatedRawScores.row(0) = equatedRawScores;
        results.equatingStandardErrors.row(0) = equatingStandardErrors;

        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Implementation::Utilities::momentsFromScoreFrequencies(
          equatedRawScores,
          pData.scoreFrequenciesX
        );

        results.equatedRawScoreMoments.row(0) = moments.momentValues;
      }

    private:
      boost::math::normal_distribution<> normalDist;

      /*--------------------------------------------------------------------------
      KernelContinuPdf
      
      functionality

      Computes kernel smoohxing function from von Davier, Holland, & Thayer 
      (2004). The Kernel Method of Test Equating.
      
      author: Tianyou Wang 10/29/2004.
      
      input
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        fd          vector containing the relative frequency distribution
        hx          bandwidhx for the kernel smoohxing
        x           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        The function returns the kernel pdf
    --------------------------------------------------------------------------*/
      double kernelPdf(const size_t& numberOfScores,
                       const Eigen::VectorXd& scores,
                       const Eigen::VectorXd& scoreRelativeFreqDist,
                       const double& bandwithX,
                       const double& scoreValue) {
        double zscore; /*the standardized normal random variable */

        // double ax, sumsq = 0.0; /*ax, and second moment */
        // double temp1;           /*temporary variable */
        double pdf = 0;

        /*the mean and standard deviation of the unsmoothed distribution */
        double mu = scores.dot(scoreRelativeFreqDist);
        double sigma = (scores.cwiseProduct(scores)).dot(scoreRelativeFreqDist) - std::pow(mu, 2);

        double ax = std::sqrt(sigma / (sigma + std::pow(bandwithX, 2)));

        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          double zscore = (scoreValue - ax * scores(scoreIndex) - (1.0 - ax) * mu) / ax / bandwithX;
          double temp1 = stdNormalPdf(zscore);
          pdf += scoreRelativeFreqDist(scoreIndex) * temp1 / ax / bandwithX;
        }

        return pdf;
      }

      /*--------------------------------------------------------------------------
      KernelContinuCdf
      
      functionality:

      Computes kernel smoohxing function from von Davier, Holland, & Thayer 
      (2004). The Kernel Method of Test Equating.
      
      author: Tianyou Wang 10/29/2004.
      
      input:
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        fd          vector containing the relative frequency distribution
        hx          bandwidhx for the kernel smoohxing
        x           a particular score for which the pdf and cdf is 
                    generated            

    
      output:
        The function returns the kernel cdf
    --------------------------------------------------------------------------*/
      double kernelCdf(const size_t& numberOfScores,
                       const Eigen::VectorXd& scores,
                       const Eigen::VectorXd& scoreRelativeFreqDist,
                       const double& bandwithX,
                       const double& scoreValue) {
        double mu = scores.dot(scoreRelativeFreqDist);
        double sigma = (scores.cwiseProduct(scores)).dot(scoreRelativeFreqDist) - std::pow(mu, 2);

        double ax = sqrt(sigma / (sigma + std::pow(bandwithX, 2)));

        double cdf = 0;

        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          double zscore = (scoreValue - ax * scores(scoreIndex) - (1.0 - ax) * mu) / ax / bandwithX;
          double temp1 = stdNormalCdf(zscore);
          cdf += scoreRelativeFreqDist(scoreIndex) * temp1;
        }

        return cdf;
      }

      /*--------------------------------------------------------------------------
      StdNormalCdf
      
      functionality

      Computes the cdf for a z score from a standard normal distribution
      
      author: Tianyou Wang 10/29/2004.
      
      input
        x           a z score        

    
      output
        the function returns the normal cdf for x.
    --------------------------------------------------------------------------*/
      double stdNormalCdf(const double& x) {
        double result = boost::math::cdf(normalDist, x);
        return result;
      }

      /*--------------------------------------------------------------------------
      StdNormalPdf
      
      functionality

      Computes the pdf for a z score from a standard normal distribution
      
      author: Tianyou Wang 10/29/2004.
      
      input
        x           a z score        

    
      output
        the function returns the normal cdf for x.
    --------------------------------------------------------------------------*/
      double stdNormalPdf(const double& x) {
        double result = boost::math::pdf(normalDist, x);
        return result;
      }

      /*--------------------------------------------------------------------------
      NormalPdf
      
      functionality

      Computes the pdf for a z score from a normal distribution
      
      author: Tianyou Wang 10/29/2004.
      
      input
        x           a z score        

    
      output
        the function returns the normal cdf for x.
    --------------------------------------------------------------------------*/
      double normalPdf(const Eigen::VectorXd& para,
                       const double& x) {
        double mu;
        double sigma;
        if (para.size() == 2) {
          mu = para(0);
          sigma = para(1);
        } else {
          mu = 0.0;
          sigma = 1.0;
        }

        boost::math::normal_distribution<> unstdNormalDist(mu, sigma);
        double result = boost::math::pdf(unstdNormalDist, x);

        return result;
      }

      /*------------------------------------------------------------------------------
      CalcKernelMoments	
      
      functionality

      calculates mean, sd, skewness and kurtosis for a continuous distribution.

      author: Tianyou Wang 10/29/2004.

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
      void calcKernelMoments(const double& a,
                             const double& b,
                             const size_t& numberOfScores,
                             const Eigen::VectorXd& score,
                             const Eigen::VectorXd& scoreRelativeFreqDist,
                             const double& bandwithX,
                             Eigen::VectorXd& moments) {
        size_t j;
        double xr, xm, dx, ss, s1, s2, s3, s4;
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

        xm = 0.5 * (b + a);
        xr = 0.5 * (b - a);

        ss = 0;
        for (j = 0; j < 32; j++) {
          dx = xr * x[j];
          ss += w[j] * (kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm + dx)) +
                        kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm - dx)));
        }
        ss *= xr;

        s1 = 0;
        for (j = 0; j < 32; j++) {
          dx = xr * x[j];
          s1 += w[j] * (kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm + dx)) / ss * (xm + dx) +
                        kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm - dx)) / ss * (xm - dx));
        }
        s1 *= xr;

        s2 = 0;
        for (j = 0; j < 32; j++) {
          dx = xr * x[j];
          s2 += w[j] * (kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm + dx)) / ss * pow(xm + dx - s1, 2) +
                        kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm - dx)) / ss * pow(xm - dx - s1, 2));
        }
        s2 *= xr;

        s3 = 0;
        for (j = 0; j < 32; j++) {
          dx = xr * x[j];
          s3 += w[j] * (kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm + dx)) / ss * pow(xm + dx - s1, 3) +
                        kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm - dx)) / ss * pow(xm - dx - s1, 3));
        }
        s3 *= xr / pow(s2, 1.5);

        s4 = 0;
        for (j = 0; j < 32; j++) {
          dx = xr * x[j];
          s4 += w[j] * (kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm + dx)) / ss * pow(xm + dx - s1, 4) +
                        kernelPdf(numberOfScores, score, scoreRelativeFreqDist, bandwithX, (xm - dx)) / ss * pow(xm - dx - s1, 4));
        }
        s4 *= xr / pow(s2, 2);

        moments(0) = s1;
        moments(1) = sqrt(s2);
        moments(2) = s3;
        moments(3) = s4;
      }

      /*****************************************************************************
      Pen1

      Functionality:

      Compute the Penality function PEN1 of von Davier, Holland, & Thayer 
      (2004, p. 62).

      author: Tianyou Wang 1/5/2005.

      Input:
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        fd          vector containing the frequency distribution
        hx          bandwidhx for the kernel smoohxing

      Output:
            Pen1        The panelty function PEN1 value
    *****************************************************************************/
      double pen1(const size_t& numberOfScores,
                  const Eigen::VectorXd& scores,
                  const Eigen::VectorXd& scoreRelativeFreqDist,
                  const double& bandwithX) {
        double sum = 0;

        for (size_t i = 0; i < numberOfScores; i++) {
          double pdf = kernelPdf(numberOfScores, scores, scoreRelativeFreqDist, bandwithX, scores(i));
          sum += (pdf - scoreRelativeFreqDist(i)) * (pdf - scoreRelativeFreqDist(i));
        }

        return sum;
      }

      /*****************************************************************************
      Pen2

      Functionality:

      Compute the Penality function PEN2 of von Davier, Holland, & Thayer 
      (2004, p. 63).

      author: Tianyou Wang 1/5/2005.

      Input:
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        fd          vector containing the frequency distribution
        hx          bandwidhx for the kernel smoohxing

      Output:
            Pen1        The panelty function PEN1 value
    *****************************************************************************/
      double pen2(const size_t& numberOfScores,
                  const Eigen::VectorXd& scores,
                  const Eigen::VectorXd& scoreRelativeFreqDist,
                  const double& bandwithX) {
        double sum = 0;

        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          double left = scores(scoreIndex) - 0.25;
          double right = scores(scoreIndex) + 0.25;
          double pdfPL = kernelPdfDerivative(numberOfScores, scores, scoreRelativeFreqDist, bandwithX, left);
          double pdfPR = kernelPdfDerivative(numberOfScores, scores, scoreRelativeFreqDist, bandwithX, right);
          double A = 0;
          double B = 1;

          if (pdfPL < 0) {
            A = 1;
          }

          if (pdfPR > 0) {
            B = 0;
          }

          sum += A * (1 - B);
        }

        return sum;
      }

      /*--------------------------------------------------------------------------
      KernelPdfDerivative
      
      functionality

      Computes the derivative of kernel pdf function from von Davier, Holland, & 
      Thayer (2004, p. 63). The Kernel Method of Test Equating.
      
      author: Tianyou Wang 1/5/2005.
      
      input
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        inc         Increment between consecutive raw scores for old form.
        fd          vector containing the relative frequency distribution
        hx          bandwidhx for the kernel smoohxing
        x           a particular score for which the pdf and cdf is 
                    generated            

    
      output
        The function returns the kernel pdf
    --------------------------------------------------------------------------*/
      double kernelPdfDerivative(const size_t& numberOfScores,
                                 const Eigen::VectorXd& scores,
                                 const Eigen::VectorXd& scoreRelativeFreqDist,
                                 const double& bandwithX,
                                 const double& scoreValue) {
        double mu = scores.dot(scoreRelativeFreqDist);
        double sigma = (scores.cwiseProduct(scores)).dot(scoreRelativeFreqDist) - std::pow(mu, 2);
        double ax = sqrt(sigma / (sigma + std::pow(bandwithX, 2)));

        double pdfDer = 0;

        for (size_t scoreIndex = 0; scoreIndex < numberOfScores; scoreIndex++) {
          double zscore = (scoreValue - ax * scores(scoreIndex) - (1 - ax) * mu) / ax / bandwithX;
          double temp1 = stdNormalPdf(zscore);
          pdfDer += scoreRelativeFreqDist(scoreIndex) * temp1 / ax / ax / bandwithX / bandwithX * zscore;
        }

        pdfDer *= -1.0;

        return pdfDer;
      }

      /*****************************************************************************
      Pen

      Functionality:

      Compute the combined Penality function of von Davier, Holland, & Thayer 
      (2004, p. 63).

      author: Tianyou Wang 1/5/2005.

      Input:
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        fd          vector containing the frequency distribution
        hx          bandwidhx for the kernel smoohxing
        K           The weight for Pen2

      Output:
            return the combined panelty function value
    *****************************************************************************/
      double Pen(const size_t& numberOfScores,
                 const Eigen::VectorXd& scores,
                 const Eigen::VectorXd& scoreRelativeFreqDist,
                 const double& bandwithX,
                 const double& pen2Weight) {
        double PEN1 = Pen1(numberOfScores, scores, scoreRelativeFreqDist, bandwithX);
        double PEN2 = Pen2(numberOfScores, scores, scoreRelativeFreqDist, bandwithX);
        double PEN = PEN1 + pen2Weight * PEN2;

        return PEN;
      }

      /*****************************************************************************
      Optimalh

      Functionality:
            Find the optimal bandwidhx parameter for kernel continuization hx based on 
            von Davier, Holland, & Thayer (2004, p. 63). The Kernel Method of Test Equating.

      author: Tianyou Wang 1/5/2005.

      Input:
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        inc         Increment between consecutive raw scores for old form.
        fd          vector containing the frequency distribution
        K           The weight for PEN2
      
      Output:
            return the optimal hx
    *****************************************************************************/
      double optimalBandwithXParameter(const size_t& numberOfScores,
                                       const Eigen::VectorXd& scores,
                                       const Eigen::VectorXd& scoreRelativeFreqDist,
                                       const double& pen2Weight) {
        double eps = .0001;
        double hxl = .0001;
        double hxu = 3;
        double hxlplus = .0002;
        double hxuminxs = 2.9;
        double hxb = 1.1;
        double hxx;
        double optimhx;

        double pl = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxl, pen2Weight);
        double pu = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxu, pen2Weight);
        double plplus = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxlplus, pen2Weight);
        double puminxs = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxuminxs, pen2Weight);
        double pb = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxb, pen2Weight);

        if (pl < pb && pb < pu && pl < plplus) {
          optimhx = hxl;
        } else if (pu < pb && pb < pl && pu < puminxs) {
          optimhx = hxu;
        } else {
          for (size_t iter = 1; iter <= 200; iter++) {
            hxb = .38197 * hxu + .61803 * hxl;
            hxx = .61803 * hxu + .38197 * hxl;
            double absdif = std::abs(hxu - hxl);
            if (absdif < eps) {
              break;
            }
            pb = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxb, pen2Weight);
            double px = Pen(numberOfScores, scores, scoreRelativeFreqDist, hxx, pen2Weight);

            if (px <= pb) {
              hxl = hxb;
            }

            if (px > pb) {
              hxu = hxx;
            }
          }

          optimhx = 0.5 * (hxb + hxx);
        }

        return optimhx;
      }

      /*--------------------------------------------------------------------------
      KernelInverseCdf
      
      functionality:

      Computes the inverse of the cdf in von Davier, Holland, & Thayer 
      (2004). The Kernel Method of Test Equating.
      
      author: Tianyou Wang 1/5/2005.
      
      input:
        ncat    Number of discrete score categories
        scores      vector containing the discrete scores
        fd          vector containing the frequency distribution
        h           bandwidhx for the kernel smoohxing
        cdf         a particular cdf for which the score is found
    
      output:
        The function returns the inverse of cdf
    --------------------------------------------------------------------------*/
      double kernelInverseCdf(const size_t& ncat,
                              const Eigen::VectorXd& scores,
                              const Eigen::VectorXd& fd,
                              const double& h,
                              const double& cdf) {
        double absdif;
        double eps = .000001;
        double half;

        double lb = scores(0) - 5.0;
        double ub = scores(ncat - 1) + 5.0;
        double cdfl = kernelCdf(ncat, scores, fd, h, lb);
        double cdfu = kernelCdf(ncat, scores, fd, h, ub);

        if (cdf < cdfl) {
          half = scores(0);
        } else if (cdf > cdfu) {
          half = scores(ncat - 1);
        } else {
          for (size_t iter = 1; iter <= 200; iter++) {
            half = 0.5 * (lb + ub);
            double cdfhalf = kernelCdf(ncat, scores, fd, h, half);
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
        KernelEquate
        
        functionality:

        Computes kernel equating function based on continuized cdf in von Davier, 
        Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
        
        author: Tianyou Wang 1/5/2005.
        
        input:
          ncatx   Number of discrete score categories for the new form
          scoresx     vector containing the discrete scores for the new form
          fdx         vector containing the relative frequency distribution  for the new form
          hx          bandwidhx for the kernel smoohxing  for the new form
          ncaty   Number of discrete score categories for the old form
          scoresy     vector containing the discrete scores for the old form
          fdy         vector containing the relative frequency distribution  for the old form
          hy          bandwidhx for the kernel smoohxing  for the old form

      
        output:
          Equatedx   a vector containing the equated score 
      --------------------------------------------------------------------------*/
      void kernelEquate(const size_t& numberOfScoresX,
                        const Eigen::VectorXd& scoresX,
                        const Eigen::VectorXd& scoreRelativeFreqDistX,
                        const double& bandwithX,
                        const size_t& numberOfScoresY,
                        const Eigen::VectorXd& scoresY,
                        const Eigen::VectorXd& scoreRelativeFreqDistY,
                        const double& bandwithY,
                        Eigen::VectorXd& equatedX) {
        equatedX.resize(numberOfScoresX);

        for (size_t scoreIndex = 0; scoreIndex < numberOfScoresX; scoreIndex++) {
          double cdfx = kernelCdf(numberOfScoresX, scoresX, scoreRelativeFreqDistX, bandwithX, scoresX(scoreIndex));
          equatedX(scoreIndex) = kernelInverseCdf(numberOfScoresY, scoresY, scoreRelativeFreqDistY, bandwithY, cdfx);
        }
      }

      /*--------------------------------------------------------------------------
      ComputeCmatrix
      
      functionality: 

      Computes the C matrix in von Davier, Holland, & Thayer 
      (2004, Equation 3.10). The Kernel Method of Test Equating.
      call the numerical recipe function qrdcmp to do QR decomposition.
      
      author: Tianyou Wang 3/11/2005.
      
      input:
        ncat    Number of discrete score categories
        degree      The highest polynomial degree of the log-linear model
        np          sample size
        B           B matrix (design matrix)
        fd          vector containing the frequency distribution
    
      output:
        Cr          C matrix in von Davier et. al. Equation 3.10
    --------------------------------------------------------------------------*/
      void computeCmatrix(const size_t& numberOfScores,
                          const size_t& numberOfDegreesSmoothing,
                          const size_t& numberOfExaminees,
                          const Eigen::VectorXd& designMatrix,
                          const Eigen::VectorXd& scoreRelativeFreqDist,
                          Eigen::MatrixXd& cMatrix) {
        // int i, j, k, l;
        // double *D, **A, **a, **qt, **q, **r, rrij;
        // double *c, *d;
        // double con;

        Eigen::VectorXd D = Eigen::VectorXd::Zero(numberOfScores);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(numberOfScores, numberOfScores);
        Eigen::MatrixXd a = Eigen::MatrixXd::Zero(1, numberOfScores + 1);
        Eigen::MatrixXd qt = Eigen::MatrixXd::Zero(numberOfScores, numberOfScores);
        Eigen::MatrixXd q = Eigen::MatrixXd::Zero(numberOfScores, numberOfScores);
        Eigen::MatrixXd r = Eigen::MatrixXd::Zero(numberOfScores, numberOfScores);
        Eigen::VectorXd c = Eigen::VectorXd::Zero(numberOfScores + 1);
        Eigen::VectorXd d = Eigen::VectorXd::Zero(numberOfScores + 1);

        cMatrix.setZero(numberOfScores, numberOfDegreesSmoothing);

        D = scoreRelativeFreqDist.cwiseSqrt();

        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfScores; j++) {
            double rrij = std::sqrt(scoreRelativeFreqDist(i)) * scoreRelativeFreqDist(j);
            if (i == j) {
              A(i, j) = D(i) - rrij;
            } else {
              A(i, j) = -1.0 * rrij;
            }
          }
        }

        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfScores; j++) {
            for (size_t k = 0; k < numberOfScores; k++) {
              a(i + 1, j + 1) += A(i, k) * designMatrix(j, k);
            }
          }
        }

        Eigen::HouseholderQR qrDecomp = a.householderQr();
        Eigen::MatrixXd qrMatrix = qrDecomp.matrixQR();

        // er_qrdcmp(a, ncat, ncat, c, d);

        /* compute the Q and R matrices */
        for (size_t k = 0; k < numberOfScores; k++) {
          for (size_t l = 0; l < numberOfScores; l++) {
            if (l > k) {
              r(k, l) = qrMatrix(k + 1, l + 1);
              q(k, l) = 0.0;
            } else if (l < k) {
              r(k, l) = 0.0;
              q(k, l) = 0.0;
            } else {
              r(k, l) = d(k + 1);
              q(k, l) = 1.0;
            }
          }
        }

        for (size_t i = numberOfScores - 2; i >= 0; i--) {
          double con = 0.0;

          for (size_t k = i; k < numberOfScores; k++) {
            con += std::pow(a(k + 1, i + 1), 2);
          }

          con /= 2.0;

          for (size_t k = i; k < numberOfScores; k++) {
            for (size_t l = i; l < numberOfDegreesSmoothing; l++) {
              qt(k, l) = 0.0;
              for (size_t j = i; j < numberOfScores; j++) {
                qt(k, l) += q(j, l) * a(k + 1, i + 1) * a(j + 1, i + 1) / con;
              }
            }
          }

          for (size_t k = i; k < numberOfScores; k++) {
            for (size_t l = i; l < numberOfDegreesSmoothing; l++) {
              q(k, l) -= qt(k, l);
            }
          }
        }

        /* compute the Cr matrix */
        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfDegreesSmoothing; j++) {
            cMatrix(i, j) += (1.0 / sqrt(static_cast<double>(numberOfExaminees))) * D(i) * q(i, j);
          }
        }
      }

      /*--------------------------------------------------------------------------
      ComputeCmatrixGen
      
      functionality: 

      Computes the C matrix in von Davier, Holland, & Thayer 
      (2004, Equation 3.10) for a general case where non-square matrixes 
      are involved. The Kernel Method of Test Equating.
      call the numerical recipe function qrdcmp to do QR decomposition.
      
      author: Tianyou Wang 3/11/2005.
      
      input:
        ncat    Number of discrete score categories
        degree      The highest polynomial degree of the log-linear model
        np          sample size
        B           B matrix (design matrix)
        fd          vector containing the frequency distribution
    
      output:
        Cr          C matrix in von Davier et. al. Equation 3.10
    --------------------------------------------------------------------------*/
      void computeCmatrixGen(const size_t& numberOfScores,
                             const size_t& numberOfDegreesSmoothing,
                             const size_t& numberOfExaminees,
                             const Eigen::MatrixXd& designMatrix,
                             const Eigen::VectorXd& scoreRelativeFreqDist,
                             Eigen::MatrixXd& cMatrix) {
        // int i, j, k, l;
        // double *D, **A, **a, **b, **qt, **q, **r, rrij, sum;
        // double *c, *d;
        // double con;
        // FILE* outf;
        // char outfname[] = "Bmatrix.out";
        // static int ii = 0;

        Eigen::VectorXd D = Eigen::VectorXd::Zero(numberOfScores);                                   // dvector(0, ncat - 1);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(numberOfScores, numberOfScores);                   // dmatrix(0, ncat - 1, 0, ncat - 1);
        Eigen::MatrixXd a = Eigen::MatrixXd::Zero(numberOfScores + 1, numberOfDegreesSmoothing + 1); // dmatrix(1, ncat, 1, degree);
        Eigen::MatrixXd b = Eigen::MatrixXd::Zero(numberOfScores + 1, numberOfDegreesSmoothing + 1); // dmatrix(1, ncat, 1, degree);
        Eigen::MatrixXd qt = Eigen::MatrixXd::Zero(numberOfScores, numberOfDegreesSmoothing);        // dmatrix(0, ncat - 1, 0, degree - 1);
        Eigen::MatrixXd q = Eigen::MatrixXd::Zero(numberOfScores, numberOfDegreesSmoothing);         // dmatrix(0, ncat - 1, 0, degree - 1);
        Eigen::MatrixXd r = Eigen::MatrixXd::Zero(numberOfScores, numberOfDegreesSmoothing);         // dmatrix(0, ncat - 1, 0, degree - 1);
        Eigen::VectorXd c = Eigen::VectorXd::Zero(numberOfScores + 1);                               // dvector(1, ncat);
        Eigen::VectorXd d = Eigen::VectorXd::Zero(numberOfScores + 1);                               // dvector(1, ncat);

        D = scoreRelativeFreqDist.cwiseSqrt();

        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfScores; j++) {
            double rrij = std::sqrt(scoreRelativeFreqDist(i)) * scoreRelativeFreqDist(j);
            if (i == j) {
              A(i, j) = D(i) - rrij;
            } else {
              A(i, j) = -1.0 * rrij;
            }
          }
        }

        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfDegreesSmoothing; j++) {
            for (size_t k = 0; k < numberOfScores; k++) {
              a(i + 1, j + 1) += A(i, k) * designMatrix(j, k);
            }

            b(i + 1, j + 1) = a(i + 1, j + 1);
          }
        }

        Eigen::HouseholderQR qrDecomp = a.householderQr();
        Eigen::MatrixXd qrMatrix = qrDecomp.matrixQR();

        /* compute the Q and R matrices */
        for (size_t k = 0; k < numberOfScores; k++) {
          for (size_t l = 0; l < numberOfDegreesSmoothing; l++) {
            if (l > k) {
              r(k, l) = qrMatrix(k + 1, l + 1);
              q(k, l) = 0.0;
            } else if (l < k) {
              r(k, l) = 0.0;
              q(k, l) = 0.0;
            } else {
              r(k, l) = d(k + 1);
              q(k, l) = 1.0;
            }
          }
        }

        /* in the next line, original i=ncat-2, now i=degree-1. i=degree-2 does not work */
        for (size_t i = numberOfDegreesSmoothing - 1; i >= 0; i--) {
          double con = 0.0;

          for (size_t k = i; k < numberOfScores; k++) {
            con += std::pow(a(k + 1, i + 1), 2);
          }

          con /= 2.0;

          for (size_t k = i; k < numberOfScores; k++) {
            for (size_t l = i; l < numberOfDegreesSmoothing; l++) {
              qt(k, l) = 0.0;

              for (size_t j = i; j < numberOfScores; j++) {
                qt(k, l) += q(j, l) * a(k + 1, i + 1) * a(j + 1, i + 1) / con;
              }
            }
          }

          for (size_t k = i; k < numberOfScores; k++) {
            for (size_t l = i; l < numberOfDegreesSmoothing; l++) {
              q(k, l) *= -1.0;
            }
          }
        }

        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfDegreesSmoothing; j++) {
            a(i + 1, j + 1) = 0;

            for (size_t k = 0; k < numberOfDegreesSmoothing; k++) {
              a(i + 1, j + 1) += q(i, k) * r(k, j);
            }
          }
        }

        for (size_t i = 0; i < numberOfDegreesSmoothing; i++) {
          for (size_t j = i; j < numberOfDegreesSmoothing; j++) {
            double sum = 0;

            for (size_t k = 0; k < numberOfScores; k++) {
              sum += q(k, i) * q(k, j);
            }
          }
        }

        /* compute the Cr matrix */
        for (size_t i = 0; i < numberOfScores; i++) {
          for (size_t j = 0; j < numberOfDegreesSmoothing; j++) {
            cMatrix(i, j) += (1.0 / std::sqrt(static_cast<double>(numberOfExaminees))) * D(i) * q(i, j);
          }
        }
      }

      /*--------------------------------------------------------------------------
      PartialFPartialr
      
      functionality: 

      Computes the partial derivative of F to r in von Davier, Holland, & Thayer 
      (2004, Equation 5.21). The Kernel Method of Test Equating.
      call the numerical recipe function qrdcmp to do QR decomposition.
      
      author: Tianyou Wang 3/15/2005.
      
      input:
        ncat    Number of discrete score categories
        degree      The highest polynomial degree of the log-linear model
        np          sample size
        scores      vector containing the discrete scores
        fd          vector containing the frequency distribution
        hx          kernel continuizing parameter
        x           X score values at which the partial derivatives are evaluated.
                    x can be non-integer values.
    
      output:
        Fr          vector containing partial derivatives of F with respect to r
                    in von Davier et. al. Equation 5.12
    --------------------------------------------------------------------------*/
      void partialFPartialr(const size_t& numberOfScores,
                            const Eigen::VectorXd& scores,
                            const Eigen::VectorXd& scoreRelativeFreqDist,
                            const double& bandwithX,
                            Eigen::VectorXd& Fr,
                            const double& x) {
        double mu = scores.dot(scoreRelativeFreqDist);
        double sigma = (scores.cwiseProduct(scores)).dot(scoreRelativeFreqDist) - std::pow(mu, 2);

        double ax = std::sqrt(sigma / (sigma + std::pow(bandwithX, 2)));

        Fr.resize(numberOfScores);

        for (size_t j = 0; j < numberOfScores; j++) {
          double Rjx = (x - ax * scores(j) - (1 - ax) * mu) / ax / bandwithX;
          double temp2 = stdNormalCdf(Rjx);
          double zjx = (scores(j) - mu) / sqrt(sigma);
          double Mjx = .5 * (x - mu) * (1 - ax * ax) * zjx * zjx + (1 - ax) * scores[j];
          double pdf = kernelPdf(numberOfScores, scores, scoreRelativeFreqDist, bandwithX, x);
          Fr(j) = temp2 - Mjx * pdf;
        }
      }

      /*--------------------------------------------------------------------------
      FrCrSqNorm
      
      functionality: 

      Computes the sqaured norm of Fr multiplied by C matrix in von Davier, Holland, & Thayer 
      (2004, Equation 7.5). The Kernel Method of Test Equating.
      call the numerical recipe function qrdcmp to do QR decomposition.
      
      author: Tianyou Wang 3/15/2005.
      
      input:
        ncat    Number of discrete score categories
        degree      The highest polynomial degree of the log-linear model
        np          sample size
        scores      vector containing the discrete scores
        fd          vector containing the frequency distribution
        Fr          vector containing partial derivatives of F with respect to r
    
      output:
        return the square of the norm in von Davier, Holland, & Thayer 
          (2004, Equation 7.5)
    --------------------------------------------------------------------------*/
      double frCrSqNorm(const size_t& numberOfScores,
                        const size_t& numberofDegreesSmoothing,
                        const Eigen::VectorXd& Fr,
                        const Eigen::MatrixXd& cMatrix) {
        double norm = 0;

        for (size_t j = 0; j < numberofDegreesSmoothing; j++) {
          double temp = 0;
          for (size_t i = 0; i < numberOfScores; i++) {
            temp += Fr(i) * cMatrix(i, j);
          }

          norm += temp * temp;
        }

        return norm;
      }

      /*--------------------------------------------------------------------------
      vPMN
      
      functionality:

      Computes vectorized P and M and N matrix from the bivariate distribution 
      defined in von Davier, 	Holland, & Thayer (2004, Equations 2.8, 2.9, 2.10). 
      The Kernel Method of Test Equating.
      
      author: Tianyou Wang 1/16/2005.
      
      input:
        bdist       bivariate fitted distribution
        ncatx       Number of discrete score categories for the new form
        ncaty       Number of discrete score categories for the old form

    
      output:
          vP          vectorized P 
        M           M matrix
        N           N matrix
    --------------------------------------------------------------------------*/
      void vPMN(const size_t& numberOfScoresX,
                const size_t& numberOfScoresY,
                const Eigen::MatrixXd& bivariateFittedDistribution,
                Eigen::VectorXd& vP,
                Eigen::MatrixXd& M,
                Eigen::MatrixXd& N) {
        size_t totalNumberOfScores = numberOfScoresX * numberOfScoresY;

        for (size_t i = 0; i < numberOfScoresY; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            vP(i * numberOfScoresX + j) = bivariateFittedDistribution(j, i);
          }
        }

        M.setZero();

        for (size_t i = 0; i < numberOfScoresX; i++) {
          for (size_t j = 0; j < numberOfScoresY; j++) {
            M(i, j * numberOfScoresX + i) = 1;
          }
        }

        N.setZero();

        for (size_t i = 0; i < numberOfScoresY; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            N(i, i * numberOfScoresX + j) = 1;
          }
        }
      }

      /*--------------------------------------------------------------------------
      vPP
      
      functionality:

      Computes vectorized P and M and N matrix from the bivariate distribution 
      defined in von Davier, 	Holland, & Thayer (2004, Equations 2.8, 2.9, 2.10). 
      The Kernel Method of Test Equating.
      
      author: Tianyou Wang 1/16/2005.
      
      input:
        bdist       bivariate fitted distribution
        ncatx       Number of discrete score categories for the new form
        ncaty       Number of discrete score categories for the old form

    
      output:
          vPP         vectorized P 
    --------------------------------------------------------------------------*/
      void vPT(const size_t& numberOfScoresX,
               const size_t& numberOfScoresY,
               const Eigen::MatrixXd& bivariateFittedDistribution,
               Eigen::VectorXd& vPP) {
        for (size_t i = 0; i < numberOfScoresX; i++) {
          for (size_t j = 0; j < numberOfScoresY; j++) {
            vPP(i * numberOfScoresY + j) = bivariateFittedDistribution(i, j);
          }
        }
      }

      /*-----------------------------------------------------------------------------------------
      KernelEquateSEERG
      
      functionality: 

      Computes the standard error of equating (SEE) for the Random Groups Design
      (called EG design in von Davier, Holland, & Thayer 
      (2004, Equation 7.5)). The Kernel Method of Test Equating.
      call the numerical recipe function qrdcmp to do QR decomposition.
      
      author: Tianyou Wang 3/15/2005.
      
      input:
        ncatx        Number of discrete score categories for the new form
        degreex      The highest polynomial degree of the log-linear model for the new form
        npx          sample size for the new form
        scoresx      vector containing the discrete scores for the new form
        fdx          vector containing the relative frequency distribution for the new form  
        hx           kernel continuizing parameter for the new form
        Equatedx     the equated X scores
        ncaty        Number of discrete score categories for the old form
        degreey      The highest polynomial degree of the log-linear model for the old form
        npy          sample size for the old form
        scoresy      vector containing the discrete scores for the old form
        fdy          vector containing the relative frequency distribution for the old form  
        hy           kernel continuizing parameter for the old form
    
      output:
          Equatedx     vector containing equated scores
        SEE          vector containing the standard error of equating 
    -----------------------------------------------------------------------------------------*/
      void kernelEquateSEERG(const size_t& numberOfScoresX,
                             const size_t& numberOfDegreesSmoothingX,
                             const size_t& numberOfExamineesX,
                             const Eigen::VectorXd& scoresX,
                             const Eigen::VectorXd& scoreRelativeFreqDistX,
                             const size_t& numberOfScoresY,
                             const size_t& numberOfDegreesSmoothingY,
                             const size_t& numberOfExamineesY,
                             const Eigen::VectorXd& scoresY,
                             const Eigen::VectorXd& scoreRelativeFreqDistY,
                             Eigen::VectorXd& equatedX,
                             Eigen::VectorXd& see) {
        size_t nparax = numberOfDegreesSmoothingX;
        size_t nparay = numberOfDegreesSmoothingY;
        Eigen::MatrixXd Crx(numberOfScoresX, numberOfScoresX);
        Eigen::MatrixXd Cry(numberOfScoresY, numberOfScoresY);
        Eigen::MatrixXd Bx(nparax, numberOfScoresX);
        Eigen::MatrixXd By(nparay, numberOfScoresY);
        Eigen::VectorXd Fr(numberOfScoresX);
        Eigen::VectorXd Gs(numberOfScoresY);

        for (size_t i = 0; i < nparax; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            Bx(i, j) = 2.0 * std::pow(scoresX(j), static_cast<double>(i + 1));
          }
        }

        for (size_t i = 0; i < nparay; i++) {
          for (size_t j = 0; j < numberOfScoresY; j++) {
            By(i, j) = 2.0 * std::pow(scoresY(j), static_cast<double>(i + 1));
          }
        }

        computeCmatrixGen(numberOfScoresX, nparax, numberOfExamineesX, Bx, scoreRelativeFreqDistX, Crx);
        computeCmatrixGen(numberOfScoresY, nparay, numberOfExamineesY, By, scoreRelativeFreqDistY, Cry);

        double bandwithX = optimalBandwithXParameter(numberOfScoresX, scoresX, scoreRelativeFreqDistX, 1);
        double bandwithY = optimalBandwithXParameter(numberOfScoresY, scoresY, scoreRelativeFreqDistY, 1);

        kernelEquate(numberOfScoresX, scoresX, scoreRelativeFreqDistX, bandwithX, numberOfScoresY, scoresY, scoreRelativeFreqDistY, bandwithY, equatedX);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          partialFPartialr(numberOfScoresX, scoresX, scoreRelativeFreqDistX, bandwithX, Fr, scoresX(i));
          partialFPartialr(numberOfScoresY, scoresY, scoreRelativeFreqDistY, bandwithY, Gs, equatedX(i));
          double Gp = kernelPdf(numberOfScoresY, scoresY, scoreRelativeFreqDistY, bandwithY, equatedX(i));
          double SqNormx = FrCrSqNorm(numberOfScoresX, nparax, Fr, Crx);
          double SqNormy = FrCrSqNorm(numberOfScoresY, nparay, Gs, Cry);
          see(i) = std::sqrt(SqNormx + SqNormy) / Gp;
        }
      }

      /*--------------------------------------------------------------------------
      KernelEquateSEESG
      
      functionality:

      Computes kernel equating function and SEE for the Single Group (SG) design
      based on continuized cdf in von Davier, 
      Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
      
      author: Tianyou Wang 3/26/2005.
      
      input:
        bivar       information about bivariate distribution of X and Y
    
      output:
        Equatedx    a vector containing the equated score 
        SEE         a vector containing standard error of equating 
    --------------------------------------------------------------------------*/
      void kernelEquateSEESG(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingUV,
                             Eigen::VectorXd& equatedX,
                             Eigen::VectorXd& see) {
        size_t numberOfCrossProductMoments = bivariateLogLinearSmoothingUV.numberOfCrossProductMoments;
        size_t numberOfDegreesOfSmoothingU = bivariateLogLinearSmoothingUV.numberOfDegreesOfSmoothingU;
        size_t numberOfDegreesOfSmoothingV = bivariateLogLinearSmoothingUV.numberOfDegreesOfSmoothingV;

        Eigen::MatrixXd cpm = bivariateLogLinearSmoothingUV.crossProductMoments(Eigen::seq(0, 1),
                                                                                Eigen::seq(0, numberOfCrossProductMoments - 1))
                                  .cast<double>();

        size_t numberOfScoresX = bivariateLogLinearSmoothingUV.numberOfScoresU;
        size_t numberOfScoresY = bivariateLogLinearSmoothingUV.numberOfScoresV;
        size_t npara = bivariateLogLinearSmoothingUV.numberOfDegreesOfSmoothingU +
                       bivariateLogLinearSmoothingUV.numberOfDegreesOfSmoothingV +
                       bivariateLogLinearSmoothingUV.numberOfCrossProductMoments;

        Eigen::VectorXd interx = cpm.row(0);
        Eigen::VectorXd intery = cpm.row(1);
        size_t numberOfExaminees = bivariateLogLinearSmoothingUV.numberOfExaminees;

        Eigen::VectorXd scoresx(numberOfScoresX);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          scoresx(i) = bivariateLogLinearSmoothingUV.minimumRawScoreX +
                       static_cast<double>(i) * bivariateLogLinearSmoothingUV.scoreIncrementX;
        }

        Eigen::VectorXd scoresy(numberOfScoresY);
        for (size_t i = 0; i < numberOfScoresY; i++) {
          scoresy(i) = bivariateLogLinearSmoothingUV.minimumRawScoreV +
                       static_cast<double>(i) * bivariateLogLinearSmoothingUV.scoreIncrementV;
        }

        size_t totalNumberOfScores = numberOfScoresX * numberOfScoresY;

        Eigen::MatrixXd Cr(totalNumberOfScores, npara);
        Eigen::MatrixXd B(npara, totalNumberOfScores);
        Eigen::MatrixXd U(numberOfScoresX, npara);
        Eigen::MatrixXd V(numberOfScoresY, npara);
        Eigen::VectorXd Fr(numberOfScoresX);
        Eigen::VectorXd Gs(numberOfScoresY);
        Eigen::VectorXd FrU(npara);
        Eigen::VectorXd GsV(npara);
        Eigen::MatrixXd fitbdist(numberOfScoresX, numberOfScoresY);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          for (size_t j = 0; j < numberOfScoresY; j++) {
            fitbdist(i, j) = bivariateLogLinearSmoothingUV.fittedBivariateFreqDist(i, j) /
                             static_cast<double>(numberOfExaminees);
          }
        }

        B.setZero();

        /* The following code set a natrual design matrix corresponding to a natrual basis of 
					polynomials. The B matrix is a transpose of the design matrix in loglinear model
				  (without the first column)
					First, assign natrual design matrix first for rows corresponding to x */

        for (size_t i = 0; i < numberOfDegreesOfSmoothingU; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            for (size_t k = 0; k < numberOfScoresY; k++) {
              B(i, k * numberOfScoresY + j) = std::pow(static_cast<double>(j), i + 1);
            }
          }
        }

        /* then assign values to rows corresponding to y */
        for (size_t i = numberOfDegreesOfSmoothingU; i < numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            for (size_t k = 0; k < numberOfScoresY; k++) {
              B(i, k * numberOfScoresX + j) = std::pow(static_cast<double>(k), (i - numberOfDegreesOfSmoothingU + 1));
            }
          }
        }

        /* assign value to the last rows corresponding to the interaction terms */
        for (size_t i = 0; i < numberOfCrossProductMoments; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            for (size_t k = 0; k < numberOfScoresY; k++) {
              B(i + numberOfDegreesOfSmoothingU + numberOfDegreesOfSmoothingV, k * numberOfScoresX + j) =
                  std::pow(static_cast<double>(j), interx(i)) *
                  std::pow(static_cast<double>(k), intery(i));
            }
          }
        }

        Eigen::VectorXd vP(totalNumberOfScores);
        Eigen::VectorXd vPP(totalNumberOfScores);
        Eigen::MatrixXd M(numberOfScoresX, numberOfScoresY);
        Eigen::MatrixXd N(numberOfScoresY, totalNumberOfScores);
        Eigen::VectorXd r(numberOfScoresX);
        Eigen::VectorXd s(numberOfScoresY);

        vPMN(numberOfScoresX, numberOfScoresY, fitbdist, vP, M, N);
        vPT(numberOfScoresX, numberOfScoresY, fitbdist, vPP);

        r = M * vP;
        s = N * vP;

        computeCmatrixGen(totalNumberOfScores, npara, numberOfExaminees, B, vP, Cr);

        U = M * Cr;
        V = N * Cr;

        double bandwithX = optimalBandwithXParameter(numberOfScoresX, scoresx, r, 1);
        double bandwithY = optimalBandwithXParameter(numberOfScoresY, scoresy, s, 1);

        kernelEquate(numberOfScoresX, scoresx, r, bandwithX, numberOfScoresY, scoresy, s, bandwithY, equatedX);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          partialFPartialr(numberOfScoresX, scoresx, r, bandwithX, Fr, scoresx(i));
          partialFPartialr(numberOfScoresY, scoresy, s, bandwithY, Gs, equatedX(i));
          double Gp = kernelPdf(numberOfScoresY, scoresy, s, bandwithY, equatedX(i));

          FrU = Fr * U;
          GsV = Gs * V;

          for (size_t j = 0; j < npara; j++) {
            FrU[j] -= GsV[j];
          }

          see(i) = FrU.norm() / Gp;
        }
      }

      /*--------------------------------------------------------------------------
			KernelEquateSG
			
			functionality:

			Computes kernel equating function and SEE for the Single Group (SG) design
			based on continuized cdf in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			author: Tianyou Wang 3/26/2005.
			
			input:
				bivar       information about bivariate distribution of X and Y
		
			output:
				Equatedx    a vector containing the equated score 
		--------------------------------------------------------------------------*/
      void kernelEquateSG(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothingUV,
                          Eigen::VectorXd& equatedX) {
        size_t numberOfScoresX = bivariateLogLinearSmoothingUV.numberOfScoresX;
        size_t numberOfScoresY = bivariateLogLinearSmoothingUV.numberOfScoresV;
        size_t numberOfExaminees = bivariateLogLinearSmoothingUV.numberOfExaminees;

        Eigen::VectorXd scoresX(numberOfScoresX);
        Eigen::VectorXd scoresY(numberOfScoresY);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          scoresX(i) = bivariateLogLinearSmoothingUV.minimumRawScoreX + i * bivariateLogLinearSmoothingUV.scoreIncrementX;
        }

        for (size_t i = 0; i < numberOfScoresY; i++) {
          scoresY(i) = bivariateLogLinearSmoothingUV.minimumRawScoreV + i * bivariateLogLinearSmoothingUV.scoreIncrementV;
        }

        size_t totalNumberOfScores = numberOfScoresX * numberOfScoresY;

        Eigen::MatrixXd fitbdist(numberOfScoresX, numberOfScoresY);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          for (size_t j = 0; j < numberOfScoresY; j++) {
            fitbdist(i, j) = bivariateLogLinearSmoothingUV.fittedBivariateFreqDist(i, j) /
                             static_cast<double>(numberOfExaminees);
          }
        }

        Eigen::VectorXd vP(totalNumberOfScores);
        Eigen::MatrixXd M(numberOfScoresX, totalNumberOfScores);
        Eigen::MatrixXd N(numberOfScoresY, totalNumberOfScores);
        Eigen::VectorXd r(numberOfScoresX);
        Eigen::VectorXd s(numberOfScoresY);

        vPMN(numberOfScoresX, numberOfScoresY, fitbdist, vP, M, N);

        r = M * vP;
        s = N * vP;

        double bandwithX = optimalBandwithXParameter(numberOfScoresX, scoresX, r, 1);
        double bandwithY = optimalBandwithXParameter(numberOfScoresY, scoresY, s, 1);

        kernelEquate(numberOfScoresX, scoresX, r, bandwithX, numberOfScoresY, scoresY, s, bandwithY, equatedX);
      }

      /*--------------------------------------------------------------------------
			KernelEquateSEECB
			
			functionality:

			Computes kernel equating function and SEE for the counter balance (CB) 
			design based on continuized cdf in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			author: Tianyou Wang 3/27/2005.
			
			input:
				bivar1      information about bivariate distribution for form X
				bivar2      information about bivariate distribution for form Y
				wts         vector containing weights for the form taken first.
										wts[0] is for form X, wts[1] is for form Y

		
			output:
				Equatedx    a vector containing the equated score 
				SEE         a vector containing standard error of equating 
		--------------------------------------------------------------------------*/
      void kernelEquateSEECB(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing1,
                             const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing2,
                             const Eigen::VectorXd& wts,
                             Eigen::VectorXd& equatedX,
                             Eigen::VectorXd& see) {
        size_t numberOfCrossProductMoments = bivariateLogLinearSmoothing1.numberOfCrossProductMoments;
        size_t numberOfDegreesSmoothingX = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingU;
        size_t numberOfDegreesSmoothingY = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingV;
        Eigen::MatrixXd crossProductMoments1(2, numberOfCrossProductMoments + 1);

        for (size_t i = 0; i < numberOfCrossProductMoments; i++) {
          for (size_t j = 0; j < 2; j++) {
            crossProductMoments1(j, i) = static_cast<double>(bivariateLogLinearSmoothing1.crossProductMoments(i, j));
          }
        }

        size_t numberOfScoresX = bivariateLogLinearSmoothing1.numberOfScoresX;
        size_t numberOfScoresY = bivariateLogLinearSmoothing1.numberOfScoresV;
        size_t npara = numberOfDegreesSmoothingX + numberOfDegreesSmoothingY + numberOfCrossProductMoments;

        Eigen::VectorXd interx = crossProductMoments1.row(0);
        Eigen::VectorXd intery = crossProductMoments1.row(1);

        size_t numberOfExaminees1 = bivariateLogLinearSmoothing1.numberOfExaminees;
        size_t numberOfExaminees2 = bivariateLogLinearSmoothing2.numberOfExaminees;

        Eigen::VectorXd scoresX(numberOfScoresX);
        Eigen::VectorXd scoresY(numberOfScoresY);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          scoresX(i) = bivariateLogLinearSmoothing1.minimumRawScoreX +
                       static_cast<double>(i) * bivariateLogLinearSmoothing1.scoreIncrementX;
        }

        for (size_t i = 0; i < numberOfScoresY; i++) {
          scoresY(i) = bivariateLogLinearSmoothing1.minimumRawScoreV +
                       static_cast<double>(i) * bivariateLogLinearSmoothing1.scoreIncrementV;
        }

        size_t totalNumberOfScores = numberOfScoresX * numberOfScoresY;

        Eigen::VectorXd vP12(totalNumberOfScores);                   // = dvector(0, ncat - 1);
        Eigen::VectorXd vP21(totalNumberOfScores);                   // = dvector(0, ncat - 1);
        Eigen::VectorXd Fr(numberOfScoresX);                         // = dvector(0, ncatx - 1);
        Eigen::VectorXd r(numberOfScoresX);                          // = dvector(0, ncatx - 1);
        Eigen::VectorXd r1(numberOfScoresX);                         // = dvector(0, ncatx - 1);
        Eigen::VectorXd r2(numberOfScoresX);                         // = dvector(0, ncatx - 1);
        Eigen::VectorXd Gs(numberOfScoresY);                         // = dvector(0, ncaty - 1);
        Eigen::VectorXd s(numberOfScoresY);                          // = dvector(0, ncaty - 1);
        Eigen::VectorXd s1(numberOfScoresY);                         // = dvector(0, ncaty - 1);
        Eigen::VectorXd s2(numberOfScoresY);                         // = dvector(0, ncaty - 1);
        Eigen::VectorXd FrU12(npara);                                // = dvector(0, npara - 1);
        Eigen::VectorXd FrU21(npara);                                // = dvector(0, npara - 1);
        Eigen::VectorXd GsV12(npara);                                // = dvector(0, npara - 1);
        Eigen::VectorXd GsV21(npara);                                // = dvector(0, npara - 1);
        Eigen::MatrixXd Cr12(totalNumberOfScores, npara);            // = dmatrix(0, ncat - 1, 0, npara - 1);
        Eigen::MatrixXd Cr21(totalNumberOfScores, npara);            // = dmatrix(0, ncat - 1, 0, npara - 1);
        Eigen::MatrixXd M(numberOfScoresX, totalNumberOfScores);     // = dmatrix(0, ncatx - 1, 0, ncat - 1);
        Eigen::MatrixXd fitbdist1(numberOfScoresX, numberOfScoresY); // = dmatrix(0, ncatx - 1, 0, ncaty - 1);
        Eigen::MatrixXd fitbdist2(numberOfScoresX, numberOfScoresY); // = dmatrix(0, ncatx - 1, 0, ncaty - 1);
        Eigen::MatrixXd U12(numberOfScoresX, npara);                 // = dmatrix(0, ncatx - 1, 0, npara - 1);
        Eigen::MatrixXd U21(numberOfScoresX, npara);                 // = dmatrix(0, ncatx - 1, 0, npara - 1);
        Eigen::MatrixXd N(numberOfScoresY, totalNumberOfScores);     // = dmatrix(0, ncaty - 1, 0, ncat - 1);
        Eigen::MatrixXd V12(numberOfScoresY, npara);                 // = dmatrix(0, ncaty - 1, 0, npara - 1);
        Eigen::MatrixXd V21(numberOfScoresY, npara);                 // = dmatrix(0, ncaty - 1, 0, npara - 1);
        Eigen::MatrixXd B(npara, totalNumberOfScores);               // = dmatrix(0, npara - 1, 0, ncat - 1);

        /*The following code set a natural design matrix corresponding to a natural basis of 
					polynomials. note that this B is the transpose of the design matrix in the loglinear 
				model. */

        B.setZero();

        /*First, assign natrual design matrix first for rows corresponding to x */
        for (size_t i = 0; i < numberOfDegreesSmoothingX; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            for (size_t k = 0; k < numberOfScoresY; k++) {
              B(i, k * numberOfScoresX + j) = std::pow(static_cast<double>(j), i + 1);
            }
          }
        }

        /* then assign values to rows corresponding to y */
        for (size_t i = numberOfDegreesSmoothingX; i < numberOfDegreesSmoothingX + numberOfDegreesSmoothingY; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            for (size_t k = 0; k < numberOfScoresY; k++) {
              B(i, k * numberOfScoresX + j) = std::pow(static_cast<double>(k), i - numberOfDegreesSmoothingX + 1);
            }
          }
        }

        /* assign value to the last rows corresponding to the interaction terms */
        for (size_t i = 0; i < numberOfDegreesSmoothingX; i++) {
          for (size_t j = 0; j < numberOfScoresX; j++) {
            for (size_t k = 0; k < numberOfScoresY; k++) {
              B(i + numberOfDegreesSmoothingX + numberOfDegreesSmoothingY,
                k * numberOfScoresX + j) = std::pow(static_cast<double>(j), interx(i)) *
                                           std::pow(static_cast<double>(k), intery(i));
            }
          }
        }

        for (size_t i = 0; i < numberOfScoresX; i++) {
          for (size_t j = 0; j < numberOfScoresY; j++) {
            fitbdist1(i, j) = bivariateLogLinearSmoothing1.fittedBivariateFreqDist(i, j) / static_cast<double>(numberOfExaminees1); // bivar1->bfd[i][j] / np1;
            fitbdist2(i, j) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(i, j) / static_cast<double>(numberOfExaminees2); // bivar2->bfd[i][j] / np2;
          }
        }

        vPMN(numberOfScoresX, numberOfScoresY, fitbdist1, vP12, M, N);
        vPMN(numberOfScoresX, numberOfScoresY, fitbdist2, vP21, M, N);

        r1 = M * vP12;
        s2 = N * vP12;
        r2 = M * vP21;
        s1 = N * vP21;

        for (size_t i = 0; i < numberOfScoresX; i++) {
          r(i) = wts(0) * r1(i) + (1.0 - wts(0)) * r2(i);
        }

        for (size_t i = 0; i < numberOfScoresY; i++) {
          s(i) = wts(1) * s1(i) + (1.0 - wts(1)) * s2(i);
        }

        computeCmatrixGen(totalNumberOfScores, npara, numberOfExaminees1, B, vP12, Cr12);
        computeCmatrixGen(totalNumberOfScores, npara, numberOfExaminees2, B, vP21, Cr21);

        U12 = M * Cr12; // MatrixMultiMatrix(ncatx, ncat, npara, M, Cr12, U12);
        V12 = N * Cr12; // MatrixMultiMatrix(ncaty, ncat, npara, N, Cr12, V12);
        U21 = M * Cr21; // MatrixMultiMatrix(ncatx, ncat, npara, M, Cr21, U21);
        V21 = N * Cr21; // MatrixMultiMatrix(ncaty, ncat, npara, N, Cr21, V21);

        double bandwithX = optimalBandwithXParameter(numberOfScoresX, scoresX, r, 1);
        double bandwithY = optimalBandwithXParameter(numberOfScoresY, scoresY, s, 1);

        kernelEquate(numberOfScoresX, scoresX, r, bandwithX, numberOfScoresY, scoresY, s, bandwithY, equatedX);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          partialFPartialr(numberOfScoresX, scoresX, r, bandwithX, Fr, scoresX(i));
          partialFPartialr(numberOfScoresY, scoresY, s, bandwithY, Gs, equatedX(i));

          double Gp = kernelPdf(numberOfScoresY, scoresY, s, bandwithY, equatedX(i));

          FrU12 = Fr * U12; // VectorMultiMatrix(ncatx, npara, Fr, U12, FrU12);
          GsV12 = Gs * V12; // VectorMultiMatrix(ncaty, npara, Gs, V12, GsV12);
          FrU21 = Fr * U21; // VectorMultiMatrix(ncatx, npara, Fr, U21, FrU21);
          GsV21 = Gs * V21; // VectorMultiMatrix(ncaty, npara, Gs, V21, GsV21);

          for (size_t j = 0; j < npara; j++) {
            FrU12(j) *= wts[0];
            FrU12(j) -= (1 - wts(1)) * GsV12(j);
            FrU21(j) *= (1 - wts(0));
            FrU21(j) -= wts(1) * GsV21(j);
          }

          see(i) = std::sqrt((FrU12.squaredNorm() + FrU21.squaredNorm()) / Gp); // (VectorNormSq(npara, FrU12) + VectorNormSq(npara, FrU21)) / Gp;
        }
      }

      /*--------------------------------------------------------------------------
			KernelEquateNEATPS
			
			functionality:

			Computes kernel equating function for the NEAT design 
			with post stratification method based on procedures in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			author: Tianyou Wang 3/29/2005.
			
			input:
				bivar1      information about bivariate distribution for form X
				bivar2      information about bivariate distribution for form Y
				wts         wts for group taking form X.
		
			output:
				Equatedx    a vector containing the equated score 
		--------------------------------------------------------------------------*/
      void kernelEquateNEATPS(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing1,
                              const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing2,
                              const double& wts,
                              Eigen::VectorXd& equatedX) {
        size_t numberOfScoresX = bivariateLogLinearSmoothing1.numberOfScoresX;
        size_t numberOfScoresY = bivariateLogLinearSmoothing2.numberOfScoresX;
        size_t numberOfScoresV = bivariateLogLinearSmoothing1.numberOfScoresV;
        size_t numberOfExaminees1 = bivariateLogLinearSmoothing1.numberOfExaminees;
        size_t numberOFExaminees2 = bivariateLogLinearSmoothing2.numberOfExaminees;

        Eigen::VectorXd scoresX(numberOfScoresX);
        Eigen::VectorXd scoresY(numberOfScoresY);

        for (size_t i = 0; i < numberOfScoresX; i++) {
          scoresX(i) = bivariateLogLinearSmoothing1.minimumRawScoreX +
                       static_cast<double>(i) * bivariateLogLinearSmoothing1.scoreIncrementX;
        }

        for (size_t i = 0; i < numberOfScoresY; i++) {
          scoresY(i) = bivariateLogLinearSmoothing2.numberOfScoresX +
                       static_cast<double>(i) * bivariateLogLinearSmoothing2.scoreIncrementX;
        }

        size_t numberOfScores1 = numberOfScoresX * numberOfScoresV;
        size_t numberOfScores2 = numberOfScoresY * numberOfScoresV;

        Eigen::VectorXd vP1(numberOfScores1);   // = dvector(0, ncat1 - 1);
        Eigen::VectorXd vP2(numberOfScores2);   // = dvector(0, ncat2 - 1);
        Eigen::VectorXd vP121(numberOfScores1); // = dvector(0, ncat1 - 1);
        Eigen::VectorXd vP212(numberOfScores2); // = dvector(0, ncat2 - 1);
        Eigen::VectorXd r1(numberOfScoresX);    // = dvector(0, ncatx - 1);
        Eigen::VectorXd s1(numberOfScoresY);    // = dvector(0, ncaty - 1);
        Eigen::VectorXd r2(numberOfScoresX);    // = dvector(0, ncatx - 1);
        Eigen::VectorXd s2(numberOfScoresY);    // = dvector(0, ncaty - 1);
        Eigen::VectorXd r(numberOfScoresX);     // = dvector(0, ncatx - 1);
        Eigen::VectorXd s(numberOfScoresY);     // = dvector(0, ncaty - 1);
        Eigen::VectorXd v1(numberOfScoresV);    // = dvector(0, ncatv - 1);
        Eigen::VectorXd v2(numberOfScoresV);    // = dvector(0, ncatv - 1);

        Eigen::MatrixXd M1(numberOfScoresX, numberOfScores1);        // = dmatrix(0, ncatx - 1, 0, ncat1 - 1);
        Eigen::MatrixXd N1(numberOfScoresV, numberOfScores1);        // = dmatrix(0, ncatv - 1, 0, ncat1 - 1);
        Eigen::MatrixXd M2(numberOfScoresY, numberOfScores2);        // = dmatrix(0, ncaty - 1, 0, ncat2 - 1);
        Eigen::MatrixXd N2(numberOfScoresV, numberOfScores2);        // = dmatrix(0, ncatv - 1, 0, ncat2 - 1);
        Eigen::MatrixXd fitbdist1(numberOfScoresX, numberOfScoresV); // = dmatrix(0, ncatx - 1, 0, ncatv - 1);
        Eigen::MatrixXd fitbdist2(numberOfScoresY, numberOfScoresV); // = dmatrix(0, ncaty - 1, 0, ncatv - 1);

        for (size_t j = 0; j < numberOfScoresV; j++) {
          for (size_t i = 0; i < numberOfScoresX; i++) {
            fitbdist1(i, j) = bivariateLogLinearSmoothing1.fittedBivariateFreqDist(i, j) /
                              static_cast<double>(numberOfExaminees1); // bivar1->bfd[i][j] / np1;
          }

          for (size_t i = 0; i < numberOfScoresY; i++) {
            fitbdist2(i, j) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(i, j) /
                              static_cast<double>(numberOFExaminees2); // fitbdist2[i][j] = bivar2->bfd[i][j] / np2;
          }
        }

        vPMN(numberOfScoresX, numberOfScoresV, fitbdist1, vP1, M1, N1);
        vPMN(numberOfScoresY, numberOfScoresV, fitbdist2, vP2, M2, N2);

        r1 = M1 * vP1; // MatrixMultiVector(numberOfScoresX, numberOfScores1, M1, vP1, r1);
        v1 = N1 * vP1; // MatrixMultiVector(numberOfScoresV, numberOfScores1, N1, vP1, v1);
        s2 = M2 * vP2; // MatrixMultiVector(numberOfScoresY, numberOfScores2, M2, vP2, s2);
        v2 = N2 * vP2; // MatrixMultiVector(ncatv, ncat2, N2, vP2, v2);

        for (size_t j = 0; j < numberOfScoresV; j++) {
          for (size_t i = 0; i < numberOfScoresX; i++) {
            vP121(i + j * numberOfScoresX) = vP1(i + j * numberOfScoresX) * v2(j) / v1(j);
          }
        }

        for (size_t j = 0; j < numberOfScoresV; j++) {
          for (size_t i = 0; i < numberOfScoresY; i++) {
            vP212(i + j + numberOfScoresY) = vP2(i + j * numberOfScoresY) * v1(j) / v2(j);
          }
        }

        r2 = M1 * vP121;
        s1 = M2 * vP212;

        for (size_t i = 0; i < numberOfScoresX; i++) {
          r(i) = wts * r1(i) + (1.0 - wts) * r2(i);
        }

        for (size_t i = 0; i < numberOfScoresY; i++) {
          s[i] = wts * s1[i] + (1 - wts) * s2[i];
        }

        double bandwithX = optimalBandwithXParameter(numberOfScoresX, scoresX, r, 1);
        double bandwithY = optimalBandwithXParameter(numberOfScoresY, scoresY, s, 1);

        kernelEquate(numberOfScoresX, scoresX, r, bandwithX, numberOfScoresY, scoresY, s, bandwithY, equatedX);
      }

      /*--------------------------------------------------------------------------
			KernelEquateSEENEATPS3
			
			functionality:

			Computes kernel equating function and SEE for the NEAT design 
			with post stratification method based on procedures in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			author: Tianyou Wang 3/29/2005.
			
			input:
				bivar1      information about bivariate distribution for form X
				bivar2      information about bivariate distribution for form Y
				wts         wts for group taking form X.
		
			output:
				Equatedx    a vector containing the equated score 
				SEE         a vector containing standard error of equating 
		--------------------------------------------------------------------------*/
      void kernelEquateSEENEATPS3(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing1,
                                  const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing2,
                                  const double& wts,
                                  Eigen::VectorXd& equatedX,
                                  Eigen::VectorXd& see) {
        // int i, j, k, ncat1, ncat2, ncatx, ncaty, ncatv, npara1, npara2;
        // int cu1, cv1, cuv1, cu2, cv2, cuv2;
        // long np1, np2;
        // double *vP1, *vP2, *vPP1, *vPP2, **M1, **N1, **M2, **N2, *r1, *s1, *r2, *s2, *r, *s, *v1, *v2;
        // /* r1 is for X taken first(first group), r2 is for X taken second (second group),
        // s1 is for Y taken first(second group), s2 is for Y taken second (first group)
        // vP1 is for first group, vP2 is for second group */
        // double *vP121, *vP212;
        // double **UR, **US, **VR, **VS, **UP, **UPStar, **UQ, **UQStar, **Utemp, **Vtemp;
        // double *pl, *vPl, *ql, *vQl;
        // double hx, hy;
        // double **Cr1, **Cr1Star, **Cr2, **Cr2Star, **B1, **B2;
        // double *Fr, *Gs, Gp, *FrUR, *GsVR, *FrUS, *GsVS;
        // double *scoresx, *scoresy;
        // int *interx, *interv1, *interv2, *intery;
        // double **fitbdist1, **fitbdist2, sum1 = 0, sum2 = 0;
        // int currentrow1, currentrow2;
        // int **cpm1, **cpm2;

        size_t cuv1 = bivariateLogLinearSmoothing1.numberOfCrossProductMoments;
        size_t cu1 = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingU;
        size_t cv1 = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingV;
        size_t cuv2 = bivariateLogLinearSmoothing2.numberOfCrossProductMoments;
        size_t cu2 = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingU;
        size_t cv2 = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingV;

        Eigen::MatrixXi cpm1(2, cuv1);

        for (size_t i = 0; i < cuv1; i++) {
          for (size_t j = 0; j < 2; j++) {
            cpm1(j, i) = bivariateLogLinearSmoothing1.crossProductMoments(i, j);
          }
        }

        Eigen::MatrixXi cpm2(2, cuv2);

        for (size_t i = 0; i < cuv2; i++) {
          for (size_t j = 0; j < 2; j++) {
            cpm2(j, i) = bivariateLogLinearSmoothing2.crossProductMoments(i, j);
          }
        }

        size_t ncatx = bivariateLogLinearSmoothing1.numberOfScoresX;
        size_t ncaty = bivariateLogLinearSmoothing2.numberOfScoresX;
        size_t ncatv = bivariateLogLinearSmoothing1.numberOfScoresV;
        size_t npara1 = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingU +
                        bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingV +
                        bivariateLogLinearSmoothing1.numberOfCrossProductMoments + 10;
        size_t npara2 = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingU +
                        bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingV +
                        bivariateLogLinearSmoothing2.numberOfCrossProductMoments + 10;
        Eigen::VectorXd interx = cpm1.row(0).cast<double>();
        Eigen::VectorXd interv1 = cpm1.row(1).cast<double>();
        Eigen::VectorXd intery = cpm2.row(0).cast<double>();
        Eigen::VectorXd interv2 = cpm2.row(1).cast<double>();

        size_t np1 = bivariateLogLinearSmoothing1.numberOfExaminees;
        size_t np2 = bivariateLogLinearSmoothing2.numberOfExaminees;
        Eigen::VectorXd scoresX(ncatx);
        for (size_t i = 0; i < ncatx; i++) {
          scoresX(i) = bivariateLogLinearSmoothing1.minimumRawScoreX +
                       static_cast<double>(i) * bivariateLogLinearSmoothing1.scoreIncrementX;
        }

        Eigen::VectorXd scoresY(ncaty);

        for (size_t i = 0; i < ncaty; i++) {
          scoresY(i) = bivariateLogLinearSmoothing2.minimumRawScoreX +
                       static_cast<double>(i) * bivariateLogLinearSmoothing2.scoreIncrementX;
        }

        size_t ncat1 = ncatx * ncatv;
        size_t ncat2 = ncaty * ncatv;

        Eigen::VectorXd vP1(ncat1);   // vP1 =
        Eigen::VectorXd vP121(ncat1); // vP121 =
        Eigen::VectorXd vPP1(ncat1);  // vPP1 =
        Eigen::VectorXd vP2(ncat2);   // vP2 =
        Eigen::VectorXd vP212(ncat2); // vP212 =
        Eigen::VectorXd vPP2(ncat2);  // vPP2 =
        Eigen::VectorXd v1(ncatv);    // v1 =
        Eigen::VectorXd v2(ncatv);    // v2 =
        Eigen::VectorXd Fr(ncatx);    // Fr =
        Eigen::VectorXd pl(ncatx);    // pl =
        Eigen::VectorXd r(ncatx);     // r =
        Eigen::VectorXd r1(ncatx);    // r1 =
        Eigen::VectorXd r2(ncatx);    // r2 =
        Eigen::VectorXd Gs(ncaty);    // Gs =
        Eigen::VectorXd ql(ncaty);    // ql =
        Eigen::VectorXd s(ncaty);     // s =
        Eigen::VectorXd s1(ncaty);    // s1 =
        Eigen::VectorXd s2(ncaty);    // s2 =
        Eigen::VectorXd FrUR(npara1); // FrUR =
        Eigen::VectorXd GsVR(npara1); // GsVR =
        Eigen::VectorXd vPl(npara1);  // vPl =
        Eigen::VectorXd FrUS(npara2); // FrUS =
        Eigen::VectorXd GsVS(npara2); // GsVS =
        Eigen::VectorXd vQl(npara2);  // vQl =

        Eigen::MatrixXd Cr1(ncat1, npara1);      // Cr1 =
        Eigen::MatrixXd Cr1Star(ncat1, npara1);  // Cr1Star =
        Eigen::MatrixXd Cr2(ncat2, npara2);      // Cr2 =
        Eigen::MatrixXd Cr2Star(ncat2, npara2);  // Cr2Star =
        Eigen::MatrixXd N1(ncatv, ncat1);        // N1 =
        Eigen::MatrixXd N2(ncatv, ncat2);        // N2 =
        Eigen::MatrixXd M1(ncatx, ncat1);        // M1 =
        Eigen::MatrixXd fitbdist1(ncatx, ncatv); // fitbdist1 =
        Eigen::MatrixXd UP(ncatx, npara1);       // UP =
        Eigen::MatrixXd UPStar(ncatx, npara1);   // UPStar =
        Eigen::MatrixXd UR(ncatx, npara1);       // UR =
        Eigen::MatrixXd Utemp(ncatx, npara1);    // Utemp =
        Eigen::MatrixXd US(ncatx, npara2);       // US =
        Eigen::MatrixXd M2(ncaty, ncat2);        // M2 =
        Eigen::MatrixXd fitbdist2(ncaty, ncatv); // fitbdist2 =
        Eigen::MatrixXd VR(ncaty, npara1);       // VR =
        Eigen::MatrixXd UQ(ncaty, npara2);       // UQ =
        Eigen::MatrixXd UQStar(ncaty, npara2);   // UQStar =
        Eigen::MatrixXd VS(ncaty, npara2);       // VS =
        Eigen::MatrixXd Vtemp(ncaty, npara2);    // Vtemp =
        Eigen::MatrixXd B1(npara1, ncat1);       // B1 =
        Eigen::MatrixXd B2(npara2, ncat2);       // B2 =

        /*The following code set a natural design matrix corresponding to a natural basis of 
					polynomials. note that this B is the transpose of the design matrix in the loglinear 
				  model. */

        /*First, assign natural design matrix first for rows corresponding to x */
        for (size_t i = 0; i < cu1; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B1(i, k + ncatx + j) = std::pow(static_cast<double>(j), i + 1);
            }
          }
        }

        /* then assign values to rows corresponding to v */
        size_t currentrow1 = cu1;
        for (size_t i = 0; i < cv1; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B1(currentrow1 + i, k * ncatx + j) = std::pow(static_cast<double>(k), (i + 1));
            }
          }
        }

        /* assign value to the last rows corresponding to the interaction terms */
        currentrow1 += cv1;
        for (size_t i = 0; i < cuv1; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B1(currentrow1 + i, k * ncatx + j) = std::pow(static_cast<double>(j), interx(i)) *
                                                   std::pow(static_cast<double>(k), interv1(i));
            }
          }
        }

        /* assign values to the rows for the 0 score on form X */
        currentrow1 += cuv1;
        for (size_t j = 0; j < ncatx; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (j == 0) {
              B1(currentrow1, k * ncatx + j) = 1;
            } else {
              B1(currentrow1, k * ncatx + j) = 0;
            }
          }
        }

        /* assign values to the rows for the 0 score on common set v */
        currentrow1++;
        for (size_t j = 0; j < ncatx; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (k == 0) {
              B1(currentrow1, k * ncatx + j) = 1;
            } else {
              B1(currentrow1, k * ncatx + j) = 0;
            }
          }
        }

        /*assign value to the rows for score gaps in x 	*/
        currentrow1++;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if ((j % 5) == 0 && j != 0) {
                B1(currentrow1 + i, k * ncatx + j) =
                    std::pow(static_cast<double>(j), (double)i);
              } else {
                B1(currentrow1 + i, k * ncatx + j) = 0;
              }
            }
          }
        }

        /*assign value to the rows for score gaps in v */
        currentrow1 += 4;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if (((k) % 5) == 0 && k != 0)
                B1(currentrow1 + i, k * ncatx + j) = std::pow(static_cast<double>(k), i);
              else
                B1(currentrow1 + i, k * ncatx + j) = 0;
            }
          }
        }

        /*The following code set a natrual design matrix corresponding to a natrual basis of 
					polynomials. note that this B2 is the transpose of the design matrix in the loglinear 
				model. */

        /*First, assign natrual design matrix first for rows corresponding to y */
        for (size_t i = 0; i < cu2; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B2(i, k * ncaty + j) = std::pow(static_cast<double>(j), (i + 1));
            }
          }
        }

        /* then assign values to rows corresponding to v */
        size_t currentrow2 = cu2;
        for (size_t i = 0; i < cv2; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B2(currentrow2 + 1, k * ncaty + j) = std::pow(static_cast<double>(k), (i + 1));
            }
          }
        }

        /* assign value to the last rows corresponding to the interaction terms */
        currentrow2 += cv2;
        for (size_t i = 0; i < cuv2; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B2(currentrow2 + i, k * ncaty + j) = std::pow(static_cast<double>(j), intery(i)) *
                                                   std::pow(static_cast<double>(k), interv2(i));
            }
          }
        }

        /* assign values to the rows for the 0 score on form Y  */
        currentrow2 += cuv2;
        for (size_t j = 0; j < ncaty; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (j == 0) {
              B2(currentrow2, k * ncaty + j) = 1;
            } else {
              B2(currentrow2, k * ncaty + j) = 0;
            }
          }
        }

        /* assign values to the rows for the 0 score on common set V  */
        currentrow2++;
        for (size_t j = 0; j < ncaty; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (k == 0) {
              B2(currentrow2, k * ncaty + j) = 1;
            } else {
              B2(currentrow2, k * ncaty + j) = 0;
            }
          }
        }

        /*assign value to the rows for score gaps in y 	*/
        currentrow2++;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if ((j % 5) == 0 && j != 0) {
                B2(currentrow2 + i, k * ncaty + j) = std::pow(static_cast<double>(j), i);
              }

              else
                B2(currentrow2 + i, k * ncaty + j) = 0;
            }
          }
        }

        /*assign value to the rows for score gaps in v */
        currentrow2 += 4;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if (((k) % 5) == 0 && k != 0) {
                B2(currentrow2 + i, k * ncaty + j) = std::pow(static_cast<double>(k), i);
              } else {
                B2(currentrow2 + i, k * ncaty + j) = 0;
              }
            }
          }
        }

        double sum1 = 0;
        double sum2 = 0;
        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncatx; i++) {
            fitbdist1(i, j) = bivariateLogLinearSmoothing1.fittedBivariateFreqDist(i, j) /
                              static_cast<double>(np1);
            sum1 += fitbdist1(i, j);
          }

          for (size_t i = 0; i < ncaty; i++) {
            fitbdist2(i, j) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(i, j) /
                              static_cast<double>(np2);
            sum2 += fitbdist2(i, j);
          }
        }

        vPMN(ncatx, ncatv, fitbdist1, vP1, M1, N1);
        vPMN(ncaty, ncatv, fitbdist2, vP2, M2, N2);
        vPT(ncatx, ncatv, fitbdist1, vPP1);
        vPT(ncaty, ncatv, fitbdist2, vPP2);
        r1 = M1 * vP1; // MatrixMultiVector(ncatx, ncat1, M1, vP1, r1);
        v1 = N1 * vP1; // MatrixMultiVector(ncatv, ncat1, N1, vP1, v1);
        s2 = M2 * vP2; // MatrixMultiVector(ncaty, ncat2, M2, vP2, s2);
        v2 = N2 * vP2; // MatrixMultiVector(ncatv, ncat2, N2, vP2, v2);

        computeCmatrixGen(ncat1, npara1, np1, B1, vP1, Cr1);
        computeCmatrixGen(ncat2, npara2, np2, B2, vP2, Cr2);

        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncatx; i++) {
            vP121(i + j * ncatx) = vP1(i + j * ncatx) * v2(j) / v1(j);
            for (size_t k = 0; k < npara1; k++)
              Cr1Star(i + j * ncatx, k) = Cr1(i + j * ncatx, k) * v2[j] / v1[j];
          }
        }

        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncaty; i++) {
            vP212(i + j * ncaty) = vP2(i + j * ncaty) * v1(j) / v2(j);
            for (size_t k = 0; k < npara2; k++)
              Cr2Star(i + j * ncatx, k) = Cr2(i + j * ncatx, k) * v1(j) / v2(j);
          }
        }

        for (size_t i = 0; i < ncatx; i++) {
          for (size_t k = 0; k < npara1; k++)
            Utemp(i, k) = 0;
          for (size_t k = 0; k < npara2; k++)
            US(i, k) = 0;
        }

        for (size_t i = 0; i < ncaty; i++) {
          for (size_t k = 0; k < npara2; k++)
            Vtemp(i, k) = 0;
          for (size_t k = 0; k < npara1; k++)
            VR(i, k) = 0;
        }

        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncatx; i++)
            pl(i) = fitbdist1(i, j);
          for (size_t k = 0; k < npara1; k++) {
            vPl(k) = 0;

            for (size_t i = 0; i < ncatx; i++)
              vPl(k) += Cr1(j * ncatx + i, k);
          }
          for (size_t i = 0; i < ncatx; i++)
            for (size_t k = 0; k < npara1; k++)
              Utemp(i, k) += v2(j) / v1(j) / v1(j) * pl(i) * vPl(k);

          for (size_t i = 0; i < ncaty; i++)
            ql(i) = fitbdist2(i, j);
          for (size_t k = 0; k < npara2; k++) {
            vQl(k) = 0;

            for (size_t i = 0; i < ncaty; i++)
              vQl(k) += Cr2(j * ncaty + i, k);
          }
          for (size_t i = 0; i < ncaty; i++)
            for (size_t k = 0; k < npara2; k++)
              Vtemp(i, k) += v1(j) / v2(j) / v2(j) * ql(i) * vQl(k);

          for (size_t i = 0; i < ncaty; i++)
            for (size_t k = 0; k < npara1; k++)
              VR(i, k) += wts / v2(j) * ql(i) * vPl(k);

          for (size_t i = 0; i < ncatx; i++)
            for (size_t k = 0; k < npara2; k++)
              US(i, k) += (1 - wts) / v1(j) * pl(i) * vQl(k);
        }

        UP = M1 * Cr1;         // MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1, UP);
        UQ = M2 * Cr2;         // MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2, UQ);
        UPStar = M1 * Cr1Star; // MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1Star, UPStar);
        UQStar = M2 * Cr2Star; // MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2Star, UQStar);
        r2 = M1 * vP121;       // MatrixMultiVector(ncatx, ncat1, M1, vP121, r2);
        s1 = M2 * vP212;       // MatrixMultiVector(ncatx, ncat2, M2, vP212, s1);

        for (size_t i = 0; i < ncatx; i++)
          r(i) = wts * r1(i) + (1 - wts) * r2(i);
        for (size_t i = 0; i < ncaty; i++)
          s(i) = wts * s1(i) + (1 - wts) * s2(i);

        double hx = optimalBandwithXParameter(ncatx, scoresX, r, 1);
        double hy = optimalBandwithXParameter(ncaty, scoresY, s, 1);
        hx = 1.9243;
        hy = 2.0056;

        kernelEquate(ncatx, scoresX, r, hx, ncaty, scoresY, s, hy, equatedX);

        for (size_t i = 0; i < ncatx; i++) {
          for (size_t k = 0; k < npara1; k++)
            UR(i, k) = wts * UP(i, k) + (1 - wts) * (UPStar(i, k) - Utemp(i, k));
        }

        for (size_t i = 0; i < ncaty; i++) {
          for (size_t k = 0; k < npara2; k++)
            VS(i, k) = (1 - wts) * UQ(i, k) + wts * (UQStar(i, k) - Vtemp(i, k));
        }

        for (size_t i = 0; i < ncatx; i++) {
          partialFPartialr(ncatx, scoresX, r, hx, Fr, scoresX(i));
          partialFPartialr(ncaty, scoresY, s, hy, Gs, equatedX(i));
          double Gp = kernelPdf(ncaty, scoresY, s, hy, equatedX(i));
          FrUR = Fr * UR; // VectorMultiMatrix(ncatx, npara1, Fr, UR, FrUR);
          GsVR = Gs * VR; // VectorMultiMatrix(ncaty, npara1, Gs, VR, GsVR);
          FrUS = Fr * US; // VectorMultiMatrix(ncatx, npara2, Fr, US, FrUS);
          GsVS = Gs * VS; // VectorMultiMatrix(ncaty, npara2, Gs, VS, GsVS);
          for (size_t j = 0; j < npara1; j++)
            FrUR(j) -= GsVR(j);
          for (size_t j = 0; j < npara2; j++)
            FrUS(j) -= GsVS(j);
          see(i) = std::sqrt(FrUR.squaredNorm() + FrUS.squaredNorm()) / Gp;
        }
      }

      /*--------------------------------------------------------------------------
			KernelEquateSEENEATPS
			
			functionality:

			Computes kernel equating function and SEE for the NEAT design 
			with post stratification method based on procedures in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			Author: Tianyou Wang 3/29/2005.
			
			input:
				bivar1      information about bivariate distribution for form X
				bivar2      information about bivariate distribution for form Y
				wts         wts for group taking form X.
		
			output:
				Equatedx    a vector containing the equated score 
				SEE         a vector containing standard error of equating 
		--------------------------------------------------------------------------*/
      void kernelEquateSEENEATPS(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing1,
                                 const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing2,
                                 const double& wts,
                                 Eigen::VectorXd& equatedX,
                                 Eigen::VectorXd& see) {
        // int i, j, k, ncat1, ncat2, ncatx, ncaty, ncatv, npara1, npara2;
        // int cu1, cv1, cuv1, cu2, cv2, cuv2;
        // long np1, np2;
        // double *vP1, *vP2, *vPP1, *vPP2, **M1, **N1, **M2, **N2, *r1, *s1, *r2, *s2, *r, *s, *v1, *v2;
        // /* r1 is for X taken first(first group), r2 is for X taken second (second group),
        // s1 is for Y taken first(second group), s2 is for Y taken second (first group)
        // vP1 is for first group, vP2 is for second group */
        // double *vP121, *vP212;
        // double **UR, **US, **VR, **VS, **UP, **UPStar, **UQ, **UQStar, **Utemp, **Vtemp;
        // double *pl, *vPl, *ql, *vQl;
        // double hx, hy;
        // double **Cr1, **Cr1Star, **Cr2, **Cr2Star, **B1, **B2;
        // double *Fr, *Gs, Gp, *FrUR, *GsVR, *FrUS, *GsVS;
        // double *scoresx, *scoresy;
        // int *interx, *interv1, *interv2, *intery;
        // double **fitbdist1, **fitbdist2, sum1 = 0, sum2 = 0;
        // int currentrow1, currentrow2;
        // int **cpm1, **cpm2;

        size_t cuv1 = bivariateLogLinearSmoothing1.numberOfCrossProductMoments;
        size_t cu1 = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingU;
        size_t cv1 = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingV;
        size_t cuv2 = bivariateLogLinearSmoothing2.numberOfCrossProductMoments;
        size_t cu2 = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingU;
        size_t cv2 = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingV;
        Eigen::MatrixXi cpm1(2, cuv1);

        for (size_t i = 0; i < cuv1; i++)
          for (size_t j = 0; j < 2; j++)
            cpm1(j, i) = bivariateLogLinearSmoothing1.crossProductMoments(i, j);

        Eigen::MatrixXi cpm2(2, cuv2);
        for (size_t i = 0; i < cuv2; i++)
          for (size_t j = 0; j < 2; j++)
            cpm2(j, i) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(i, j);

        size_t ncatx = bivariateLogLinearSmoothing1.numberOfScoresX;
        size_t ncaty = bivariateLogLinearSmoothing2.numberOfScoresX;
        size_t ncatv = bivariateLogLinearSmoothing1.numberOfScoresV;
        size_t npara1 = bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingU +
                        bivariateLogLinearSmoothing1.numberOfDegreesOfSmoothingV +
                        bivariateLogLinearSmoothing1.numberOfCrossProductMoments + 10;
        size_t npara2 = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingU +
                        bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingV +
                        bivariateLogLinearSmoothing2.numberOfCrossProductMoments + 10;

        Eigen::VectorXd interx = cpm1.row(0).cast<double>();
        Eigen::VectorXd interv1 = cpm1.row(1).cast<double>();
        Eigen::VectorXd intery = cpm2.row(0).cast<double>();
        Eigen::VectorXd interv2 = cpm2.row(1).cast<double>();

        size_t np1 = bivariateLogLinearSmoothing1.numberOfExaminees;
        size_t np2 = bivariateLogLinearSmoothing2.numberOfExaminees;
        Eigen::VectorXd scoresx(ncatx);

        for (size_t i = 0; i < ncatx; i++)
          scoresx(i) = bivariateLogLinearSmoothing1.minimumRawScoreX +
                       i * bivariateLogLinearSmoothing1.scoreIncrementX;

        Eigen::VectorXd scoresy(ncaty);
        for (size_t i = 0; i < ncaty; i++)
          scoresy(i) = bivariateLogLinearSmoothing2.minimumRawScoreX +
                       i * bivariateLogLinearSmoothing2.scoreIncrementX;

        size_t ncat1 = ncatx * ncatv;
        size_t ncat2 = ncaty * ncatv;

        Eigen::VectorXd vP1(ncat1);   // vP1 =
        Eigen::VectorXd vP121(ncat1); // vP121 =
        Eigen::VectorXd vPP1(ncat1);  // vPP1 =
        Eigen::VectorXd vP2(ncat2);   // vP2 =
        Eigen::VectorXd vP212(ncat2); // vP212 =
        Eigen::VectorXd vPP2(ncat2);  // vPP2 =
        Eigen::VectorXd v1(ncatv);    // v1 =
        Eigen::VectorXd v2(ncatv);    // v2 =
        Eigen::VectorXd Fr(ncatx);    // Fr =
        Eigen::VectorXd pl(ncatx);    // pl =
        Eigen::VectorXd r(ncatx);     // r =
        Eigen::VectorXd r1(ncatx);    // r1 =
        Eigen::VectorXd r2(ncatx);    // r2 =
        Eigen::VectorXd Gs(ncaty);    // Gs =
        Eigen::VectorXd ql(ncaty);    // ql =
        Eigen::VectorXd s(ncaty);     // s =
        Eigen::VectorXd s1(ncaty);    // s1 =
        Eigen::VectorXd s2(ncaty);    // s2 =
        Eigen::VectorXd FrUR(npara1); // FrUR =
        Eigen::VectorXd GsVR(npara1); // GsVR =
        Eigen::VectorXd vPl(npara1);  // vPl =
        Eigen::VectorXd FrUS(npara2); // FrUS =
        Eigen::VectorXd GsVS(npara2); // GsVS =
        Eigen::VectorXd vQl(npara2);  // vQl =

        Eigen::MatrixXd Cr1(ncat1, npara1);      // Cr1 =
        Eigen::MatrixXd Cr1Star(ncat1, npara1);  // Cr1Star =
        Eigen::MatrixXd Cr2(ncat2, npara2);      // Cr2 =
        Eigen::MatrixXd Cr2Star(ncat2, npara2);  // Cr2Star =
        Eigen::MatrixXd N1(ncatv, ncat1);        // N1 =
        Eigen::MatrixXd N2(ncatv, ncat2);        // N2 =
        Eigen::MatrixXd M1(ncatx, ncat1);        // M1 =
        Eigen::MatrixXd fitbdist1(ncatx, ncatv); // fitbdist1 =
        Eigen::MatrixXd UP(ncatx, npara1);       // UP =
        Eigen::MatrixXd UPStar(ncatx, npara1);   // UPStar =
        Eigen::MatrixXd UR(ncatx, npara1);       // UR =
        Eigen::MatrixXd Utemp(ncatx, npara1);    // Utemp =
        Eigen::MatrixXd US(ncatx, npara2);       // US =
        Eigen::MatrixXd M2(ncaty, ncat2);        // M2 =
        Eigen::MatrixXd fitbdist2(ncaty, ncatv); // fitbdist2 =
        Eigen::MatrixXd VR(ncaty, npara1);       // VR =
        Eigen::MatrixXd UQ(ncaty, npara2);       // UQ =
        Eigen::MatrixXd UQStar(ncaty, npara2);   // UQStar =
        Eigen::MatrixXd VS(ncaty, npara2);       // VS =
        Eigen::MatrixXd Vtemp(ncaty, npara2);    // Vtemp =
        Eigen::MatrixXd B1(npara1, ncat1);       // B1 =
        Eigen::MatrixXd B2(npara2, ncat2);       // B2 =

        /*The following code set a natrual design matrix corresponding to a natrual basis of 
					polynomials. note that this B is the transpose of the design matrix in the loglinear 
				model. */

        /*First, assign natrual design matrix first for rows corresponding to x */
        for (size_t i = 0; i < cu1; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B1(i, k + ncatx + j) = std::pow(static_cast<double>(j), i + 1);
            }
          }
        }

        /* then assign values to rows corresponding to v */
        size_t currentrow1 = cu1;
        for (size_t i = 0; i < cv1; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B1(currentrow1 + i, k * ncatx + j) = std::pow(static_cast<double>(k), (i + 1));
            }
          }
        }

        /* assign value to the last rows corresponding to the interaction terms */
        currentrow1 += cv1;
        for (size_t i = 0; i < cuv1; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B1(currentrow1 + i, k * ncatx + j) = std::pow(static_cast<double>(j), interx(i)) *
                                                   std::pow(static_cast<double>(k), interv1(i));
            }
          }
        }

        /* assign values to the rows for the 0 score on form X */
        currentrow1 += cuv1;
        for (size_t j = 0; j < ncatx; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (j == 0) {
              B1(currentrow1, k * ncatx + j) = 1;
            } else {
              B1(currentrow1, k * ncatx + j) = 0;
            }
          }
        }

        /* assign values to the rows for the 0 score on common set v */
        currentrow1++;
        for (size_t j = 0; j < ncatx; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (k == 0) {
              B1(currentrow1, k * ncatx + j) = 1;
            } else {
              B1(currentrow1, k * ncatx + j) = 0;
            }
          }
        }

        /*assign value to the rows for score gaps in x 	*/
        currentrow1++;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if ((j % 5) == 0 && j != 0) {
                B1(currentrow1 + i, k * ncatx + j) =
                    std::pow(static_cast<double>(j), (double)i);
              } else {
                B1(currentrow1 + i, k * ncatx + j) = 0;
              }
            }
          }
        }

        /*assign value to the rows for score gaps in v */
        currentrow1 += 4;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncatx; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if (((k) % 5) == 0 && k != 0)
                B1(currentrow1 + i, k * ncatx + j) = std::pow(static_cast<double>(k), i);
              else
                B1(currentrow1 + i, k * ncatx + j) = 0;
            }
          }
        }

        /*The following code set a natrual design matrix corresponding to a natrual basis of 
					polynomials. note that this B2 is the transpose of the design matrix in the loglinear 
				model. */

        /*First, assign natrual design matrix first for rows corresponding to y */
        for (size_t i = 0; i < cu2; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B2(i, k * ncaty + j) = std::pow(static_cast<double>(j), (i + 1));
            }
          }
        }

        /* then assign values to rows corresponding to v */
        size_t currentrow2 = cu2;
        for (size_t i = 0; i < cv2; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B2(currentrow2 + 1, k * ncaty + j) = std::pow(static_cast<double>(k), (i + 1));
            }
          }
        }

        /* assign value to the last rows corresponding to the interaction terms */
        currentrow2 += cv2;
        for (size_t i = 0; i < cuv2; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              B2(currentrow2 + i, k * ncaty + j) = std::pow(static_cast<double>(j), intery(i)) *
                                                   std::pow(static_cast<double>(k), interv2(i));
            }
          }
        }

        /* assign values to the rows for the 0 score on form Y  */
        currentrow2 += cuv2;
        for (size_t j = 0; j < ncaty; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (j == 0) {
              B2(currentrow2, k * ncaty + j) = 1;
            } else {
              B2(currentrow2, k * ncaty + j) = 0;
            }
          }
        }

        /* assign values to the rows for the 0 score on common set V  */
        currentrow2++;
        for (size_t j = 0; j < ncaty; j++) {
          for (size_t k = 0; k < ncatv; k++) {
            if (k == 0) {
              B2(currentrow2, k * ncaty + j) = 1;
            } else {
              B2(currentrow2, k * ncaty + j) = 0;
            }
          }
        }

        /*assign value to the rows for score gaps in y 	*/
        currentrow2++;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if ((j % 5) == 0 && j != 0) {
                B2(currentrow2 + i, k * ncaty + j) = std::pow(static_cast<double>(j), i);
              }

              else
                B2(currentrow2 + i, k * ncaty + j) = 0;
            }
          }
        }

        /*assign value to the rows for score gaps in v */
        currentrow2 += 4;
        for (size_t i = 0; i <= 3; i++) {
          for (size_t j = 0; j < ncaty; j++) {
            for (size_t k = 0; k < ncatv; k++) {
              if (((k) % 5) == 0 && k != 0) {
                B2(currentrow2 + i, k * ncaty + j) = std::pow(static_cast<double>(k), i);
              } else {
                B2(currentrow2 + i, k * ncaty + j) = 0;
              }
            }
          }
        }

        double sum1 = 0;
        double sum2 = 0;
        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncatx; i++) {
            fitbdist1(i, j) = bivariateLogLinearSmoothing1.fittedBivariateFreqDist(i, j) /
                              static_cast<double>(np1);
            sum1 += fitbdist1(i, j);
          }

          for (size_t i = 0; i < ncaty; i++) {
            fitbdist2(i, j) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(i, j) /
                              static_cast<double>(np2);
            sum2 += fitbdist2(i, j);
          }
        }

        vPMN(ncatx, ncatv, fitbdist1, vP1, M1, N1);
        vPMN(ncaty, ncatv, fitbdist2, vP2, M2, N2);
        vPT(ncatx, ncatv, fitbdist1, vPP1);
        vPT(ncaty, ncatv, fitbdist2, vPP2);
        r1 = M1 * vP1; // MatrixMultiVector(ncatx, ncat1, M1, vP1, r1);
        v1 = N1 * vP1; // MatrixMultiVector(ncatv, ncat1, N1, vP1, v1);
        s2 = M2 * vP2; // MatrixMultiVector(ncaty, ncat2, M2, vP2, s2);
        v2 = N2 * vP2; // MatrixMultiVector(ncatv, ncat2, N2, vP2, v2);

        computeCmatrixGen(ncat1, npara1, np1, B1, vP1, Cr1);
        computeCmatrixGen(ncat2, npara2, np2, B2, vP2, Cr2);

        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncatx; i++) {
            vP121(i + j * ncatx) = vP1(i + j * ncatx) * v2(j) / v1(j);
            for (size_t k = 0; k < npara1; k++)
              Cr1Star(i + j * ncatx, k) = Cr1(i + j * ncatx, k) * v2(j) / v1(j);
          }
        }

        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncaty; i++) {
            vP212(i + j * ncaty) = vP2(i + j * ncaty) * v1(j) / v2(j);
            for (size_t k = 0; k < npara2; k++)
              Cr2Star(i + j * ncatx, k) = Cr2(i + j * ncatx, k) * v1(j) / v2(j);
          }
        }
        for (size_t i = 0; i < ncatx; i++) {
          for (size_t k = 0; k < npara1; k++)
            Utemp(i, k) = 0;
          for (size_t k = 0; k < npara2; k++)
            US(i, k) = 0;
        }

        for (size_t i = 0; i < ncaty; i++) {
          for (size_t k = 0; k < npara2; k++)
            Vtemp(i, k) = 0;
          for (size_t k = 0; k < npara1; k++)
            VR(i, k) = 0;
        }

        for (size_t j = 0; j < ncatv; j++) {
          for (size_t i = 0; i < ncatx; i++)
            pl(i) = fitbdist1(i, j);
          for (size_t k = 0; k < npara1; k++) {
            vPl(k) = 0;
            for (size_t i = 0; i < ncatx; i++)
              vPl(k) += Cr1(j * ncatx + i, k);
          }
          for (size_t i = 0; i < ncatx; i++)
            for (size_t k = 0; k < npara1; k++)
              Utemp(i, k) += v2(j) / v1(j) / v1(j) * pl(i) * vPl(k);

          for (size_t i = 0; i < ncaty; i++)
            ql(i) = fitbdist2(i, j);
          for (size_t k = 0; k < npara2; k++) {
            vQl(k) = 0;
            for (size_t i = 0; i < ncaty; i++)
              vQl(k) += Cr2(j * ncaty + i, k);
          }
          for (size_t i = 0; i < ncaty; i++)
            for (size_t k = 0; k < npara2; k++)
              Vtemp(i, k) += v1(j) / v2(j) / v2(j) * ql(i) * vQl(k);

          for (size_t i = 0; i < ncaty; i++)
            for (size_t k = 0; k < npara1; k++)
              VR(i, k) += wts / v2(j) * ql(i) * vPl(k);

          for (size_t i = 0; i < ncatx; i++)
            for (size_t k = 0; k < npara2; k++)
              US(i, k) += (1 - wts) / v1(j) * pl(i) * vQl(k);
        }

        UP = M1 * Cr1;         // MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1, UP);
        UQ = M2 * Cr2;         // MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2, UQ);
        UPStar = M1 * Cr1Star; // MatrixMultiMatrix(ncatx, ncat1, npara1, M1, Cr1Star, UPStar);
        UQStar = M2 * Cr2Star; // MatrixMultiMatrix(ncaty, ncat2, npara2, M2, Cr2Star, UQStar);
        r2 = M1 * vP121;       // MatrixMultiVector(ncatx, ncat1, M1, vP121, r2);
        s1 = M2 * vP212;       // MatrixMultiVector(ncatx, ncat2, M2, vP212, s1);

        for (size_t i = 0; i < ncatx; i++)
          r(i) = wts * r1(i) + (1 - wts) * r2(i);
        for (size_t i = 0; i < ncaty; i++)
          s(i) = wts * s1(i) + (1 - wts) * s2(i);

        double hx = optimalBandwithXParameter(ncatx, scoresx, r, 1);
        double hy = optimalBandwithXParameter(ncaty, scoresy, s, 1);

        kernelEquate(ncatx, scoresx, r, hx, ncaty, scoresy, s, hy, equatedX);

        for (size_t i = 0; i < ncatx; i++) {
          for (size_t k = 0; k < npara1; k++)
            UR(i, k) = wts * UP(i, k) + (1 - wts) * (UPStar(i, k) - Utemp(i, k));
        }

        for (size_t i = 0; i < ncaty; i++) {
          for (size_t k = 0; k < npara2; k++)
            VS(i, k) = (1 - wts) * UQ(i, k) + wts * (UQStar(i, k) - Vtemp(i, k));
        }

        for (size_t i = 0; i < ncatx; i++) {
          partialFPartialr(ncatx, scoresx, r, hx, Fr, scoresx(i));
          partialFPartialr(ncaty, scoresy, s, hy, Gs, equatedX(i));
          double Gp = kernelPdf(ncaty, scoresy, s, hy, equatedX(i));
          FrUR = Fr * UR; // VectorMultiMatrix(ncatx, npara1, Fr, UR, FrUR);
          GsVR = Gs * VR; // VectorMultiMatrix(ncaty, npara1, Gs, VR, GsVR);
          FrUS = Fr * US; // VectorMultiMatrix(ncatx, npara2, Fr, US, FrUS);
          GsVS = Gs * VS; // VectorMultiMatrix(ncaty, npara2, Gs, VS, GsVS);
          for (size_t j = 0; j < npara1; j++)
            FrUR(j) -= GsVR(j);

          for (size_t j = 0; j < npara2; j++)
            FrUS(j) -= GsVS(j);

          see(i) = std::sqrt(FrUR.squaredNorm() + FrUS.squaredNorm()) / Gp;
        }
      }

      /*--------------------------------------------------------------------------
			KernelEquateSEENEATChn
			
			functionality:

			Computes kernel equating function and SEE for the NEAT design 
			with chained equating method based on procedures in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			author: Tianyou Wang 4/1/2005.
			
			input:
				bivar1      information about bivariate distribution for form X
				bivar2      information about bivariate distribution for form Y
				wts         wts for group taking form X.
		
			output:
				Equatedx    a vector containing the equated score 
				SEE         a vector containing standard error of equating 
		--------------------------------------------------------------------------*/
      void kernelEquateSEENEATChn(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing1,
                                  const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing2,
                                  Eigen::VectorXd& equatedX,
                                  Eigen::VectorXd& see) {
        // int i, j, k, ncat, ncatx, ncaty, ncatv, npara, cu, cv, cuv;
        // long np;
        // double *SEEA, *EquatedxA, *SEEY;
        // double** fitbdist;
        // double *vP, **M, **N, *r, *s, **U, **V;
        // double hv, hy;
        // double **Cr, **B;
        // double *Fr, *Gs, GpQ, HpQ, *FrU, *GsV;
        // double cdfv;
        // double *scoresv, *scoresy;
        // int *intery, *interv;
        // int** cpm2;

        size_t cuv = bivariateLogLinearSmoothing2.numberOfCrossProductMoments;
        size_t cu = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingU;
        size_t cv = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingV;

        Eigen::MatrixXi cpm2(2, cuv);
        for (size_t i = 0; i < cuv; i++)
          for (size_t j = 0; j < 2; j++)
            cpm2(j, i) = bivariateLogLinearSmoothing2.crossProductMoments(i, j);

        size_t ncatx = bivariateLogLinearSmoothing1.numberOfScoresX;
        Eigen::VectorXd seeA(ncatx);
        Eigen::VectorXd equatedxA(ncatx);
        Eigen::VectorXd seeY(ncatx);

        kernelEquateSEESG(bivariateLogLinearSmoothing1, equatedxA, seeA);

        size_t ncaty = bivariateLogLinearSmoothing2.numberOfScoresX;
        size_t ncatv = bivariateLogLinearSmoothing2.numberOfScoresV;
        size_t npara = bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingU +
                       bivariateLogLinearSmoothing2.numberOfDegreesOfSmoothingV +
                       bivariateLogLinearSmoothing2.numberOfCrossProductMoments;
        Eigen::VectorXd intery = cpm2.row(0).cast<double>();
        Eigen::VectorXd interv = cpm2.row(1).cast<double>();
        size_t np = bivariateLogLinearSmoothing2.numberOfExaminees;
        Eigen::VectorXd scoresy(ncaty);
        for (size_t i = 0; i < ncaty; i++)
          scoresy(i) = bivariateLogLinearSmoothing2.minimumRawScoreX +
                       i * bivariateLogLinearSmoothing2.scoreIncrementX;

        Eigen::VectorXd scoresv(ncatv);
        for (size_t i = 0; i < ncatv; i++)
          scoresv(i) = bivariateLogLinearSmoothing2.minimumRawScoreV +
                       i * bivariateLogLinearSmoothing2.scoreIncrementV;
        size_t ncat = ncatv * ncaty;

        Eigen::VectorXd Fr(ncatv);
        Eigen::VectorXd Gs(ncaty);
        Eigen::VectorXd FrU(npara);
        Eigen::VectorXd GsV(npara);

        Eigen::MatrixXd Cr(ncat, npara);
        Eigen::MatrixXd B(npara, ncat);
        Eigen::MatrixXd U(ncatv, npara);
        Eigen::MatrixXd V(ncaty, npara);
        Eigen::MatrixXd fitbdist(ncatv, ncaty);

        for (size_t i = 0; i < ncatv; i++)
          for (size_t j = 0; j < ncaty; j++)
            fitbdist(i, j) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(j, i) / np;

        /*The following code set a natrual design matrix corresponding to a natrual basis of 
					polynomials. This B matrix is the transpose of the design matrix in log-linear model
					First, assign natrual design matrix first for rows corresponding to x */

        for (size_t i = 0; i < cu; i++) {
          for (size_t j = 0; j < ncatv; j++) {
            for (size_t k = 0; k < ncaty; k++) {
              B(i, k * (ncatv) + j) = std::pow(static_cast<double>(j), (i + 1));
            }
          }
        }

        /* then assign values to rows corresponding to y */
        for (size_t i = cu; i < cu + cv; i++) {
          for (size_t j = 0; j < ncatv; j++) {
            for (size_t k = 0; k < ncaty; k++) {
              B(i, k * ncatv + j) = std::pow(static_cast<double>(k), (i - cu + 1));
            }
          }
        }

        /* assign value to the last columns corresponding to the interaction terms */
        for (size_t i = 0; i < cuv; i++) {
          for (size_t j = 0; j < ncatv; j++) {
            for (size_t k = 0; k < ncaty; k++) {
              B(i + cu + cv, k * (ncatv) + j) = std::pow(static_cast<double>(j), interv(i)) *
                                                std::pow(static_cast<double>(k), intery(i));
            }
          }
        }

        Eigen::VectorXd vP(ncat);
        // Eigen::VectorXd vP(ncat);
        Eigen::VectorXd r(ncatv);
        Eigen::VectorXd s(ncaty);
        Eigen::MatrixXd M(ncatv, ncat);
        Eigen::MatrixXd N(ncaty, ncat);

        vPMN(ncatv, ncaty, fitbdist, vP, M, N);
        r = M * vP; // MatrixMultiVector(ncatv, ncat, M, vP, r);
        s = N * vP; // MatrixMultiVector(ncaty, ncat, N, vP, s);

        computeCmatrixGen(ncat, npara, np, B, vP, Cr);
        U = M * Cr; // MatrixMultiMatrix(ncatv, ncat, npara, M, Cr, U);
        V = N * Cr; // MatrixMultiMatrix(ncaty, ncat, npara, N, Cr, V);

        double hv = optimalBandwithXParameter(ncatv, scoresv, r, 1);
        double hy = optimalBandwithXParameter(ncaty, scoresy, s, 1);

        for (size_t i = 0; i < ncatx; i++) {
          double cdfv = kernelCdf(ncatv, scoresv, r, hv, equatedxA(i));
          equatedX(i) = kernelInverseCdf(ncaty, scoresy, s, hy, cdfv);
        }

        for (size_t i = 0; i < ncatx; i++) {
          partialFPartialr(ncatv, scoresv, r, hv, Fr, equatedxA(i));
          partialFPartialr(ncaty, scoresy, s, hy, Gs, equatedX(i));
          double HpQ = kernelPdf(ncatv, scoresv, r, hv, equatedxA(i));
          double GpQ = kernelPdf(ncaty, scoresy, s, hy, equatedX(i));
          FrU = Fr * U; //VectorMultiMatrix(ncatv, npara, Fr, U, FrU);
          GsV = Gs * V; // VectorMultiMatrix(ncaty, npara, Gs, V, GsV);
          for (size_t j = 0; j < npara; j++)
            FrU(j) -= GsV(j);

          seeY(i) = FrU.norm() / GpQ;
          see(i) = std::sqrt(DSQR(HpQ / GpQ * seeA(i)) + DSQR(seeY(i)));
        }
      }

      double DSQR(double a) {
        double result;
        if (a == 0.0) {
          result = 0.0;
        } else {
          result = std::pow(a, 2);
        }

        return result;
      }

      /*--------------------------------------------------------------------------
			KernelEquateNEATChn
			
			functionality:

			Computes kernel equating function and for the NEAT design 
			with chained equating method based on procedures in von Davier, 
			Holland, & Thayer 	(2004). The Kernel Method of Test Equating.
			
			author: Tianyou Wang 4/1/2005.
			
			input:
				bivar1      information about bivariate distribution for form X
				bivar2      information about bivariate distribution for form Y
				wts         wts for group taking form X.
		
			output:
				Equatedx    a vector containing the equated score 
		--------------------------------------------------------------------------*/
      void kernelEquateNEATChn(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing1,
                               const EquatingRecipes::Structures::BivariateLogLinearSmoothing& bivariateLogLinearSmoothing2,
                               Eigen::VectorXd& equatedX) {
        // int i, j, ncat, ncatx, ncaty, ncatv;
        // long np;
        // double* EquatedxA;
        // double** fitbdist;
        // double *vP, **M, **N, *r, *s;
        // double hv, hy;
        // double cdfv;
        // double *scoresv, *scoresy;

        size_t ncatx = bivariateLogLinearSmoothing1.numberOfScoresX;
        Eigen::VectorXd equatedxA(ncatx);

        kernelEquateSG(bivariateLogLinearSmoothing1, equatedxA);

        size_t ncaty = bivariateLogLinearSmoothing2.numberOfScoresX;
        size_t ncatv = bivariateLogLinearSmoothing2.numberOfScoresV;
        size_t np = bivariateLogLinearSmoothing2.numberOfExaminees;

        Eigen::VectorXd scoresy(ncaty);
        for (size_t i = 0; i < ncaty; i++)
          scoresy(i) = bivariateLogLinearSmoothing2.minimumRawScoreX + i * bivariateLogLinearSmoothing2.scoreIncrementX;

        Eigen::VectorXd scoresv(ncatv);
        for (size_t i = 0; i < ncatv; i++)
          scoresv(i) = bivariateLogLinearSmoothing2.minimumRawScoreV + i * bivariateLogLinearSmoothing2.scoreIncrementV;

        size_t ncat = ncatv * ncaty;

        Eigen::MatrixXd fitbdist(ncatv, ncaty);

        for (size_t i = 0; i < ncatv; i++)
          for (size_t j = 0; j < ncaty; j++)
            fitbdist(i, j) = bivariateLogLinearSmoothing2.fittedBivariateFreqDist(j, i) / np;

        Eigen::VectorXd vP(ncat);
        Eigen::VectorXd r(ncatv);
        Eigen::VectorXd s(ncaty);
        Eigen::MatrixXd M(ncatv, ncat);
        Eigen::MatrixXd N(ncaty, ncat);

        vPMN(ncatv, ncaty, fitbdist, vP, M, N);
        r = M * vP; // MatrixMultiVector(ncatv, ncat, M, vP, r);
        s = N * vP; // MatrixMultiVector(ncaty, ncat, N, vP, s);

        double hv = optimalBandwithXParameter(ncatv, scoresv, r, 1);
        double hy = optimalBandwithXParameter(ncaty, scoresy, s, 1);

        for (size_t i = 0; i < ncatx; i++) {
          double cdfv = kernelCdf(ncatv, scoresv, r, hv, equatedxA(i));
          equatedX(i) = kernelCdf(ncaty, scoresy, s, hy, cdfv);
        }
      }
    };
  } // namespace Implementation
} // namespace EquatingRecipes

#endif