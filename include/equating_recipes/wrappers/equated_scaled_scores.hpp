/* 
  From Source: ERutilities.h, ERutilities.c
  Original Struct: 
  Description: 
*/

#ifndef WRAPPERS_EQUATED_SCALED_SCORES_HPP
#define WRAPPERS_EQUATED_SCALED_SCORES_HPP

#include <algorithm>
#include <limits>
#include <map>
#include <string>
#include <Eigen/Core>

#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/utilities.hpp>

/*
    Input:
    
      uses the following from struct PDATA inall
        nm        = number of methods 
        names[][] = names of methods
        min       = min raw score for x 
        max       = max raw score for x
        inc       = increment for x         
        fdx[]     = fd for new form x 
        rep       = replication number for bootstrap; = 0 if no bootstrap
        
      uses the following from struct ERAW_RESULTS r
        eraw[][] = equated raw scores
        
      minp    = minimum raw score for base form y conversion yct
      maxp    = maximum raw score for base form y conversion yct
	  incp    = increment in raw scores in base form y conversion yct
      nameyct = name of file with raw-to-scale score conversion for Y
      round: if round=i --> round to i-th decimal place
      lprss = lowest possible rounded scale score
      hprss = highest possible rounded scale score
       
    Output: 
      
      s = struct ESS_RESULTS:
            essu[][] = unrounded equated scale scores
            essr[][] = rounded equated scale scores  
            mtsu[][] = moments for equated unrounded scale scores 
            mtsr[][] = moments for equated rounded scale scores 
        
      NOTE: The following are added to struct PDATA inall: minp, maxp, nameyct,
            round, lprss, hprss, and **yct
        
      NOTE: struct ESS_RESULTS must be different struct's for actual 
            equating and bootstrap use of Wrapper_ESS()
        
      NOTE:       
      yct[][] = conversion table for Y; first dimension ranges from
	            0 to (nscores(maxp,minp,inc)+1
				(see comments for ReadSSConvTableForY() for details)

  Function calls other than C or NR utilities: 
    nscores()
    ReadSSConvTableForY()
    Equated_ss()
    MomentsFromFD()
                                               
  R. L. Brennan

  Date of last revision: 6/30/08                                     
*/

namespace EquatingRecipes {
  namespace Wrappers {
    struct EquatedScaledScores {
      // struct PDATA *inall, struct ERAW_RESULTS *r, double minp,
      //           double maxp, double incp, char *nameyct, int round,
      //           int lprss, int hprss, struct ESS_RESULTS *s
      EquatingRecipes::Structures::EquatedScaledScoresResults run(EquatingRecipes::Structures::PData& pData,
                                                                  const EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults,
                                                                  const double& minimumRawScoreYct,
                                                                  const double& maximumRawScoreYct,
                                                                  const double& scoreIncrementYct,
                                                                  const size_t& roundToNumberOfDecimalPlaces,
                                                                  const int& lowestObservableRoundedScaledScore,
                                                                  const int& highestObservableRoundedScaledScore,
                                                                  const EquatingRecipes::Structures::RawToScaledScoreTable& rawToScaledScoreTable) {
        EquatingRecipes::Structures::EquatedScaledScoresResults results;

        size_t numberOfScores = EquatingRecipes::Utilities::numberOfScores(equatedScaledScoresInput.minimumRawScoreX,
                                                                           equatedScaledScoresInput.maximumRawScoreX,
                                                                           equatedScaledScoresInput.scoreIncrementX);

        pData.minimumRawScoreYct = minimumRawScoreYct;
        pData.maximumRawScoreYct = maximumRawScoreYct;
        pData.scoreIncrementYct = scoreIncrementYct;
        pData.roundToNumberOfDecimalPlaces = roundToNumberOfDecimalPlaces;
        pData.lowestObservableRoundedScaledScore = lowestObservableRoundedScaledScore;
        pData.highestObservableRoundedScaledScore = highestObservableRoundedScaledScore;
        pData.rawToScaledScoreTable = rawToScaledScoreTable;
        
        if (pData.bootstrapReplicationNumber <= 1) {
          // actual equating or 1st bootstrap replication
          results.unroundedEquatedScaledScores.setZero(pData.methods.size());
          results.roundedEquatedScaledScores.setZero(pData.methods.size());

          results.unroundedEquatedScaledScoreMoments.setZero(pData.methods.size(), 4);
          results.roundedEquatedScaledScoreMoments.setZero(pData.methods.size(), 4);
        }

        for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
          
        }
            // for(i=0;i<=inall->nm-1;i++)
            //   Equated_ss(inall->min,inall->max,inall->inc,minp,maxp,incp,r->eraw[i],
            //             inall->yct,inall->round,inall->lprss,inall->hprss,s->essu[i],s->essr[i]);

            // /* compute moments:  Note that when inc==1, essu[*][min-minp+1] is the
            //   unrounded scale score associated with fdx[0], where fdx[0]
            //   is associated with scores ranging from min to max.
            //   Recall that essu[*][0] is the unrounded scale score
            //   associated with minp-inc/2 = minp-.5 when inc = 1.  In this example,
            //   min-minp+1 = loc(min,minp,inc) + 1

            // */
            // if (inall->fdx != NULL) {
            // for(i=0;i<=inall->nm-1;i++){
            //   MomentsFromFD(inall->min, inall->max, inall->inc,
            //                 s->essu[i],inall->fdx, s->mtsu[i]);
            //   MomentsFromFD(inall->min, inall->max, inall->inc,
            //                 s->essr[i],inall->fdx, s->mtsr[i]);
            // }
            // }

            return results;
      }
    };
  } // namespace Wrappers
} // namespace EquatingRecipes

#endif