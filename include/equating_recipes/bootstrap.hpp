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

#include <random>
#include <Eigen/Core>

#include <equating_recipes/beta_binomial.hpp>
#include <equating_recipes/cg_equipercentile_equating.hpp>
#include <equating_recipes/cg_no_smoothing.hpp>
#include <equating_recipes/log_linear_equating.hpp>
#include <equating_recipes/rg_and_sg_equating.hpp>
#include <equating_recipes/score_statistics.hpp>
#include <equating_recipes/structures/beta_binomial_smoothing.hpp>
#include <equating_recipes/structures/bivariate_log_linear_smoothing.hpp>
#include <equating_recipes/structures/bivariate_statistics.hpp>
#include <equating_recipes/structures/bootstrap_equated_raw_score_results.hpp>
#include <equating_recipes/structures/bootstrap_equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/equated_scaled_scores_results.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/structures/univariate_statistics.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class Bootstrap {
  public:
    void runBootstrap(EquatingRecipes::Structures::PData& pData,
                      const size_t& numberOfReplications,
                      long& idum,
                      EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& bootstrapEquatedRawScoreResults,
                      EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& bootstrapEquatedScaledScoresResults) {}

  private:
    void initializeEquatedRawScoreResults(EquatingRecipes::Structures::PData& pData,
                                          EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& bootstrapEquatedRawScoreResults) {}

    void bootstrapUnivariateStatistics(EquatingRecipes::Structures::UnivariateStatistics& xUnivariateStatistics,
                                       long& idum,
                                       const size_t& replicationNumber,
                                       EquatingRecipes::Structures::UnivariateStatistics& xBootstrapUnivariateStatistics) {}

    void bootstrapBivaraiteStatistics(EquatingRecipes::Structures::BivariateStatistics& xvBivariateStatistics,
                                      long& idum,
                                      const size_t& replicationNumber,
                                      EquatingRecipes::Structures::BivariateStatistics& s) {}

    void bootstrapAccumulateEquatedRawScores(EquatingRecipes::Structures::PData& pData,
                                             EquatingRecipes::Structures::EquatedRawScoreResults& b,
                                             EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& t) {}

    void bootstrapStandardErrorsEquatedRawScores(EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::BootstrapEquatedRawScoreResults& t) {}

    // void Print_Boot_se_eraw(FILE *fp, char tt[], struct PDATA *inall,
    //                         struct ERAW_RESULTS *r,
    //                         struct BOOT_ERAW_RESULTS *t, int mdiff);

    void bootstraptInitializeEquatedScaledScores(EquatingRecipes::Structures::PData& pData,
                                                 EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& u) {}

    void bootstrapAccumulateEquatedScaledScores(EquatingRecipes::Structures::PData& pData,
                                                EquatingRecipes::Structures::EquatedScaledScoresResults& s,
                                                EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& u) {}

    void bootstrapStandardErrorsEquatedScaledScores(EquatingRecipes::Structures::PData& pData,
                                                    EquatingRecipes::Structures::BootstrapEquatedScaledScoresResults& u) {}

    // void Print_Boot_se_ess(FILE *fp, char tt[], struct PDATA *inall,
    //                        struct ESS_RESULTS *s,
    //                        struct BOOT_ESS_RESULTS *u, int mdiff);

    void parametricBootstraptUnivariateBetaBinomial(EquatingRecipes::Structures::BetaBinomialSmoothing& x,
                                                    long& idum,
                                                    size_t& replicationNumber,
                                                    EquatingRecipes::Structures::BetaBinomialSmoothing& btx) {}

    void parametricBootstraptUnivariateLogLinear(EquatingRecipes::Structures::UnivariateLogLinearSmoothing& x,
                                                 long& idum,
                                                 size_t& replicationNumber,
                                                 EquatingRecipes::Structures::UnivariateLogLinearSmoothing& btx) {}

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
    void parametricBootstraptBivariate(const EquatingRecipes::Structures::BivariateLogLinearSmoothing& xv,
                                       long& idum,
                                       size_t& replicationNumber,
                                       EquatingRecipes::Structures::BivariateLogLinearSmoothing& bootstrapXV) {
      // int i, j,
      //     nsx = xv->nsx,
      //     nsv = xv->nsv,
      //     ns = nsx * nsv, /* total number of score categories */
      //     np = xv->num_persons;
      // float r;
      // double* ptr; /* pointer variable */

      size_t totalNumberOfScores = xv.numberOfScoresX * xv.numberOfScoresV;

      if (replicationNumber == 1) {
        bootstrapXV.isInternalAnchor = xv.isInternalAnchor;
        bootstrapXV.numberOfExaminees = xv.numberOfExaminees;
        bootstrapXV.numberOfScoresX = xv.numberOfScoresX;
        bootstrapXV.minimumRawScoreX = xv.minimumRawScoreX;
        bootstrapXV.scoreIncrementX = xv.scoreIncrementX;
        bootstrapXV.numberOfScoresV = xv.numberOfScoresV;
        bootstrapXV.minimumRawScoreV = xv.minimumRawScoreV;
        bootstrapXV.scoreIncrementV = xv.scoreIncrementV;
        bootstrapXV.totalNumberOfScores = xv.totalNumberOfScores;

        bootstrapXV.fittedBivariateFreqDist.resize(xv.numberOfScoresX, xv.numberOfScoresV);

        bootstrapXV.fittedFrequencesX.resize(xv.numberOfScoresX);
        bootstrapXV.fittedRawScoreDensityX.resize(xv.numberOfScoresX);
        bootstrapXV.fittedRawScoreCumulativeRelativeFreqDistX.resize(xv.numberOfScoresX);
        bootstrapXV.fittedRawScorePercentileRankDistX.resize(xv.numberOfScoresX);

        bootstrapXV.fittedFrequencesV.resize(xv.numberOfScoresV);
        bootstrapXV.fittedRawScoreDensityV.resize(xv.numberOfScoresV);
        bootstrapXV.fittedRawScoreCumulativeRelativeFreqDistV.resize(xv.numberOfScoresV);
        bootstrapXV.fittedRawScorePercentileRankDistV.resize(xv.numberOfScoresV);

        bootstrapXV.cumulativeRelativeFreqDistRowMajorVector.resize(totalNumberOfScores);
        bootstrapXV.fittedBivariateRelativeFreqDistXV.resize(xv.numberOfScoresX, xv.numberOfScoresV);
      }

      bootstrapXV.cumulativeRelativeFreqDistRowMajorVector.setZero();

      /* In the following code, xv->crfd_vector_bfd[] is the crfd  for the 
        fitted actual data, while btxv->crfd_vector_bfd[] is the bootstrap fd 
      at end of for loop */
      std::mt19937_64 gen = EquatingRecipes::Utilities::getSeedEngine();
      std::uniform_real_distribution<> distrib = EquatingRecipes::Utilities::getUniformRealDistribution();

      for (size_t examineeIndex = 1; examineeIndex <= xv.numberOfExaminees; examineeIndex++) {
        double randomValue = EquatingRecipes::Utilities::getUniformDoubleRandomNumber(distrib,
                                                                                      gen);

        for (size_t index = 0; index < totalNumberOfScores; index++) {
          if (randomValue < xv.cumulativeRelativeFreqDistRowMajorVector(index)) {
            bootstrapXV.cumulativeRelativeFreqDistRowMajorVector(index)++;
            break;
          }
        }
      }

      /* Here btxv->crfd_vector_bfd[] is bootstrap frequencies in vector format; 
      next statements put it in matrix format and store it in btxv->bfd[][] */

      ptr = btxv->crfd_vector_bfd;
      for (i = 0; i < nsx; i++)
        for (j = 0; j < nsv; j++)
          btxv->bfd[i][j] = *ptr++;

      /* get row and col marginal fd, rfd, crfd, and PR */

      for (j = 0; j < nsv; j++)
        btxv->fd_v[j] = 0.; /* initialization */

      for (i = 0; i < nsx; i++) {
        btxv->fd_x[i] = 0.; /* initialization */
        for (j = 0; j < nsv; j++) {
          btxv->fd_x[i] += btxv->bfd[i][j]; /* row fd */
          btxv->fd_v[j] += btxv->bfd[i][j]; /* col fd */
        }
        btxv->density_x[i] = btxv->fd_x[i] / np; /* row density */
      }
      for (j = 0; j < nsv; j++)
        btxv->density_v[j] = btxv->fd_v[j] / np; /* col density */

      cum_rel_freqs(0, nsx - 1, 1, btxv->density_x, btxv->crfd_x); /* row crfd */
      for (i = 0; i < nsx; i++)
        btxv->prd_x[i] = perc_rank(0, nsx - 1, 1, btxv->crfd_x, i); /* row prd */

      cum_rel_freqs(0, nsv - 1, 1, btxv->density_v, btxv->crfd_v); /* col crfd */
      for (j = 0; j < nsv; j++)
        btxv->prd_v[j] = perc_rank(0, nsv - 1, 1, btxv->crfd_v, j); /* col prd */

      /* following code gets rel fd version of btxv->bfd[]; i.e., 
     btxv->brfd[][] is the bootstrap rel freq biv dist for x by v; 
	 needed for FEorMFE_EE() */

      for (i = 0; i < nsx; i++)
        for (j = 0; j < nsv; j++)
          btxv->brfd[i][j] = btxv->bfd[i][j] / np;
    }
  };
} // namespace EquatingRecipes

#endif