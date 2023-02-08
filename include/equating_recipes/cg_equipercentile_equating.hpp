/* 

CG_EquiEquate.c - contains functions for computing equating results for 
          (a) frequency estimation (FE) and Braun-Holland 
              under FE (if requested);
          (b) modified frequency estimation (MFE) (Wang & Brennan, 2006) 
              and Braun-Holland under MFE (if requested); and
          (c) chained equipercentile equating

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

#ifndef CG_EQUIV_EQUATING_HPP
#define CG_EQUIV_EQUATING_HPP

#include <string>
#include <vector>
#include <Eigen/Core>

#include <equating_recipes/utilities.hpp>
#include <equating_recipes/structures/all_structures.hpp>

namespace EquatingRecipes {
  struct CGEquipercentileEquating {
    /*
      Computes results for common-item (CI) equipercentile equating for EITHER
      (a) frequency estimation (FE) and Braun-Holland under FE (if requested); OR
      (b) modified frequency estimation (MFE) (Wang & Brennan, 2006) and 
          Braun-Holland under MFE (if requested).

      Input
        w1        weight for pop 1
        internal  internal anchor if non-zero; external anchor if zero
      nsv       number of score categories for v
        minv      minimum score for v
        maxv      maximum score for v      
        nsx       number of score categories for x
        minx      minimum score for x
        maxx      maximum score for x
        nsy       number of score categories for y
        miny      minimum score for y
        maxy      maximum score for y
        inc       increment--assumed to be the same for x, y, and v 
        bxvin[][] biv rel freq dist for x (rows) and v -- input 
        byvin[][] biv rel freq dist for y (rows) and v -- input 
        rv1       reliability for common items in pop 1 (used for MFE)
        rv2       reliability for common items in pop 2 (used for MFE)

      Output
        fxs[]      rel freq dist for x in syn pop  
        gys[]      rel freq dist for y in syn pop  
        eraw[]     equipercentile equated raw scores 
        a[]        slope for Braun-Holland Method 
        b[]        intercept for Braun-Holland Method 
        erawBH[]   Braun-Holland linear equated raw scores 

      NOTES: (1) Results are for FE if rv1 == 0 or rv2 == 0;
                  otherwise results are for MFE
              (2) Assumes space allocated for fxs, gys, and eraw
                  before function called
              (3) Braun-Holland results (a, b, erawBH) provided
                  if a != NULL and b !=NULL and erawBH != NULL
              (4) If Braun-Holland results are desired, then before
                  function is called, space must be allocated for erawBH;
                  also, a and b variables must be declared and their
                  addresses passed
              (5) Braun-Holland results are provided along with 
                  equipercentile results; Braun Holland results not
                  provided in isolation (refer to 2, above)

      Function calls other than C or NR utilities:
      SyntheticDensities()
        cum_rel_freqs()
        perc_rank()
        EquiEquate()
        BH_LinEq()
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */

    // w1        weight for pop 1
    //   internal  internal anchor if non-zero; external anchor if zero
    // nsv       number of score categories for v
    //   minv      minimum score for v
    //   maxv      maximum score for v
    //   nsx       number of score categories for x
    //   minx      minimum score for x
    //   maxx      maximum score for x
    //   nsy       number of score categories for y
    //   miny      minimum score for y
    //   maxy      maximum score for y
    //   inc       increment--assumed to be the same for x, y, and v
    //   bxvin[][] biv rel freq dist for x (rows) and v -- input
    //   byvin[][] biv rel freq dist for y (rows) and v -- input
    //   rv1       reliability for common items in pop 1 (used for MFE)
    //   rv2       reliability for common items in pop 2 (used for MFE)

    // Output
    //   fxs[]      rel freq dist for x in syn pop
    //   gys[]      rel freq dist for y in syn pop
    //   eraw[]     equipercentile equated raw scores
    //   a[]        slope for Braun-Holland Method
    //   b[]        intercept for Braun-Holland Method
    //   erawBH[]   Braun-Holland linear equated raw scores

    EquatingRecipes::Structures::CGEquipercentileEquatingResults equipercentileEquating(const double& population1Weight,
                                                                                        const bool& isInternalAnchor,
                                                                                        const size_t& numberOfScoresV,
                                                                                        const double& mininumScoreV,
                                                                                        const double& maximumScoreV,
                                                                                        const size_t& numberOfScoresX,
                                                                                        const double& mininumScoreX,
                                                                                        const double& maximumScoreX,
                                                                                        const size_t& numberOfScoresY,
                                                                                        const double& mininumScoreY,
                                                                                        const double& maximumScoreY,
                                                                                        const double& scoreIncrement,
                                                                                        const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsXV,
                                                                                        const EquatingRecipes::Structures::BivariateStatistics& bivariateStatisticsYV,
                                                                                        const double& reliabilityCommonItemsPopulation1 = 0.0,
                                                                                        const double& reliabilityCommonItemsPopulation2 = 0.0) {
      EquatingRecipes::Structures::CGEquipercentileEquatingResults results;

      return results;
    }

    void syntheticDensities(const double& population1Weight,
                            const bool& isInternalAnchor,
                            const size_t& numberOfScoresV,
                            const size_t& numberOfScoresX,
                            Eigen::MatrixXd& bivariateRelativeFreqDistXV,
                            const size_t& numberOfScoresY,
                            Eigen::MatrixXd& bivariateRelativeFreqDistYV,
                            const double& reliabilityCommonItemsPopulation1,
                            const double& reliabilityCommonItemsPopulation2,
                            Eigen::VectorXd& syntheticPopulationRelativeFreqDistX,
                            Eigen::VectorXd& syntheticPopulationRelativeFreqDistY) {}

    void mixSmooth(const size_t& numberOfScoresV,
                   const size_t& numberOfScoresX,
                   const double& mixingProportion,
                   const Eigen::MatrixXd& observedBivariateRelativeFreqDistXV,
                   Eigen::MatrixXd& smoothedBivariateProbabilitiesXV,
                   Eigen::VectorXd& smoothedMarginalProbabilitiesX,
                   Eigen::VectorXd& smoothedMarginalProbabilitiesV) {}

    /*
      Conditional bivariate distribution of x given v. (Obviously applies
      as well to conditional bivariate distribution of y given v.)

      Input
        nsv    number of categoies for v
        nsx    number of categories for x 
        bxv    bivariate relative freq distribution of x (rows) and v (cols)
        hv     mariginal rel freq distribution of v
    
      Output
        bxv	   conditional distribution of x given v (or y given v)
                  where v is rows

      Function calls other than C or NR utilities: None
                                                    
      B. A. Hanson with updates by R. L. Brennan

      Date of last revision: 6/30/08 
    */
    Eigen::MatrixXd conditionBivariateDistribution(const size_t& numberOfScoresV,
                                                   const size_t& numberOfScoresX,
                                                   const Eigen::MatrixXd& bivariateRelativeFreqDistXV,
                                                   const Eigen::VectorXd& marginalRelativeFreqDistV) {
      Eigen::MatrixXd condBivDist;
      return condBivDist;
    }

    /*
      Braun-Holland linear equating (based on FE or MFE systhetic densities

      Input
        minx  min score for x
        maxx  max score for x
        miny  min score for y
        maxy  max score for y
        inc   increment for both x and y
        fxs   synthetic density for x
        gsy   synthetic density for y

      Output
        a    slope
        b    intercept

      Function calls other than C or NR utilities:
        MomentsFromRFD()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void braunHollandLinearEquating(const size_t& minimumScoreX,
                                    const size_t& maximumScoreX,
                                    const size_t& minimumScoreY,
                                    const size_t& maximumScoreY,
                                    const double& scoreIncrement,
                                    const Eigen::VectorXd& syntheticDensityX,
                                    const Eigen::VectorXd& syntheticDensityY,
                                    double& slope,
                                    double& intercept) {}

    /*
      For MFE method, modify conditional bivariate distribution of
      x given v with x as rows.  Description here worded 
      in terms of x given v, and code written in same manner, but 
      obviously function applies to y given v, as well.

      Input
        internal internal anchor if non-zero; external anchor if zero
        nsv      number of score categories for v
        nsx      number of score categories for x
        rv1      rel of common items in pop 1
        rv2      rel of common items in pop 2
        muv1     mean of v in pop 1
        muv2     mean of v in pop 2
        bxv      bivariate distribution of x given v in pop 1 with x as rows

      Output
        bxv      modified biv dist of x given v in pop 2 with x as rows

      NOTE:  There is a complexity here for an internal anchor. Let 
              nsu = nsx - nsv = # categories for non-equating items. Recall that
              x = 0,...,nsx-1 and v = 0,..., nsv-1.  In
              the bxv matrix there is a structural 0 whenever
              (x<v || x>nsu-1+v), which creates linear interpolation problems
              that are circumvented here by collapsing the xv matrix to a 
              uv matrix, doing the linear interpolation, and then expanding
              the uv back to an xv matrix. (Note that there is a one-to-one 
              relationship between elements of original xv matrix and 
              collapsed uv matrix.)

      Function calls other than C or NR utilities: 
        interpolate()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    void modifiedConditionalBivariateDistribution(const bool& isInternalAnchor,
                                                  const size_t& numberOfScoresV,
                                                  const size_t& numberOfScoresX,
                                                  const double& reliabilityCommonItemsPopulation1,
                                                  const double& reliabilityCommonItemsPopulation2,
                                                  const double& population1MeanV,
                                                  const double& population2MeanV,
                                                  const Eigen::MatrixXd& bivariateDistributionXVPopulation1,
                                                  Eigen::MatrixXd& bivariateDistributionXVPopulation2) {
    }

    /*
      Chained equipercentile equating using the composed function
      discussed in Kolen and Brennan (2004, pp. 145-147)

      Input

        nsx    = number of score categories for X in pop 1
        prx1[] = PR for X in pop 1
        minv   = min score for V in both pop 1 and pop 2
        maxv   = max score for V in both pop 1 and pop 2
        incv   = increment for V in both pop 1 and pop 2
        nsv    = # of score categories for V in both pops
        Gv1[]  = crfd for V in pop 1

        miny   = min score for Y in pop 2
        incy   = increment for Y in pop 2
        nsy    = number of score categories for Y in pop 2
        Gy2[]  = crfd for Y in pop 2
        Gv2[]  = crfd for V in pop 2

      Output

        eraw[] = chained equipercentile Y-equivalents of X
        
      Notes.  (a) It is assumed that space for eraw[] allocated
                  prior to function call
              (b) minv, maxv, incv, and nsv are 
                  necessarily the same for both pops

      Function calls other than C or NR utilities:
        EquiEquate()
        perc_rank()
                                                    
      R. L. Brennan

      Date of last revision: 6/30/08 
    */
    Eigen::VectorXd ChainedEquipercentileEquating(const size_t& numberOfScoresXPopulation1,
                                                  const Eigen::VectorXd& percentileRanksXPopulation1,
                                                  const double& minimumScoreV,
                                                  const double& maximumScoreV,
                                                  const double& scoreIncrement,
                                                  const size_t& numberOfScoresXPopulation2,
                                                  const Eigen::VectorXd& cumlativeRelativeFreqDistYPopulation2,
                                                  const Eigen::VectorXd& cumlativeRelativeFreqDistVPopulation2) {
      Eigen::VectorXd equatedRawScores;
      return equatedRawScores;
    }

    // void Print_SynDens(FILE *fp, char tt[], struct PDATA *inall, struct ERAW_RESULTS *r);
  };
} // namespace EquatingRecipes

#endif