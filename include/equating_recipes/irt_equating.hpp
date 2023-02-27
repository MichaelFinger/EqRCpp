/*
  IRT True Score Equating for Mixed-Format Tests

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

  For details about partial derivatives, refer to the following report:

  Kim, S., & Kolen, M.J. (2005). Methods for obtaining a common scale under
      unidimensional IRT models: A technical review and further extensions
      (Iowa Testing Programs Occasional Paper, No. 52). The University of
      Iowa.

  Note: All pointers are of 0-offset unless otherwise indicated.
*/

#ifndef IRT_EQUATING_HPP
#define IRT_EQUATING_HPP

#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <fmt/core.h>

#include <equating_recipes/irt_scale_transformation.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/equated_raw_score_results.hpp>
#include <equating_recipes/structures/irt_equating_results.hpp>
#include <equating_recipes/structures/irt_fitted_distribution.hpp>
#include <equating_recipes/structures/irt_input.hpp>
#include <equating_recipes/structures/irt_method.hpp>
#include <equating_recipes/structures/irt_model.hpp>
#include <equating_recipes/structures/irt_scale_transformation_control.hpp>
#include <equating_recipes/structures/item_specification.hpp>
#include <equating_recipes/structures/method.hpp>
#include <equating_recipes/structures/p_data.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  class IRTEquating {
  public:
    /*
      Wrapper for IRT equating, including both true score equating and observed score equating
      
      Input:
      
        design:         'R' (Random groups)
                        'S' (Single group)
                        'C' (CINEG)
        method:         'T' = IRT true score equating
                        'O' = IRT observed score equating
                        'A' = true score + observed score
        w1:             weight for X in the synthetic population
        ItemNewFile[]:  name of file containing item parameters for new form X converted
                        to old-form Y scale
        ItemOldFile[]:  name of file containing item parameters for old form Y
        DistNewFile[]:  name of file containing quadrature ability scale for new form X
                        converted to scale of old-form Y
        DistOldFile[]:  name of file containing quadrature ability scale for old form Y 
        NewItems:       struct for parameters for X items (on scale of Y)
        OldItems:       struct for parameters for Y items 
        NewForm:        struct for IRT fitted dist for X (on scale of Y)
        OldForm:        struct for IRT fitted dist for Y 
        RawEq:          struct for IRT true and observed equivalents
        StInfo:         struct containing quadrature distributions 
        NewFD:          new form actual relative frequency distribution;
                        if NULL then no moments provided based on new form
                        actual freq dist.

      Note: if only IRT true score equating is not being requested,
            DistNewFile[] and DistOldFile[] are ignored; in this case, it is
            best to set them to NULL can be set to NULL.  Also, w1 is
            ignored.   
      
      Output:
        
        irtall:      IRT_INPUT structure that contains all IRT results
        pinall:      PDATA structure that contains "all" input    
        r:           ERAW_RESULTS structure that contains raw equivalents and moments     

                    if Method =='T' eraw[0] and mts[0] contain IRT true score eq results

                    if Method =='O' eraw[0] and mts[0] contain IRT obs score eq results

                    if Method =='A' eraw[0] and mts[0] contain IRT true score eq results;
                                    eraw[1] and mts[1] contain IRT obs score eq results
                                  

      NOTE:  No capability built in for bootstrapping
                                                
      Function calls other than C or NR utilities:
        RawFitMem()
        RawEqResultsMem()
        trueScoreEq()     
        IRTmixObsEq()    
        MomentsFromFD()
        MomentsFromRFD()
                                                  
      Author's: R. L. Brennan and T. D. Wang
      Date of last revision 9/15/08
    */
    void runIRTEquating(const EquatingRecipes::Structures::Design& design,
                        const EquatingRecipes::Structures::IRTMethod& irtMethod,
                        const double& w1,
                        const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                        const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems,
                        EquatingRecipes::Structures::IRTFittedDistribution& newForm,
                        EquatingRecipes::Structures::IRTFittedDistribution& oldForm,
                        EquatingRecipes::Structures::IRTEquatingResults& irtEquatingResults,
                        EquatingRecipes::Structures::IRTScaleTransformationData& stInfo,
                        const Eigen::VectorXd& newFormFrequencyDistribution,
                        EquatingRecipes::Structures::IRTInput& irtall,
                        EquatingRecipes::Structures::PData& pData,
                        EquatingRecipes::Structures::EquatedRawScoreResults& equatedRawScoreResults) {
      std::vector<std::string> names {"IRT Tr Scr", "IRT Ob Scr"}; /* method names */

      /* The following code segment reads quadrature distributions
      only if IRT observed score equating is requested, since
      quadrature distributions are not required for IRT
      true score equating */

      if (irtMethod == EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE ||
          irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE) {
        // TODO: import quadrature points and weights
      }

      /* Allocate memory */
      initializeRawEqResults(newItems,
                             stInfo,
                             irtEquatingResults);

      initializeRawFitMem(oldItems,
                          false,
                          stInfo,
                          oldForm);

      initializeRawFitMem(newItems,
                          true,
                          stInfo,
                          newForm);

      irtall.method = irtMethod;
      irtall.irtScaleTransformationData = stInfo;
      irtall.newItems = newItems;
      irtall.oldItems = oldItems;
      irtall.irtEquatingResults = irtEquatingResults;
      irtall.newFormIRTFittedDistribution = newForm;
      irtall.oldFormIRTFittedDistribution = oldForm;

      pData.bootstrapReplicationNumber = 0;
      pData.design = design;
      pData.weightSyntheticPopulation1 = w1;
      pData.mininumScoreX = stInfo.minimumRawScoreNewForm;
      pData.maximumScoreX = stInfo.maximumRawScoreNewForm;
      pData.scoreIncrementX = stInfo.rawScoreIncrementNewForm;

      switch (irtMethod) {
        case EquatingRecipes::Structures::IRTMethod::TRUE_SCORE:
          pData.methods.push_back(names[0]);
          break;

        case EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE:
          pData.methods.push_back(names[1]);
          break;

        case EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE:
          pData.methods.push_back(names[0]);
          pData.methods.push_back(names[1]);
          break;

        default:
          break;
      }

      /* allocation and assignments for r (equivalents and moments) */
      equatedRawScoreResults.equatedRawScores.resize(pData.methods.size(),
                                                     newForm.numberOfRawScoreCategories);

      equatedRawScoreResults.equatedRawScoreMoments.resize(pData.methods.size(),
                                                           4);

      if (irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_SCORE) {
        trueScoreEquating(stInfo,
                          newItems,
                          oldItems,
                          newForm.numberOfRawScoreCategories,
                          newForm.rawScores,
                          irtEquatingResults.unroundedEquatedTrueScore,
                          irtEquatingResults.thetaEquivalentFormXScore,
                          irtEquatingResults.minimumTrueScoreNewForm,
                          irtEquatingResults.minimumTrueScoreOldForm);

        equatedRawScoreResults.equatedRawScores.row(0) = irtEquatingResults.unroundedEquatedTrueScore;
      }

      /* IRT observed score equating.  First code segment declares an error and exits
      if IRT observed score equating cannot be done because 
      w1<0 || w1>1 || DistNewFile == NULL || DistOldFile == NULL*/

      if (irtMethod == EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE ||
          irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE) {
        if (w1 < 0 || w1 > 1) {
          throw std::runtime_error("\nIRT observed-score equating cannot be performed because\nw1<0 || w1>1 || DistNewFile == NULL || DistOldFile == NULL \n");
        }
      }

      if (irtMethod == EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE) {
        irtMixObsEq(stInfo,
                    newItems,
                    oldItems,
                    w1,
                    1.0 - w1,
                    newForm,
                    oldForm,
                    irtEquatingResults);

        equatedRawScoreResults.equatedRawScores.row(0) = irtEquatingResults.unroundedEquatedObservedScore;
      }

      /* Both IRT true-score and observed-score equating */
      if (irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE) {
        trueScoreEquating(stInfo,
                          newItems,
                          oldItems,
                          newForm.numberOfRawScoreCategories,
                          newForm.rawScores,
                          irtEquatingResults.unroundedEquatedTrueScore,
                          irtEquatingResults.thetaEquivalentFormXScore,
                          irtEquatingResults.minimumTrueScoreNewForm,
                          irtEquatingResults.minimumTrueScoreOldForm);

        irtMixObsEq(stInfo,
                    newItems,
                    oldItems,
                    w1,
                    1.0 - w1,
                    newForm,
                    oldForm,
                    irtEquatingResults);

        equatedRawScoreResults.equatedRawScores.row(0) = irtEquatingResults.unroundedEquatedTrueScore;
        equatedRawScoreResults.equatedRawScores.row(1) = irtEquatingResults.unroundedEquatedObservedScore;
      }

      /* get moments based on actual frequencies for group that took new form*/

      if (newFormFrequencyDistribution.size() >= 1) {
        irtall.newFormFrequencyDistribution = newFormFrequencyDistribution;
        pData.scoreFrequenciesX = newFormFrequencyDistribution;

        for (size_t methodIndex = 0; methodIndex < pData.methods.size(); methodIndex++) {
          EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(equatedRawScoreResults.equatedRawScores.row(methodIndex),
                                                                                                                 newFormFrequencyDistribution);

          equatedRawScoreResults.equatedRawScoreMoments.row(methodIndex) = moments.momentValues;
        }
      }

      /* calculate raw score moments for both new form and old form
        using IRT fitted distributions for new, old, and synthetic gps
        and the quadrature distributions for the new and old groups.
        To get these results, the quadrature distributions must be
        provided, as well as w1.  These results have an ambiguous status 
        in the context of IRT true-score equating which does not 
        involve synthetic groups. 

        From here to the end of the function there are 8 calls to MomentsFromRFD().
        The first three parameters of each call were revised on 3/8/09  */
      if (irtMethod == EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE ||
          irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE) {
        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(newForm.rawScores,
                                                                                                               newForm.fittedDistributionNewGroup);
        newForm.momentsFittedDistributionNewGroup = moments.momentValues;

        moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(newForm.rawScores,
                                                                          newForm.fittedDistributionOldGroup);
        newForm.momentsFittedDistributionOldGroup = moments.momentValues;

        moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(newForm.rawScores,
                                                                          newForm.fittedDistributionSyntheticGroup);
        newForm.momentsFittedDistributionSyntheticGroup = moments.momentValues;

        moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(oldForm.rawScores,
                                                                          oldForm.fittedDistributionNewGroup);
        oldForm.momentsFittedDistributionNewGroup = moments.momentValues;

        moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(oldForm.rawScores,
                                                                          oldForm.fittedDistributionOldGroup);
        oldForm.momentsFittedDistributionOldGroup = moments.momentValues;

        moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(oldForm.rawScores,
                                                                          oldForm.fittedDistributionSyntheticGroup);
        oldForm.momentsFittedDistributionSyntheticGroup = moments.momentValues;
      }

      /* moments for true-score and observed-score equivalents 
        using quadrature distribution for group that took new form.
        These differ from the GroupFormMoments() in that these
        employ the equivalents. Caveat: the results are somewhat 
        ambiguous for IRT true-score equating which does not i
        involve synthetic groups */
      if (irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_SCORE ||
          irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE) {
        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(irtEquatingResults.unroundedEquatedTrueScore,
                                                                                                               newForm.fittedDistributionNewGroup);

        irtEquatingResults.momentsEquatedTrueScores = moments.momentValues;
      }

      if (irtMethod == EquatingRecipes::Structures::IRTMethod::OBSERVED_SCORE ||
          irtMethod == EquatingRecipes::Structures::IRTMethod::TRUE_AND_OBSERVED_SCORE) {
        EquatingRecipes::Structures::Moments moments = EquatingRecipes::Utilities::momentsFromScoreFrequencies(irtEquatingResults.unroundedEquatedObservedScore,
                                                                                                               newForm.fittedDistributionNewGroup);

        irtEquatingResults.momentsEquatedObservedScores = moments.momentValues;
      }

      pData.irtInput = irtall;
    }

    /*---------------------------------------------------------------------------
      Functionality:
        Performs IRT true score equating with mixed-format tests.
      
      Input:
        Handle      A pointer to control the computing environment for equating
                    The same type of structure as used for IRT scale transformation
                    is used.
        NewItems    A pointer to designate items on the new form (0-offset)
        OldItems    A pointer to designate items on the old form (0-offset)
        nScores     Nnumber of new form scores for which old form equivalents
                    are to be computed
        newScores   New form scores for which old form equivalents are to be
                    computed. These scores are assumed to be equally spaced;
                    the difference between consecutive elements is assumed to be
                    constant (0-offset).

      Output:
        eqvOld      vector of old form equivalents of new form scores; 0-offset
        theta       vector of theta values which produce the vector of integer new
                    form true scores.
                    
        *NewMin     the lowest possible score on the new form.
        *OldMin     the lowest possible score on the old form.

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    bool trueScoreEquating(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                           const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                           const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems,
                           const size_t& numberOfScores,
                           const Eigen::VectorXd& newScores,
                           Eigen::VectorXd& oldFormEquivalents,
                           Eigen::VectorXd& theta,
                           double& lowestObservableScoreNewForm,
                           double& lowestObservableScoreOldForm) {
      double trueMinNew = 0.0;
      double trueMinOld = 0.0;
      double oldScoreMax = 0.0; /* Maximum old form score */
      double newTestMin = 0;    /* Minimum new form test score */
      double oldTestMin = 0;    /* Minimum old form test score */

      this->controlHandle = handle;
      this->controlNewItems = newItems;

      std::for_each(newItems.begin(),
                    newItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& newItem) {
                      if (newItem.irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC) {
                        trueMinNew += newItem.scoringFunctionValues(0) * (1.0 - newItem.c(1) +
                                                                          newItem.scoringFunctionValues(1) * newItem.c(1));
                      } else {
                        trueMinNew += newItem.scoringFunctionValues(0);
                      }

                      newTestMin += newItem.scoringFunctionValues(0);
                    });

      std::for_each(oldItems.begin(),
                    oldItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& oldItem) {
                      if (oldItem.irtModel == EquatingRecipes::Structures::IRTModel::THREE_PARAMETER_LOGISTIC) {
                        trueMinOld += (oldItem.scoringFunctionValues(0) + (oldItem.scoringFunctionValues(1)) * oldItem.c(1));
                      } else {
                        trueMinOld += oldItem.scoringFunctionValues(0);
                      }

                      size_t maximumCategory = oldItem.numberOfCategories;
                      oldScoreMax += oldItem.scoringFunctionValues(maximumCategory - 1);
                      oldTestMin += oldItem.scoringFunctionValues(0);
                    });

      /* Find index of largest score below sum of c's */
      size_t chance;
      for (chance = 0; chance < numberOfScores - 1; chance++) {
        if (newScores(chance) > trueMinNew) {
          break;
        }
      }

      /* Compute slope and intercept of interpolating function for 
        new form scores below chance */
      double slope = trueMinNew - newTestMin;

      if (slope > 0) {
        slope = (trueMinOld - oldTestMin) / slope;
      }

      double intercept = oldTestMin - slope * newTestMin;

      double xh = 99;
      double xl = -99;

      this->trueS = newScores(chance); /* Check if true score at theta=-99 is greater      */
      if (f_mix(xl) < 0.0) {           /* than computed largest score less than the        */
        ++chance;                      /* the sum of the c's, if so increase chance index. */
      }

      /* Assign score equivalents */
      for (size_t scoreLocation = 0; scoreLocation <= chance; scoreLocation++) {
        oldFormEquivalents(scoreLocation) = intercept + newScores(scoreLocation) * slope;
        theta(scoreLocation) = xl;
      }

      double x1 = -99.0;
      double x2 = 99.0;

      for (size_t scoreLocation = chance + 1; scoreLocation < numberOfScores - 1; scoreLocation++) {
        this->trueS = newScores(scoreLocation);

        for (size_t r = 0; r <= 10; r++) {
          theta(scoreLocation) = er_rtsafe(x1, x2, 0.00001);

          if (theta(scoreLocation) == -9999.0) {
            theta(scoreLocation) = -99.0; /* tianyou added this line 7/18/08 */
          }

          if (theta(scoreLocation) != -9999.0 &&
              theta(scoreLocation) != -99.0 &&
              theta(scoreLocation) != 99.0) {
            break;
          } else {
            x1 = -1.0 * std::pow(2.0, static_cast<double>(r));
            x2 = std::pow(2.0, static_cast<double>(r));
          }
        }

        oldFormEquivalents(scoreLocation) = trueScore(oldItems,
                                                      theta(scoreLocation));
      }

      /* Convert maximum score on new form to maximum score on old form */
      oldFormEquivalents(numberOfScores - 1) = oldScoreMax;
      theta(numberOfScores - 1) = xh;

      lowestObservableScoreNewForm = trueMinNew;
      lowestObservableScoreOldForm = trueMinOld;

      return 0;
    }

  private:
    EquatingRecipes::Structures::IRTScaleTransformationControl controlHandle;
    std::vector<EquatingRecipes::Structures::ItemSpecification> controlNewItems;
    EquatingRecipes::Structures::Symmetry controlSymmetry;
    bool controlFunctionStandardization;
    double trueS;

    /* IRT observed score equating functions */

    /*------------------------------------------------------------------------------
      Functionality:
        Performs calculation of IRT observed score equivalents of new form scores.

      Input:
        Handle      A pointer to control the computing environment for equating;
                    The same type of structure as used for IRT scale transformation.
        NewItems    A pointer to designate items on the new form (0-offset)
        OldItems    A pointer to designate items on the old form (0-offset)
        wNew        New group weight for a synthetic group
        wOld        Old group weight for a synthetic group; wNew + wOld = 1.0
        newForm     A pointer to an object of the RawFitDist structure for the new
                    form
        oldForm     A pointer to an object of the RawFitDist structure for the old
                    form
        RawEq       A pointer to an object of the RawTruObsEquiv structure to save
                    the equating results for raw scores.

      Output:
        newForm    Fitted distributions for the new form with the new, old, and
                  synthetic groups
        oldForm    Fitted distributions for the old form with the new, old, and
                  synthetic groups
        RawEq      Old form raw score equivalents of new form scores;
                  The results are saved into RawEq->unroundedEqObs.


      Author: Seonghoon Kim
      Date of last revision 9/25/08

    ------------------------------------------------------------------------------*/
    void irtMixObsEq(const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                     const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                     const std::vector<EquatingRecipes::Structures::ItemSpecification>& oldItems,
                     const double& wNew,
                     const double& wOld,
                     EquatingRecipes::Structures::IRTFittedDistribution& newForm,
                     EquatingRecipes::Structures::IRTFittedDistribution& oldForm,
                     EquatingRecipes::Structures::IRTEquatingResults& irtEquatingResults) {
      size_t newFormMaximumScore;
      size_t oldFormMaximumScore;
      Eigen::VectorXd newFormCumulativeFreqDist;
      Eigen::VectorXd oldFormCumulativeFreqDist;
      Eigen::VectorXd newFormPercentileRankDist;

      newFormMaximumScore = newForm.numberOfRawScoreCategories;
      oldFormMaximumScore = oldForm.numberOfRawScoreCategories;

      /* fitted distribution for the new form with the new group */
      irtMixObsDist(newItems,
                    handle.numberOfItemsNewForm,
                    newFormMaximumScore,
                    handle.quadratureNewForm.thetaValues,
                    handle.quadratureNewForm.thetaWeights,
                    newForm.numberOfRawScoreCategories,
                    newForm.rawScores,
                    newForm.fittedDistributionNewGroup);

      /* fitted distribution for the new form with the old group */
      irtMixObsDist(newItems,
                    handle.numberOfItemsNewForm,
                    newFormMaximumScore,
                    handle.quadratureOldForm.thetaValues,
                    handle.quadratureOldForm.thetaWeights,
                    newForm.numberOfRawScoreCategories,
                    newForm.rawScores,
                    newForm.fittedDistributionOldGroup);

      /* fitted distribution for the new form with the synthetic group */
      for (size_t scoreLocation = 0; scoreLocation < newForm.numberOfRawScoreCategories; scoreLocation++) {
        newForm.fittedDistributionSyntheticGroup(scoreLocation) = wNew * newForm.fittedDistributionNewGroup(scoreLocation) +
                                                                  wOld * newForm.fittedDistributionOldGroup(scoreLocation);
      }

      /* fitted distribution for the old form with the new group */
      irtMixObsDist(oldItems,
                    handle.numberOfItemsOldForm,
                    oldFormMaximumScore,
                    handle.quadratureNewForm.thetaValues,
                    handle.quadratureNewForm.thetaWeights,
                    oldForm.numberOfRawScoreCategories,
                    oldForm.rawScores,
                    oldForm.fittedDistributionNewGroup);

      if (oldForm.numberOfRawScoreCategories != oldFormMaximumScore) {
        throw std::runtime_error("\nPossibly wrong execution in IRTmixObsEq()\n");
      }

      /* fitted distribution for the old form with the old group */
      irtMixObsDist(oldItems,
                    handle.numberOfItemsOldForm,
                    oldFormMaximumScore,
                    handle.quadratureOldForm.thetaValues,
                    handle.quadratureOldForm.thetaWeights,
                    oldForm.numberOfRawScoreCategories,
                    oldForm.rawScores,
                    oldForm.fittedDistributionOldGroup);

      if (oldForm.numberOfRawScoreCategories != oldFormMaximumScore) {
        throw std::runtime_error("\nPossibly wrong execution in IRTmixObsEq()\n");
      }

      /* fitted distribution for the old form with the synthetic group */
      for (size_t scoreLocation = 0; scoreLocation < oldForm.numberOfRawScoreCategories; scoreLocation++) {
        oldForm.fittedDistributionSyntheticGroup(scoreLocation) = wNew * oldForm.fittedDistributionNewGroup(scoreLocation) +
                                                                  wOld * (oldForm.fittedDistributionOldGroup(scoreLocation));
      }

      newFormCumulativeFreqDist.resize(newForm.numberOfRawScoreCategories);

      newFormCumulativeFreqDist(0) = newForm.fittedDistributionSyntheticGroup(0);
      for (size_t scoreLocation = 1; scoreLocation < newForm.numberOfRawScoreCategories; scoreLocation++) {
        newFormCumulativeFreqDist(scoreLocation) = newFormCumulativeFreqDist(scoreLocation - 1) + newForm.fittedDistributionSyntheticGroup(scoreLocation);
      }

      oldFormCumulativeFreqDist.resize(oldForm.numberOfRawScoreCategories);

      oldFormCumulativeFreqDist(0) = oldForm.fittedDistributionSyntheticGroup(0);
      for (size_t scoreLocation = 1; scoreLocation < oldForm.numberOfRawScoreCategories; scoreLocation++) {
        oldFormCumulativeFreqDist(scoreLocation) = oldFormCumulativeFreqDist(scoreLocation - 1) + oldForm.fittedDistributionSyntheticGroup(scoreLocation);
      }

      /* conduct equipercentile equating with synthetic distributions */
      double oldFormMinimumScore = oldForm.rawScores(0);
      double oldFormScoreIncrement = oldForm.rawScores(1) - oldFormMinimumScore;

      /* compute the percentile rank function --Tianyou added 8/18/08 */
      newFormPercentileRankDist = EquatingRecipes::Utilities::percentileRanks(0,
                                                                              newForm.numberOfRawScoreCategories - 1,
                                                                              1,
                                                                              newFormCumulativeFreqDist);

      /* calling ERutilities equipercentile equating function (TW 8/18/08) */
      irtEquatingResults.unroundedEquatedObservedScore = EquatingRecipes::Utilities::getEquipercentileEquivalents(oldForm.numberOfRawScoreCategories,
                                                                                                                  oldFormMinimumScore,
                                                                                                                  oldFormScoreIncrement,
                                                                                                                  oldFormCumulativeFreqDist,
                                                                                                                  newForm.numberOfRawScoreCategories,
                                                                                                                  newFormPercentileRankDist);
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Calculates the marginal distribution of total score for IRT models. The IRT
        models include the 3PL, LGR, GPC, and NR models.

      Input:
        Items:   pointer to designate items on a test form (0-offset)
        n :      number of items on a test form
        nq:      number of quadrature points
        xqpts:   vector for quadrature points (0-offset)
        xqwts:   vector for quadrature weights (0-offset)

      Output:
        nscr:    number of score categories for observed score distribution
        xscr:    vector of scores associated with each score category
                (are consecutive integers); 0-offset
        xmarg:   vector of marginal probabilities associated with each score
                category; 0-offset 

      Author: Seonghoon Kim
      Date of last revision 9/25/08

    ------------------------------------------------------------------------------*/
    void irtMixObsDist(const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                       const size_t& numberOfItemsOnTestForm,
                       const size_t& maximumScorePoint,
                       const Eigen::VectorXd& quadraturePoints,
                       const Eigen::VectorXd& quadratureWeights,
                       size_t& numberOfScores,
                       Eigen::VectorXd& scores,
                       Eigen::VectorXd& marginalResponseProbabilities) {
      // int i, j, k, MaxCat;
      // 	double theta, xx, *xnew;

      /* finds the maximum number of categories across items */
      size_t maximumNumberOfCategories = 2;
      std::for_each(items.begin(),
                    items.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      maximumNumberOfCategories = std::max(maximumNumberOfCategories, item.numberOfCategories);
                    });

      Eigen::VectorXd xnew(maximumNumberOfCategories);

      marginalResponseProbabilities.setZero();

      for (size_t quadraturePointIndex = 0; quadraturePointIndex < quadraturePoints.size(); quadraturePointIndex++) {
        double theta = quadraturePoints(quadraturePointIndex);

        ObsDistGivenTheta(theta,
                          items,
                          numberOfItemsOnTestForm,
                          maximumNumberOfCategories,
                          maximumScorePoint,
                          numberOfScores,
                          scores,
                          xnew);

        marginalResponseProbabilities += quadratureWeights(quadraturePointIndex) * xnew;
      }

      marginalResponseProbabilities /= marginalResponseProbabilities.sum();
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Calculates the conditional distribution of total score given theta
        for IRT models. The IRT models include the 3PL, LGR, GPC, and NR models.
        Uses Hanson's (1994) generalization of the Lord-Wingersky (1982) recursive
        algorithm.

      Input:
        theta:   examinee ability
        Items:   pointer to designate items on a test form (0-offset)
        n :      number of items on a test form
        MaxCat:  the maximum number of categories across items
        MaxScrP: maximum number of score points

      Output:
        nscr:    number of score categories for observed score distribution
        xscr:    vector of scores associated with each score category
                (are consecutive integers); 0-offset
        xnew:    vector of probabilities associated with each score category;
                0-offset

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void ObsDistGivenTheta(const double& theta,
                           const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                           const size_t& numberOfItemsOnForm,
                           const size_t& maximumCategoryIndex,
                           const size_t& maximumScorePoint,
                           size_t& numberOfScores,
                           Eigen::VectorXd& scores,
                           Eigen::VectorXd& xnew) {
      // int i, j, k, index;
      // int mino, maxo, minn, maxn;
      Eigen::VectorXd xitem(maximumCategoryIndex + 1); /* zero-offset, but not use xitem[0] */
      Eigen::VectorXd xold(maximumScorePoint);         /* zero-offset */

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      /* calculates probabilities for Item 1 */
      for (size_t categoryIndex = 0; categoryIndex < items[0].numberOfCategories; categoryIndex++) {
        xitem(categoryIndex) = irtModelFunctions.itemResponseFunction(items[0], categoryIndex, theta);
      }

      double mino = items[0].scoringFunctionValues(0);
      double maxo = items[0].scoringFunctionValues(items[0].numberOfCategories - 1);
      double minn = mino;
      double maxn = maxo;

      xold.setZero();

      for (size_t categoryIndex = 0; categoryIndex < items[0].numberOfCategories; categoryIndex++) {
        size_t index = static_cast<size_t>(items[0].scoringFunctionValues(categoryIndex) - minn);
        xold(index) = xitem(categoryIndex); /* mino associated with index of 0 */
      }                                     /* mino does vary; see below      */

      xnew = xold;

      if (numberOfItemsOnForm == 1) {
        size_t maxMinusMin = static_cast<size_t>(maxn - minn);

        for (size_t score = 0; score <= maxMinusMin; score++) {
          scores(score) = static_cast<double>(score) + minn;
        }

        numberOfScores = static_cast<size_t>(maxn - minn + 1);

        return;
      }

      /* updates distribution for items 2 through nitems */
      for (size_t itemIndex = 1; itemIndex < numberOfItemsOnForm; itemIndex++) {
        for (size_t categoryIndex = 0; categoryIndex < items[itemIndex].numberOfCategories; categoryIndex++) {
          xitem(categoryIndex) = irtModelFunctions.itemResponseFunction(items[itemIndex],
                                                                        categoryIndex,
                                                                        theta);
        }

        recurs(mino,
               maxo,
               xold,
               items[itemIndex].numberOfCategories,
               items[itemIndex].scoringFunctionValues,
               xitem,
               minn,
               maxn,
               xnew);

        mino = minn;
        maxo = maxn;

        size_t maxIndex = static_cast<size_t>(maxn - minn);

        xold(Eigen::seq(0, maxIndex)) = xnew(Eigen::seq(0, maxIndex));
      }

      size_t maxIndex = static_cast<size_t>(maxn - minn);
      for (size_t index = 0; index <= maxIndex; index++) {
        scores(index) = static_cast<double>(index) + minn;
      }

      numberOfScores = static_cast<size_t>(maxn - minn + 1);
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Updates a distribution of scores using Hanson's (1994) generalization of
        the Lord-Wingersky (1982) formula.

        Assumes that test scores are consecutive integers
        Assumes that item scores are consecutive integers
        Assumes that test and item scores are sorted from low to high.

      Input (r-1 added items):
        mino : minimum integer score for old distribution
        maxo : maximum integer score for old distribution
        xold : probability array for each score point from mino to maxo
              (f_r-1 (x|theta)); 0-offset
        mitem: number of distinct score points for the added (new) item
        iitem: values of the scoring function for the added item; 0-offset
        xitem: values of category response functions for the added item; 0-offset

      Output (r added items):
        minn : minimum integer score for new (updated) distribution
        maxn : maximum integer score for new (updated) distribution
        xnew : probability array for each score point from minn to maxn
              (f_r(x|theta)); 0-offset

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void recurs(const int& mino,
                const int& maxo,
                const Eigen::VectorXd& xold,
                const size_t& numberOfCategories,
                const Eigen::VectorXd& iitem,
                const Eigen::VectorXd& xitem,
                double& minn,
                double& maxn,
                Eigen::VectorXd& xnew) {
      minn = mino + iitem(0);
      maxn = maxo + iitem(numberOfCategories - 1);

      for (double i = minn; i <= maxn; i += 1.0) {
        double in = i - minn;

        xnew(static_cast<size_t>(in)) = 0.0;

        for (size_t j = 0; j < numberOfCategories; j++) {
          double io = i - iitem(j) - mino;

          if (io >= 0 && io <= maxo - mino) {
            xnew(static_cast<size_t>(in)) += xold(static_cast<size_t>(io)) * xitem(j);
          }
        }
      }
    }

    /*
      Purpose:  
          This function implements bisection method to find x such that

              f(x)==0.0.

      Input:             
          funcd   function pointer whose root we are trying to find.
        x0      left end of the function domain 
        x1      right end of the function domain

      Output:                
          x that satisfies f(x)==0.0
      
      Precondition:                                                     
        f(x) is defined on [x0,x1] with f(x0)*f(x1) < = 0.
    
      Author: Jaehoon Seol
      Date: August 20, 2009
      Version: 1.0
      References : 
          J. Stoer & R. Bulirsch, Introduction to Numerical analysis, 2nd edition 
        David S. Watkins, Fundamentals of Matrix Computations,2002
      Comments:
    */
    double er_rtsafe(double x0, double x1, double error) {
      double left = f_mix(x0);
      double right = f_mix(x1);

      /* assert that there is a root between x1 and x2 */
      if (left * right > 0) { /* function values have the same sign */
        std::string msg = "Equating Recipes error occured\n";
        msg.append("Source: er_find_root, Error: no root found exists on the interval\n");
        throw std::runtime_error(msg);
      }

      double diff = std::abs(right - left);
      double mid = (left + right) / 2.0;

      while (diff > error) {
        mid = (left + right) / 2.0;
        double side1 = f_mix(mid);
        double side2 = f_mix(right);

        if (side1 * side2 <= 0) {
          left = mid;
        } else {
          right = mid;
        }

        diff = std::abs(right - left);
      }

      return mid;
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Allocate memory to the pointer members of the RawFitDist structure.
        
        Input:
          Items       A pointer to an array of the ItemSpec structure
          oldOrnew    A string indicating that the array is about either the old or
                      new form (either "old" or "new" must be used as characters)
          Handle      A pointer to a variable of the IRTstControl structure 
          Form        A pointer to an object of the RawFitDist structure

      Author: Seonghoon Kim (with some modifications by Tianyou Wang and R. L. Brennan)
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void initializeRawFitMem(const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                             const bool& isNewForm,
                             EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                             EquatingRecipes::Structures::IRTFittedDistribution& irtFittedDistribution) {
      int minimumNumberOfScoreCategories;
      int maximumNumberOfScoreCategories;
      double minimumScore;
      double maximumScore;

      minMaxNumberOfScoreCategories(items,
                                    minimumNumberOfScoreCategories,
                                    maximumNumberOfScoreCategories,
                                    minimumScore,
                                    maximumScore);

      if (isNewForm) {
        handle.minimumRawScoreNewForm = static_cast<double>(minimumNumberOfScoreCategories);
        handle.maximumRawScoreNewForm = maximumScore;
        handle.rawScoreIncrementNewForm = 1;
      } else {
        handle.minimumRawScoreOldForm = static_cast<double>(minimumNumberOfScoreCategories);
        handle.maximumRawScoreOldForm = maximumScore;
        handle.rawScoreIncrementOldForm = 1;
      }

      irtFittedDistribution.numberOfRawScoreCategories = static_cast<size_t>(maximumNumberOfScoreCategories);
      irtFittedDistribution.rawScores.resize(irtFittedDistribution.numberOfRawScoreCategories);
      irtFittedDistribution.fittedDistributionNewGroup.resize(irtFittedDistribution.numberOfRawScoreCategories);
      irtFittedDistribution.fittedDistributionOldGroup.resize(irtFittedDistribution.numberOfRawScoreCategories);
      irtFittedDistribution.fittedDistributionSyntheticGroup.resize(irtFittedDistribution.numberOfRawScoreCategories);

      for (size_t scoreLocation = 0; scoreLocation < irtFittedDistribution.numberOfRawScoreCategories; scoreLocation++) {
        irtFittedDistribution.rawScores(scoreLocation) = static_cast<double>(minimumNumberOfScoreCategories) +
                                                         static_cast<double>(scoreLocation);
      }
    }

    /*------------------------------------------------------------------------------
      Functionality:
        Allocate memory to the pointer members of the RawTruObsEquiv structure.
        
        Input:
          NewItems    A pointer to an array of the ItemSpec structure for the new
                      form 
          Handle      A pointer to a variable of the IRTstControl structure 
          RawEq       A pointer to an object of the RawTruObsEquiv structure

      Author: Seonghoon Kim
      Date of last revision 9/25/08
    ------------------------------------------------------------------------------*/
    void initializeRawEqResults(const std::vector<EquatingRecipes::Structures::ItemSpecification>& newItems,
                                const EquatingRecipes::Structures::IRTScaleTransformationControl& handle,
                                EquatingRecipes::Structures::IRTEquatingResults& irtEquatingResults) {
      int minimumNumberOfScoreCategories;
      int maximumNumberOfScoreCategories;
      double maximumScore;
      double minimumScore;

      minMaxNumberOfScoreCategories(newItems,
                                    minimumNumberOfScoreCategories,
                                    maximumNumberOfScoreCategories,
                                    minimumScore,
                                    maximumScore);

      irtEquatingResults.numberOfRawScoreCategoriesNewForm = static_cast<size_t>(maximumNumberOfScoreCategories);
      irtEquatingResults.thetaEquivalentFormXScore.resize(static_cast<size_t>(maximumNumberOfScoreCategories));
      irtEquatingResults.unroundedEquatedTrueScore.resize(static_cast<size_t>(maximumNumberOfScoreCategories));
      irtEquatingResults.unroundedEquatedObservedScore.resize(static_cast<size_t>(maximumNumberOfScoreCategories));
      irtEquatingResults.roundedEquatedTrueScore.resize(static_cast<size_t>(maximumNumberOfScoreCategories));
      irtEquatingResults.roundedEquatedObservedScore.resize(static_cast<size_t>(maximumNumberOfScoreCategories));
    }

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double f_mix(const double& theta) {
      double expectedRawScore = trueScore(this->controlNewItems,
                                          theta);

      return (this->trueS - expectedRawScore);
    }

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double f_mixDer(const double& theta) {
      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      double v = 0.0;
      std::for_each(this->controlNewItems.begin(),
                    this->controlNewItems.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      for (size_t categoryIndex = 0; categoryIndex < item.scoringFunctionValues.size(); categoryIndex++) {
                        double pd = irtModelFunctions.itemResponseFunctionDerivative(item, categoryIndex, theta);

                        v += pd * item.scoringFunctionValues(categoryIndex);
                      }
                    });

      v *= -1.0;

      return v;
    }

    /*
      Author: Seonghoon Kim
      Date of last revision 9/25/08
    */
    double trueScore(const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                     const double& theta) {
      double expectedRawScore = 0.0;

      EquatingRecipes::IRTModelFunctions irtModelFunctions;

      std::for_each(items.begin(),
                    items.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      for (size_t categoryIndex = 0; categoryIndex < item.scoringFunctionValues.size(); categoryIndex++) {
                        double probResponse = irtModelFunctions.itemResponseFunction(item,
                                                                   categoryIndex,
                                                                   theta);

                        expectedRawScore += probResponse * item.scoringFunctionValues(categoryIndex);
                      }
                    });

      return expectedRawScore;
    }

    void minMaxNumberOfScoreCategories(const std::vector<EquatingRecipes::Structures::ItemSpecification>& items,
                                       int& minimumNumberOfScoreCategories,
                                       int& maximumNumberOfScoreCategories,
                                       double& minimumScore,
                                       double& maximumScore) {
      /* Calculate the largest possible number of observed score categories */
      maximumScore = 0;
      minimumScore = 0;

      std::for_each(items.begin(),
                    items.end(),
                    [&](const EquatingRecipes::Structures::ItemSpecification& item) {
                      minimumScore += item.scoringFunctionValues(0);
                      maximumScore += item.scoringFunctionValues(item.scoringFunctionValues.size() - 1);
                    });

      if (minimumScore < 0.0) {
        minimumNumberOfScoreCategories = static_cast<int>(minimumScore - 0.5);
      } else {
        minimumNumberOfScoreCategories = static_cast<int>(minimumScore + 0.5);
      }

      maximumNumberOfScoreCategories = static_cast<int>(maximumScore) - minimumNumberOfScoreCategories + 1;
    }
  };
} // namespace EquatingRecipes

#endif