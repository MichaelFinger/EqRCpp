/* 
  From Source: ERutilities.h, ERutilities.c
  Original Struct: 
  Description: 
*/

/*
  ERutilities.c 

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

  NOTES:

       Unless otherwise noted, raw scores can differ by any constant 
       positive real number, and raw scores can be negative. (An obvious
       exception is the beta-binomial model.)  In general, the user specifies
       the lowest raw score (min), the highest raw score (max), and the 
       increment (inc) between any two adjacent raw scores. Freqencies and 
       related statistics for raw scores are stored in zero-offset vectors 
       or matrices. So, if fd[] is a frequency distribution, then fd[0] is
       associated with min, and fd[(int) ((max-min)/inc + .5)] is associated
       with max.
 
       In these utilities, mind and maxd are the minimum and maximum scores,
       respectively, in the data, whereas min and max are the user-specified
       minimum and maximum scores. It must always be true that min <= mind 
       and max >= maxd.  Often, directly or indirectly, the user will set 
       min = mind and max = maxd, especially for number correct scores.  
       This distinction allows zero frequencies at the low and high ends of
       distributions.       
       
       minp and maxp are the minimum and maximum possible raw scores.
       It must be true that minp<=min and maxp>=max.  minp and maxp are
       used exclusively with scale scales; whereas min and max are used 
       exclusively with raw scores.  Distinguishing between min and minp
       (as well as between max and maxp) can be necessary (or at least
       desirable) in cases such as the following:
         (a) a particular new form has fewer items than the old form,
             but the scale score conversion for the new form must be
             provided for all of the possible old-form scores;
         (b) min is set equal to mind and max is set equal to maxd 
       For more information about minp and maxp, see comments in 
       ReadSSConvTableForY().
       
       No function that involves scale scores is
       design/methodology/smoothing dependent. 

       The single-group design (S or SG) requires one struct BSTATS.
       The random-groups design (R or RG) requires two struct USTATS.
       The common-item non-equivalent groups design (C or CG)
         requires two struct BSTATS. 

  R. L. Brennan

  March 18, 2009                           
*/

#ifndef STRUCTURES_UTILITIES_HPP
#define STRUCTURES_UTILITIES_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <string>

#include <Eigen/Dense>
#include <fmt/core.h>

#include <equating_recipes/structures/raw_to_scaled_score_table.hpp>
#include <equating_recipes/structures/design.hpp>
#include <equating_recipes/structures/method.hpp>

namespace EquatingRecipes {
  struct Utilities {
    /*
        Location of score x in zero-offset vector; 

        Input:
          x = score
          min = min score
          inc = score increment

          NOTE:  best if x, min, and inc are specified very precisely 
                (say at least eight digits)

        Returns location of score x in zero-offset vector
        
        Function calls other than C or NR utilities: None
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08         
      */
    static size_t getScoreLocation(const double& score,
                                   const double& minimumScore,
                                   const double& scoreIncrement) {
      size_t location = static_cast<size_t>((score - minimumScore) / scoreIncrement + 0.5);

      return location;
    }

    /*
        Number of scores (or categories in zero-offset vector; 

        Input:
          max = max score
          min = min score
          inc = score increment

          NOTE:  best if max, min, and inc are specified very precisely 
                (say at least eight digits)

        Returns number of scores (or categories) in zero-offset vector
        with (implicit) scores of min(inc)max;   

        Function calls other than C or NR utilities: None
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08         
      */
    static size_t numberOfScores(const double& minimumScore,
                                 const double& maximumScore,
                                 const double& scoreIncrement) {
      size_t nScores = getScoreLocation(maximumScore, minimumScore, scoreIncrement) + 1;

      return nScores;
    }

    /*
        Score associated with location in zero-offset vector 

        Input:
          loc = location in zero-offset vector
          min = min score
          inc = score increment

          NOTE:  best if min and inc are specified very precisely 
                (say at least eight digits)

        Returns score associated with location in zero-offset vector

        Function calls other than C or NR utilities: None
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08  
      */
    static double getScore(const size_t& scoreLocation,
                           const double& minimumScore,
                           const double& scoreIncrement) {
      double scoreValue = minimumScore + static_cast<double>(scoreLocation) * scoreIncrement;

      return scoreValue;
    }

    /*
        Compute cumulative relative frequencies from relative frequencies

        Input
          min = min score
          max = max score
          inc = increment
          rfd[] = rel freq dist

        Output
          crfd[] = cum rel freq dist (space already allocated)

        Function calls other than C or NR utilities:
          loc()
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08  
      */
    static Eigen::VectorXd cumulativeRelativeFreqDist(const double& minimumScore,
                                                      const double& maximumScore,
                                                      const double& scoreIncrement,
                                                      const Eigen::VectorXd& relativeFreqDist) {
      size_t maximumScoreLocation = Utilities::getScoreLocation(maximumScore,
                                                                minimumScore,
                                                                scoreIncrement);

      Eigen::VectorXd cumRelFreqDist = Eigen::VectorXd::Zero(maximumScoreLocation + 1);

      cumRelFreqDist(0) = relativeFreqDist(0);

      for (size_t freqDistIndex = 1; freqDistIndex < maximumScoreLocation; ++freqDistIndex) {
        cumRelFreqDist(freqDistIndex) = cumRelFreqDist(freqDistIndex - 1) + relativeFreqDist(freqDistIndex);
      }

      return cumRelFreqDist;
    }

    /*
        Compute percentile rank (pr) given cumulative relative frequencies crfd[]

        pr can be computed for any score x, not simply those scores that are 
          actually achieved (i.e, associated with a distribution). Code would be
          simpler for only "achievable" scores, and simpler still for only 
          integer scores.  

        Formula used is the analogue of Equation 2.14 in Kolen nd Brennan (2004) 
          for the more general case of x scores being any real numbers 
          differing by a constant amount (inc) and ranging from min to max

        In effect, the inverse of this function is the perc_point() function.

        Input
          min = min score
          max = max score
          inc = increment
          crfd[] = cum rel freq dist
          x = score for which pr is desired 

        Returns percentile rank

        Function calls other than C or NR utilities:
          loc()
          score()
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08  
      */
    static double getPercentileRank(const double& minimumScore,
                                    const double& maximumScore,
                                    const double& scoreIncrement,
                                    const Eigen::VectorXd& cumulativeRelativeFreqDist,
                                    const double& score) {
      double percRank;

      if (score < minimumScore - scoreIncrement / 2.0) {
        percRank = 0.0;
      } else if (score < minimumScore + scoreIncrement / 2.0) {
        percRank = 100 * ((score - (minimumScore - scoreIncrement / 2.0)) / scoreIncrement) * cumulativeRelativeFreqDist[0];
      } else if (score >= maximumScore + scoreIncrement / 2.0) {
        percRank = 100.0;
      } else {
        size_t maxScoreLocation = Utilities::getScoreLocation(maximumScore, minimumScore, scoreIncrement);

        size_t scoreLocation;
        double xStar;

        for (scoreLocation = 1; scoreLocation <= maxScoreLocation; ++scoreLocation) {
          xStar = Utilities::getScore(scoreLocation, minimumScore, scoreIncrement);
          if (score < xStar + scoreIncrement / 2.0) {
            break;
          }
        }

        percRank = 100 * (cumulativeRelativeFreqDist(scoreLocation - 1) +
                          ((score - (xStar - scoreIncrement / 2.0)) / scoreIncrement) *
                              (cumulativeRelativeFreqDist(scoreLocation) - cumulativeRelativeFreqDist(scoreLocation - 1)));
      }

      return percRank;
    }

    /* 
        Interpolated value of f at x assuming there are ns score
        categories.  Function treats first category as having 
        a score of 0; last category has score of ns-1. 
        For x < 0 or x > ns-1, extrapolated value is provided.

        Input
          x = score
          ns = number of score categories
          f = f[0]...f[ns-1]: usually (relative) frequencies

        Returned: interpolated value of f at score of x

        Function calls other than C or NR utilities: None
                                                    
        R. L. Brennan

        Date of last revision: 6/30/08  
      */
    static double interpolate(const double& score,
                              const int& numberOfScoreCategories,
                              const Eigen::VectorXd& frequencies) {
      int truncatedScore = static_cast<int>(score);
      size_t maximumScore = numberOfScoreCategories - 1;
      double value;

      if (score <= 0.0) { /* extrapolate low */
        value = frequencies(0) + score * (frequencies(1) - frequencies(0));
      } else if (score >= static_cast<double>(maximumScore)) { /* extrapolate high */
        value = frequencies(maximumScore) + (score - maximumScore) * (frequencies(maximumScore) - frequencies(maximumScore - 1));
      } else { /* interpolate */
        value = frequencies(truncatedScore) + (score - truncatedScore) * (frequencies(truncatedScore + 1) - frequencies(truncatedScore));
      }

      /* an extrapolated value can be less than 0.  Hence ... */
      double result = (value > 0) ? value : 1.0e-10;

      return result;
    }

    // perc_point, ERUtilities.h, ERUtilities.c
    static double percentilePoint(const size_t& numberOfScores,
                                  const double& minimumScore,
                                  const double& scoreIncrement,
                                  const Eigen::VectorXd& cumulativeRelativeFreqDist,
                                  const double& percentileRank) {
      double percentileRankProportion = percentileRank / 100.0;
      double percentilePointLowerValue;
      double percentilePointUpperValue;

      double percPoint;

      if (percentileRankProportion <= 1.0e-8) {
        /* Special case: PR=0 and 0 freqs at bottom of rfd[].
        First line of code means that prp <= 1.0e-8 
        is an operational definition of prp==0.  Remaining code is
        to handle possibility of 0 freq's at bottom of rfd[].
        Note that for loop starts at 0 */

        percentilePointLowerValue = -0.5;

        size_t scoreLocation;
        for (scoreLocation = 0; scoreLocation < numberOfScores - 1; scoreLocation++) {
          if (cumulativeRelativeFreqDist(scoreLocation) > 1.0e-8) {
            break;
          }
        }

        percentilePointUpperValue = static_cast<double>(scoreLocation) - 0.5;

        percPoint = minimumScore + scoreIncrement * ((percentilePointUpperValue + percentilePointLowerValue) / 2.0);
      } else if (percentileRankProportion >= 1.0 - 1.0e-8) {
        /* Special case: PR=1 and 0 freqs at top of rfd[].
        First line of code means that 1 - 1.0e-8 is 
        an operational definition of prp==1.  Remaining code is
        to handle possibility of 0 freq's at top of rfd[].
        Not that for loop starts at ns-1 */

        percentilePointUpperValue = static_cast<double>(numberOfScores) - 0.5;

        size_t scoreLocation;
        for (scoreLocation = numberOfScores - 1; scoreLocation >= 0; scoreLocation--) {
          if (cumulativeRelativeFreqDist(scoreLocation) < 1.0 - 1.0e-8) {
            break;
          }
        }

        percentilePointLowerValue = 1 + scoreLocation + 0.5;

        percPoint = minimumScore + scoreIncrement * ((percentilePointUpperValue + percentilePointLowerValue) / 2.0);
      } else if (cumulativeRelativeFreqDist(0) > percentileRankProportion) {
        /* Special case:  crfd[0] > prp can happen in equipercentile
          equating. If this occurs, then we have the following anomalous
          circumstances: 
            (a) x*_U = 0 by defn on p. 45, in which case
                crfd[x*_U - 1] = crfd[-1] which is undefined; and 
            (b) x*_L is non-existent by defn on p. 45.
          Next two lines of code are consistent with graphical 
          procedures in Kolen and Brennan (sect. 2.5.1). 
        */
        percentilePointUpperValue = 2.0 * percentileRankProportion / cumulativeRelativeFreqDist(0);
        percentilePointLowerValue = -1.0;

        percPoint = minimumScore + scoreIncrement * ((percentileRankProportion / cumulativeRelativeFreqDist(0)) - 0.5);
      } else {
        /* upper pp -- get x*_U */

        size_t scoreLocation;
        for (scoreLocation = 1; scoreLocation < numberOfScores; scoreLocation++) {
          if (cumulativeRelativeFreqDist(scoreLocation) > percentileRankProportion) {
            break;
          }
        }

        if (cumulativeRelativeFreqDist(scoreLocation) != cumulativeRelativeFreqDist(scoreLocation - 1)) {
          percentilePointUpperValue = ((percentileRankProportion - cumulativeRelativeFreqDist(scoreLocation - 1)) /
                                       (cumulativeRelativeFreqDist(scoreLocation) - cumulativeRelativeFreqDist(scoreLocation - 1))) +
                                      (static_cast<double>(scoreLocation) - 0.5);
        } else {
          percentilePointUpperValue = static_cast<double>(scoreLocation) - 0.5;
        }

        /* lower pp -- get x*_L */

        for (scoreLocation = numberOfScores - 2; scoreLocation >= 0; scoreLocation--) {
          if (cumulativeRelativeFreqDist(scoreLocation) < percentileRankProportion) {
            break;
          }
        }

        if (cumulativeRelativeFreqDist(scoreLocation + 1) != cumulativeRelativeFreqDist(scoreLocation)) {
          percentilePointLowerValue = ((percentileRankProportion - cumulativeRelativeFreqDist(scoreLocation)) /
                                       (cumulativeRelativeFreqDist(scoreLocation + 1) - cumulativeRelativeFreqDist(scoreLocation))) +
                                      (static_cast<double>(scoreLocation) + 0.5);
        } else {
          percentilePointUpperValue = scoreLocation + 0.5;
        }

        /* return area */

        percPoint = minimumScore + scoreIncrement * ((percentilePointUpperValue + percentilePointLowerValue) / 2.0);
      }
      return percPoint;
    }

    /*
      EquiEquate
  
      Computes equipercentile equivalents on scale of y for percentile
      ranks on scale of x. See comments in perc_point() for details.
    */
    static Eigen::VectorXd getEquipercentileEquivalents(const size_t& numberOfRawScoreCategoriesY,
                                                        const double& minimumRawScoreY,
                                                        const double& rawScoreIncrementY,
                                                        const Eigen::VectorXd& cumulativeRelativeFreqDistY,
                                                        size_t numberOfRawScoreCategoriesX,
                                                        const Eigen::VectorXd& percentileRankDistX) {
      Eigen::VectorXd equipercentileEquivalents(numberOfRawScoreCategoriesX);

      // for (i = 0; i <= nsx - 1; i++)
      //   eraw[i] = perc_point(nsy, miny, incy, crfdy, prdx[i]);
      for (size_t scoreLocation = 0; scoreLocation < numberOfRawScoreCategoriesX; scoreLocation++) {
        equipercentileEquivalents(scoreLocation) = Utilities::percentilePoint(numberOfRawScoreCategoriesY,
                                                                              minimumRawScoreY,
                                                                              rawScoreIncrementY,
                                                                              cumulativeRelativeFreqDistY,
                                                                              percentileRankDistX(scoreLocation));
      }

      return equipercentileEquivalents;
    }

    /* 
      From Source: ERutilities.h 
      Original Struct: Equated_ss
      Description: 
        Get SS conversion given raw score equating results (eraw) and raw-to-ss
        conversion Y (yct).  This function is for one method only, and it 
        does not make any assumptions about the equating method used.  For the 4 
        linear results this function can be called 4 times using eraw[0]...eraw[3] 

        NOTE: It is assumed that increment for x is the same as inc in the 
              conversion table for the base form Y

        Note that:
          (a) eraw[] goes from eraw[0] to
              eraw[loc(max,min,inc)] = eraw[nscores(max,min.inc)-1]
          (b) yct[*][] goes from yct[*][0] to 
              yct[*][loc(maxp,minp,inc)+2] = yct[*][nscores(maxp,minp,inc)+1]
          (c) raw score associated with eraw[i] is score(i,min,inc) 
    */

    static void getEquatedScaledScores(const double& minimumRawScoreX,
                                       const double& maximumRawScoreX,
                                       const double& rawScoreIncrement,
                                       const double& minimumRawScoreYct,
                                       const double& maximumRawScoreYct,
                                       const double& scoreIncrementYct,
                                       const Eigen::VectorXd& equatedRawScores,
                                       const EquatingRecipes::Structures::RawToScaledScoreTable& rawToScaledScoreTable,
                                       const size_t& roundToNumberOfDecimalPlaces,
                                       const int lowestObservableRoundedScaledScore,
                                       const int highestObservableRoundedScaledScore,
                                       Eigen::VectorXd& unroundedEquatedScaledScores,
                                       Eigen::VectorXd& roundedEquatedScaledScores) {
      size_t numberOfRawScoresYct = Utilities::numberOfScores(maximumRawScoreYct,
                                                              minimumRawScoreYct,
                                                              scoreIncrementYct);

      size_t numberOfEquatedRawScores = Utilities::numberOfScores(maximumRawScoreX,
                                                                  minimumRawScoreX,
                                                                  rawScoreIncrement);

      EquatingRecipes::Structures::RawToScaledScoreTable::Entry firstRawToScaledScoreEntry = rawToScaledScoreTable.getFirstEntry();
      EquatingRecipes::Structures::RawToScaledScoreTable::Entry lastRawToScaledScoreEntry = rawToScaledScoreTable.getLastEntry();

      size_t lastScoreLocation = 0;

      for (size_t scoreLocation = 0; scoreLocation < numberOfEquatedRawScores; scoreLocation++) {
        double equatedRawScore = equatedRawScores(scoreLocation);

        if (equatedRawScore <= firstRawToScaledScoreEntry.rawScore) {
          unroundedEquatedScaledScores(scoreLocation) = firstRawToScaledScoreEntry.scaledScore;
        } else if (equatedRawScore >= lastRawToScaledScoreEntry.rawScore) {
          unroundedEquatedScaledScores(scoreLocation) = lastRawToScaledScoreEntry.scaledScore;
        } else {
          /* find j such that yct[0][j] > y */
          size_t innerScoreLocation;

          for (innerScoreLocation = lastScoreLocation; innerScoreLocation < numberOfRawScoresYct; innerScoreLocation++) {
            EquatingRecipes::Structures::RawToScaledScoreTable::Entry entry = rawToScaledScoreTable.getEntry(innerScoreLocation);

            if (entry.rawScore > equatedRawScore) {
              break;
            }
          }

          EquatingRecipes::Structures::RawToScaledScoreTable::Entry entryLocationMinus1 = rawToScaledScoreTable.getEntry(innerScoreLocation - 1);
          EquatingRecipes::Structures::RawToScaledScoreTable::Entry entryLocation = rawToScaledScoreTable.getEntry(innerScoreLocation);
          unroundedEquatedScaledScores(scoreLocation) = entryLocationMinus1.scaledScore +
                                                        (entryLocation.scaledScore - entryLocationMinus1.scaledScore) *
                                                            (equatedRawScore - entryLocationMinus1.rawScore) / (entryLocation.rawScore - entryLocationMinus1.rawScore);

          lastScoreLocation = innerScoreLocation;
        }
      }

      for (size_t scoreLocation = 0; scoreLocation < numberOfEquatedRawScores; scoreLocation++) {
        if (roundToNumberOfDecimalPlaces >= 1) {
          roundedEquatedScaledScores(scoreLocation) = std::pow(static_cast<double>(10.0), static_cast<double>(roundToNumberOfDecimalPlaces - 1)) *
                                                      std::trunc(unroundedEquatedScaledScores(scoreLocation) /
                                                                     std::pow(static_cast<double>(10.0), static_cast<double>(roundToNumberOfDecimalPlaces - 1)) +
                                                                 0.5);

          std::clamp(roundedEquatedScaledScores(scoreLocation),
                     static_cast<double>(lowestObservableRoundedScaledScore),
                     static_cast<double>(highestObservableRoundedScaledScore));
        } else {
          roundedEquatedScaledScores(scoreLocation) = unroundedEquatedScaledScores(scoreLocation);
        }
      }
    }

    //----------------------------------------------------------------------------------------------------
    // Custom Function Written for EqRCpp
    //----------------------------------------------------------------------------------------------------
    static Eigen::VectorXd getRawScoreFrequencyDistribution(const Eigen::VectorXd& rawScores,
                                                            const double& minimumScore,
                                                            const double& maximumScore,
                                                            const double& scoreIncrement = 1,
                                                            const bool& includeRawScoresWithZeroFrequency = true) {
      size_t minimumScoreLocation = Utilities::getScoreLocation(minimumScore, minimumScore, scoreIncrement);
      size_t maximumScoreLocation = Utilities::getScoreLocation(maximumScore, minimumScore, scoreIncrement);
      size_t numberOfScores = Utilities::numberOfScores(minimumScore, maximumScore, scoreIncrement);

      Eigen::VectorXd freqDist = Eigen::VectorXd::Zero(numberOfScores);

      for (size_t scoreLocation = 0; scoreLocation < numberOfScores; ++scoreLocation) {
        double score = Utilities::getScore(scoreLocation,
                                           minimumScore,
                                           scoreIncrement);

        int scoreFreq = static_cast<int>(rawScores.cwiseEqual(score).count());

        if (scoreFreq >= 1 || (scoreFreq == 0 && includeRawScoresWithZeroFrequency)) {
          freqDist(scoreLocation) = scoreFreq;
        }
      }

      return freqDist;
    }

    static std::string vectorXiToString(const Eigen::VectorXi& vec, const bool& asRow) {
      std::string value = "";

      if (vec.size() >= 1) {
        std::string sep = asRow ? ", " : "\n";

        value.append(fmt::format("{}", vec(0)));

        for (Eigen::Index index = 1; index < vec.size(); ++index) {
          value.append(fmt::format("{}{}", sep, vec(index)));
        }
      }

      return value;
    }

    static std::string vectorXdToString(const Eigen::VectorXd& vec, const bool& asRow) {
      std::string value = "";

      if (vec.size() >= 1) {
        std::string sep = asRow ? ", " : "\n";

        value.append(fmt::format("{}", vec(0)));

        for (Eigen::Index index = 1; index < vec.size(); ++index) {
          value.append(fmt::format("{}{}", sep, vec(index)));
        }
      }

      return value;
    }

    static std::string matrixXiToString(const Eigen::MatrixXi& mat) {
      std::string value = "";

      if (mat.rows() >= 1 && mat.cols() >= 1) {
        for (Eigen::Index rowIndex = 0; rowIndex < mat.rows(); ++rowIndex) {
          std::string rowStr = fmt::format("{}", mat(rowIndex, 0));

          for (Eigen::Index columnIndex = 1; columnIndex < mat.cols(); ++columnIndex) {
            rowStr.append(fmt::format(", {}", mat(rowIndex, columnIndex)));
          }

          value.append(fmt::format("{}\n", rowStr));
        }
      }
      return value;
    }

    static std::string matrixXdToString(const Eigen::MatrixXd& mat) {
      std::string value = "";

      if (mat.rows() >= 1 && mat.cols() >= 1) {
        for (Eigen::Index rowIndex = 0; rowIndex < mat.rows(); ++rowIndex) {
          std::string rowStr = fmt::format("{}", mat(rowIndex, 0));

          for (Eigen::Index columnIndex = 1; columnIndex < mat.cols(); ++columnIndex) {
            rowStr.append(fmt::format(", {}", mat(rowIndex, columnIndex)));
          }

          value.append(fmt::format("{}\n", rowStr));
        }
      }

      return value;
    }

    static Eigen::VectorXd percentileRanks(const double& minimumScore,
                                           const double& maximumScore,
                                           const double& scoreIncrement,
                                           const Eigen::VectorXd& cumulativeRelativeFreqDist) {
      size_t numberOfScores = Utilities::numberOfScores(minimumScore,
                                                        maximumScore,
                                                        scoreIncrement);

      Eigen::VectorXd percRanks = Eigen::VectorXd::Zero(numberOfScores);

      for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
        double score = Utilities::getScore(scoreLocation,
                                           minimumScore,
                                           scoreIncrement);

        percRanks(scoreLocation) = Utilities::getPercentileRank(minimumScore,
                                                                maximumScore,
                                                                scoreIncrement,
                                                                cumulativeRelativeFreqDist,
                                                                score);
      }
      return percRanks;
    }

    static double getFirstObservedScore(const Eigen::VectorXd& scoreFrequencies,
                                        const double& minimumScore,
                                        const double& maximumScore,
                                        const double& scoreIncrement,
                                        const bool& searchForward = true) {
      double firstObservedScore = std::numeric_limits<double>::quiet_NaN();
      size_t scoreLocation = scoreFrequencies.size();

      if (searchForward) {
        auto iter = std::find_if(scoreFrequencies.begin(),
                                 scoreFrequencies.end(),
                                 [](const double& scoreFrequency) {
                                   return scoreFrequency > 0.0;
                                 });

        scoreLocation = std::distance(scoreFrequencies.begin(), iter);
      } else {
        for (size_t index = scoreFrequencies.size() - 1; 0; index--) {
          if (scoreFrequencies(index) > 0.0) {
            scoreLocation = index;
          }
        }
      }

      if (scoreLocation < scoreFrequencies.size()) {
        firstObservedScore = Utilities::getScore(scoreLocation,
                                                 minimumScore,
                                                 scoreIncrement);
      }

      return firstObservedScore;
    }

    static Eigen::FullPivLU<Eigen::MatrixXd> getFullPivotLUDecomposition(const Eigen::MatrixXd& mat) {
      Eigen::FullPivLU<Eigen::MatrixXd> lu(mat);

      return lu;
    }

    static Eigen::PartialPivLU<Eigen::MatrixXd> getPartialPivotLUDecomposition(const Eigen::MatrixXd& mat) {
      Eigen::PartialPivLU<Eigen::MatrixXd> lu(mat);

      return lu;
    }

    static EquatingRecipes::Structures::Method getMethod(const EquatingRecipes::Structures::Design& design,
                                                         const std::string& methodCode) {
      //   'M' = mean; 'L' = lin; 'E' equi (for RG and CG designs);
      // 	 'For CG design,
      // 		  'E' = Freq esti FE) with Braun-Holland (BH-FE)
      //      'F' = Modified freq est (MFE) with Braun-Holland (BH-MFE)
      //      'G' = FE + BH-FE + MFE + BH-MFE
      //      'C' = Chained
      // 		  'H' = FE + BH-FE + Chained
      //      'A' = FE + BH-FE + MFE + BH-MFE + Chained
      if (methodCode == "M") {
        return EquatingRecipes::Structures::Method::MEAN;

      } else if (methodCode == "L") {
        return EquatingRecipes::Structures::Method::LINEAR;

      } else if (methodCode == "E" &&
                 design == EquatingRecipes::Structures::Design::RandomGroups) {
        return EquatingRecipes::Structures::Method::EQUIPERCENTILE;

      } else if (design == EquatingRecipes::Structures::Design::CommonItenNonEquivalentGroups) {
        if (methodCode == "E") {
          return EquatingRecipes::Structures::Method::FE_BH;

        } else if (methodCode == "F") {
          return EquatingRecipes::Structures::Method::MFE_BH;

        } else if (methodCode == "G") {
          return EquatingRecipes::Structures::Method::FE_BH_MFE_BH;

        } else if (methodCode == "C") {
          return EquatingRecipes::Structures::Method::CHAINED;

        } else if (methodCode == "H") {
          return EquatingRecipes::Structures::Method::FE_BH_CHAINED;

        } else if (methodCode == "A") {
          return EquatingRecipes::Structures::Method::FE_BH_MFE_BH_CHAINED;

        } else {
          throw std::runtime_error("Invalid method code for CINEG design.");
        }
      } else {
        throw std::runtime_error("Invalid method code.");
      }
    }

    static std::string getMethodCode(const EquatingRecipes::Structures::Method& method) {
      switch (method) {
        case EquatingRecipes::Structures::Method::MEAN:
          return "M";
        case EquatingRecipes::Structures::Method::LINEAR:
          return "L";
        case EquatingRecipes::Structures::Method::EQUIPERCENTILE:
          return "E";
        case EquatingRecipes::Structures::Method::FE_BH:
          return "E";
        case EquatingRecipes::Structures::Method::MFE_BH:
          return "F";
        case EquatingRecipes::Structures::Method::FE_BH_MFE_BH:
          return "G";
        case EquatingRecipes::Structures::Method::CHAINED:
          return "C";
        case EquatingRecipes::Structures::Method::FE_BH_CHAINED:
          return "H";
        case EquatingRecipes::Structures::Method::FE_BH_MFE_BH_CHAINED:
          return "A";
        default:
          throw std::runtime_error("Invalid method.");
      }
    }

    static std::mt19937_64 getSeedEngine() {
      std::random_device rd;  // Will be used to obtain a seed for the random number engine
      std::mt19937_64 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

      return gen;
    }

    static std::uniform_int_distribution<> getUniformIntegerDistribution(const int& lowerBound,
                                                                         const int& upperBound) {
      std::uniform_int_distribution<> distrib(lowerBound,
                                              upperBound);

      return distrib;
    }

    static std::uniform_real_distribution<> getUniformRealDistribution() {
      std::uniform_real_distribution<> distrib(0.0, 1.0);

      return distrib;
    }

    static int getUniformIntegerRandomNumber(std::uniform_int_distribution<>& distrib,
                                             std::mt19937_64& gen) {
      int rndNumber = distrib(gen);
      return rndNumber;
    }

    static double getUniformDoubleRandomNumber(std::uniform_real_distribution<>& distrib,
                                               std::mt19937_64& gen) {
      double rndNumber = distrib(gen);
      return rndNumber;
    }
  };
} // namespace EquatingRecipes

#endif