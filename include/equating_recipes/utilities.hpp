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

#include <map>
#include <string>
#include <Eigen/Core>

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
    static double percentileRank(const double& minimumScore,
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

    //----------------------------------------------------------------------------------------------------
    // Custom Function Written for EqRCpp
    //----------------------------------------------------------------------------------------------------
    static Eigen::VectorXi getRawScoreFrequencyDistribution(const Eigen::VectorXd& rawScores,
                                                            const double& minimumScore,
                                                            const double& maximumScore,
                                                            const double& scoreIncrement = 1,
                                                            const bool& includeRawScoresWithZeroFrequency = true) {
      size_t minimumScoreLocation = Utilities::getScoreLocation(minimumScore, minimumScore, scoreIncrement);
      size_t maximumScoreLocation = Utilities::getScoreLocation(maximumScore, minimumScore, scoreIncrement);
      size_t numberOfScores = Utilities::numberOfScores(minimumScore, maximumScore, scoreIncrement);

      Eigen::VectorXi freqDist = Eigen::VectorXi::Zero(numberOfScores);

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
      std::string sep = asRow ? ", " : "\n";

      value.append(fmt::format("{}", vec(0)));

      for (Eigen::Index index = 1; index < vec.size(); ++index) {
        value.append(fmt::format("{}{}", sep, vec(index)));
      }

      return value;
    }

    static std::string vectorXdToString(const Eigen::VectorXd& vec, const bool& asRow) {
      std::string value = "";
      std::string sep = asRow ? ", " : "\n";

      value.append(fmt::format("{}", vec(0)));

      for (Eigen::Index index = 1; index < vec.size(); ++index) {
        value.append(fmt::format("{}{}", sep, vec(index)));
      }

      return value;
    }
    static std::string matrixXiToString(const Eigen::MatrixXi& mat) {
      std::string value = "";

      for (Eigen::Index rowIndex = 0; rowIndex < mat.rows(); ++rowIndex) {
        std::string rowStr = fmt::format("{}", mat(rowIndex, 0));

        for (Eigen::Index columnIndex = 1; columnIndex < mat.cols(); ++columnIndex) {
          rowStr.append(fmt::format(", {}", mat(rowIndex, columnIndex)));
        }

        value.append(fmt::format("{}\n", rowStr));
      }

      return value;
    }

    static std::string matrixXdToString(const Eigen::MatrixXd& mat) {
      std::string value = "";

      for (Eigen::Index rowIndex = 0; rowIndex < mat.rows(); ++rowIndex) {
        std::string rowStr = fmt::format("{}", mat(rowIndex, 0));

        for (Eigen::Index columnIndex = 1; columnIndex < mat.cols(); ++columnIndex) {
          rowStr.append(fmt::format(", {}", mat(rowIndex, columnIndex)));
        }

        value.append(fmt::format("{}\n", rowStr));
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

        percRanks(scoreLocation) = Utilities::percentileRank(minimumScore,
                                                             maximumScore,
                                                             scoreIncrement,
                                                             cumulativeRelativeFreqDist,
                                                             score);
      }
      return percRanks;
    }
  };
} // namespace EquatingRecipes

#endif