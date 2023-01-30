/* 
  From Source: ERutilities.h, ERutilities.c 
  Original Method: ReadRawGet_moments 
  Description: Compute moments from raw scores in file;
    scores need not be integers or positive;
    frequency distribution NOT computed or output
    assumes space allocated for moments[4];
    assumes data are whitespaced delimited
  
    Input
      fname[] = name of file for input
      scol =  column for scores (whitespace delimited columns)
    
    Output
      moments[] = mean, sd, skew, kurt (space already allocated)
      mind = lowest score (double) in file
      maxd = highest score (double) in file
      
    Returns n = number of persons

    Function calls other than C or NR utilities:
      flushline()
                                                
    R. L. Brennan

    Date of last revision: 6/30/08  
*/

#ifndef STRUCTURES_MOMENTS_HPP
#define STRUCTURES_MOMENTS_HPP

#include <cmath>
#include <limits>
#include <map>

#include <Eigen/Core>

#include <equating_recipes/structures/moments.hpp>
#include <equating_recipes/utilities.hpp>

namespace EquatingRecipes {
  namespace Structures {
    struct Moments {
      Eigen::VectorXd momentValues;
      double minimumObservedScore;
      double maximumObservedScore;
      size_t numberOfExaminees;

      // template <typename Derived>
      // void print_size(const Eigen::EigenBase<Derived>& b)
      // {
      //   std::cout << "size (rows, cols): " << b.size() << " (" << b.rows()
      //             << ", " << b.cols() << ")" << std::endl;
      // }

      static Moments fromScores(const Eigen::VectorXd& scores) {
        Moments scoreMoments;

        scoreMoments.numberOfExaminees = scores.size();

        scoreMoments.minimumObservedScore = scores.minCoeff();
        scoreMoments.maximumObservedScore = scores.maxCoeff();

        scoreMoments.momentValues.setZero(4);
        scoreMoments.momentValues(0) = scores.mean();

        Eigen::VectorXd meanVector = Eigen::VectorXd::Constant(scores.size(), scoreMoments.momentValues(0));
        Eigen::VectorXd deviations = scores - meanVector;

        double variance = deviations.array().pow(2).mean();
        double skewness = deviations.array().pow(3).mean();
        double kurtosis = deviations.array().pow(4).mean();

        scoreMoments.momentValues(1) = std::sqrt(variance);
        scoreMoments.momentValues(2) = skewness / std::pow(scoreMoments.momentValues(1), 3);
        scoreMoments.momentValues(3) = kurtosis / std::pow(scoreMoments.momentValues(1), 4);

        return scoreMoments;
      }

      static Moments fromScoreFrequencies(const Eigen::VectorXi& scoreFrequencies,
                                          const double& minimumScore,
                                          const double& maximumScore,
                                          const double& scoreIncrement) {
        Moments scoreMoments;

        size_t numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore,
                                                                           maximumScore,
                                                                           scoreIncrement);

        size_t maximumScoreLocation = EquatingRecipes::Utilities::getScoreLocation(maximumScore,
                                                                                   minimumScore,
                                                                                   scoreIncrement);

        scoreMoments.minimumObservedScore = std::numeric_limits<double>::max();
        scoreMoments.maximumObservedScore = std::numeric_limits<double>::min();
        scoreMoments.numberOfExaminees = scoreFrequencies.sum();
        scoreMoments.momentValues.setZero(4);

        Eigen::VectorXd deviations(numberOfScores);

        for (size_t scoreLocation = 0; scoreLocation <= maximumScoreLocation; scoreLocation++) {
          double score = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                              minimumScore,
                                                              scoreIncrement);

          int scoreFreq = scoreFrequencies(scoreLocation);

          deviations(scoreLocation) = score;

          if (scoreFreq >= 1) {
            scoreMoments.minimumObservedScore = std::min(score, scoreMoments.minimumObservedScore);
            scoreMoments.maximumObservedScore = std::max(score, scoreMoments.maximumObservedScore);
          }

          scoreMoments.momentValues(0) += score * static_cast<double>(scoreFreq);
        }

        scoreMoments.momentValues(0) /= static_cast<double>(scoreMoments.numberOfExaminees);

        deviations = deviations - Eigen::VectorXd::Constant(deviations.size(), scoreMoments.momentValues(0));

        for (size_t powCoeff = 2; powCoeff <= 4; powCoeff++) {
          scoreMoments.momentValues(powCoeff - 1) = (deviations.array().pow(powCoeff) * scoreFrequencies.array().cast<double>()).sum() /
                                                    static_cast<double>(scoreMoments.numberOfExaminees);
        }

        scoreMoments.momentValues(1) = std::sqrt(scoreMoments.momentValues(1));
        scoreMoments.momentValues(2) /= std::pow(scoreMoments.momentValues(1), 3);
        scoreMoments.momentValues(3) /= std::pow(scoreMoments.momentValues(1), 4);

        return scoreMoments;
      }

      // Based on MomentsFromRFD in ERUtilities.h, ERUtilities.c
      static Moments fromScoreRelativeFrequencies(const Eigen::VectorXd& relativeFrequencies,
                                                  const double& minimumScore,
                                                  const double& maximumScore,
                                                  const double& scoreIncrement,
                                                  const Eigen::VectorXd& scores) {
        Moments moments;

        /*
      int i,
      ns = nscores(max,min,inc);
  double mean=0., var=0., skew=0., kurt=0., dev, dev2;
  double *dscores;
  
  if(scores==NULL){
    dscores = dvector(0,ns-1);
    for(i=0;i<=ns-1;i++) dscores[i] = score(i,min,inc);
  }
  else
    dscores = scores;
 
  for(i=0;i<=ns-1;i++) {  
    mean += dscores[i]*rfd[i];  
  }

  moments[0] = mean;
  
  for(i=0;i<=ns-1;i++) { 
    dev = dscores[i] - mean;
    dev2 = dev*dev;
    var += dev2*rfd[i];
    dev *= dev2*rfd[i];
    skew += dev;
    dev2 = dev2*dev2*rfd[i];
    kurt += dev2; 
  }
  
  moments[1] = sqrt(var);
  var *= moments[1];
  moments[2] = skew / var;
  var *= moments[1];
  moments[3] = kurt / var;
  
      */

        // size_t numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore);

        return moments;
      }
      
      static Moments fromScoreRelativeFrequencies(const Eigen::VectorXd& relativeFrequencies,
                                                  const double& minimumScore,
                                                  const double& maximumScore,
                                                  const double& scoreIncrement) {
        size_t numberOfScores = EquatingRecipes::Utilities::numberOfScores(minimumScore,
                                                                           maximumScore,
                                                                           scoreIncrement);

        Eigen::VectorXd scores(numberOfScores);

        for (size_t scoreLocation = 0; scoreLocation < numberOfScores; scoreLocation++) {
          scores(scoreLocation) = EquatingRecipes::Utilities::getScore(scoreLocation,
                                                                       minimumScore,
                                                                       maximumScore);
        }

        Moments moments = Moments::fromScoreRelativeFrequencies(relativeFrequencies,
                                                                minimumScore,
                                                                maximumScore,
                                                                scoreIncrement,
                                                                scores);

        return moments;
      }
    };
  } // namespace Structures
} // namespace EquatingRecipes

#endif