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

#include <map>
#include <Eigen/Core>

namespace EquatingRecipes {
  namespace Structures {
    struct Moments {
      Eigen::VectorXd momentValues;
      double minimumObservedScore;
      double maximumObservedScore;
      size_t numberOfPersons;

      static Moments getScoreMoments(const Eigen::VectorXd& scores);
      static Moments getScoreMoments(const std::map<double, int>& scoreFreqDist);
    };
  }
}

#endif