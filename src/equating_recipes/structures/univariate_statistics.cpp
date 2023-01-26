#include <equating_recipes/structures/univariate_statistics.hpp>

namespace EquatingRecipes {
  namespace Structures {
    UnivariateStatistics UnivariateStatistics::create(std::map<double, int> scoreFreqDist,
                                         double minimumScore,
                                         double maximumScore,
                                         double scoreIncrement,
                                         std::string id) {
// s->id = id;
// s->min = min;
// s->max = max; 
// s->inc = inc;
// s->ns = nscores(max,min,inc);

// s->fd = ivector(0,s->ns-1);                            /* allocations */
// s->dbl_fd = dvector(0,s->ns-1);
// s->cfd = ivector(0,s->ns-1);
// s->rfd = dvector(0,s->ns-1); 
// s->crfd = dvector(0,s->ns-1); 
// s->prd = dvector(0,s->ns-1);   

// fp = fopen(fname,"r");
// if (fp==NULL) runerror("No file open\n");

// for(i=0;i<=s->ns-1;i++) s->fd[i] = 0;                     /* initialize */
// for(;;) {
//   for(i=1;i<=scol-1;i++) 
//     if(fscanf(fp,"%*s") == EOF) break;                     /* skip to x */
//   if(fscanf(fp,"%lf",&x) == EOF) break;                       /* read x */
  
//   for(i=1;i<=fcol-scol-1;i++) 
//     if(fscanf(fp,"%*s") == EOF) break;                     /* skip to f */
//   if(fscanf(fp,"%d",&f) == EOF) break;                        /* read f */

//   if(x<min) {
//     printf("In file %s, there is at least one ",fname);
//     runerror("score < min");
//   }
//   if(x>max){
//     printf("In file %s, there is at least one ",fname);
//     runerror("score > max");
//   }
  
//   s->fd[loc(x,min,inc)] = f;                    /* freq for score of x */
//   n += f;
//   flushline(fp);  
// }
// for(i=0;i<=s->ns-1;i++) if(s->fd[i]) break; 
// s->mind = score(i,min,inc);
// for(i=s->ns-1;i>=0;i--) if(s->fd[i]) break; 
// s->maxd = score(i,min,inc);
// if((s->n = MomentsFromFD(min, max, inc, NULL, s->fd, s->mts)) != n)
//   runerror("\nError somewhere in ReadFdGet_USTATS()");

// s->cfd[0] = s->fd[0];
// for(i=1;i<=loc(max,min,inc);i++) s->cfd[i] = s->cfd[i-1] + s->fd[i];
// for(i=0;i<=s->ns-1;i++) s->rfd[i] = s->fd[i]/((double) n);
// cum_rel_freqs(min,max,inc,s->rfd,s->crfd);

// for(i=0;i<=loc(max,min,inc);i++){
//   s->prd[i] = perc_rank(min,max,inc,s->crfd,score(i,min,inc));
//   s->dbl_fd[i] = s->fd[i];
// }
    }
  }
}