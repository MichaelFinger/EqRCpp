/*

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

#include "ERutilities.h"
#include "ER.h"

int main(void)
{
  struct USTATS x;
  struct BSTATS xv;
  int fieldsACT[2][3] = {{0,37,38},{0,39,40}};
  FILE *outf;                                         
                  
  outf = fopen("Chap 2 out","w");

  /*************************************************************/
  
  /* Random Groups Design: Kolen and Brennan (2004)
     Chapter 2 example (see pp. 50-52) */

  ReadFdGet_USTATS("actmathfreq.dat",1,2,0,40,1,'X',&x);
  Print_USTATS(outf,"ACT Math X",&x);  

  /* Common-items Nonequivalent Groups Design: 
     Kolen and Brennan (2004) Chapter 4 example (see page 123)*/ 
 
  convertFtoW("mondatx.dat",2,fieldsACT,"mondatx-temp");
  ReadRawGet_BSTATS("mondatx-temp",1,2,0,36,1,0,12,1,'X','V',&xv);
  Print_BSTATS(outf,"xv",&xv,1);
 
  fclose(outf);
  return 0;
}
