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

