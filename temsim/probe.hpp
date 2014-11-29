/*      probe.hpp

------------------------------------------------------------------------
Copyright 1998-2013 Earl J. Kirkland

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 

------------------------------------------------------------------------

  header file for probe.cpp
  class with member subroutines to
  calculate STEM coherent probe wavefunction

  put the partial cross section (integrated over the ADF detector
  angles) at a single pixels (position of the corresponding atom)
  and convolve with the point spread function (focused probe intensity)

  reference:

  [1] E. Kirkland, "Advanced Computing in Electron Microscopy",
        Plenum 1998, 2nd edit. Springer 2010
  
  this file is formatted for a tab size of 4 char

   start conversion to separate class file 1-jul-2013 ejk
   move invert2D() to here (was in incostem) 7-jul-2013 ejk
 */

#ifndef PROBE_HPP   // only include this file if its not already

#define PROBE_HPP   // remember that this has been included


#include <cstdio>  /* standard ANSI libraries */
#include <cstdlib>
#include <cmath>

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "slicelib.hpp"    /* misc. routines for multislice */

//------------------------------------------------------------------
class probe{

public:
    
    probe( );         // constructor functions
    
    ~probe();        //  destructor function

    int makeProbe(  cfpix &cpix, int nx, int ny, double xp, double yp,
    float p[], double wavlen, double k2max, double pixel, int multiMode,
    int ismoth, float kx[], float kx2[], float ky[], float ky2[]);

    double prbSize( float** pixsq, int nx, int ny,
        double xp, double yp, double ax, double by );

    void invert2D( float** pix, long nx, long ny );

private:

    void message(  char msg[],  int level = 0 );  // common error message handler

}; // end incostem::

#endif
