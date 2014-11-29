/*
      *** autosstem.hpp ***

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

  ANSI C and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o autoslic autoslic.cpp autosliccmd.cpp slicelib.o
                       tiffsubs.o  cfpix.o -lfftw3f

  Transmit an electron wave through a specimen using the
  multislce method with automatic slicing.  Read in the (x,y,z)
  coordinates of the whole specimen and break into slices
  on-the-fly.

   ANSI-C and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)
  see:   www.fftw.org

  to output the position averaged CBED pattern include the filename
  after the program name on the command line (in 2D image mode only)

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o autostem autostem.cpp slicelob.o
                      floatTIFF.o  cfpix.o -lfftw3f

  Calculate STEM images and line scans from a non-periodic
  distribution of (x,y,z) atomic coordinates using repetitive multislice
  
  The transmission function for each slice can take a lot of
  computer time.  This program propagates the STEM probe for many
  beam position through the specimen at the same time to avoid recalculating
  the specimen transmission function for each probe position. This
  requires a huge amount of memory.  In 1D line scan mode a 2D probe
  wave function is stored for all positions.  In 2D image mode the 2D
  probe wave functions for a whole line are stored simulataneously.

  this file is formatted for a tab size of 4 characters

  multithreaded code using openMP may be turned on/off with
  symbol USE_OPENMP (ignore undefined pragma warning if off)

  Ref:
  [1] Zhiheng Yu, Philip Batson, John Silcox, "Artifacts in Aberration
      Corrected ADF-STEM Imaging", Ultramicroscopy 96 (2003) p.275-284
  [2] J. M. LeBeau et al, "Position averaged convergent beam electron
      diff.: Theory and Applications", Ultramic. 110 (2010) p.118-125

  started from stemslic.c and autoslic.c 19-jul-1998 E. Kirkland
  convert to simultaneous transmission of many probes at the same
    time to reuse transmission functions  28-jul-1998 ejk
  finished 2D image mode (1D mode works OK) 29-jul-1998 ejk
  last modified 29-jul-1998 ejk
  add Mote-Carlo integration of source size 24-aug-1998 ejk
  fixed typo in random dither in y direction of probe position
    1-feb-1999 ejk
  fixed error in nprobes in image mode (should be nyout
      but it was nxout- at top question) 3-july-2001 ejk
  update memory allocation routines,
     change void main() to int main() for better portability,
     and add 5th order spherical aberration  3-mar-2004 ejk
  start modification to keep multiple thickness's 21-jul-2005 ejk
      - in working form 26-jul-2005 ejk
  put in faster sorting routine sortByZ() 6-feb-2006 ejk
  add some openMP multithreading 23-may-2006 ejk
  move openMP stuff into conditional compilation
     so I only need one version of the code,
     and move sortbyz() to sliclib 23-aug-2006 ejk
  add periodic() to put probe position back into the supercell
     and recal. x,y position without source size wobble
     for output in 1D mode 31-aug-2006 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better 14-sep-2006 ejk
  change range to start at -0.25*deltaz 26-sep-2006 ejk
  minor cosmetic changes 11-jul-2007 ejk
  minor fixes for openMP stuff 27-feb-2008 ejk
  move vzatomLUT() to slicelib.c and reformat for TAB size of 4 char
               11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  fix bug in multithreading (add more private var.) 
      5-nov-2008 ejk
  add confocal mode 31-jul-2009 E. Kirkland
  offset confocal collector by zslice (same ref as obj) 4-aug-2009 ejk
  fix param listing in confocal mode 3,6-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4 
         10-feb-2010 ejk
  start conversion to FFTW  10-feb-2010 ejk
     in working form 16-feb-2010 ejk
  move vz LUT init (for openMP) to top 22-feb-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add multipole aberrations to probe (but not collector yet)
       9-may-2011 ejk
  change detector max test to < from <= so adding many
     ADF detectors together will work 6-jul-2011 ejk
  add trap for undersampling the probe+aperture 3-mar-2012 ejk
  start conversion to floatTIFF and C++ 24-may-2012 ejk
  add option to save position averaged CBED 30-jun-2012 ejk
  convert to cfpix/fftw class from raw fftw 07-nov-2012 to 10-nov-2012 ejk
  fix problem with different probe size (in pixels) that was
     created when multipole aberr added 7-apr-2013 ejk
  start conversion to autostem class with separate cmd line front end
       27-jul-2013 ejk
  add CountBeams() 17-aug-2013 ejk
  change obj. apert. to radians in param[] to be consistent with other programs
      -was in mrad  22-aug-2013 ejk
  add string message 23-aug-2013 ejk
  add confocal parameters 24-aug-2013 ejk
  change probe to new and add delete and return status value in calculate
      14-sep-2013 ejk
  change RNG seed argument to reference so it get updated for successive calls
       an lverbose 21-sep-2013 ejk
  fix confocal detection 5-oct-2013 ejk
  remove tString() from here and use slicelib 5-dec-2013 ejk

  this file is formatted for a TAB size of 8 characters 
  
*/

#ifndef AUTOSTEM_HPP   // only include this file if its not already
#define AUTOSTEM_HPP   // remember that this has been included


#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <string>	// STD string class

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "slicelib.hpp"    /* misc. routines for multislice */

enum{ ADF=0, CONFOCAL=1};        // modes of collector

//------------------------------------------------------------------
class autostem{

public:
    
    autostem( );         // constructor functions
    
    ~autostem();        //  destructor function

    // (input) control flags- should be set externally
    //   things that don't fit in param[]
    int lwobble, l1d, lxzimage, lpacbed, lverbose;

    //  misc info that may be used in calling program
    long nbeamt;
    double totmin, totmax, xmin, ymin, xmax, ymax;

    void CountBeams( float param[], int &nbeamp, int &nbeampo, float &res, float &almax );

    //  main calculation
    int calculate( float param[], int multiMode, int natom, unsigned long *iseed,
        int Znum[], float xa[], float ya[], float za[], float occ[], float wobble[],
        double xi, double xf, double yi, double yf, int nxout, int nyout,
        double ThickSave[], int nThick,
        double almin[], double almax[], int collectorMode[], int ndetect,
        float ***pixr, float  **rmin, float **rmax,
        float **pacbedPix );

private:

        int NZMAX, NSMAX, NCMAX, NRMAX;

        int nx, ny, nxprobe, nyprobe, nslice;
        int natom;

        //  TRUE and FALSE seemed to be predefined in the GUI so change name slightly
        enum{ xFALSE=0, xTRUE=1};

        float BW;

        double ABERR, RMIN, RMAX;

        std::string sbuffer;

        void messageAST( std::string &smsg, int level = 0 );  // common error message handler

        cfpix *probe, trans;     // complex probes and transmission functions

        float zmin, zmax;
        float *kx, *ky, *kx2, *ky2, *kxp, *kyp, *kxp2, *kyp2;
        float *propxr, *propxi;   /* complex propagator vs x */
        float *propyr, *propyi;   /* complex propagator vs y */
        float *xp, *yp;
        float *xa, *ya, *za, *occ, *wobble;
        float *xa2, *ya2, *za2, *occ2;
        int *Znum, *Znum2;
        double wavlen, k2maxp, Cs3,Cs5, df,apert1, apert2, pi, keV;
        float ax, by, cz;                   /*  specimen dimensions */
        double *k2max, *k2min, deltaz;

        double periodic( double pos, double size );
        void STEMsignals( double x[], double y[], int npos, float p[],
                 int multiMode, double ***detect, int ndetect,
                 double ThickSave[], int nThick, double sum[], int collectorMode[] );
        void invert2D( float** pix, long nx, long ny );    /*   for CBED pix */

        void trlayer( const float x[], const float y[], const float occ[],
            const int Znum[], const int natom, 
            const float ax, const float by, const float kev,
            cfpix &trans, const long nx, const long ny,
            double *phirms, long *nbeams, const float k2max );

        /* extra for confocal mode */
        int doConfocal;
        float dfa2C, dfa2phiC, dfa3C, dfa3phiC;   /* astigmatism parameters */
        double *collectMin, *collectMax;
        double Cs3C, Cs5C, dfC, apert1C, apert2C; /* aberrations of collector lens */

}; // end autostem::

#endif