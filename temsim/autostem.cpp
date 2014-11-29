/*      *** autostem.cpp ***

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

-----------------------------------------------------------------------------
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
  change RNG seed argument to referenece so it get updated for successive calls
       an lverbose 21-sep-2013 ejk
  fix confocal detection 5-oct-2013 ejk
  convert detector angles to radians 19-oct-2013 ejk
  remove tString() from here and use slicelib 5-dec-2013 ejk

    this file is formatted for a TAB size of 8 characters 
*/

#include <ctime>  /* ANSI C libraries used */

#include "autostem.hpp"   // header for this class

#include <sstream>	// string streams

//
#define USE_OPENMP      /* define to use openMP */

#define MANY_ABERR      //  define to include many aberrations 

//=============================================================
//---------------  creator and destructor --------------

autostem::autostem()
{
        BW= (2.0F/3.0F);  // bandwidth limit
        ABERR= 1.0e-4;    // max error for a,b

        NSMAX= 1000;   // max number of slices
        NCMAX= 1024;   // max characters in file names
        NZMAX= 103;    // max atomic number Z

        NRMAX=   100;    // number of in look-up-table in vzatomLUT
        RMIN=    0.01;   // r (in Ang) range of LUT for vzatomLUT()
        RMAX=    5.0;

        pi = 4.0 * atan( 1.0 );

        // init control flags
        lwobble = l1d = lxzimage = lpacbed = 0;
        doConfocal = 0;

        return;

};   //  end autostem::autostem()

autostem::~autostem()
{
}

//=============================================================
/*  calculate()

  input:
        param[] = image parameters (most will not be changed here)
        multimode = flag controlling multipole aberrations
        natomin = number of atoms
        x[],y[],z[] = atomic coord
        Znum[] atomic number of each atom
        occ[] = occupancy of each atomic site
        wobble[] = thermal oscillation amplitude
        iseed = random number seed - will be changed here

        xi,xf,yi,yf = range of output
        nxout, nyout = size of output in pixels 
                (nxout must be =1 for 1D mode)

        ThickSave[] = thickness levels to save
        nThick = number of thickness levels to save

        almin[], almax[] = detector angles in radians or Angst.
        collectorMode[] = type of detectors (ADF vs. confocal)
        ndetect = number of detector geometries

  mode flags:  lwobble, l1d, lxzimage, lpacbed

  output:
        pacbedPix[][] = (real) 2D image with position averaged CBED (nxprobe x nyprobe)
                                if lpacbed == 1

        pixr[][][] = (real) array of 2D output images
                        size; (ndetect*nThick) x nxout x nyout
                        indexed as pixr[id+it*ndetect][ix][iy]
                    for 1D line scan nxout=1 and nyout is length of line in pixels

        rmin[][], rmax[][] = range of output pixr indexed as [it][id]
                                size; nThick x ndetect

  return value is +1 for successs and negative for failure 
*/
int autostem::calculate( float param[], int multiMode, int natomin, unsigned long *iseed,
        int Znum[], float xa[], float ya[], float za[], float occ[], float wobble[],
        double xi, double xf, double yi, double yf, int nxout, int nyout,
        double ThickSave[], int nThick,
        double almin[], double almax[], int collectorMode[], int ndetect,
        float ***pixr, float  **rmin, float **rmax,
        float **pacbedPix )
{
    int ix, iy, i, idetect, iwobble, nwobble,
        nprobes, ip, it, nbeamp, nbeampo, ix2, iy2;

    float prr, pri, temp, temperature;

    double dxp, dyp;

    double scale, *x, *y, sum, *sums, w, ***detect,
       tctx, tcty, dx, dy, ctiltx, ctilty, sourcesize, sourceFWHM,
       vz, rsq, k2maxa, k2maxb, k2;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    cz = param[ pCZ ];
    nx = ToInt( param[ pNX ] );		// size of transmission function (in pixels)
    ny = ToInt( param[ pNY ] );
    keV = param[ pENERGY ];		// electron beam energy in keV
    df = param[pDEFOCUS];		// defocus
    ctiltx = param[ pXCTILT ];          // crystal tilt
    ctilty = param[ pYCTILT ];
    apert1 = param[ pOAPMIN ];		// obj. apert size in radians (min, max for annular apert)
    apert2 = param[ pOAPERT ];
    temperature = param[ pTEMPER ];		// temperature
    nwobble = ToInt( param[ pNWOBBLE ] );       // number config. to average
    deltaz = param[ pDELTAZ ];			// slice thickness
    sourceFWHM = param[ pSOURCE ];		// source size  - doesn't work so not really used

    nxprobe = ToInt( param[ pNXPRB ] );  // probe size in pixels
    nyprobe = ToInt( param[ pNYPRB ] );
    
    //  confocal collector lens parameters
    dfC = param[ pCDF ];		// collector defocus
    dfa2C = param[ pCDFA2 ];		// collector astig 2nd order
    dfa2phiC = param[ pCDFA2PHI ];
    dfa3C = param[ pCDFA3 ];		// collector astig 2nd order
    dfa3phiC = param[ pCDFA3PHI ];
    Cs3C = param[ pCCS3 ];		// collector spherical aberr.
    Cs5C =  param[ pCCS5 ];
    apert1C = param[ pCCAPMIN ];	//  collector apert. in radians
    apert2C = param[ pCCAPMAX ];

    natom = natomin;
    wavlen = wavelength( keV );

    //  this code now works more consistently for both 1D and 2D
    nprobes = nyout;
    if( nxout < 1 ) {
        sbuffer = "nxout must be > 1 in autostem but it is "+toString(nxout);
        messageAST( sbuffer, 2 );
        return( -1 );
    }
    if( nyout < 1 ) {
        sbuffer = "nyout must be > 1 in autostem but it is "+toString(nyout);
        messageAST( sbuffer, 2 );
        return( -2 );
    }

    //  remember that nwobble must be at least one even if there is no phonon wobble
    if( nwobble < 1 ) nwobble = 1;

    /* convert FWHM to standard deviation 
            by dividing by 2*sqrt(2*ln(2)) 
        --  monte carlo of poistion does NOT converge very well for source size so
            do NOT use this, but leave the code here in case I ever figure out a better way */
    sourcesize = sourceFWHM / 2.354820045;

    //  see if confocal needed
    doConfocal = xFALSE;
    for( i=0; i<ndetect; i++)
            if( CONFOCAL == collectorMode[i]  ) doConfocal = xTRUE;

#ifdef USE_OPENMP
    /*  force LUT init. to avoid redundant init in parallel form */ 
    rsq = 0.5;  /* arbitrary position */   
    for( i=0; i<natom; i++) vz =  vzatomLUT( Znum[i], rsq );
#endif

    if( lwobble == 0 ) {
        sbuffer = "Sorting atoms by depth...";
        messageAST( sbuffer, 0 );
        sortByZ( xa, ya, za, occ, Znum, natom );
    }
    /* to add random offsets */
    xa2 = (float*) malloc1D( natom, sizeof(float),  "xa2" );
    ya2 = (float*) malloc1D( natom, sizeof(float), "ya2" );
    za2 = (float*) malloc1D( natom, sizeof(float), "za2" );
    Znum2 = (int*) malloc1D( natom, sizeof(int), "Znum2" );
    occ2 = (float*) malloc1D( natom, sizeof(float), "occ2" );

#ifdef NOT_HERE
    /*  calculate the total specimen volume and echo */
    xmin = xmax = xa[0];
    ymin = ymax = ya[0];
    zmin = zmax = za[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( xa[i] < xmin ) xmin = xa[i];
        if( xa[i] > xmax ) xmax = xa[i];
        if( ya[i] < ymin ) ymin = ya[i];
        if( ya[i] > ymax ) ymax = ya[i];
        if( za[i] < zmin ) zmin = za[i];
        if( za[i] > zmax ) zmax = za[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    printf("Total specimen range is\n %g to %g in x\n"
           " %g to %g in y\n %g to %g in z\n", xmin, xmax,
           ymin, ymax, zmin, zmax );
    if( lwobble == 1 )
        printf("Range of thermal rms displacements (300K) = %g to %g\n",
            wmin, wmax );
    /*  check for valid scan coordinates  */

    if( (xi < 0.0) || (xi > ax) ||
        (xf < 0.0) || (xf > ax) ||
        (yi < 0.0) || (yi > by) ||
        (yf < 0.0) || (yf > by) ) {
            sbuffer = "WARNING: Coordinates out of range; will be made periodic.\n xi,xf,yi,yf= "
              + toString(xi)+", "+toString(xf)+", "+ toString(yi)+", "+ toString(yf);
            messageAST( sbuffer, 0 );
    }
#endif

    /*  check that requested probe size is not bigger 
        than transmission function size (or too small)
    */
    if( (nxprobe > nx) || (nxprobe < 2) ) {
        nxprobe = nx;
        sbuffer = "Probe size reset to nx= " + toString( nxprobe);
        messageAST( sbuffer, 0 );
    }

    if( (nyprobe > ny) || (nyprobe < 2) ) {
        nyprobe = ny;
        sbuffer = "probe size reset to ny= " + toString( nyprobe );
        messageAST( sbuffer, 0 );
    }

    /*  calculate spatial frequencies for future use
        (one set for transmission function and one for probe
        wavefunction)
    NOTE: zero freg is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for FFT - don't waste time rearranging

    remember : the x values are the same for both sets
    
    x2, y2 are used for confocal
    
    */

    kx  = (float*) malloc1D( nx, sizeof(float), "kx" );
    ky  = (float*) malloc1D( ny, sizeof(float), "ky" );
    kx2 = (float*) malloc1D( nx, sizeof(float), "kx2" );
    ky2 = (float*) malloc1D( ny, sizeof(float), "ky2" );
    xp  = (float*) malloc1D( nx, sizeof(float), "x2" );
    yp  = (float*) malloc1D( ny, sizeof(float), "y2" );

    freqn( kx, kx2, xp, nx, ax );
    freqn( ky, ky2, yp, ny, by );

    kxp  = (float*) malloc1D( nxprobe, sizeof(float), "kxp" );
    kyp  = (float*) malloc1D( nyprobe, sizeof(float), "kyp" );
    kxp2 = (float*) malloc1D( nxprobe, sizeof(float), "kxp2" );
    kyp2 = (float*) malloc1D( nyprobe, sizeof(float), "kyp2" );
    
    freqn( kxp, kxp2, xp, nxprobe, dxp = ax*((double)nxprobe)/nx );
    freqn( kyp, kyp2, yp, nyprobe, dyp = by*((double)nyprobe)/ny );

    /* impose anti-aliasing bandwidth limit on transmission functions */

    sum = ((double)nx)/(2.0*ax);
    k2maxp = ((double)ny)/(2.0*by);
    if( sum < k2maxp ) k2maxp = sum;
    k2maxp= BW * k2maxp;
    //printf("Bandwidth limited to a real space resolution of %f Angstroms\n",   //???
    //                 1.0F/k2maxp);
    //printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
    //     wavlen*k2maxp*1000.0F);
    k2maxp = k2maxp * k2maxp;

    /*  allocate some more arrays and initialize propagator */

    propxr = (float*) malloc1D( nxprobe, sizeof(float), "propxr" );
    propxi = (float*) malloc1D( nxprobe, sizeof(float), "propxi" );
    propyr = (float*) malloc1D( nyprobe, sizeof(float), "propyr" );
    propyi = (float*) malloc1D( nyprobe, sizeof(float), "propyi" );

    /* calculate propagator functions with probe sample size
        impose anti-aliasing bandwidth limit  */
    tctx = 2.0 * tan(ctiltx);
    tcty = 2.0 * tan(ctilty);

    scale = pi * deltaz;
    for( ix=0; ix<nxprobe; ix++) {
        w = scale * ( kxp2[ix] * wavlen - kxp[ix]*tctx );
        propxr[ix]= (float)  cos(w);
        propxi[ix]= (float) -sin(w);
    }

    for( iy=0; iy<nyprobe; iy++) {
        w = scale * ( kyp2[iy] * wavlen - kyp[iy]*tcty );
        propyr[iy]= (float)  cos(w);
        propyi[iy]= (float) -sin(w);
    }

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 /wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavlen;
    k2maxb = k2maxb * k2maxb;
    
    nbeamp = nbeampo = 0;
    for( iy=0; iy<nyprobe; iy++)
    for( ix=0; ix<nxprobe; ix++) {
        k2 = kyp2[iy] + kxp2[ix];
        if( k2 < k2maxp ) nbeamp++;
        if( (k2 >= k2maxa) && (k2 <= k2maxb) ) nbeampo++;
    }

    //  output this in outer calling program if wanted - not every time here
    //sbuffer = "Number of symmetrical anti-aliasing beams in probe= "
    //       + toString(nbeamp);
    //messageAST( sbuffer, 0 );
    //sbuffer = "Number of beams in probe aperture= " + toString(nbeampo);
    //messageAST( sbuffer, 0 );

    if( nbeamp < 200 ) {
        sbuffer = "WARNING: the probe is under sampled, this is a bad idea...";
        messageAST( sbuffer, 1 );
    }
    if( nbeampo < 100 ) {
        sbuffer = "WARNING: the probe aperture is under sampled, this is a bad idea...";
        messageAST( sbuffer, 2 );
        exit( EXIT_FAILURE );
    }

    /*  convert aperture dimensions */

    k2min = (double*) malloc1D( ndetect, sizeof(double), "k2min" );
    k2max = (double*) malloc1D( ndetect, sizeof(double), "k2max" );

    for( idetect=0; idetect<ndetect; idetect++) {
        if( ADF == collectorMode[idetect] ) {
            k2max[idetect] = almax[idetect]/wavlen;
            k2max[idetect] = k2max[idetect] * k2max[idetect];
            k2min[idetect] = almin[idetect]/wavlen;
            k2min[idetect] = k2min[idetect] * k2min[idetect];
        } else if( CONFOCAL == collectorMode[idetect] ) {
            k2max[idetect] = almax[idetect] * almax[idetect];
            k2min[idetect] = almin[idetect] * almin[idetect];
        }
    }

    /*  init the min/max record of total integrated intensity */

    totmin =  10.0;
    totmax = -10.0;
    detect  = (double***) malloc3D( nThick, ndetect, nprobes,
        sizeof(double), "detect" );
    sums = (double*) malloc1D( nprobes, sizeof(double), "sums" ); 

    // allocate probe wave function and transmission function arrays

    probe = new cfpix[ nprobes ];
    if( NULL == probe ) {
        sbuffer = "Cannot allocate probe array";
        messageAST( sbuffer, 2 );
        exit( EXIT_FAILURE );
    }
    ix = probe[0].resize(nxprobe, nyprobe );
    if( ix < 0 ) {
        sbuffer = "Cannot allocate probe array storage";
        messageAST( sbuffer, 2 );
        exit( EXIT_FAILURE );
    }
    probe[0].init();
    if( nprobes > 1 ) for( ip=1; ip<nprobes; ip++){
        ix = probe[ip].resize(nxprobe, nyprobe );
        if( ix < 0 ) {
            sbuffer = "Cannot allocate probe array storage";
            messageAST( sbuffer, 2 );
            exit( EXIT_FAILURE );    //  should do something better here ??
        }
        probe[ip].copyInit( probe[0] );
    }

    trans.resize( nx, ny );
    trans.init();

    if( lpacbed == xTRUE ) {
        for( ix=0; ix<nxprobe; ix++) for( iy=0; iy<nyprobe; iy++)
                pacbedPix[ix][iy] = 0;
    }

/* ------------- start here for a full image output -------------- */
/*
  do one whole line at once NOT the whole image (which may be huge)
*/
    if( l1d == 0 ) {
       sbuffer = "output file size in pixels is " + toString(nxout) +" x " 
               + toString(nyout);
       messageAST( sbuffer, 0 );
       if( nprobes != nyout ) {
           sbuffer = "Error, nprobes= "+ toString(nprobes)+" must be the same as nyout= "
               +toString(nyout)+", in image mode.";
           messageAST( sbuffer, 2 );
           exit( 0 );
       }
       /* double up first index to mimic a 4D array */
       for( i=0; i<(nThick*ndetect); i++) {
           for( ix=0; ix<nxout; ix++)
           for( iy=0; iy<nyout; iy++)
            pixr[i][ix][iy] = 0.0F;
       }

       /*  iterate the multislice algorithm proper for each
           position of the focused probe */

       if( nxout > 1 ) dx = (xf-xi)/((double)(nxout-1));
       else dx = 1.0;
       if( nyout > 1 ) dy = (yf-yi)/((double)(nyout-1));
       else dy = 1.0;
       x = (double*) malloc1D( nprobes, sizeof(double), "x" );
       y = (double*) malloc1D( nprobes, sizeof(double), "y" );

        /*  add random thermal displacements 
               scaled by temperature if requested 
            remember that initial wobble is at 300K for
               each direction */

        for( iwobble=0; iwobble<nwobble; iwobble++) {
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    za2[i] = za[i] + 
                            (float)(wobble[i]*rangauss(iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                sbuffer = "configuration # " + toString( iwobble+1 );
                messageAST( sbuffer, 0 );
                sbuffer = "The new range of z is "
                    + toString(za2[0]) + " to " + toString( za2[natom-1] );
                messageAST( sbuffer, 0 );
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];  /* reset zmin/max after wobble */
            zmax = za2[natom-1];
    
            for( ix=0; ix<nxout; ix++) {
    
                for( iy=0; iy<nyout; iy++) {
                    x[iy] = xi + dx * ((double) ix);
                        //  + sourcesize * rangauss(iseed);  - does not converge well
                    y[iy] = yi + dy * ((double) iy);
                        //  + sourcesize * rangauss(iseed);  - does not converge well
                    x[iy] = periodic( x[iy], ax );   /* put back in supercell */
                    y[iy] = periodic( y[iy], by );   /* if necessary */
                }

                STEMsignals( x, y, nyout, param, multiMode, detect, ndetect, 
                    ThickSave, nThick, sums, collectorMode );
                for( iy=0; iy<nyout; iy++) {
                    if( sums[iy] < totmin ) totmin = sums[iy];
                    if( sums[iy] > totmax ) totmax = sums[iy];
                    for( it=0; it<nThick; it++){
                        for( idetect=0; idetect<ndetect; idetect++)
                        pixr[idetect + it*ndetect][ix][iy] += (float)
                            (detect[it][idetect][iy]/((double)nwobble));
                    }
                    if( sums[iy] < 0.9) {
                        sbuffer = "Warning integrated intensity too small, = "
                           + toString(sums[iy])+" at "+toString(x[iy])+", "+toString(y[iy]);
                        messageAST( sbuffer, 0 );
                    }
                    if( sums[iy] > 1.1) {
                        sbuffer =  "Warning integrated intensity too large, = "
                           + toString(sums[iy])+" at "+toString(x[iy])+", "+toString(y[iy]);
                        messageAST( sbuffer, 0 );
                    }
                }

                /*   sum position averaged CBED if requested 
                     - assume probe still left from stemsignal()  */
                if( lpacbed == xTRUE ) {
                    for( iy=0; iy<nyout; iy++) {
                        for( ix2=0; ix2<nxprobe; ix2++)
                        for( iy2=0; iy2<nyprobe; iy2++) {
                            prr = probe[iy].re(ix2,iy2);
                            pri = probe[iy].im(ix2,iy2);
                            pacbedPix[ix2][iy2] += (prr*prr + pri*pri);
                       }
                    }
                }   /*  end if( lpacbed.... */
            
            } /* end for(ix...) */
    
        } /* end for(iwobble... ) */

        /*  find range to output data files  */
        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            rmin[it][i] = rmax[it][i] = pixr[i+it*ndetect][0][0];
            for( ix=0; ix<nxout; ix++)
            for( iy=0; iy<nyout; iy++) {
                temp = pixr[i+it*ndetect][ix][iy];
                if( temp < rmin[it][i] )rmin[it][i] = (float) temp;
                if( temp > rmax[it][i] )rmax[it][i] = (float) temp;
            }
        }
        if( lpacbed == xTRUE ) {
            invert2D( pacbedPix, nxprobe, nyprobe );  /*  put zero in middle */
         }

    /* ------------- start here for 1d line scan ---------------- */

    } else { // the only other posibility is; if ( l1d == 1 ) {

       if( lpacbed == xTRUE ) {
          sbuffer = "warning: cannot do pos. aver. CBED in 1d";
          messageAST( sbuffer, 0 );
       }
       if( nxout > 1 ) {
               sbuffer="nxout must be 1 in 1D mode but is "+toString(nxout);
               messageAST( sbuffer, 0 );
       }
       if( nyout > 1 ) dx = (xf-xi)/((double)(nyout-1));
       else dx = 1.0;
       if( nyout > 1 ) dy = (yf-yi)/((double)(nyout-1));
       else dy = 1.0;
       x = (double*) malloc1D( nprobes, sizeof(double), "x" );
       y = (double*) malloc1D( nprobes, sizeof(double), "y" );
       for( ip=0; ip<nyout; ip++) {
            for( it=0; it<nThick; it++)
            for( idetect=0; idetect<ndetect; idetect++)
                    pixr[idetect+it*ndetect][0][ip] = 0.0F;
       }

       /*  add random thermal displacements scaled by temperature
            if requested 
        remember that initial wobble is at 300K for each direction */

       for( iwobble=0; iwobble<nwobble; iwobble++) {
    
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    za2[i] = za[i] + 
                        (float)(wobble[i]*rangauss(iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                sbuffer = "configuration # " + toString( iwobble+1 );
                messageAST( sbuffer, 0 );
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                sbuffer = "The new range of z is "
                    + toString(za2[0]) + " to " + toString( za2[natom-1] );
                messageAST( sbuffer, 0 );
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];      /* reset zmin/max after wobble */
            zmax = za2[natom-1];
            for( ip=0; ip<nyout; ip++) {
                x[ip] = xi + dx * ((double)ip);
                            //  + sourcesize * rangauss(iseed);  - does not converge well
                y[ip] = yi + dy * ((double)ip);
                            //  + sourcesize * rangauss(iseed);  - does not converge well
                x[ip] = periodic( x[ip], ax );   // put back in supercell
                y[ip] = periodic( y[ip], by );   // if necessary
            }
         
            STEMsignals( x, y, nprobes, param, multiMode, detect, ndetect, 
                ThickSave, nThick, sums, collectorMode );
            for( ip=0; ip<nprobes; ip++) {
                if( sums[ip] < totmin ) totmin = sums[ip];
                if( sums[ip] > totmax ) totmax = sums[ip];
                for( it=0; it<nThick; it++){
                for( idetect=0; idetect<ndetect; idetect++)
                   //  nwobble should be small so its prob. safe to sum into single prec. var.
                   pixr[idetect+it*ndetect][0][ip] += 
                        (float) ( detect[it][idetect][ip]/((double)nwobble) );
                }
                if( sums[ip] < 0.9) {
                    sbuffer = "Warning integrated intensity too small, = "
                           + toString(sums[ip])+" at "+toString(x[ip])+", "+toString(y[ip]);
                    messageAST( sbuffer, 0 );
                }
                if( sums[ip] > 1.1) {
                    sbuffer =  "Warning integrated intensity too large, = "
                           + toString(sums[ip])+" at "+toString(x[ip])+", "+toString(y[ip]);
                    messageAST( sbuffer, 0 );
                }
            }

       }  /* end for(iwobble... */


    } /* end if( l1d.. ) */

    //----------- end:  free scratch arrays and exit --------------------

    free( xa2 );
    free( ya2 );
    free( za2 );
    free( Znum2 );
    free( occ2 );

    free( kx );
    free( kx2 );
    free( xp );

    free( ky );
    free( ky2 );
    free( yp );

    free( kxp );
    free( kxp2 );

    free( kyp );
    free( kyp2 );

    free( propxr );
    free( propxi );
    free( propyr );
    free( propyi );

    free( k2min );
    free( k2max );

    free( x );
    free( y );

    free ( sums );

    free3D( (void***) detect, nThick, ndetect );

    delete [] probe;

    return( + 1 );
         
}  // end autostem::calculate()

/* -------------------  CountBeams() -------------------

   count the number of beams (fourier coefficients) in the
   transmission function and the probe for informational purposes
   - to communicate to the main calling program
   -  must be recalculated here (should be identical to that
      used in calculate()
  input:
        param[] = array of parameters

  output:
        nbeamp  = number of symmetrical beams in the probe
        nbeampo = number of beams in the probe aperture
        res  = bandwidth limited resolution in Angstroms
        almax = alpha max in radians
*/
void autostem::CountBeams( float param[], int &nbeamp, int &nbeampo, float &res, float &almax )
{
    //  use all local array variables so I don't accidentally
    //     disturb the main calculation
    int ix, iy, nx1, ny1, nxprobe1, nyprobe1;
    float *xp1, *yp1, *kxp1, *kyp1, *kxp21, *kyp21;
    float ax, by, wavl, apert1, apert2;
    double sum, k2maxp, k2maxa, k2maxb, k2;

    ax = param[ pAX ];			// supercell size in Ang.
    by = param[ pBY ];

    nx1 = ToInt( param[ pNX ] );	// size of transmission function (in pixels)
    ny1 = ToInt( param[ pNY ] );

    nxprobe1 = ToInt( param[ pNXPRB ] );	// probe size in pixels
    nyprobe1 = ToInt( param[ pNYPRB ] );

    wavl = (float) wavelength( param[ pENERGY ] );	//  electron wavelength in Ang.

    apert1 = param[ pOAPMIN ];	// obj. apert size in radians (min, max for annular apert)
    apert2 = param[ pOAPERT ];

    //  transmission function sampling
    // impose anti-aliasing bandwidth limit on transmission functions
    sum = ((double)nx1)/(2.0*ax);
    k2maxp = ((double)ny1)/(2.0*by);
    if( sum < k2maxp ) k2maxp = sum;
    k2maxp= BW * k2maxp;
    res = (float) ( 1.0/k2maxp);
    almax = (float) (wavl*k2maxp);
    k2maxp = k2maxp * k2maxp;

    //  probe sampling
    //  the only way to do this is to do part of the calculation and
    //    throw away the results (kx,ky etc.) - but only a small bit wasted
    xp1  = (float*) malloc1D( nxprobe1, sizeof(float), "xp1" );
    yp1  = (float*) malloc1D( nyprobe1, sizeof(float), "yp1" );
    kxp1  = (float*) malloc1D( nxprobe1, sizeof(float), "kxp1" );
    kyp1  = (float*) malloc1D( nyprobe1, sizeof(float), "kyp1" );
    kxp21 = (float*) malloc1D( nxprobe1, sizeof(float), "kxp21" );
    kyp21 = (float*) malloc1D( nyprobe1, sizeof(float), "kyp21" );
    
    freqn( kxp1, kxp21, xp1, nxprobe1, ax*((double)nxprobe1)/nx1 );
    freqn( kyp1, kyp21, yp1, nyprobe1, by*((double)nyprobe1)/ny1 );

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 /wavl;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavl;
    k2maxb = k2maxb * k2maxb;
    
    nbeamp = nbeampo = 0;
    for( iy=0; iy<nyprobe1; iy++)
    for( ix=0; ix<nxprobe1; ix++) {
        k2 = kyp21[iy] + kxp21[ix];
        if( k2 < k2maxp ) nbeamp++;
        if( (k2 >= k2maxa) && (k2 <= k2maxb) ) nbeampo++;
    }

    free( xp1 );
    free( kxp1 );
    free( kxp21 );

    free( yp1 );
    free( kyp1 );
    free( kyp21 );

};  // end autostem::CountBeams()


/* -------------------  messageAST() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void autostem::messageAST( std::string &smsg,  int level )
{
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

}  // end autoslic::messageAST()

/*------------------------ periodic() ---------------------*/
/*
     make probe positions periodic in the supercell
     in case some wobble off the edge with source size of user excess

    pos = input position (x or y);
    size = supercell size ( 0 to size)

    return positive value  0 <= x < size
*/
double autostem::periodic( double pos, double size )
{
    double x=pos;
    while( x < 0 ) x += size;
    x = fmod( x, size );
    return( x );
}

/*------------------------ STEMsignals() ---------------------*/
/*

  NOTE: this is NOT the same as STEMsignal() in stemslice.c

  subroutine to calculate the stem signal arising from a given
  probe position

  iterate the multislice algorithm proper for each position of
  the focused probe

  This version uses massive amounts of memory to avoid
  recalculating the transmission functions more than necessary

   note zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for fft - don't waste time rearranging

     change to propagate thru whole unit cell not just
     range of atoms on 14-sep-2006 ejk
     add multipole aberrations 9-may-2011 ejk
     change to cfpix for probe and trans 10-nov-2012 ejk

  x[],y[]     = real positions of the incident probe
  npos        = int number of positions
  param[]     = parameters of probe
  multiMode   = flag to add multipole aberrations
  detect[][][]= real array to get signal into each detector
            for each probe position and thickness
  ndetect     = number of detector geometries
  ThickSave[] = thicknesses at which to save data (other than the last)
  nThick      = number of thickness levels (including the last)
  sum         = real total integrated intensity
  
  the assumed global variables are:
  
  nxprobe,nyprobe   = int size of probe wavefunction in pixels
  nx,ny         = int size of transmission function in pixels
  trans                  = float complex transmission function
  propxr[][], propxi[][] = float real,imag propagator vs x
  propyr[][], propyi[][] = float real,imag propagator vs y
  ax,by,cz      = float unit cell size in Angs
  kxp[], kyp[]      = float spatial frequencies vs x, y
  kxp2[], kyp2[]    = float square of kxp[], kyp[]
  xp[], yp[]        = float real space positions in probe (confocal)
  apert1, apert2    = double min,max objective aperture (in radians)
  k2maxp            = double max spatial freq of probe squared
  pi                = double constant PI
  wavlen            = double electron wavelength in Angs
  df                = double defocus (in Ang)
  Cs3,Cs5           = double spherical aberration (in Ang)

  xa[],ya[],za[]    = atom coordinates
  occ[]         = atomic occupancy
  Znum[]        = atomic numbers
  natom         = number of atoms
  deltaz        = slice thickness
  v0            = beam energy
  nbeamt        = number of beams in transmission function
  zmin, zmax    = range of z coord. of the atoms
  nslice        = number of slices
  doConfocal    = flag indicating confocal is needed
  
    NOTE:  too many thing come in as globals, but...

*/

void autostem::STEMsignals( double x[], double y[], int npos, float p[],
         int multiMode, double ***detect, int ndetect,
         double ThickSave[], int nThick, double sum[], int collectorMode[] )
{
    int ix, iy, ixt, iyt, idetect, *ixoff, *iyoff, ixmid, iymid;
    int istart, na, ip, i, it;

    long nxprobel, nyprobel, nxl, nyl;

    float scale, prr, pri, tr, ti;

    double *xoff, *yoff, chi0, chi1, k2maxa, k2maxb,
        w, k2, phi, phirms, alx, aly;
    double sum0, sum1, delta, zslice, totalz;
    
    /* extra for confocal */
    float hr, hi;
    double chi2C, chi3C, k2maxaC, k2maxbC, r2, rx2;
    cfpix cpix;            /* complex confocal image */

    /* ------ make sure x,y are ok ------ */

    for( ip=0; ip<npos; ip++) {
        if( (x[ip] < 0.0) || (x[ip] > ax) ||
            (y[ip] < 0.0) || (y[ip] > by) ) {
            sum[ip] = -1.2345;
            sbuffer = "bad x=%f,y=%f in STEMsignals()\n" + toString(x[ip]) + toString(y[ip]);
            messageAST( sbuffer, 0 );
            return;
        }
    }

    /*  generate a probe at position x,y which must be inside the
        0->ax and 0->by region
        normalize to have unit integrated intensity
        probe position (x,y) = (0,0) = lower left corner

    NOTE: The probe wavefunction is nxprobe x nyprobe which may
        be smaller than the transmission function.
        Position the probe near the center of the nxprobe x nyprobe
        region because it does not wrap around properly
        (i.e. manually wrap around in the nx x ny region if
         the probe is near a boundary of the nx x ny region))
    */

    ixmid = nxprobe/2;
    iymid = nyprobe/2;
    chi1 = pi * wavlen;
    k2maxa = apert1 /wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 /wavlen;
    k2maxb = k2maxb * k2maxb;
    
    /* extra for confocal */
    chi2C = 0.5 * Cs3C *wavlen*wavlen;
    chi3C = Cs5C * wavlen*wavlen*wavlen*wavlen /3.0;
    k2maxaC = apert1C /wavlen;
    k2maxaC = k2maxaC *k2maxaC;
    k2maxbC = apert2C /wavlen;
    k2maxbC = k2maxbC * k2maxbC;

    ixoff = (int*) malloc1D( npos, sizeof(int), "ixoff" );
    iyoff = (int*) malloc1D( npos, sizeof(int), "iyoff" );
    xoff = (double*) malloc1D( npos, sizeof(double), "xoff" );
    yoff = (double*) malloc1D( npos, sizeof(double), "yoff" );

    /* ------- calculate all of the probe wave functions at once ------
        to reuse the transmission functions which takes a long
        time to calculate*/

/*  paralleling this loop has little effect */
#pragma omp parallel for private(ix,iy,sum0,k2,w,chi0,scale,tr,ti,alx,aly) 
    for( ip=0; ip<npos; ip++) {
        ixoff[ip] = (int) floor( x[ip]*((double)nx) / ax ) - ixmid;
        xoff[ip]  = x[ip] - ax*((double)ixoff[ip])/((double)nx);

        iyoff[ip] = (int) floor( y[ip]*((double)ny) / by ) - iymid;
        yoff[ip]  = y[ip] - by*((double)iyoff[ip])/((double)ny);

        sum0 = 0.0;
        for( ix=0; ix<nxprobe; ix++) {
            alx = wavlen * kxp[ix];  /* x component of angle alpha */
            for( iy=0; iy<nyprobe; iy++) {
                aly = wavlen * kyp[iy];  /* y component of angle alpha */
                k2 = kxp2[ix] + kyp2[iy];
                if( (k2 >= k2maxa) && (k2 <= k2maxb) ) {
                    w = 2.*pi* ( xoff[ip]*kxp[ix] + yoff[ip]*kyp[iy] );
                    chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode );
                    chi0= - chi0 + w;
                    probe[ip].re(ix,iy) = tr = (float) cos( chi0 );
                    probe[ip].im(ix,iy) = ti = (float) sin( chi0 );
                    sum0 += (double) (tr*tr + ti*ti);
                } else {
                    probe[ip].re(ix,iy) = 0.0F;
                    probe[ip].im(ix,iy) = 0.0F;
                }
            }
        }  /* end for( ix... */

        scale = (float) ( 1.0/sqrt(sum0) );  
        probe[ip] *= scale;

    }  /* end for( ip...) */

    /* -------- transmit thru nslice layers ------------------------
        beware ixoff,iyoff must be with one nx,ny
            of the array bounds
    */          

    nxprobel = (long) nxprobe;
    nyprobel = (long) nyprobe;

    nxl = (long) nx;
    nyl = (long) ny;
    
    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    nslice = 0;
    it = 0;     /* thickness level index */          

    if( zmax > cz ) totalz = zmax;
        else totalz = cz;
    sbuffer= "specimen range is 0 to "+ toString(totalz) + " Ang.";
    messageAST( sbuffer, 0 );
    
    /* range of unit cell */
    while(  (zslice < (totalz+0.25*deltaz)) || (istart<natom) ) {

       /* find range of atoms for current slice */
       na = 0;
       for(i=istart; i<natom; i++)
            if( za2[i] < zslice ) na++; else break;

       if( 0 != lverbose ) {
           sbuffer= "slice ending at z= "+toString(zslice)+" Ang. with "+toString(na)+" atoms";
           messageAST( sbuffer, 0 );
           //ss.str("");
           //ss << "slice ending at z= " << zslice << " Ang. with " << na << " atoms";
           //messageAST( ss.str(), 0 );
       }

       /* calculate transmission function and bandwidth limit */
       if( na > 0 ) trlayer( &xa2[istart], &ya2[istart], &occ2[istart],
            &Znum2[istart], na, (float)ax, (float)by, (float)keV,
            trans, nxl, nyl, &phirms, &nbeamt, (float) k2maxp );

            /*----- one multislice trans/prop cycle for all probes ---- */
#pragma omp parallel for private(ix,iy,ixt,iyt,prr,pri)
       for( ip=0; ip<npos; ip++) {
           /* apply transmission function if there are atoms in this slice */
           if( na > 0 ) {
                probe[ip].ifft();
                for( ix=0; ix<nxprobe; ix++) {
                    ixt = ix + ixoff[ip];
                    if( ixt >= nx ) ixt = ixt - nx;
                    else if( ixt < 0 ) ixt = ixt + nx;
                    for( iy=0; iy<nyprobe; iy++) {
                        iyt = iy + iyoff[ip];
                        if( iyt >= ny ) iyt = iyt - ny;
                        else if( iyt < 0 ) iyt = iyt + ny;
                        prr = probe[ip].re(ix,iy);
                        pri = probe[ip].im(ix,iy);
                        probe[ip].re(ix,iy) =  prr*trans.re(ixt,iyt)
                                              - pri*trans.im(ixt,iyt);  // real
                        probe[ip].im(ix,iy) =  prr*trans.im(ixt,iyt)
                                              + pri*trans.re(ixt,iyt);  // imag
                    } /* end for(iy...) */
                }  /* end for(ix...) */
                probe[ip].fft();
           }
    
            /*  multiplied by the propagator function */
            propagate( probe[ip], propxr, propxi, propyr, propyi,
                kxp2, kyp2, (float)k2maxp, nxprobe, nyprobe );

        }  /* end  for( ip=... */

        /*  if this is a good thickness level then save the ADF or confocal signals
           - remember that the last level may be off by one layer with
              thermal displacements so special case it
        */
       
        /*  look at all values because they may not be in order */
        for( it = 0; it<nThick; it++ ) 
        if( fabs(ThickSave[it]-zslice)<fabs(0.5*deltaz)) {
            
            sbuffer= "save ADF/confocal signals, thickness level " + toString(it);  // diagnostic
            messageAST( sbuffer, 0 );
    
            /*  loop over all probes again */
#pragma omp parallel for private(ix,iy,idetect,prr,pri,delta,k2,cpix,phi,chi0,hr,hi,sum0,sum1,rx2,r2)
            for( ip=0; ip<npos; ip++) {

                /*  zero sum count */
                sum[ip] = 0.0;
                for(ix=0; ix<ndetect; ix++) detect[it][ix][ip] = 0.0;
   
                /*  sum intensity incident on the ADF detector
                        and calculate total integrated intensity
                    - changed detector limits to >= min and < max
                       so many concentric ADF detectors sum correctly
                       7-jul-2011
                */
                for( ix=0; ix<nxprobe; ix++) {
                    for( iy=0; iy<nyprobe; iy++) {

                        prr = probe[ip].re(ix,iy);
                        pri = probe[ip].im(ix,iy);
                        delta = prr*prr + pri*pri;
                        sum[ip] += delta;
                        k2 = kxp2[ix] + kyp2[iy];
                        for( idetect=0; idetect<ndetect; idetect++) {
                            if( ADF == collectorMode[idetect] ) {
                                if( (k2 >= k2min[idetect] ) &&
                                (k2 < k2max[idetect] ) )
                                detect[it][idetect][ip] += delta;
                            }
                        }
                    } /* end for(iy..) */
                }  /* end for(ix...) */

                /*  transform back if confocal needed 
                    - use copy of probe so original can continue in use  */
                if( doConfocal == xTRUE ) {
                    /*  allocate/deallocate here so openMP will work 
                        otherwise have to allocate nxprobe cpix arrays 
                        - a littel slow but a lot less memory */
                    cpix.resize(nxprobe,nyprobe);
                    cpix.copyInit(probe[0]);
                    sum0 = 0;
                    for( ix=0; ix<nxprobe; ix++) {
                        for( iy=0; iy<nyprobe; iy++) {
                            k2 = kxp2[ix] + kyp2[iy];
                            if( (k2 >= k2maxaC) && (k2 <= k2maxbC) ) {
                                phi = atan2( ky[iy], kx[ix] );
                                /*  offset defocus by zslice so both lens referenced to 
                                   entrance surface of specimen  */
                                chi0 = chi1*k2* ( (chi2C + chi3C*k2)*k2 - dfC + zslice
                                    + dfa2C*sin( 2.0*(phi-dfa2phiC) ) 
                                    + 2.0F*dfa3C*wavlen*sqrt(k2)*
                                    sin( 3.0*(phi-dfa3phiC) )/3.0 );
                                chi0= - chi0;
                                hr = (float) cos( chi0 );
                                hi = (float) sin( chi0 );
                                prr = probe[ip].re(ix,iy);  // real
                                pri = probe[ip].im(ix,iy);  // imag
                                cpix.re(ix,iy) = prr*hr -pri*hi;
                                cpix.im(ix,iy) = prr*hi +pri*hr;
                                sum0 += prr*prr + pri*pri;
                            } else {
                                cpix.re(ix,iy) = 0.0F;
                                cpix.im(ix,iy) = 0.0F;
                            }
                        }  /*  end for( iy... )  */
                    }  /*  end for( ix... )  */
            
                    cpix.ifft();
                    
                    /* find normalization constant 
                      i.e. correct for constants in the FFT */
                    sum1 = 0.0;
                    for( ix=0; ix<nxprobe; ix++) {
                         for( iy=0; iy<nyprobe; iy++) {
                            prr = cpix.re(ix,iy);
                            pri = cpix.im(ix,iy);
                            sum1 += prr*prr + pri*pri;
                        }
                   }
            
                   /* integrate over real space detector and normalize */
                   for( ix=0; ix<nxprobe; ix++) {
                           rx2 = xoff[ip] - xp[ix];
                           rx2 = rx2*rx2;
                           for( iy=0; iy<nyprobe; iy++) {
                               r2 = yoff[ip] - yp[iy];
                               r2 = rx2 + r2*r2;
                               prr = cpix.re(ix,iy);
                               pri = cpix.im(ix,iy);
                               delta = prr*prr + pri*pri;
                               for( idetect=0; idetect<ndetect; idetect++) {
                                   if( CONFOCAL == collectorMode[idetect] ) {
                                     if( (r2 >= k2min[idetect] ) &&
                                               (r2 < k2max[idetect] ) )
                                         detect[it][idetect][ip] += delta*(sum0/sum1);
                                   }
                               }
                           }  /* end for( iy... )*/
                   }  /*  end for( ix....) */

                }  /* end if( doConfocal==TRUE) */
            
            }  /* end for( ip.. */

        }  /* end if( ((it...*/

        sbuffer = "xxx bad";
       nslice++;
       zslice += deltaz;
       istart += na;

    }  /* end while( istart...) */

    free( ixoff );
    free( iyoff );
    free( xoff );
    free( yoff );

    return;

}/* end autostem::STEMsignals() */


/*--------------------- trlayer() -----------------------*/
/*  same subroutine in autoslic.c and autostem.c

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  trans   = 2D array to get complex specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

  convert to cfpix class for trans 10-nov-2012 ejk

*/
void autostem::trlayer( const float x[], const float y[], const float occ[],
    const int Znum[], const int natom, 
    const float ax, const float by, const float kev,
    cfpix &trans, const long nx, const long ny,
    double *phirms, long *nbeams, const float k2max )
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;

    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rmin2 = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++)
           trans.re(ix,iy) = 0.0F;
    }
    
/*  paralleling this loop has little effect   */
/*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,r,ixw,iyw,vz,rsq)*/
#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq)
    for( i=0; i<natom; i++) {
        ixo = (int) ( x[i]/scalex );
        iyo = (int) ( y[i]/scaley );
        nx1 = ixo - idx;
        nx2 = ixo + idx;
        ny1 = iyo - idy;
        ny2 = iyo + idy;

    /* add proj. atomic potential at a local region near its center
       taking advantage of small range of atomic potential */

        for( ix=nx1; ix<=nx2; ix++) {
            rx2 = x[i] - ((double)ix)*scalex;
            rx2 = rx2 * rx2;
            ixw = ix;
            while( ixw < 0 ) ixw = ixw + nx;
            ixw = ixw % nx;
            for( iy=ny1; iy<=ny2; iy++) {
                rsq = y[i] - ((double)iy)*scaley;
                rsq = rx2 + rsq*rsq;
                if( rsq <= rmax2 ) {
                  iyw = iy;
                  while( iyw < 0 ) iyw = iyw + ny;
                  iyw = iyw % ny;
                  if( rsq < rmin2 ) rsq = rmin2;
                  /*r = sqrt( rx2 + r*r );
                  vz = occ[i] * vzatom( Znum[i], r ); slow */
                  vz = occ[i] * vzatomLUT( Znum[i], rsq );
                  trans.re(ixw,iyw) += (float) vz;
                }
            } /* end for(iy... */
       }  /* end for(ix... */

    } /* end for(i=0... */

    /* convert phase to a complex transmission function */
    sum = 0;
    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++) {
            vz = scale * trans.re(ix,iy);
            sum += vz;
            trans.re(ix,iy) = (float) cos( vz );
            trans.im(ix,iy) = (float) sin( vz );
        }
    }

    *phirms = sum / ( ((double)nx)*((double)ny) );

    /* bandwidth limit the transmission function */
    *nbeams = 0;
    trans.fft();
    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++) {
            k2 = ky2[iy] + kx2[ix];
            if (k2 < k2max) *nbeams += 1;
            else trans.re(ix,iy) = trans.im(ix,iy) = 0.0F;
        }
    }
    trans.ifft();
    
    return;

 };  /* end autostem::trlayer() */

/*------------------------- invert2D() ----------------------*/
/*
        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void autostem::invert2D( float** pix, long nx, long ny )
{
#define SWAP(a,b)       {t=a; a=b; b=t;}

        long ix, iy, ixmid, iymid;
        float t;

        ixmid = nx/2;
        iymid = ny/2;

        for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<iymid; iy++)
                SWAP( pix[ix][iy], pix[ix][iy+iymid] );

        for( ix=0; ix<ixmid; ix++) 
        for( iy=0; iy<ny; iy++)
                SWAP( pix[ix][iy], pix[ix+ixmid][iy] );

#undef SWAP
};  // end autostem::invert2D()
