/*
      *** autoslice.cpp ***

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

  started 24-july-1996 E. Kirkland
  working 19feb-1997 ejk
  added look-up-table vzatomLUT() for 3X-4X increase 
        in speed 23-may-1997 ejk
  put bandwith limit inside trlayer() 1-oct-1997 ejk
  added Gaussian thermal displacements 1-oct-1997 ejk
  removed /sqrt(3) from Thermal rms displacements 
    to be consistent with Int'l X-ray tables 22-dec-1997 ejk
  corrected zmin/max error with thermal displac. 24-dec-1997 ejk
  fixed small aliasing problem 5-jan-1998 ejk
  added unit cell replication option and moved ReadXYZcoord()
    into slicelib.c  11-jan-1998 ejk
  added astigmatism and modify to use different set of
    random offsets on each illum. angle with partial coherence
         5-feb-1998 ejk
  fix typo in z range message with partial coherence and
    thermal vibrations 9-july-1998 ejk
  update memory allocation routines 19-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  fixed bug in zmin/zmax calculation in coherent mode
     (move to after sortByZ() - it was before ) 8-jan-2002 ejk
  add cross section option (in non-partial coherence mode only)
        27-may-2005 ejk
  convet to faster sortByZ() 8-feb-2006 ejk
  move sortbyz() to slicelib.c 5-sep-2006 ejk
  add echo on y position in pixels for xz mode 4-may-2007 ejk
  update data type of nxl,nyl to be consistent with new tiffsubs
     17-jul-2007 ejk
  move xz depthpix save to be after transmit+propagate to get a
     full slice and proper anti-aliasing and also be consisten
     with what you get doing it by hand  and increase possible
     slices output (nz was off by one) 24-jan-2008 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better (consistent with autostem) 23-mar-2008 ejk
  take small things out of loop in trlayer() 14-may-2008 ejk
  parameterize vzatomLUT() vs r^2 instead of r to avoid a lot of sqrt()
      calls (a little faster)  6-jun-2008 ejk
  move vzatomLUT() to slicelib.c  11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  add Cs5 (and Cs->Cs3) 15-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4
     and convert to 4 char TAB size formatting 21-feb-2010 ejk
  add parallel computing of a few parts 21-feb-2010 ejk
  start conversion to faster FFTW 24-feb-2010 ejk
  move some things into slicelibW.c to share 6-mar-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add option to average over many frozen phonon
      configurations 3-aug-2010 ejk
  add multipole aberrations to probe 12-may-2011 ejk
  start conversion to floatTIFF.cpp and C++ 28-may-2012 ejk
  working 3-jun-2012 ejk
  convert to cfpix/fftw class from raw fftw 13-nov-2012 to 21-nov-2012 ejk
  move calculation into a class with separate command line front end
      29-may-2013 ejk
  fix typo in starting wave function 11-jun-2013 ejk
  fix bug to restore orginal df value, and add code to handle
     sigmf=0 problem 29-jun-2013 ejk
  fix minor format issue in nbeams message (%ld to %d) 19-jul-2013 ejk
  convert to string message 9-sep-2013 ejk
  change RNG seed argument to referenece so it get updated for 
      successive calls 21-sep-2013 ejk
  move toString() to slicelib from here 28-nov-2013 ejk
  add abbPhase2D() to calculate the 2D phase abb function 21-aug-2014 ejk 

  ax,by,cz  = unit cell size in x,y
  BW     = Antialiasing bandwidth limit factor
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  df0    = defocus (mean value)
  sgmaf = defocus spread (standard deviation)
  dfdelt = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 8 characters 
  
*/

#include "slicelib.hpp"    // misc. routines for multislice
#include "cfpix.hpp"       // complex image handler with FFT
#include "autoslic.hpp"    // header for this class

#include <sstream>	// string streams

//=============================================================
//---------------  creator and destructor --------------

autoslic::autoslic()
{
        BW= (2.0F/3.0F);  // bandwidth limit
        ABERR= 1.0e-4;    // max error for a,b

        NSMAX= 1000;  // max number of slices
        NCMAX= 256;   // max characters in file names
        NZMAX= 103;   // max atomic number Z

        twopi = 2.0 * (4.0 * atan( 1.0 ));

        // init control flags
        lcross = 0;
        lpartl = 0;
        lstart = 0;
        lwobble = 0;
        lbeams = 0;

        return;

};   //  end autoslic::autoslic()

autoslic::~autoslic()
{
}
//=============================================================
/*  abbPhase2D()

  calculate the phase of the abberation function in 2D
  in the objective aperture plane 
  (but out to the max angle allowed by sampling)
  - mainly just to look at

  added 21-aug-2014 ejk

  ab2D() = will get the 2D image of the abb. phase
  param[] = holds image parameters
  multiMode = flag, if not 0 then include all multipole abberations

*/
void autoslic::abbPhase2D( cfpix &ab2D, float param[], int multiMode )
{
    int ix, iy, ixmid, iymid, nx, ny;
    float k2, k2max, v0, wavlen, ax, by, pi, t;
    float *kx, *ky, *xpos, *ypos, *kx2, *ky2;
    double chi0, alx, aly;

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV

    wavlen = (float) wavelength( v0 );

    pi = (float) (4.0 * atan( 1.0 ));

    //----- calculate spatial frequencies and positions for future use 

    kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
    kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
    xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
    freqn( kx, kx2, xpos, nx, ax );

    ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
    ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
    ypos = (float*) malloc1D( ny, sizeof(float), "ypos" );
    freqn( ky, ky2, ypos, ny, by );

    //  rearrange frequencies before calculation so we don't
    //  have to rearrange the whole image
    ixmid = nx/2;
    iymid = ny/2;
    for( ix=0; ix<ixmid; ix++) {
        t=kx[ix];  kx[ix] = kx[ix+ixmid];  kx[ix+ixmid] =t;
        t=kx2[ix]; kx2[ix]= kx2[ix+ixmid]; kx2[ix+ixmid]=t;
    }
    for( iy=0; iy<iymid; iy++) {
        t=ky[iy];  ky[iy] = ky[iy+iymid];  ky[iy+iymid] =t;
        t=ky2[iy]; ky2[iy]= ky2[iy+iymid]; ky2[iy+iymid]=t;
    }

    //-----  make array right size if needed
    ab2D.resize( nx, ny );
    ab2D = 0.0F;

    //-----  calculate max sampling angles
    k2max = nx/(2.0F*ax);
    t = ny/(2.0F*by);
    if( t < k2max ) k2max = t;
    k2max = BW * k2max;
    k2max = k2max*k2max;
    
    //------- calcualte the phase and file the array
    for( ix=0; ix<nx; ix++) {
        alx = wavlen * kx[ix];  // x component of angle alpha
        for( iy=0; iy<ny; iy++) {
            aly = wavlen * ky[iy];  // y component of angle alpha
            k2 = kx2[ix] + ky2[iy];
            if( k2 <= k2max ) {
                chi0 = (2.0*pi/wavlen) * chi( param, 
                               alx, aly, multiMode );
                //  make phase modulo 2pi - should be a better way to do this (?)
                //   remember that % operator only works with int
                while( chi0 < -pi) chi0 += twopi;
                while( chi0 >  pi) chi0 -= twopi;

                ab2D.re(ix, iy ) = (float) chi0;
            }
        } // end for( iy=...
    }  // end for(ix=....

    //----------- end:  free scratch arrays and exit --------------------

    free( kx );
    free( kx2 );
    free( xpos );

    free( ky );
    free( ky2 );
    free( ypos );

    return;
 
}  // end abbPhase2D()

//=============================================================
/*  calculate()

  input:
        wave0 = complex image with starting image
                (ignored if lstart = 0 or partial coherence calculated)
        param[] = image parameters (most will not be changed here)
        multimode = flag controlling multipole aberrations
        natom = number of atoms
        x[],y[],z[] = atomic coord
        Znum[] atomic number of each atom
        occ[] = occupancy of each atomic site
        wobble[] = thermal oscillation amplitude 
        hb[],kb[] = (h,k) indexes of beams to monitor during
                          propagation
                (ignore if lbeams = 0)
        nbout  = number of beams to record
        ycross = y position to save xz cross section

  mode flags:  lbeams, lcross, lpartl, lstart, lwobble

  output:
        pix = complex image to get results, orig. data lost and may be resized
                (complex for coherent calc. and real for partial coherence)
        beams = track specified beams as wave propagates
                (ignored if lbeams = 0)
*       depthpix = will be resized and get xz cross section image
                save intensity in real part (imag part not used)
                (ignored if lcross = 0)
*/
void autoslic::calculate(cfpix &pix, cfpix &wave0, cfpix &depthpix,
        float param[], int multiMode, int natom, unsigned long *iseed,
        int Znum[], float x[], float y[], float z[], float occ[], float wobble[],
        cfpix &beams, int hb[], int kb[], int nbout, float ycross, float dfdelt )
{
    int i, ix, iy, iz, ixmid, iymid, nx, ny, nz, iycross, istart, nwobble, nbeams,
        nacx,nacy, iqx, iqy, iwobble, ndf, idf, ib, na, islice, nzout, nzbeams,
        n1, n2;
    int *Znum2, *hbeam, *kbeam;

    float wmin, wmax, xmin,xmax, ymin, ymax, zmin, zmax;
    float k2, k2max, scale, v0, mm0, wavlen, rx, ry, rx2,ry2,
        ax, by, pi, rmin, rmax, aimin, aimax,
        ctiltx, ctilty, tctx, tcty, acmin, acmax, df, df0, sigmaf,
        aobj, qx, qy, qy2, q2, q2min, q2max, sumdf, pdf, k2maxo,
        temperature;
    float tr, ti, wr, wi;

    float *kx, *ky, *xpos, *ypos, *kx2, *ky2;
    float *x2, *y2, *z2, *occ2;
    float *propxr, *propxi, *propyr, *propyi;

    double sum, xdf, chi0, t, zslice, deltaz, phirms, rsq, vz, alx, aly;

    cfpix wave;            // complex probe wave functions
    cfpix trans;           // complex transmission functions
    cfpix temp ;           // complex scratch wavefunction

    // ---- get setup parameters from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );
    v0 = param[pENERGY];                // electron beam energy in keV
    df0 = param[pDEFOCUS];              // defocus
    sigmaf = param[pDDF];
    ctiltx = param[ pXCTILT ];          // crystal tilt
    ctilty = param[ pYCTILT ];
    acmax = param[pCAPERT];             // condencer angles
    acmin = param[pCAPERTMIN];
    aobj = param[ pOAPERT ];            // objective aperture
    temperature = param[ pTEMPER ];     // temperature
    nwobble = ToInt( param[ pNWOBBLE ] );       //  number config. to average
    deltaz = param[ pDELTAZ ];          // slice thickness

    if( nwobble < 1 ) nwobble = 1;      //  has to be at least one configuration

    pi = (float) (4.0 * atan( 1.0 ));

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        sbuffer="bad size parameters in autoslic::calculate()";
        messageAS( sbuffer );
        exit( 0 );
    }
            
    /*  calculate relativistic factor and electron wavelength */
    mm0 = 1.0F + v0/511.0F;
    wavlen = (float) wavelength( v0 );
    //printf("electron wavelength = %g Angstroms\n", wavlen);

    /*  calculate the total specimen volume and echo */
    xmin = xmax = x[0];
    ymin = ymax = y[0];
    zmin = zmax = z[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( x[i] < xmin ) xmin = x[i];
        if( x[i] > xmax ) xmax = x[i];
        if( y[i] < ymin ) ymin = y[i];
        if( y[i] > ymax ) ymax = y[i];
        if( z[i] < zmin ) zmin = z[i];
        if( z[i] > zmax ) zmax = z[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    // --- leave this in main calling program
    //sprintf(stemp, "Total specimen range is\n %g to %g in x\n"
    //       " %g to %g in y\n %g to %g in z",
    //         xmin, xmax, ymin, ymax, zmin, zmax );
    //messageAS( stemp );

    // --- leave this in main calling program
    //if( lwobble == 1 ) {
    //    sprintf(stemp, "Range of thermal rms displacements (300K) = %g to %g\n",
    //        wmin, wmax );
    //    messageAS( stemp );
    //}
    
#ifdef USE_OPENMP
    /*  force LUT init. to avoid redundant init in parallel form */ 
    rsq = 0.5;  /* arbitrary position */   
    for( i=0; i<natom; i++) vz =  vzatomLUT( Znum[i], rsq );
#endif

/*  calculate spatial frequencies and positions for future use */

    rx = 1.0F/ax;
    rx2= rx*rx;
    ry = 1.0F/by;
    ry2= ry*ry;
    ixmid = nx/2;
    iymid = ny/2;

    kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
    kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
    xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
    freqn( kx, kx2, xpos, nx, ax );

    ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
    ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
    ypos = (float*) malloc1D( ny, sizeof(float), "ypos" );
    freqn( ky, ky2, ypos, ny, by );

    /*---- allocate some more arrays and initialize wavefunction ----*/

    if( (nbout > 0) && (lbeams ==1) ) {
        hbeam = (int*) malloc1D( nbout, sizeof(int), "hbeam" );
        kbeam = (int*) malloc1D( nbout, sizeof(int), "kbeam" );
    }

    trans.resize( nx, ny );
    wave.resize( nx, ny );

    if( (lstart == 0) || (nx!= wave0.nx()) || (ny!=wave0.ny()) ) {
             wave = 1.0F;  
    } else { wave = wave0; }

    trans.init();
    wave.copyInit( trans );   //  must be after "wave = wave0"

    if( lcross == 1 ) {
        /* nz may be too small with thermal vibrations so add a few extra */
        nz = (int) ( (zmax-zmin)/ deltaz + 3.5);
        depthpix.resize( nx, nz );
        for( ix=0; ix<nx; ix++)
        for( iz=0; iz<nz; iz++)  depthpix.re(ix,iz) = depthpix.im(ix,iz) = 0.0F;
        iycross = (int) ( 0.5 + (ny * ycross / by));
        while( iycross < 0 ) iycross += ny;
        iycross = iycross%ny;  /* make periodic in ny */
        sbuffer= "save xz cross section at iy= "+toString(iycross)+" pixels";
        messageAS( sbuffer );
    }

 /*  calculate propagator function  */
 
    k2max = nx/(2.0F*ax);
    tctx = ny/(2.0F*by);
    if( tctx < k2max ) k2max = tctx;
    k2max = BW * k2max;
    sbuffer= "Bandwidth limited to a real space resolution of "+toString(1.0F/k2max)
            +" Angstroms";
    messageAS( sbuffer );
    sbuffer= "   (= " + toString(wavlen*k2max*1000.0F)
            + " mrad)  for symmetrical anti-aliasing.";
    messageAS( sbuffer );
    k2max = k2max*k2max;

    tctx = (float) (2.0 * tan(ctiltx));
    tcty = (float) (2.0 * tan(ctilty));
    
    propxr = (float*) malloc1D( nx, sizeof(float), "propxr" );
    propxi = (float*) malloc1D( nx, sizeof(float), "propxi" );
    propyr = (float*) malloc1D( ny, sizeof(float), "propyr" );
    propyi = (float*) malloc1D( ny, sizeof(float), "propyi" );

    scale = pi * ((float)deltaz);

    for( ix=0; ix<nx; ix++) {
        t = scale * ( kx2[ix]*wavlen - kx[ix]*tctx );
        propxr[ix] = (float)  cos(t);
        propxi[ix] = (float) -sin(t);
    }
    for( iy=0; iy<ny; iy++) {
        t = scale * ( ky2[iy]*wavlen - ky[iy]*tcty );
        propyr[iy] = (float)  cos(t);
        propyi[iy] = (float) -sin(t);
    }
 
/*  iterate the multislice algorithm proper

   NOTE: zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

  partial coherence method
   force the integrals to include the origin and to be symmetric
   about the origin and to have the same periodic boundary
   conditions as the sampling grid
*/
    if( lpartl == 1 ) {
        
        //sprintf(stemp,"Illumination angle sampling (in mrad) = %f, %f\n",
        //    1000.*rx*wavlen, 1000.*ry*wavlen);
        //messageAS( stemp );

        pix.resize(nx,ny);
        pix = 0.0F;     // start with zero and sum into this pix
        
        temp.resize( nx, ny );
        temp.copyInit( trans );

        if( fabs( (double) dfdelt ) < 1.0 ) ndf = 1; 
        else ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );

        nacx = (int) ( ( acmax / ( wavlen * rx ) ) + 1.5F );
        nacy = (int) ( ( acmax / ( wavlen * ry ) ) + 1.5F );

        q2max = acmax / wavlen;
        q2max = q2max*q2max;

        q2min = acmin / wavlen;
        q2min = q2min*q2min;

        k2maxo = aobj / wavlen;
        k2maxo = k2maxo*k2maxo;

        nillum = 0;

        /* for Monte Carlo stuff */
        x2      = (float*) malloc1D( natom, sizeof(float), "x2" );
        y2      = (float*) malloc1D( natom, sizeof(float), "y2" );
        z2      = (float*) malloc1D( natom, sizeof(float), "z2" );
        occ2    = (float*) malloc1D( natom, sizeof(float), "occ2" );
        Znum2   =   (int*) malloc1D( natom, sizeof(int), "Znum2" );

        if( lwobble == 0 ) sortByZ( x, y, z, occ, Znum, natom );

        /*  integrate over the illumination angles */

        for( iwobble=0; iwobble<nwobble; iwobble++) {
            if( lwobble == 1 ) {
                sbuffer= "configuration # " + toString( iwobble+1 );
                messageAS( sbuffer );
            }
            for( iqy= -nacy; iqy<=nacy; iqy++) {
                qy = iqy * ry;
                qy2 = qy * qy;
        
                for( iqx= -nacx; iqx<=nacx; iqx++) {
                    qx = iqx * rx;
                    q2 = qx*qx + qy2;
        
                    if( (q2 <= q2max) && (q2 >= q2min) ) {
                        nillum += 1;
                        for( ix=0; ix<nx; ix++) {
                            for( iy=0; iy<ny; iy++) {
                                t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
                                wave.re(ix,iy) = (float) cos(t);  /* real */
                                wave.im(ix,iy) = (float) sin(t);  /* imag */
                            }
                        }
                        /*  add random thermal displacements scaled by temperature
                                if requested 
                            remember that initial wobble is at 300K for each direction */
                        if( lwobble == 1 ){
                            scale = (float) sqrt(temperature/300.0) ;
                            for( i=0; i<natom; i++) {
                                x2[i] = x[i] + (float)(wobble[i]*rangauss(iseed)*scale);
                                y2[i] = y[i] + (float)(wobble[i]*rangauss(iseed)*scale);
                                z2[i] = z[i] + (float)(wobble[i]*rangauss(iseed)*scale);
                                occ2[i] = occ[i];
                                Znum2[i] = Znum[i];
                            }
                            sbuffer= "Sorting atoms by depth...";
                            messageAS( sbuffer );
                            sortByZ( x2, y2, z2, occ2, Znum2, natom );
                            zmin = z2[0];       /* reset zmin/max after wobble */
                            zmax = z2[natom-1];
                            sbuffer= "Thickness range with thermal displacements is "
                                   + toString(zmin) + " to " + toString(zmax)+" (in z)";
                            messageAS( sbuffer );
                        } else for( i=0; i<natom; i++) {
                            x2[i] = x[i];
                            y2[i] = y[i];
                            z2[i] = z[i];
                            occ2[i] = occ[i];
                            Znum2[i] = Znum[i];
                        }
            
                        zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
                        istart = 0;
            
                        while( istart < natom ) {
            
                            /* find range of atoms for current slice */
                            na = 0;
                            for(i=istart; i<natom; i++) 
                            if( z2[i] < zslice ) na++; else break;
            
                            /* calculate transmission function, skip if layer empty */
                            if( na > 0 ) {
                                trlayer( &x2[istart], &y2[istart], &occ2[istart],
                                    &Znum2[istart],na, ax, by, v0, 
                                    trans,  nx, ny, kx2, ky2, &phirms, &nbeams, k2max );
                
                                wave *= trans;  // transmit
                            }
            
                            /* remember: prop needed here to get anti-aliasing
                                    right */
                            wave.fft();
                            propagate( wave, propxr, propxi, propyr, propyi,
                                kx2,  ky2,  k2max, nx, ny );
                            wave.ifft();
            
                            zslice += deltaz;
                            istart += na;
            
                        } /* end while(zslice<=..) */
          
                        scale = 1.0F / ( ((float)nx) * ((float)ny) );
                        sum = 0.0;
                        for( ix=0; ix<nx; ix++) {
                            for( iy=0; iy<ny; iy++)
                                sum += wave.re(ix,iy)*wave.re(ix,iy)
                                    + wave.im(ix,iy)*wave.im(ix,iy);
                        }
                        sum = sum * scale;
             
                        sbuffer=  "Illum. angle = " + toString(1000.*qx*wavlen) +
                                ", "+toString(1000.*qy*wavlen) +
                                " mrad, integ. intensity= "+toString(sum);
                        messageAS( sbuffer );
            
                        /*-------- integrate over +/- 2.5 sigma of defocus ------------ */
                        //   should convert to Gauss-Hermite quadrature sometime
                        wave.fft();
                        sumdf = 0.0F;

                        if( fabs( (double) sigmaf ) < 1.0 ) n1 = n2 = 0; 
                        else {
                                n1 = -ndf;
                                n2 = ndf;
                        }

                        for( idf= n1; idf<=n2; idf++) {
                            param[pDEFOCUS] = df = df0 + idf*dfdelt;
            
                            for( ix=0; ix<nx; ix++) {
                                alx = wavlen * kx[ix];  /* x component of angle alpha */
                                for( iy=0; iy<ny; iy++) {
                                    aly = wavlen * ky[iy];  /* y component of angle alpha */
                                    k2 = kx2[ix] + ky2[iy];
                                    if( k2 <= k2maxo ) {
                                        chi0 = (2.0*pi/wavlen) * chi( param, 
                                                alx, aly, multiMode );
                                        tr = (float)  cos(chi0);
                                        ti = (float) -sin(chi0);
                                        wr = wave.re(ix,iy);
                                        wi = wave.im(ix,iy);
                                        temp.re(ix,iy) = wr*tr - wi*ti;
                                        temp.im(ix,iy) = wr*ti + wi*tr;
                                    } else {
                                        temp.re(ix,iy) = 0.0F;  /* real */
                                        temp.im(ix,iy) = 0.0F;  /* imag */
                                    }
                                }  /*  end for( iy=0... ) */
                            }   /*  end for( ix=0... ) */

                            temp.ifft();
            
                            if( (0==n1) && (0==n2) ) pdf = 1;
                            else {
                                xdf = (double) ( (df - df0) /sigmaf );
                                pdf = (float) exp( -0.5F * xdf*xdf );
                            }
                            sumdf += pdf;
            
                            for( ix=0; ix<nx; ix++) {
                                for( iy=0; iy<ny; iy++) {
                                    wr = temp.re(ix,iy);
                                    wi = temp.im(ix,iy);
                                    pix.re(ix,iy) += pdf* ( wr*wr + wi*wi );
                                }
                            }
            
                        }/* end for(idf..) */

                        param[ pDEFOCUS ] = df0;  // return to original value

                    }/* end if( q2...) */
        
                } /* end for( iqx..) */
            } /* end for( iqy..) */
        } /* end for( iwobble...) */

        //----  put these in main calling program if neede
        //sprintf(stemp, "Total number of illumination angle = %ld",
        //        nillum);
        //message ( stemp );
        //sprintf(stemp, "Total number of defocus values = %d", 2*ndf+1);
        //messageAS( stemp );

        /*  remember that nillum already includes nwobble so don't
             divide by nwobble! */
        scale = 1.0F / ( ((float)nillum) * sumdf );
        rmin  = pix.re(0,0) * scale;
        rmax  = rmin;
        aimin = 0.0F;
        aimax = 0.0F;

        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pix.re(ix,iy) *= scale;
            if( pix.re(ix,iy) < rmin ) rmin = pix.re(ix,iy);
            if( pix.re(ix,iy) > rmax ) rmax = pix.re(ix,iy);
        }

/* ---- start coherent method below ----------------
        (remember that waver,i[][] was initialize above) */

    } else {

        if( lbeams ==1 ) {
            nzbeams = (int) ( (zmax-zmin)/ deltaz + 3.5);
            beams.resize(nbout, nzbeams );   // to save values
            for(ib=0; ib<nbout; ib++) {
                    hbeam[ib] = hb[ib];
                    kbeam[ib] = kb[ib];
            }

            //  make them all positive just in case
            for( ib=0; ib<nbout; ib++) {
                if( hbeam[ib] < 0 ) hbeam[ib] = nx + hbeam[ib];
                if( kbeam[ib] < 0 ) kbeam[ib] = ny + kbeam[ib];
                if( hbeam[ib] < 0 ) hbeam[ib] = 0;
                if( kbeam[ib] < 0 ) kbeam[ib] = 0;
                if( hbeam[ib] > nx-1 ) hbeam[ib] = nx-1;
                if( kbeam[ib] > ny-1 ) kbeam[ib] = ny-1;
            }
        }

        /*  add random thermal displacements scaled by temperature if requested 
            remember that initial wobble is at 300K for each direction */
        if( lwobble == 1 ){
            scale = (float) sqrt(temperature/300.0) ;
            for( i=0; i<natom; i++) {
                x[i] += (float) (wobble[i] * rangauss( iseed ) * scale);
                y[i] += (float) (wobble[i] * rangauss( iseed ) * scale);
                z[i] += (float) (wobble[i] * rangauss( iseed ) * scale);
            }
        }

        sbuffer= "Sorting atoms by depth...";
        messageAS( sbuffer );
        sortByZ( x, y, z, occ, Znum, natom );

        if( lwobble == 1 ){
            zmin = z[0];        /* reset zmin/max after wobble */
            zmax = z[natom-1];
            sbuffer="Thickness range with thermal displacements"
                " is "+toString(zmin)+" to "+toString(zmax)+" (in z)";
            messageAS( sbuffer );
        }

        scale = 1.0F / ( ((float)nx) * ((float)ny) );

        zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
        istart = 0;
        islice = 1;

        while( (istart < natom) && ( zslice < (zmax+deltaz) ) ) {

            /* find range of atoms for current slice */
            na = 0;
            for(i=istart; i<natom; i++) 
            if( z[i] < zslice ) na++; else break;

            /* calculate transmission function, skip if layer empty */
            if( na > 0 ) {
                trlayer( &x[istart], &y[istart], &occ[istart],
                    &Znum[istart], na, ax, by, v0, trans,
                    nx, ny, kx2, ky2, &phirms, &nbeams, k2max );
    
                /*??? printf("average atompot comparison = %g\n", 
                           phirms/(wavlen*mm0) ); */
    
                wave *= trans;    //  transmit
            }

            /*  bandwidth limit */
            wave.fft();

            if( (lbeams== 1) && (islice<nzbeams) && (islice>0) )  {
                for( ib=0; ib<nbout; ib++) {
                    beams.re(ib,islice-1) = scale*wave.re(hbeam[ib],kbeam[ib] );   // real
                    beams.im(ib,islice-1) = scale*wave.im(hbeam[ib],kbeam[ib] );   // imag
                }
            }

            /* remember: prop needed here to get anti-aliasing right */
            propagate( wave, propxr, propxi,
                propyr, propyi, kx2,  ky2,  k2max, nx, ny );
            wave.ifft();

            /* save depth cross section if requested */
            if( (lcross == 1) && (islice<=nz) ) {
                for( ix=0; ix<nx; ix++) {
                    depthpix.re(ix, islice-1) = 
                        wave.re(ix,iycross)*wave.re(ix,iycross)
                           + wave.im(ix,iycross)*wave.im(ix,iycross);
                }
                nzout = islice;
            }

            sum = 0.0;
            for( ix=0; ix<nx; ix++) {
                for( iy=0; iy<ny; iy++)
                    sum += wave.re(ix,iy)*wave.re(ix,iy) +
                        wave.im(ix,iy)*wave.im(ix,iy);
            }
            sum = sum * scale;

            sbuffer= "z= " + toString(zslice)+" A, " + toString(nbeams) + " beams, "
                    + toString(na)+" coord., \n"
                    + "     aver. phase= "+toString(phirms)
                    +", total intensity = "+toString(sum) ;
            messageAS( sbuffer );

            zslice += deltaz;
            istart += na;
            islice++;

        } /* end while(istart<natom..) */

        pix.resize(nx,ny);
        pix = wave;

    } /* end else .. coherent section */

    //----------- end:  free scratch arrays and exit --------------------

    free( kx );
    free( kx2 );
    free( xpos );

    free( ky );
    free( ky2 );
    free( ypos );

    free( propxr );
    free( propxi );
    free( propyr );
    free( propyi );

    if( (nbout > 0) && (lbeams ==1) ) {
        free( hbeam );
        free( kbeam );
    }

    if( lpartl == 1 ) {
        free( x2 );
        free( y2 );
        free( z2 );
        free( occ2 );
        free( Znum2 );
    }

    return;

};   //  end autoslic::calculate()

//=============================================================
/* -------------------  messageAS() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void autoslic::messageAS( std::string &smsg,  int level )
{
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

}  // end autoslic::messageAS()

//=============================================================
/*--------------------- trlayer() -----------------------*/
/*   same subroutine in autoslic.c and autostem.c

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  transr  = 2D array to get real part of specimen
        transmission function
  transi  = 2D array to get imag. part of specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

*/
void autoslic::trlayer(  const float x[], const float y[], const float occ[],
            const int Znum[], const int natom, const float ax, const float by,
            const float kev, cfpix &trans, const int nx, const int ny,
            const float kx2[], const float ky2[],
            double *phirms, int *nbeams, const float k2max  )
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rminsq, sum, scale, scalex, scaley;

    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rminsq = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) {
        for( iy=0; iy<ny; iy++)
           trans.re(ix,iy) = 0.0F;    /* real part trans[iy + ix*ny][0] */
    }

/*  run this in parallel  */
/*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq) */
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
                if( rsq < rminsq ) rsq = rminsq;
                /* r = sqrt( r );
                vz = occ[i] * scale * vzatom( Znum[i], r ); slow */
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

 }  /* end autoslic::trlayer() */
 
