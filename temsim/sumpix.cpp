/*  ***  sumpix.cpp ***

------------------------------------------------------------------------
Copyright 1998-2014 Earl J. Kirkland

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

  read in two or more images and average
  optionally take the FFT and sum square magnitudes
  (for frozen phonon averaging of CBED patterns)

  started 7-mar-1997 E. Kirkland
  added log option 8-mar-1997 ejk
  added real space averaging and changed name from sumCBED to sumpix
        6-jan-1998 ejk
  in working form 16-jan-1998 ejk
  fixed small problem with complex images in real space 19-jan-1998 ejk
  changed greyscale scaling on log() option 25-feb-1998 ejk
  update memory allocation routines 20-nov-1999 ejk
  change void main() to int main() for better portability
         22-jan-2000 ejk
  small cosmetic changes 19-jul-2007 ejk
  convert to GPL 3-jul-2008 ejk
  get return value of scanf() to remove warnings from gcc 4.4
     and convert to 4 char TAB size formatting 23-may-2010 ejk
  convert to FFTW in a very crude way 26-may-2010 ejk
  convert to floatTIFF.cpp and C++, remove summation option for
     plain TIFF images (no longer supported in floatTIFF) 1-apr-2012 ejk
  convert to cfpix/fftw class from raw fftw 6-nov-2012 to 30-oct-2012 ejk
  convert DateTime in floatTIFF to string 12-feb-2014 ejk
  convert to streams and strings 6-may-2014 ejk
*/

#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include "cfpix.hpp"      // complex image handler with FFT
#include "slicelib.hpp"   // misc. routines for multislice 
#include "floatTIFF.hpp"  // file I/O routines in TIFF format 

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output

using namespace std;

const int integerPIX=  0;   // pix type flags 
const int floatPIX=    1;

/*   subroutines at end of file */
void fft2d( float **pixr, float **pixi, const long nx,
                   const long ny, int inverse );
void invert2D( float** pix, long nx, long ny );


int main()
{
    string datetime, fileout;

    int i, ipix, ix, iy, nx, ny, nxold, nyold, ixmid, iymid, npix, npixold,
        ninput, nsum, nh, logpix, PowerSpectra, pixtype, NPARAM;
    long *nhist;

    float scale, pixc, rmin,rmin2,rmax, aimin,aimax,tr, ti, dx, dy;
    float *param;
    float  **pixr, **pixi, **pixout;
    double sum, *hist, ax, by, rx, ry2;

    ofstream fp;

    floatTIFF myFile;

    /*--------  get input file names etc. ------------ */
    cout << "sumpix version dated 6-may-2014 ejk" << endl;
    cout << "Copyright (C) 1998-2014 Earl J. Kirkland"  << endl;
    cout << "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
         << " under the GNU general public license\n"  << endl;

    cout << "Sum multiple image or wave function files,\n"
         << "complex images will be converted to squared "
         << "magnitude before summing." << endl;
    cout << "All input images must be the same type and size.\n" << endl;
    cout << "Type number of input image files" << endl;
    cin >> ninput;
    string *filein = new string[ ninput ];
    for( ipix=0; ipix<ninput; ipix++) {
        cout << "input " << ipix << " : " ;
        cin >> filein[ipix];
    }
    cout << endl;

    cout << "Type name of output file:" << endl;
    cin >> fileout;

    logpix = askYN( "Do you want to display on log scale");

    PowerSpectra = askYN( "Do you want to convert to a power spectra");

/* get image size and type from the first input pix
    all successive images have to be the same type and size !!! 
   -remember that floatTIFF cannot handle plain integer TIFF images

 -------- read floating point images and average --------

   remember that complex images are stacked side by side
    with npix=2 and nx twice its real value 
    (real images have npix=1 and nx its normal value)
*/
    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "param" );

    for( ipix=0; ipix<ninput; ipix++) {

        for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;
        if( myFile.read( filein[ipix].c_str() ) != 1 ) {
            cout << "Cannot open file " << filein[ipix] << endl;
            exit( 0 );
        }
        datetime = myFile.getDateTime( );
        nx = (int) myFile.nx();
        ny = (int) myFile.ny();
        npix = myFile.getnpix();
        if( 0 == ipix ) {
            npixold = npix;
            nxold = nx;
            nyold = ny;
            pixr = (float**) malloc2D( 2*nx, ny, sizeof(int), "pixr-1" );  //npix ????
            pixout = (float**) malloc2D( nx, ny, sizeof(int), "pixout-1" );
            for( ix=0; ix<nx; ix++) 
                for( iy=0; iy<ny; iy++) pixout[ix][iy] = 0.0F;
            pixtype = floatPIX;
            cout << "Image size : Nx= " << nx << ", Ny= " << ny  << endl;
        } else if( (nx != nxold) || (ny != nyold) ) {
            cout << "different size in file " << filein[ipix] << ", nx= "
                 << nx << ", ny= " << ny << endl;
            exit( 0 );
        }
        if( npix != npixold ) {
            cout << "Can't mix real and complex images"
                " in file: " << filein[ipix] << endl;
            exit( 0 );
        }
        if( (npix<1) || (npix>2) ) {
            cout << "bad npix = " << npix << " in TIFF file "
                << filein[ipix] << endl;
            exit( 0 );
        }

        //  copy both real+imag back to old style array to re-use old code
        //        (not optimal but works for now)
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
                pixr[ix][iy] = myFile(ix,iy);

        nx = nx /npix;
        ax = myFile.getParam(pDX) * ((float)nx);
        by = myFile.getParam(pDY) * ((float)ny);
        rmin = myFile.getParam(pRMIN);
        rmax = myFile.getParam(pRMAX);
        aimin = myFile.getParam(pIMIN);
        aimax = myFile.getParam(pIMAX);
        if( npix == 2 ) {
            cout << "pix "<< ipix <<" created " << datetime << 
                ", range: "<< rmin << " to " << rmax << " (real),\n"
                "     and " <<aimin << " to " << aimax << " (imag)" << endl;
        } else if( npix == 1 ) {
            cout << "pix " << ipix << " created " << datetime <<
                ", range: " << rmin << " to " << rmax << " (real) " << endl;
        }
        if( PowerSpectra == 1 ) {
            if( npix == 1 ) {
                if( 0 == ipix )
                    pixi = (float**) malloc2D( nx, ny, sizeof(float),
                        "pixi-2" );
                for( ix=0; ix<nx; ix++) 
                    for( iy=0; iy<ny; iy++)   pixi[ix][iy] = 0.0F;
            } else if( (npix==2) && (ipix==0) ) pixi = pixr + nx;
            npix = 2;
            fft2d ( pixr, pixi, nx, ny, +1);
        }

        if( npix == 1 ) {       /* real pix */
            for( ix=0; ix<nx; ix++) 
            for( iy=0; iy<ny; iy++) 
                pixout[ix][iy] += pixr[ix][iy];
        } else if( npix == 2 ) {    /* complex pix */
            if( 0 == ipix ) pixi = pixr + nx;
            for( ix=0; ix<nx; ix++) 
            for( iy=0; iy<ny; iy++) {
                tr = pixr[ix][iy];
                ti = pixi[ix][iy];
                pixout[ix][iy] += ( tr*tr + ti*ti);
            }
        }

    }  // end for(ipix=... )
    

/*  Output results and find min and max to echo
     NOTE the logarithmic scaling of diffraction pattern
    is taken from Gonzalez and Wintz pg 48
    added scaling trick from showpix.f  9-aug-1995 ejk
*/
    cout << "Output pix size : Nx= " << nx << ", Ny= " << ny << endl;

    if( (PowerSpectra == 1) && ( pixtype == floatPIX ) ) {

        /* put (0,0) in the center */
        invert2D( pixout, nx, ny);

        /* histogram the azimutal average */
        hist = (double*) malloc1D( (nx+ny), sizeof(double), "hist" );
        nhist = (long*) malloc1D( (nx+ny), sizeof(long), "nhist" );
        for( ix=0; ix<(nx+ny); ix++) {
            hist[ix] = 0.0;
            nhist[ix] = 0;
        }

        scale = 1.0F / ( ((float)nx) * ((float)ny) );

        sum = 0.0;
        nsum = 0;
        nh = 0;
        ixmid = nx/2;
        iymid = ny/2;

        for( iy=0; iy<ny; iy++) {
            ry2 = (double) ( iy-iymid);
            ry2 = ry2*(ax/by);
            ry2 = ry2*ry2;
            for( ix=0; ix<nx; ix++) {
                pixc = pixout[ix][iy];
                rx = (double) (ix-ixmid);
                i = (int) ( sqrt( rx*rx + ry2 ) + 0.5);
                hist[i] += pixc;
                nhist[i]++;
                if( i > nh ) nh = i;
                if( logpix == 1 ) {
                    if( pixc > 1.e-10F)  pixc = 
                        (float) log( (double) fabs(pixc) );
                    else pixc = -23.0F;
                    pixout[ix][iy] = pixc;
                }
                if( (ix == 0) && (iy == 0) ) {
                    rmin = pixc;
                    rmax = rmin;
                } else if( (ix != ixmid) && (iy != iymid) ) {
                    if( pixc < rmin ) rmin = pixc;
                    if( pixc > rmax ) rmax = pixc;
                }
                if( (ix>(3*nx)/8) && (ix<(5*nx)/8) &&
                    (iy>(3*ny)/8) && (iy<(5*ny)/8) ) {
                    sum = sum + pixc;
                    nsum += 1;
                }

            }  /* end for ix... */
        } /* end for iy... */

        cout << "write azimuthal averaged intensity vs. \n"
            "  spatial frequency k, into file azimuth.dat" << endl;
        fp.open( "azimuth.dat" );
        if( fp.bad() ) {
            cout << "cannot open file azimuthal.dat" << endl;
            exit( 0 );
        }
        for( i=0; i<=nh; i++) {
            hist[i] = hist[i] / nhist[i];
            fp << setw(16) << ((double)i)/ax << setw(16) << hist[i] << endl;
        }
        fp.close( );

        myFile.resize( nx, ny );  // in case it was complex
        myFile.setnpix( 1 );
        myFile.setParam( pRMAX, rmax);
        myFile.setParam( pIMAX, aimax = 0.0F);
        myFile.setParam( pIMIN, aimin = 0.0F);
        myFile.setParam( pRMIN, rmin );
        myFile.setParam( pDX,  dx = 1.0F / ((float)ax) );
        myFile.setParam( pDY,  dy = 1.0F / ((float)by) );
        cout << "output image size: " << nx*dx << " to " <<  ny*dy << " /Angstroms" << endl;
        cout << "Power Spectra range " << rmin << " to " << rmax << endl;

        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
            myFile(ix,iy) = pixout[ix][iy];
        rmin2= (float) (0.05*rmin + 0.95*sum/nsum);   //  somtimes better
        myFile.write( fileout.c_str(), rmin2, rmax, aimin, aimax, dx, dy );

    } else if(  (pixtype == floatPIX) && (PowerSpectra == 0) ) {

        for( iy=0; iy<ny; iy++) {
            for( ix=0; ix<nx; ix++) {
                pixc = pixout[ix][iy];
                if( logpix == 1 ) {
                    if( pixc > 1.e-30F)  pixc = 
                        (float) log( (double) fabs(pixc) );
                    else pixc = -100.0F;
                    pixout[ix][iy] = pixc;
                }
                if( (ix == 0) && (iy == 0) ) {
                    rmin = pixc;
                    rmax = rmin;
                } else {
                    if( pixc < rmin ) rmin = pixc;
                    if( pixc > rmax ) rmax = pixc;
                }

            }  /* end for ix... */
        } /* end for iy... */

        myFile.resize( nx, ny );  // in case it was complex
        myFile.setnpix( 1 );
        myFile.setParam( pRMAX, rmax);
        myFile.setParam( pIMAX, 0.0F);
        myFile.setParam( pIMIN, 0.0F);
        myFile.setParam( pRMIN, rmin);
        cout << "Summed pix range " << rmin << " to " << rmax << endl;

        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
            myFile(ix,iy) = pixout[ix][iy];
        dx = myFile.getParam( pDX );
        dy = myFile.getParam( pDY );
        aimin = aimax = 0.0F;
        myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax, dx, dy );

    }

    return EXIT_SUCCESS;

} /* end main() */


/*  mimic old fft2d subroutines using FFTW
    this is very crude and not efficient but
    works for now
*/
/*-----------------  fft2d() ---------------- */
/*
        2D complex to complex FFT

        pixr[ix][iy] = real part of 2D pix to Fourier transform
        pixi[ix][iy] = imag part of 2D pix to Fourier transform
        nx,ny = (long int) size of array
                ix=0,1,2...(nx-1)
                iy=0,1,2...(ny-1)
        inverse = if > 0 then forward transform
                  if < 0 then forward transform

        On exit pixr and pixi will have the Fourier transform 
        and the original data will be lost.      
*/

void fft2d( float **pixr, float **pixi, const long nx,
                   const long ny, int inverse )
{
    int ix, iy;

    cfpix cpix( nx, ny );    /* complex arrays for fft */
    cpix.init( 1 );

    /*------  copy to FFTW style array -------
        not very elegant but... */
    for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            cpix.re(ix,iy) = pixr[ix][iy];  /* real */
            cpix.im(ix,iy) = pixi[ix][iy];  /* imag */
    }

    if( inverse > 0 ) cpix.fft();
    else if( inverse < 0 ) cpix.ifft();

    /*------  copy back to old style array -------
        not very elegant but... */
    for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            pixr[ix][iy] = cpix.re(ix,iy);  /* real */
            pixi[ix][iy] = cpix.im(ix,iy);  /* imag */
    }

}

/*------------------------- invert2D() ----------------------*/
/*
        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2D( float** pix, long nx, long ny )
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
}

