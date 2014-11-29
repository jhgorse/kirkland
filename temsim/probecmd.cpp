/*      *** probecmd.cpp ***

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

    front end for probe.cpp using a command line user interface

    ANSI-C++ version
    this version uses FFTW 3

    see:   www.fftw.org

    on Windows file libfftw3f-3.dll must be in the PATH

    on Linux build as:
    g++ -O -o probe probe.cpp slicelib.o cfpix.o
                       floatTIFF.cpp -lfftw3f

    Calculate a focused probe wavefunction in real space

    this file is formatted for a tab size of 4 characters

    rewritten in C 6-dec-1995 ejk
    fixed sign error in aberration function 1-mar-1997 ejk
    removed commas from keyboard input 3-oct-1997 ejk
    updated memory allocation routines 20-nov-1999 ejk
    change void main() to int main() for better portability
         22-jan-2000 ejk
    add Cs5=5th order aberration and astigmatism  19-jul-2005 ejk
    small cosmetic changes 18-jul-2007 ejk
    convert to GPL 3-jul-2008 ejk
        convert to large list aber. with coma option 23-nov-2008 ejk
    get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 10-apr-2010 ejk
    convert to FFTW 9-may-2010 ejk
    fix C34a,b terms 10-may-2010 ejk
    change astig. parameters to a,b from mag.+angle 30-may-2010 ejk
    add more aberations to 5th order 30-jun-2010 to 4-jul-2010 ejk
    add probe size calculation 5-jul-2010 ejk
    split up into subroutine 16-jul-2010 ejk
    fix a few things in prbSize() 5-sep-2010 ejk
    switch to storing multipole aberr. in param[] to save in
       output file 2-may-2011 ejk
    add multiMode in chi() to avoid extra calculations if there
       are no multipole aberrations 8-may-2011 ejk
    convert to C++ and floatTIFF.cpp  22-mar-2012 ejk
    convert to cfpix/fftw class from raw fftw 5-nov-2012 ejk
    fix small bug in aspect ratio of display image (actual
       floating point image was fine) and add probe position
       parameter symbolic constants  6-apr-2013 ejk
  convert to streams and strings 22-mar-2014 ejk
*/

#include <cstdio>  /*  ANSI-C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include <string>
#include <iostream>  //  C++ stream IO
#include <fstream>
#include <iomanip>   //  to format the output

using namespace std;

// define the following symbol to enable bound checking
//     define for debugging, and undefine for final run
//#define CFPIX_BOUNDS_CHECK
#include "cfpix.hpp"       // complex image handler with FFT
#include "slicelib.hpp"    // define parameter offsets
#include "floatTIFF.hpp"   // file I/O libraries
#include "probe.hpp"       // probe calculation

#define MANY_ABERR      /*  define to include many aberrations */

/*------------------------ main() ---------------------*/
int main()
{
    string fileout, cline;
    const char version[] = "21-apr-2014 (ejk)";
    int ix, iy, nx, ny, ixmid, iymid, i, ismoth, npixels,
        done, status, multiMode, NPARAM;
    float rmin, rmax, aimin, aimax, p2, scale;
    float *param, pixr, pixi, **pixsq;
    float *kx, *ky, *xpos, *ypos, *kx2, *ky2;
    double k2max, keV, wavlen, ax, by, rx, ry,
        rx2, ry2, pi, dx, dy, pixel, Cs3, Cs5, df, time, apert;
    double  x;

    cfpix cpix;
    probe prb;

    floatTIFF myFile;

    /*  Echo version date etc.  */
    cout << "probe version dated " << version << endl;
    cout << "Copyright (C) 1998-2014 Earl J. Kirkland" << endl;
    cout << "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        << " under the GNU general public license\n" << endl;

#ifdef MANY_ABERR
    cout << "calculate a focused probe wave function including multiple aberr.\n" << endl;
#else
    cout << "calculate a focused probe wave function\n" << endl;
#endif

    pi = 4.0 * atan( 1.0 );

    /*  memory to store parameters */
    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "probe-param" );
    for( i=0; i<NPARAM; i++) param[i] = 0.0F;

    /* ---- Get desired image size, parameters etc. ------------- */

    cout << "Name of file to get focused probe wave function:" << endl;
    cin >> fileout ;

    cout << "Desired size of output image in pixels Nx,Ny:" << endl;
    cin >> nx >> ny;

    cout << "Size of output image in Angstroms ax,by:" << endl;
    cin >>  ax >> by ;

    cout << "Probe parameters, V0(kv), Cs3(mm), Cs5(mm),"
           << " df(Angstroms), apert(mrad):" << endl;
    cin >> keV >> Cs3 >> Cs5 >> df >> apert;
    param[pDEFOCUS] = (float) df;
    param[pCS]  = (float) ( Cs3*1.0e7 );
    param[pCS5] = (float) ( Cs5*1.0e7 );

    cout << "Type 1 for smooth aperture:" << endl;
    cin >> ismoth;

    cout << "Probe position x,y in Ang.:" << endl;
    cin >> dx >> dy;

#ifdef MANY_ABERR
    /*   get higher order aberrations if necessary */
    cout << "type higher order aber. name (as C32a, etc.) followed\n"
        << " by a value in mm. (END to end)" << endl;
    done = multiMode = 0;
    do{
        cin >> cline ;
        if( cline.compare( "END" ) == 0  ) {
            done = 1;
        } else {
            cin >> x;
            status = readCnm( cline, param, x );        
            if( status < 0 ) {
                cout << "unrecognized aberration, exit..." << endl;
                exit( EXIT_SUCCESS );
            } else multiMode = 1;
        }
    } while( !done );

#endif

    /* ------- Calculate misc constants ------------ */

    time = cputim( );
    
    rx  = 1.0/ax;
    rx2 = rx * rx;
    ry  = 1.0/by;
    ry2 = ry * ry;
    
    ixmid = nx/2;
    iymid = ny/2;
    
    wavlen = wavelength( keV );
    cout << "electron wavelength = " << wavlen << " Angstroms" << endl;

    k2max = apert*0.001/wavlen;
    k2max = k2max * k2max;

    /* ------- allocate memory ------------ */

    pixsq = (float**) malloc2D( nx, ny, sizeof(float), "pixsq" );

    kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
    kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
    xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
    freqn( kx, kx2, xpos, nx, ax );

    ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
    ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
    ypos = (float*) malloc1D( ny, sizeof(float), "ypos" );
    freqn( ky, ky2, ypos, ny, by );

    cpix.resize( nx, ny );
    cpix.init( 1 );        //  only fast init and slow execution needed here

    /* --------- calculate probe wavefunction -------- */
    pixel = ( rx2 + ry2 );
    npixels = prb.makeProbe( cpix, nx, ny, dx, dy, 
        param, wavlen, k2max, pixel, multiMode, ismoth, kx, kx2, ky, ky2);
    cout << "there were " << npixels
        << " pixels inside the aperture" << endl;

    scale = (float) sqrt( ((double)nx) * ((double)ny) );
    cpix *= scale;   //   switch to normalization in reciprocal space

    /* -----  copy back for output ----- */
    myFile.resize( 2*nx, ny);
    myFile.setnpix( 2 );
    for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++) {
        myFile(ix,iy)    = pixr = cpix.re(ix,iy);  // real
        myFile(ix+nx,iy) = pixi = cpix.im(ix,iy);  // imag
        pixsq[ix][iy] = p2 = pixr*pixr + pixi*pixi;
    }

    //------- Output results and find min and max to echo ---------------

    rmin  = myFile.min(0);    // real part
    rmax  = myFile.max(0);
    aimin = myFile.min(1);   // imaginary
    aimax = myFile.max(1);

    param[pRMAX] = rmax;
    param[pIMAX] = aimax;
    param[pRMIN] = rmin;
    param[pIMIN] = aimin;
    param[pDEFOCUS]= (float) df;
    param[pDX]= (float) (ax / nx);
    param[pDY]= (float) (by / ny);
    param[pENERGY]= (float) keV;
    param[pWAVEL]= (float) ( sqrt(k2max) * wavlen);
    param[pCS]= (float) Cs3;
    param[pPPOSX]= (float) dx;
    param[pPPOSY]= (float) dy;

    param[pMODE] = 4;  // save mode = probe

    for( i=0; i<NPARAM; i++) myFile.setParam( i, param[i] );  // not very elegant

    if( myFile.write( fileout.c_str(), rmin, rmax, aimin, aimax,
            param[pDX], param[pDY] ) != 1 )
        cout << "probe cannot write an output file." << endl;

    cout << "Pix range "<< rmin << " to " << rmax << " real,\n"
        << "      and " << aimin << " to " << aimax << " imaginary\n" << endl;

    /*------- calculate probe size ---------------*/

    x = prb.prbSize( pixsq, nx, ny, dx, dy, ax, by );
    cout << "probe size (FWHM-II) = " << x << " Ang."  << endl;

    /*------- exit ---------------*/

    time = cputim() - time;
    cout << "\nCPU time = " << time << " sec." << endl;

    return EXIT_SUCCESS;

}  /* end main() */
