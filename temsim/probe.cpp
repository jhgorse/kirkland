/*      *** probe.cpp ***

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

    ANSI-C++ version

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
    fix small bug in aspect ration of display image (actual
       floating point image was fine) and add probe position
       parameter symbolic constants  6-apr-2013 ejk
   start conversion to separate class file 1-jul-2013 ejk
   move invert2D() to here (was in incostem) 7-jul-2013 ejk
*/

#include "probe.hpp"   /* file I/O libraries */

//=============================================================
//---------------  creator and destructor --------------

probe::probe()
{
}

probe::~probe()
{
}

/*------------------------- invert2D() ----------------------*/
/*----- below is from fft2dc.c ----------------------------- */
/*
        rearrange pix with corners moved to center (a la FFT's)
	prbsize() doesn't deal with wrap around properely

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)
*/
void probe::invert2D( float** pix, long nx, long ny )
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


/*   ------------------ makeProbe() ----------------------
    calculate aberration limited probe wavefunction in real space
    normalized to an integrated intensity of unity
  input:
    cpix[] = fftw array to get wave function
    nx, ny = size of wavefunciton in pixels
    xp, yp = probe position in Ang.
    p[]    = aberration values
    k2max  = max k^2
    pixel  = smoothing range
    multiMode = if not 0 then include multipole aberrations
    ismoth = flag to control smoothing at the edge

  returned value = number of pixels in the aperture

  add arguments  kx[], kx2[], ky[], ky2[] on 14-apr-2013 ejk
*/

int probe::makeProbe( cfpix &cpix, int nx, int ny, double xp, double yp,
    float p[], double wavlen, double k2max, double pixel, int multiMode,
    int ismoth, float kx[], float kx2[], float ky[], float ky2[] )
{ 
    int ix, iy, npixels;
    float scale;
    double alx, aly, k2, chi0, pi, dx2p, dy2p, sum, tr, ti;

    /*    PIXEL = diagonal width of pixel squared
        if a pixel is on the aperture boundary give it a weight
        of 1/2 otherwise 1 or 0
    */
    npixels = 0;
    pi = 4.0 * atan( 1.0 );

    dx2p = 2.0*pi*xp;
    dy2p = 2.0*pi*yp;

    for( iy=0; iy<ny; iy++) {
        aly = wavlen * ky[iy];  /* y component of angle alpha */
        for( ix=0; ix<nx; ix++) {
            k2 = kx2[ix] + ky2[iy];
            alx = wavlen * kx[ix];  /* x component of angle alpha */
            if ( ( ismoth != 0) && 
                ( fabs(k2-k2max) <= pixel) ) {
                   chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode )
                           - ( (dx2p*kx[ix]) + (dy2p*ky[iy]) );
                cpix.re(ix,iy) = (float) ( 0.5 * cos(chi0));  /* real */
                cpix.im(ix,iy) = (float) (-0.5 * sin(chi0));  /* imag */
                /* printf("smooth by 0.5 at ix=%d, iy=%d\n", ix, iy ); */
            } else if ( k2 <= k2max ) {
                   chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode )
                           - ( (dx2p*kx[ix]) + (dy2p*ky[iy]) );
                cpix.re(ix,iy) = (float)  cos(chi0);  /* real */
                cpix.im(ix,iy) = (float) -sin(chi0);  /* imag */
                npixels++;
            } else {
                cpix.re(ix,iy) = cpix.im(ix,iy) = 0.0F;
            }
        }  /* end for( ix=0... )  */
    }  /* end for( iy=0... )  */

    cpix.ifft();  //  inverse transform back to real space

    // -----  normalize probe intensity -----
    sum = 0.0;
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
            tr = cpix.re(ix,iy); 
            ti = cpix.im(ix,iy);
            sum += tr*tr + ti*ti;
    }
        
    // ----- Normalize probe intensity to unity and add weight ------------ 
    scale = (float)sqrt( 1.0 / sum ); 
    cpix *= scale;

    return( npixels );

}  // end probe::makeProbe()  

/* -------------------  message() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void probe::message( char msg[],  int level )
{
        messageSL( msg, level );  //  just call slicelib version for now

}  // end probe::message()

/* -------------------  prbsSize() -------------------
   calculate probe size 

  input:
     pixsq[][] = intensity in probe
     nx, ny = size of probe in pixels
     xp, yp = original probe position (in Ang.)
     ax, by = size of cell in Ang.

  return value:  size of probe in Ang.

*/
double probe::prbSize( float** pixsq, int nx, int ny,
    double xp, double yp, double ax, double by )
{
    int ix, iy, ir, nr;
    double *ivsr, dr, scale, rx, ry, x, y, y2;

    dr = ax/nx;
    scale = by/ny;
    if( scale < dr ) dr = scale;  //  smallest radial spacing
    dr = 0.5*dr;                  //  use sub pixel sampling
    if( ny > nx ) nr = ny;
    else nr = nx;
    nr = 2*nr;
    ivsr = (double*) malloc1D( nr, sizeof(double), "ivsr");
    for( ir=0; ir<nr; ir++) ivsr[ir] = 0.0;

    /*  find curent vs. radius = azimuthal average */
    /*    orig. probe position = (dx,dy)  */
    rx = ax/nx;
    ry = by/ny;
    for( iy=0; iy<ny; iy++) {
        y = (iy*ry) - yp;
        y2 = y*y;
        for( ix=0; ix<nx; ix++) {
            scale = pixsq[ix][iy];
            x = (ix*rx) - xp;
            ir = (int)( sqrt( x*x + y2 )/dr + 0.5);
        if( ir < 0 ) ir = 0;
        if( ir >= nr ) ir = nr-1;
            ivsr[ir] += scale;
       }
    }

    /*  integrate current */
    for( ir=1; ir<nr; ir++) ivsr[ir] += ivsr[ir-1];

    /*  find half current radius */
    scale = 0.5 * ivsr[nr-1];
    ix = 0;
    for( ir=0; ir<nr; ir++) {
        ix = ir;
        if( ivsr[ir] > scale ) break;
    }

    y = ivsr[ix] - ivsr[ix-1];   /* interpolate */
    if( fabs(y) > 0.0 ) y = (scale-ivsr[ix-1])/y;
    else y = 0.0;
    x = 2.0 * ( (ix-1)*dr + y*dr );

    free( ivsr );

    return( x );

}  //  end probe::prbSize()

