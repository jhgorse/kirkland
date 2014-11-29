/*      incostem.cpp

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

  class with member subroutines to
  calculate images in the incoherent STEM approximation

  put the partial cross section (integrated over the ADF detector
  angles) at a single pixels (position of the corresponding atom)
  and convolve with the point spread function (focused probe intensity)

  reference:

  [1] E. Kirkland, "Advanced Computing in Electron Microscopy",
        Plenum 1998, 2nd edit. Springer 2010
  
  this file is formatted for a tab size of 4 char

  started 16-mar-1998 E. Kirkland
  change float1D() etc to malloc1D() etc. and get rid of
     memory.h 18-jun-2007 ejk
  add cputim() echo 27-jun-2007 ejk
  put in adaptive quadrature for r and k integration 
        - faster and more accurate 28-jun-07 ejk
  add source size 13-jul-2008 ejk
  start adding defocus spread (Cc effects) 1-mar-2009 ejk
  add Cs5 on 7-mar-2009 ejk
  add Gauss-Hermite integration of defocus spread 15-mar-2009 ejk
  change wording on initial dialog 14-nov-2009 ejk
  fix bug in FWHM conversion (factor of two error) 13-nov-2010 ejk
  start convert to FFTW  8-feb-2011 ejk
  adding many aberration mode (C32a etc)  13-feb-2011
  get return value of scanf() to remove warnings from gcc 4.4
      13-feb-2011 ejk
  convert many aberration mode to use chi() from slicelib
      and fix a few small typos 21-may-2011 ejk
  convert &lf to &lg reading multipole aberration and fix
      small typo in comments 20-jun-2011 ejk
  fix small typo's in comments and formatting 27-jul-2011 ejk
  convert to C++ and floatTIFF.cpp  14-aug-2012 ejk
  convert to cfpix/fftw class from raw fftw 26-feb-2013 to 30-oct-2012 ejk
  fix typo in final param[] saving loop  17-mar-2013 ejk 
  fix typo in final param[] updates - keep existing aberr.
     and get rid of redundant variable wavelen (use wavl) 30-mar-2013 ejk
  move main calulation into subroutine incostemCal() so I can put this 
    whole thing somewhere else if needed 13-apr-2013 ejk
  separate calculation into a separate class for easy use
     in other places 20-apr-2013 ejk
  add message() here 21-apr-2013 ejk
  consolodate makeProbe() and prbSize() in probe.cpp and call it from
      here to avoid duplicating code 05-jul-2013 ejk
  fix a few small typo's 18-jul-2013 ejk
  convert calculate() to calculate2D() an add *occ[] 25-jul-2013 ejk
  convert message() to use string data 5-sep-2013 ejk
  fix bug in reusing xpos[] if ny>nx 4-dec-2013 ejk
*/

#include <cstdio>  /* standard ANSI libraries */
#include <cstdlib>
#include <cmath>

#include "cfpix.hpp"       /* complex image handler with FFT */
#include "probe.hpp"       //  probe calculation
#include "slicelib.hpp"    /* misc. routines for multislice */
#include "incostem.hpp"    /* file I/O routines in TIFF format */

#include <sstream>	// string streams

#define MANY_ABERR      /*  define to include many aberrations */

//=============================================================
//---------------  creator and destructor --------------

incostem::incostem()
{
        NZMAX=   103; // max atomic number Z

        FCNatomf = 0;  // integrate functions numbers
        FCNfemr  = 1;
        FCNfemi  = 2;

        twopi = 2.0 * (4.0 * atan( 1.0 ));
}

incostem::~incostem()
{
}

//=============================================================

//-------------------- atomsignal() ------------------------------
/*
    calculate single atom partial cross section 
    integrated from kmin to kmax (i.e. the ADF detector)
    (does NOT include objective aberration function)

    zatom = atomic number Z
    keV   = beam energy in keV
    thetamin  = minimum scattering angle in radians
    thetamax  = maximum scattering angle in radians

    the returned value is the scattering cross section
    in square-Angstroms integrated between thetamin and thetamax

    assumed global: twopi
*/
double incostem::atomsignal( int zatom, double keV, double thetamin, double thetamax )
{   

    double tol, t2, p[3];
    int maxsteps;
    
    tol = 1.0e-5;       /* only need approx. answer here */
    maxsteps = 10000;
    p[0] = zatom;

    /*  this requires that integrate45() be safe to call recursively */

    //t2 = twopi * integrate45( atomf, p, thetamin, thetamax, tol, maxsteps );
    t2 = twopi * integrate45( FCNatomf, p, thetamin, thetamax, tol, maxsteps );

    return( t2 );

} /* end incostem::atomsignal() */

//-------------------- atomf() ------------------------------
/*
    dummy function for atomsignal() to integrate

    t = scatt. angl. (in rad.)

    p[0] = zatom = atomic number Z

    assumed global constants: twopi, wavl
*/
double incostem::atomf( double t, double p[] )
{   
    int zatom;
    double k, rfe, ife, t1;

    zatom = (int)( p[0] + 0.1);
    k = t / wavl;

    feMoliere( k, zatom, &rfe, &ife );
    t1 = t * (rfe*rfe + ife*ife);  /* theta*fe^2 */

    return( t1 );

} /* end incostem::atomf() */


//---------------------------- BJ0 ----------------------------
/*
  real function to return zeroth order Bessel function of the argument
   using the polynomial expression given on page 369 of
       Abromowitz and Stegum
*/
double incostem::BJ0( double x )
{
    int i;
    double x3, p0, f0, y0, xa;
    
    static double p[] = {1.0, -2.2499997,  1.2656208, -0.3163866,
                   0.0444479, -0.0039444, 0.0002100 };
    static double f[] = { 0.79788456, -0.00000077, -0.00552740,
        -0.00009512,  0.00137237, -0.00072805, 0.00014476 };
    static double y[]  = {  -0.78539816,  -0.04166397, -0.00003954,
         0.00262573,  -0.00054125, -0.00029333, 0.00013558 };

    xa = fabs( x );
    if( xa < 3.0 ) {    
        x3 = xa/ 3.0;   x3 = x3 * x3;   p0 = p[6];
        for( i=5; i>=0; i=i-1)  p0 = p0 * x3 + p[i];
        return( p0 );
    } else {
        x3 = 3.0/xa;    f0 = f[6];  y0 = y[6];
        for(i=5; i>=0; i=i-1) {
            y0 = y0 * x3 + y[i];
            f0 = f0 * x3 + f[i];
        }
        y0 = xa + y0;
        return( f0*cos(y0)/sqrt( xa ) );
    } 

}  /* end incostem::BJ0 */

/*---------------------------- calculate2D() ----------------------------*/
/*
    main calculation of incostem
    collect here to facilitate moving into other programs

    started 13-apr-2013 ejk
    change name to just calculate() 21-apr-2013 ejk
    change to calculate2D() 25-jul-2013 ejk
    fix bug in reusing xpos[] if ny>nx 4-dec-2013 ejk

*/
void incostem::calculate2D( cfpix &pix, float param[], int multiMode, int natom,
    int Znum[], float x[], float y[], float occ[] )
{
    int i, ix, iy, nx, ny, ixmid, iymid, idf, ndf, ismoth, npixels;
    long nxl, nyl;

    float **pixsq, **ctf, *atoms;
    float *kx, *ky, *xpos, *kx2, *ky2;
    float ax, by, rx, ry, wr, wi, tr, ti, scale, weight;

    double sum, rx2, ry2, df, df0, ddf, ddf2, aobj, k2max, dx, dy, keV,
            xp, yp, pixel, ds, ds2, k2, thetamin, thetamax, pi, dsource;

    /*  absiccas and weights for Gauss-Hermite Quadrature 
       with exp(-x*x) weighting of integrand 
       from Abramowitz and Stegun, and Numerical Recipes */
    const int NGH=9;   /* number of Gauss-Hermete coeff. to use  */
    double xGH[]={ 3.190993201781528, 2.266580584531843, 1.468553289216668,
        0.723551018752838, 0.000000000000000, -0.723551018752838,
        -1.468553289216668,-2.266580584531843,-3.190993201781528};
    double wGH[]={3.960697726326e-005, 4.943624275537e-003 ,8.847452739438e-002,
        4.326515590026e-001, 7.202352156061e-001, 4.326515590026e-001,
        8.847452739438e-002, 4.943624275537e-003, 3.960697726326e-005};

    cfpix trans;
    probe prb;

    // ---- get setup params from param[]
    ax = param[ pAX ];
    by = param[ pBY ];
    nx = ToInt( param[ pNX ] );
    ny = ToInt( param[ pNY ] );

    if( (nx < 1) || (ny < 1) || (ax<0.0) || (by<0.0) ){
        // gcc requires this step
        sbuff = "bad parameters in incostem::calculate()";
        messageIN( sbuff );
        exit( 0 );
    }

    df0 = param[pDEFOCUS];
    ddf = param[pDDF];
    wavl = param[pWAVEL];
    aobj = param[pOAPERT];
    keV = param[pENERGY];
    thetamin = param[pMINDET];
    thetamax = param[pMAXDET];
    dsource = param[pSOURCE];

    pi = (float) (4.0 * atan( 1.0 ));

    wavl = wavelength( keV );
    sigmae = sigma( keV )/ 1000.0;

   // ------  calculate spatial frequencies and positions for future use ------ */
    rx = 1.0F/ax;
    rx2= rx*rx;
    ry = 1.0F/by;
    ry2= ry*ry;
    ixmid = nx/2;
    iymid = ny/2;
    nxl = nx;
    nyl = ny;

    //  xpos[] and ypos[]
    kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
    kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
    xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
    freqn( kx, kx2, xpos, nx, ax );

    //  reuse xpos[] if possible - its not used anywhere else in the function
    //  but remember to make it bigger if ny>nx
    if( ny > nx ){
        free( xpos );
        xpos = (float*) malloc1D( ny, sizeof(float), "xpos" );
    }

    ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
    ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
    freqn( ky, ky2, xpos, ny, by );

    /*------  make FFTW arrays and plans ------- */
    trans.resize( nx, ny );
    trans.init( 1 );

/* ------  Calculate the transfer function ------  */
    sbuff ="calculate the transfer function...";   // gcc requires this step 
    messageIN( sbuff );
    ctf = (float**) malloc2D( nx, ny, sizeof(float), "ctf" );
    pixsq = (float**) malloc2D( nx, ny, sizeof(float), "pixsq" );

    k2max = aobj / wavl;
    k2max = k2max * k2max;
    dx = ( ax/((float)nx) );
    dy = ( by/((float)ny) );

    if( ddf > 1.0 ) {       /*  df integration parameters */
        ndf = NGH;
        /* ddf2 = sqrt(log(2.0)/(ddf*ddf));   convert from HWHM */
        ddf2 = sqrt(log(2.0)/(ddf*ddf/4.0));  /* convert from FWHM */
    } else {
        ndf = 1;
        ddf2 = 0.0;
    }

    for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++) ctf[ix][iy] = 0.0;

    xp = 0.0;  /* probe position */
    yp = 0.0;
    pixel = 0.0;  /* smoothing size = not used */
    ismoth = 0;

    /*---- integrate over defocus spread if ddf is large enough ----*/
    /* use Gauss-Hermite quadrature and convert exp(-a^2df^2) to exp(-u^2)  */
    for( idf=0; idf<ndf; idf++) {

        if( ndf > 1 ){
                df = df0 + xGH[idf]/ddf2;
                weight = (float) wGH[idf];
        }  else {
                df = df0;
                weight = 1.0;
        }
        /*  diagnostic */
        sbuff = "df step "+ toString(idf)+ ", df= " +toString(df) + 
                ", weight= " + toString(weight);
        messageIN( sbuff, 0 );
            
        param[pDEFOCUS] = (float) df;
            
        /* --------- calculate probe wavefunction -------- */
        npixels = prb.makeProbe( trans, nx, ny, xp, yp, 
            param, wavl, k2max, pixel, multiMode, ismoth,
            kx, kx2, ky, ky2 );

        /* -----  normalize and save probe intensity ----- */
        sum = 0.0;
        for( ix=0; ix<nx; ix++)
        for( iy=0; iy<ny; iy++) {
            tr = trans.re(ix,iy); 
            ti = trans.im(ix,iy);
            pixsq[ix][iy] = tr*tr + ti*ti;
            sum +=  pixsq[ix][iy];
        }
        
        /* ----- Normalize probe intensity to unity and add weight ------------ */
        scale = weight * ((float)sqrt( 1.0 / sum ) ); 
        
        for( ix=0; ix<nx; ix++) 
           for( iy=0; iy<ny; iy++) {
            ctf[ix][iy] +=  scale * pixsq[ix][iy];
        }

    }  /*  end for(idf...)  */

    param[pDEFOCUS] = (float) df0;   // put back the original defocus

    sum = 0.0;
    for( ix=0; ix<nx; ix++)     /* find integrated intensity */
    for( iy=0; iy<ny; iy++)  sum += ctf[ix][iy];

    scale = (float) (1.0 / (dx*dy*sum) );
    for( ix=0; ix<nx; ix++)     /* normalize integrated intensity */
    for( iy=0; iy<ny; iy++) {
        trans.re(ix,iy) = scale * ctf[ix][iy];
        trans.im(ix,iy) = 0.0F;
    }

    /*  find probe size just for fun: 
      remember that prbSize() doesn't handle wrap-around properly so
      move probe to the center first */
    prb.invert2D( ctf, nx, ny );
    xp = 0.5*ax;
    yp = 0.5*by;
    sum = prb.prbSize( ctf, nx, ny, xp, yp, ax, by );  /* about orig center */
    sbuff = "probe size (FWHM-II) = " + toString(sum) + " Ang.";
    messageIN( sbuff, 0 );

    /* -----  make final transfer function ----- */
    trans.fft();

    /*------  allocate scratch arrays and cross section look-up-table array */

    atoms = (float*) malloc1D( NZMAX+1, sizeof(float), "atoms" );   /* cross-section LUT */
    for( i=0; i<=NZMAX; i++) atoms[i] = -1.0F;

    /* ------  Calculate 2D distribution of partial cross sections in specimen 
    use look-up-table for cross sections because it takes
    a lot of CPU cycles to calculate ------ */

    sbuff ="calculate the 2D specimen function...";  // gcc requires this step
    messageIN( sbuff );
    pix.resize( nx, ny );
    pix.copyInit( trans );

    for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++)
        pix.re(ix,iy) = pix.im(ix,iy) = 0.0F;

    for( i=0; i<natom; i++) {
        ix = (int) ( nx*x[i]/ax );      /* make sure coord. is inside image */
        while( ix < 0 ) ix = ix + nx;
        ix = ix % nx;

        iy = (int) ( ny*y[i]/by );
        while( iy < 0 ) iy = iy + ny;
        iy = iy % ny;

        if( atoms[Znum[i]] < 0.0F ) {
            atoms[Znum[i]] = (float)
            atomsignal( Znum[i], keV, thetamin, thetamax);
            sbuff = "the partial cross section for Z = " + toString(Znum[i])
                + " is " + toString(atoms[Znum[i]]) + " sq-Ang";
            messageIN( sbuff, 0 );
        }

        //pix.re(ix,iy) += atoms[Znum[i]];  //-replaced 25-jul-2013 ejk
        pix.re(ix,iy) += occ[i]*atoms[Znum[i]];
    }

/*------   Convolve specimen function with transfer function ------ */
/*  NOTE: could do this faster because both the psf and the image are real
      - use real to complex FFT 28-jun-2007 ejk */

    /* fft2d( pixr, pixi, nx, ny, +1); */
    pix.fft();

    dsource = 0.5* dsource;   /* convert diameter to radius */
    ds = pi*pi * dsource*dsource/log(2.0);  /* source size factor- convert to FWHM*/
    for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++) {
        k2 = kx2[ix] + ky2[iy];
        if( k2 <= 4.0*k2max ) {
            ds2 = exp( -ds*k2 );
            tr = trans.re(ix,iy);
            ti = trans.im(ix,iy);
            wr = pix.re(ix,iy);
            wi = pix.im(ix,iy);
            pix.re(ix,iy) = (float) ( (tr*wr - ti*wi) * ds2 );
            pix.im(ix,iy) = (float) ( (ti*wr + tr*wi) * ds2 );
        } else {
            pix.re(ix,iy) = 0;
            pix.im(ix,iy) = 0;
        }
    }
    /* fft2d( pixr, pixi, nx, ny, -1);  */
    pix.ifft();

    //--- return scaratch memory and exit

    free( kx );
    free( kx2 );
    free( xpos );

    free( ky );
    free( ky2 );

    free( atoms );

    free2D( (void**) ctf, nx );
    free2D( (void**) pixsq, nx );

    return;

}  // end incostem::calculate2D()

//-------------------- feMoliere() ------------------------------
/*
    calculate complex electron scattering factor in Moliere approx
    from the projected atomic potential Vz(r)

    as in eq. 5.18 of [1]

    k = wave vector magnitude in 1/Ang
    zatom = atomic number Z

    rfe,ife = (real,image) scattering factor in Angstroms

    assumed global constants: twopi, sigmae, wavl
*/
void incostem::feMoliere( double k, int zatom, double *rfe, double *ife )
{   
    int maxsteps;
    double x, rmin=1.0e-5, rmax=6.0, tol, p[3];

    x = twopi / wavl;

    tol = 1.0e-5;       /* only need approx. answer here */
    maxsteps = 10000;
    p[0] = k;
    p[1] = zatom;

    //*rfe = x * integrate45( femr, p, rmin, rmax, tol, maxsteps );
    //*ife = x * integrate45( femi, p, rmin, rmax, tol, maxsteps );

    *rfe = x * integrate45( FCNfemr, p, rmin, rmax, tol, maxsteps );
    *ife = x * integrate45( FCNfemi, p, rmin, rmax, tol, maxsteps );

    return;

} /* end incostem::feMoliere() */

//-------------------- femi() ------------------------------
/*
    dummy function for feMoliere() to integrate
    imag. part of integrand

    r = radial coord. (in Ang.)

    p[0] = k = wave vector magnitude in 1/Ang
    p[1] = zatom = atomic number Z

    assumed global constants: twopi, sigmae
*/
double incostem::femi( double r, double p[] )
{   
    int zatom;
    double k, t1, t2;

    k = p[0];
    zatom = (int)( p[1] + 0.1);

    t1 = r * BJ0( twopi * k * r );
    t2 = sigmae * vzatom( zatom, r );

    return( t1 * ( 1.0 - cos( t2 ) ) );

} /* end incostem::femi() */


//-------------------- femr() ------------------------------
/*
    dummy function for feMoliere() to integrate
    real part of integrand

    r = radial coord. (in Ang.)

    p[0] = k = wave vector magnitude in 1/Ang
    p[1] = zatom = atomic number Z

    assumed global constants: twopi, sigmae
*/
double incostem::femr( double r, double p[] )
{   
    int zatom;
    double k, t1, t2;

    k = p[0];
    zatom = (int)( p[1] + 0.1);

    t1 = r * BJ0( twopi * k * r );
    t2 = sigmae * vzatom( zatom, r );

    return( t1 * sin( t2 ) );

} /* end incostem::femr() */

//-------------------- fint() ------------------------------
/*
    dummy function to select integration functions
*/
double incostem::fint( int FCN, double r, double p[] )
{   
    double x;

    if( FCN == FCNatomf ) x = atomf( r, p );
    else if( FCN == FCNfemr ) x = femr( r, p );
    else if( FCN == FCNfemi ) x = femi( r, p );
    else x = 0.0;
            
    return( x );

} /* end incostem::fint() */


/*---------------------------- integrate45() ----------------------------*/
/*
    adaptive quadrature using a simplified auto step size
    5th order Runge-Kutta (i.e. the function only depends on x and not y)

    integrate f(x) from x=xmin to x=xmax by solving the ODE
    dy/dx = f(x) with initial condition y(xmin)=0

   NOTE: The value of the constant SMALL assumes that f(x) is of order
     unity (1). If this is not true then you should change the value
     of SMALL (must be positive).

    fint(x,p[]) = function to integrate, x=independent variable and
            p[] is a parameter array to pass to fint()
    p[]     = parameter array to pass to fint()
    xmin,xmax   = range of integration
    maxerror    = maximum allowable error at each step
    maxsteps    = maximum number of steps
    
    Auto Step Size 5h order Runge-Kutta function with embedded
    4th order Runge-Kutta using the 
    Cash-Karp coefficients as described in Section 16.2 of
    Numerical Recipes 2nd edit. by Press et al
    
    function pointers to member functions do NOT work like I
    would like so kluge the function pointer  21-apr-2013 ejk

    started 22-jul-2002 E. Kirkland
    fixed various special case failures 23-jul-2002 ejk 
    change to function number to work on class member functions
        21-apr-2013 ejk
*/
//double incostem::integrate45( double (*fint)(double x, double p[]), double p[],
double incostem::integrate45( int fcn, double p[],
    double xmin, double xmax, double maxerror, int maxsteps )
{
    int iter;
    double k1, k2, k3, k4, k5, k6, y2;
    double  h, x, dx, delta, yfinish, scale;
    const double SMALL=1.0e-40;

    //--- the Cash-Karp Runge-Kutta coefficients ----
    static const double a2=1.0/5.0, a3=3.0/10.0, a4=3.0/5.0, a5=1.0, a6=7.0/8.0;
    static const double c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0;
    static const double cs1=2825.0/27648.0, cs3=18575.0/48384.0,
            cs4=13525.0/55296.0, cs5=277.0/14336.0, cs6=1.0/4.0;

    iter = 0;

    h = (xmax - xmin)/5.0;      //  for lack of a better start step size
    x = xmin;

    yfinish = 0.0;

    while( (x < xmax) && (iter < maxsteps) ) {
        //k1 = fint( x, p );
        k1 = fint( fcn, x, p );
        scale = fabs(yfinish) + fabs(h*k1);
        if( scale < SMALL ) scale = SMALL; // in case of zero starting point

        do {
            dx = h;
            iter += 1;
            if( iter >= maxsteps) {
                 sbuff = "Warning maxiter exceeded in integrate45()!";
                 messageIN( sbuff );
            }

            //k2 = fint( x+a2*h, p );
            //k3 = fint( x+a3*h, p );
            //k4 = fint( x+a4*h, p );
            //k5 = fint( x+   h, p );
            //k6 = fint( x+a6*h, p );

            k2 = fint( fcn, x+a2*h, p );
            k3 = fint( fcn, x+a3*h, p );
            k4 = fint( fcn, x+a4*h, p );
            k5 = fint( fcn, x+   h, p );
            k6 = fint( fcn, x+a6*h, p );

            /* the next point */
            y2 = h * ( c1*k1 + c3*k3 + c4*k4 + c6*k6 );

            /* the error term */
            delta = y2 - h * ( cs1*k1 + cs3*k3 + cs4*k4 + cs5*k5 + cs6*k6 );

            delta = fabs( delta/scale );
            
            if( (delta < 0.5*maxerror) && ( delta> SMALL)  ) 
                h = h * pow( fabs(maxerror/delta), 0.2 );
            else if( delta > maxerror ) h = h/2.0;

        } while ( delta > maxerror );
        
        if( fabs(h) < SMALL*fabs(x) )    // guard against h->0
            h = SMALL * fabs(x) * h/fabs(h);

        x = x + dx; // large round off error here but there is no other way
        if( (x+h) > xmax ) h = xmax - x;  // don't go past the end
        yfinish = y2 + yfinish;
    }

    return( yfinish );
    
}  /* end incostem::integrate45() */


/* -------------------  messageIN() -------------------
   message output
   direct all output message here to redirect to the command line
   or a GUI status line or message box when appropriate

   msg[] = character string with message to disply
   level = level of seriousness
        0 = simple status message
        1 = significant warning
        2 = possibly fatal error
*/
void incostem::messageIN( std::string &smsg,  int level )
{
        messageSL( smsg.c_str(), level );  //  just call slicelib version for now

}  // end incostem::message()

