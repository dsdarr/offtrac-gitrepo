/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *  By Robert Hallberg, November 1998                                 *
 *                                                                    *
 *    This program contains 1 subroutine.  set_metrics calculates     *
 *  the various metric terms that are used by HIM.  This routine is   *
 *  intended to be modified by the user to enable the use of any      *
 *  general orthogonal grid.                                          *
 *                                                                    *
 *    This subroutine is also used by HIMtocdf.                       *
 *                                                                    *
 *    The metric terms have the form dzp, Idzp, or dxDYp, where z can *
 *  be x or y, and p can be q, u, v, or h.  z describes the direction *
 *  of the metric, while p describes the location.  ILzp is the       *
 *  inverse of Lzp, while dxdyp is the produce of dxp and dyp.        *
 *                                                                    *
 *    Variables written all in capital letters are defined in init.h. *
 *                                                                    *
 *     A small fragment of the C-grid is shown below:                 *
 *                                                                    *
 *    j+1  x ^ x ^ x   At x:  q, dxq, Idxq, dyq, Idyq, etc.           *
 *    j+1  > o > o >   At ^:  v, dxv, Idxv, dyv, Idyv, etc.           *
 *    j    x ^ x ^ x   At >:  u, dxu, Idxu, dyu, Idyu, etc.           *
 *    j    > o > o >   At o:  h, dxh, Idxh, dyh, Idyh, dxdyh, etc.    *
 *    j-1  x ^ x ^ x                                                  *
 *        i-1  i  i+1  At x & ^:                                      *
 *           i  i+1    At > & o:                                      *
 *                                                                    *
 *  The boundaries always run through q grid points (x).              *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "init.h"
#include "metrics.h"

double latq[NYMEM];                  /* The latitude of q points.     */
double lath[NYMEM];                  /* The latitude of h points.     */
double lonq[NXMEM];                  /* The longitude of q points.    */
double lonh[NXMEM];                  /* The longitude of h points.    */
#if defined(PARALLEL_Y) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
double lath_tot[NYTOT+Y1];           /* The latitude of all points in */
double latq_tot[NYTOT+Y1];           /* the whole domain.             */
#endif
#if defined(PARALLEL_X) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
double lonh_tot[NXTOT+X1];           /* The longitude of all points   */
double lonq_tot[NXTOT+X1];           /* in the whole domain.          */
#endif
#if defined(PARALLEL_X) || defined(PARALLEL_Y)
extern int pe_here;                  /* The current processor's label.*/
#else
#define pe_here 0                    /* Without running in parallel,  */
                                     /* the current processor is 0.   */
#endif

#ifdef PARALLEL_X
extern int nx;                       /* The number of x-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
extern int X0abs;                    /* The absolute index of X0 on   */
                                     /* the current processor, rel-   */
                                     /* ative to X0 on processor 0.   */
#else
#define X0abs 0                      /* X0abs on processor 0 is 0.    */
#endif
#ifdef PARALLEL_Y
extern int ny;                       /* The number of y-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
int Y0abs;                           /* The absolute index of Y0 on   */
                                     /* the current processor, rel-   */
                                     /* ative to Y0 on processor 0.   */
#else
#define Y0abs 0                      /* Y0abs on processor 0 is 0.    */
#endif

#ifdef ISOTROPIC
double find_root(double fnval,double y,int *ittmax);
#define INTSECY(y) (((y) >= 0.0) ? (log((1.0 + sin(y))/cos(y))) : (-log((1.0 - sin(y))/cos(y))))
#endif

extern double qlat[NYMEM],hlat[NYMEM];

void set_metrics(void)
{
  int i,j;
  extern double areagr[NXMEM][NYMEM];
  double Iareagr[NXMEM][NYMEM];
  extern double D[NXMEM][NYMEM];

/*    Calculate the values of the metric terms that might be used     */
/*  and save them in arrays.                                          */

#ifdef CARTESIAN
/*   On a cartesian grid, the various DX... and DY... macros all      */
/* point to the same scalars.                                         */
  i = 0; j = 0;
  DXh(i,j) = RE * LENLON * M_PI / (180.0 * NXTOT);
  DYh(i,j) = RE * LENLAT * M_PI / (180.0 * NYTOT);
  IDXh(i,j) = 1.0 / DXh(i,j);
  IDYh(i,j) = 1.0 / DYh(i,j);
  DXDYh(i,j) = DXh(i,j) * DYh(i,j);
  IDXDYh(i,j) = IDXh(i,j) * IDYh(i,j);
  for (j=Y0-1;j<=ny+2;j++)
    latq[j] = LOWLAT + LENLAT*(double)(j-Y0+Y0abs)/(double)NYTOT;
  for (j=Y0;j<=ny;j++)
    lath[j] = LOWLAT + LENLAT*((double)(j-Y0+Y0abs)-0.5)/(double)NYTOT;
  for (i=X0-1;i<=nx+2;i++)
    lonq[i] = WESTLON + LENLON*(double)(i-X0+X0abs)/(double)NXTOT;
  for (i=X0;i<=nx;i++)
    lonh[i] = WESTLON + LENLON*((double)(i-X0+X0abs)-0.5)/(double)NXTOT;
# if defined(PARALLEL_Y) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
  for (j=Y0;j<=NYTOT+Y0;j++)
    latq_tot[j] = LOWLAT + LENLAT*(double)(j-Y0)/(double)NYTOT;
  for (j=Y0;j<=NYTOT+Y0;j++)
    lath_tot[j] = LOWLAT + LENLAT*((double)(j-Y0)-0.5)/(double)NYTOT;
# endif
# if defined(PARALLEL_X) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
  for (i=X0;i<=NXTOT+X0;i++)
    lonq_tot[i] = WESTLON + LENLON*(double)(i-X0)/(double)NXTOT;
  for (i=X0;i<=NXTOT+X0;i++)
    lonh_tot[i] = WESTLON + LENLON*((double)(i-X0)-0.5)/(double)NXTOT;
# endif
#else

/*   All of the metric terms should be defined over the domain from   */
/* X0-1 to nx+2.  Outside of the physical domain, both the metrics    */
/* and their inverses may be set to zero.                             */
/*   Any points that are outside of the computational domain should   */
/* have their values set to zero _BEFORE_ setting the other metric    */
/* terms, because these macros may or may not expand to 2-dimensional */
/* arrays.                                                            */

/*  The metric terms within the computational domain are set here.    */

# ifdef ISOTROPIC

  {
    double C0, I_C0, yq, yh, jd;
    double fnRef, yRef; /* fnRef is the value of the integral of      */
                        /* 1/(dx*cos(lat)) from the equator to a      */
                        /* reference latitude, while yRef is the      */
                        /* j-index of that reference latitude.        */
    int itt1, itt2;

    C0 = M_PI*((double) LENLON / (double) (180*NXTOT)); I_C0 = 1.0 / C0;

/*    With an isotropic grid, the north-south extent of the domain,   */
/*  the east-west extent, and the number of grid points in each       */
/*  direction are _not_ independent.  Here the north-south extent     */
/*  will be determined to fit the east-west extent and the number of  */
/*  grid points.  The grid is perfectly isotropic.                    */
#  ifdef NORTHREFERENCE
/*  The following 2 lines fixes the refererence latitude at the       */
/*  northern edge of the domain, LOWLAT+LENLAT at j=NYTOT+Y0.         */
    yRef = Y0abs - Y0 - NYTOT;
    fnRef = I_C0 * INTSECY(((LOWLAT+LENLAT)*M_PI/180.0));
#  else
/*  The following line sets the refererence latitude LOWLAT at j=Y0.  */
    yRef = Y0abs - Y0;   fnRef = I_C0 * INTSECY((LOWLAT*M_PI/180.0));
#  endif

/*  Everything else should pretty much work.                          */
    yq = LOWLAT*M_PI/180.0;
/* If the model is in parallel in the Y-direction, do the same set of */
/* calculations which would occur on a single processor.              */
    for (j=Y0-1-Y0abs;j<Y0-1;j++) {
      jd = fnRef + (double) (j + yRef) - 0.5;
      yh = find_root(jd,yq,&itt1);

      jd = fnRef + (double) (j + yRef);
      yq = find_root(jd,yh,&itt2);
#  if defined(PARALLEL_Y) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
      latq_tot[j+Y0abs] = yq*180.0/M_PI;
      lath_tot[j+Y0abs] = yh*180.0/M_PI;
#  endif
    }
    
    for (j=Y0-1;j<=ny+2;j++) {
      jd = fnRef + (double) (j + yRef) - 0.5;
      yh = find_root(jd,yq,&itt1);

      jd = fnRef + (double) (j + yRef);
      yq = find_root(jd,yh,&itt2);

      latq[j] = yq*180.0/M_PI;
      lath[j] = yh*180.0/M_PI;
#  if defined(PARALLEL_Y) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
      latq_tot[j+Y0abs] = yq*180.0/M_PI;
      lath_tot[j+Y0abs] = yh*180.0/M_PI;
#  endif

      for (i=X0-1;i<=nx+2;i++) {
        DXq(i,j) = cos(yq) * (RE * C0);
        DYq(i,j) = DXq(i,j);

        DXv(i,j) = DXq(i,j);
        DYv(i,j) = DYq(i,j);

        DXh(i,j) = cos(yh) * (RE * C0);
        DYh(i,j) = DXh(i,j);

        DXu(i,j) = DXh(i,j);
        DYu(i,j) = DYh(i,j);
      }
    }
#  if defined(PARALLEL_Y) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
    for (j=ny+3;j<=NYTOT+Y0-Y0abs;j++) {
      jd = fnRef + (double) (j + yRef) - 0.5;
      yh = find_root(jd,yq,&itt1);

      jd = fnRef + (double) (j + yRef);
      yq = find_root(jd,yh,&itt2);
      latq_tot[j+Y0abs] = yq*180.0/M_PI;
      lath_tot[j+Y0abs] = yh*180.0/M_PI;
    }
#  endif
  }
  for (i=X0-1;i<=nx+2;i++)
    lonq[i] = WESTLON + LENLON*(double)(i-X0+X0abs)/(double)NXTOT;
  for (i=X0;i<=nx;i++)
    lonh[i] = WESTLON + LENLON*((double)(i-X0+X0abs)-0.5)/(double)NXTOT;
#  if defined(PARALLEL_X) && !defined(PARALLEL_IO) && defined(NETCDF_OUTPUT)
  for (i=X0;i<=NXTOT+X0;i++)
    lonq_tot[i] = WESTLON + LENLON*(double)(i-X0)/(double)NXTOT;
  for (i=X0;i<=NXTOT+X0;i++)
    lonh_tot[i] = WESTLON + LENLON*((double)(i-X0)-0.5)/(double)NXTOT;
#  endif

# else         /*  ISOTROPIC  */

/* This code implements latitude/longitude coordinates on a sphere.   */

//BX-a
//    for (j=0;j<=NYTOT-1;j++) {
//HF    for (j=Y0-1;j<=ny+2;j++) {
   for (j=2;j<NYTOT+2;j++) {
	latq[j] = qlat[j];
	lath[j] = hlat[j];
    }
//BX-e
/* HF  for (i=0;i<=NXTOT-1;i++) {
    for (j=0;j<=NYTOT-1;j++) {
        DXq(i+2,j+2) = dxq[j][i];
	DYq(i+2,j+2) = dyq[j][i];
	DXv(i+2,j+2) = dxv[j][i];
	DYv(i+2,j+2) = dyv[j][i];
        DXh(i+2,j+2) = dxh[j][i];
	DYh(i+2,j+2) = dyh[j][i];
	DXu(i+2,j+2) = dxu[j][i];
	DYu(i+2,j+2) = dyu[j][i];
    }
    } */

//BX-a
  for (i=X0-1;i<=nx+2;i++)
    lonq[i] = WESTLON + LENLON*(double)(i-X0+X0abs)/(double)NXTOT;

  for (i=X0;i<=nx;i++)
    lonh[i] = WESTLON + LENLON*((double)(i-X0+X0abs)-0.5)/(double)NXTOT;
//BX-e
#endif


/* The remaining code should not be changed.                         */
  for (i=X0-1;i<=nx+2;i++) {
    for (j=Y0-1;j<=ny+2;j++) {
      IDXh(i,j) = (DXh(i,j) > 0.0) ? (1.0 / DXh(i,j)) : 0.0;
      IDXu(i,j) = (DXu(i,j) > 0.0) ? (1.0 / DXu(i,j)) : 0.0;
      IDXv(i,j) = (DXv(i,j) > 0.0) ? (1.0 / DXv(i,j)) : 0.0;
      IDXq(i,j) = (DXq(i,j) > 0.0) ? (1.0 / DXq(i,j)) : 0.0;
      IDYh(i,j) = (DYh(i,j) > 0.0) ? (1.0 / DYh(i,j)) : 0.0;
      IDYu(i,j) = (DYu(i,j) > 0.0) ? (1.0 / DYu(i,j)) : 0.0;
      IDYv(i,j) = (DYv(i,j) > 0.0) ? (1.0 / DYv(i,j)) : 0.0;
      IDYq(i,j) = (DYq(i,j) > 0.0) ? (1.0 / DYq(i,j)) : 0.0;

      //HF alternate def: DXDYh(i,j) = DXh(i,j) * DYh(i,j);
      //HF alternate def: IDXDYh(i,j) = IDXh(i,j) * IDYh(i,j);
      DXDYq(i,j) = DXq(i,j) * DYq(i,j);
      IDXDYq(i,j) = IDXq(i,j) * IDYq(i,j);

      Iareagr[i][j] = (areagr[i][j] > 0.0) ? (1.0 / areagr[i][j]) : 0.0;
      DXDYh(i,j) = areagr[i][j];
      IDXDYh(i,j) = Iareagr[i][j];
    }
  }

#endif /* CARTESIAN */
}


#ifdef ISOTROPIC
double find_root(double fnval,double y,int *ittmax)
{
  double ybot, ytop, fnbot, fntop, fny, dy_dfn, dy;
  double C0, I_C0;
  int itt;

  C0 = M_PI*((double) LENLON / (double) (180*NXTOT)); I_C0 = 1.0 / C0;
/*  Bracket the root. */
  ybot = y;
  fnbot = I_C0 * INTSECY(ybot) - fnval;
  while (fnbot > 0.0) {
    ybot = ((ybot - 2.0*C0*cos(ybot)) < (0.5*ybot-0.25*M_PI)) ? 
            (ybot - 2.0*C0*cos(ybot)) : (0.5*ybot-0.25*M_PI);
    fnbot = I_C0 * INTSECY(ybot) - fnval;
  }
  ytop = ((y + 2.0*C0*cos(y)) < (0.5*y+0.25*M_PI)) ? 
          (y + 2.0*C0*cos(y)) : (0.5*y+0.25*M_PI);
  fntop = I_C0 * INTSECY(ytop) - fnval;
  while (fntop < 0.0) {
    ytop = ((ytop + 2.0*C0*cos(ytop)) < (0.5*ytop+0.25*M_PI)) ? 
            (ytop + 2.0*C0*cos(ytop)) : (0.5*ytop+0.25*M_PI);
    fntop = I_C0 * INTSECY(ytop) - fnval;
  }
/*  Bisect several times to insure that the root is within the radius */
/*  of convergence in the Newton's method polisher.                   */
  for (itt=0;itt<10;itt++) {
    y = 0.5*(ybot + ytop);
    fny = I_C0 * INTSECY(y) - fnval;
    if (fny < 0.0) { fnbot = fny; ybot = y; }
    else { fntop = fny; ytop = y; }
  }

/*    Polish the root using Newton's method.                          */
  for (itt=0;itt<10;itt++) {
    dy_dfn = C0 * cos(y);
    fny = I_C0 * INTSECY(y) - fnval;

    dy = - fny * dy_dfn; y += dy;
    if (y > ytop) y = ytop; if (y < ybot) y = ybot;
    if (fabs(dy) < (8.0e-15*fabs(y)+1.e-20)) break;
  }
  if (fabs(y) < 1e-12) y = 0.0;

   *ittmax = itt;
   return(y);
}
#endif
