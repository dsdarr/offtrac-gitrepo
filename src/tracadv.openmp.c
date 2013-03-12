/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *        3-dimensional offline tracer advection routine              *
 *            David Darr - October 2002/February 2003                 *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// begin yanxu
#include <omp.h>
// end yanxu    
#include "init.h"
#include "io.h"
#include "metrics.h"
#include "offtrac.h"
#include "alloc.h"



static double dmax1, dmax2, dmin1, dmin2;
#define D_MAX(a,b) (dmax1=(a),dmax2=(b),((dmax1) > (dmax2)) ?\
                   (dmax1) : (dmax2))
#define D_MIN(a,b) (dmin1=(a),dmin2=(b),((dmin1) < (dmin2)) ?\
                   (dmin1) : (dmin2))

extern double dt;
extern double h[NZ][NXMEM][NYMEM];	/* Layer thickness, begin of advec step in m.  */

extern double hstart[NZ][NXMEM][NYMEM];	/* Layer thickness, begin of advec step in m.  */

extern double hend[NZ][NXMEM][NYMEM];	/* Layer thickness, end of advec step in m.  */
#ifdef HTEST
extern double htest[NZ][NXMEM][NYMEM];	/* Layer thickness, used for debugging  */
#endif

extern double ****tr; /* Tracer concentration g m-3.*/
extern double           D[NXMEM][NYMEM]; /* Basin depth, in m.         */

extern double ***uhtm; /* uhtm and the vhtm are the sum */
extern double ***vhtm; /* of uh and vh between tracer   */
                    /* updates, both in m3 s-1.      */

extern double umask[NXMEM][NYMEM];   /* _mask are 1 over ocean and 0  */
extern double vmask[NXMEM][NYMEM];   /* over land on the u & v grids. */

extern double ***ea;
extern double ***eb;
extern double eaml[NXMEM][NYMEM];
extern double wd[NZ+1][NXMEM][NYMEM]; /* the diapycnal velocity that is */
                                      /* calculated at the interfaces   */
extern char directory2[95];          /* The directory to use to save  */
                                     /* output files.                 */

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
#endif
#ifdef PARALLEL_Y
extern int ny;                       /* The number of y-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
#endif

// for debugging
#ifdef AGE
  extern int mAGE;
#endif
double hvolint;
extern double trintegral[NTR];

#ifdef DIFFUSE_TRACER
static void diffuse_tracer();
#endif

static void tracer_hordiff();
void print_tr(int pstage);
void hvol_integral(double ***hvol);
void hvol_kintegral(int k, double ***hvol);
void tracer_integral(int trnum, double ***hvol);
void tracer_kintegral(int trnum, int k, double ***hvol);

void tracer(int itts)
{

/*    This subroutine time steps the tracer concentration.            */
/*  A positive definite scheme is used.                               */


  double minslope;          /* The maximum concentration slope per    */
                            /* grid point consistent with mono-       */
                            /* tonicity, in conc. (nondim.).          */

  double ***hvol; /* The cell volume of an h-element   */

  double slope[NXMEM+NYMEM][NTR]; /* The concentration slope per grid */
                        /* point in units of concentration (nondim.). */
  double fluxtr[NXMEM+NYMEM][NTR];/* The flux of tracer across a      */
                        /* boundary, in m3 * conc. (nondim.).         */


  double ***uhr; /* The remaining zonal and meridional */
  double ***vhr; /* thickness fluxes, in m3.*/

  double uhh[NXMEM];        /* uhh and vhh are the reduced fluxes     */
  double vhh[NYMEM];        /* during the current iteration, in m3.d  */

  double hup, hlos;         /* hup is the upwind volume, hlos is the  */
                            /* part of that volume that might be lost */
                            /* due to advection out the other side of */
                            /* the grid box, both in m3.              */
  double ts2;
  double landvolfill;       /* An arbitrary? nonzero cell volume, m3. */

  double ***ear;
  double ***ebr;
  double ***wdh;
  
  double bet[NXMEM];        /* bet and gam are variables used by the  */
  double gam[NZ][NXMEM];    /* tridiagonal solver.                    */
  double hnew0[NXMEM];      /* The original topmost layer thickness,  */
#if defined AGE2 || defined AGE3
  //  extern double hnew[NZ][NXMEM][NYMEM];
  extern double ***hnew;
#else
  double ***hnew;
#endif

double hlst1, Ihnew;
double hlst[NYMEM];

//  double MLMIN = EPSILON;   /* min depth for ML			      */

	double MLMIN = 4.25;
	double BLMIN = 0.20;

#ifdef ENTRAIN
  double nts = dt/DT; /* number of timesteps (#day*86400/3600seconds) */
#endif
  int i, j, k, m, ii, pstage;
  int itt;
  double fract1;
  double fract2;
# ifdef WRTTS
  double wrts;
# endif

  hvol = alloc3d(NZ,NXMEM,NYMEM);
    if(hvol == NULL) {
	fprintf(stderr,"not enough memory for hvol!\n");
    }
  uhr = alloc3d(NZ,NXMEM,NYMEM);
    if(uhr == NULL) {
	fprintf(stderr,"not enough memory for uhr!\n");
    }
  vhr = alloc3d(NZ,NXMEM,NYMEM);
    if(vhr == NULL) {
	fprintf(stderr,"not enough memory for vhr!\n");
    }
  ear = alloc3d(NZ,NXMEM,NYMEM);
    if(ear == NULL) {
	fprintf(stderr,"not enough memory for ear!\n");
    }
  ebr = alloc3d(NZ,NXMEM,NYMEM);
    if(ebr == NULL) {
	fprintf(stderr,"not enough memory for ebr!\n");
    }
  wdh = alloc3d(NZ,NXMEM,NYMEM);
    if(wdh == NULL) {
	fprintf(stderr,"not enough memory for wdh!\n");
    }
#if !defined AGE2 && !defined AGE3
  hnew = alloc3d(NZ,NXMEM,NYMEM);
    if(hnew == NULL) {
	fprintf(stderr,"not enough memory for hnew!\n");
    }
#endif

    landvolfill = EPSILON*1000000.0;    /* This is arbitrary.	*/

		/* zonal re-entrance		*/

#pragma omp parallel 
{
#pragma omp for  private(j,k)
    for (j=0;j<=NYMEM-1;j++) {
      for (k=0;k<NZ;k++) {
        uhtm[k][nx+1][j] = uhtm[k][2][j];
        uhtm[k][nx+2][j] = uhtm[k][3][j];
        uhtm[k][0][j] =   uhtm[k][nx-1][j];
        uhtm[k][1][j] =   uhtm[k][nx][j];
        vhtm[k][nx+1][j] = vhtm[k][2][j];
        vhtm[k][nx+2][j] = vhtm[k][3][j];
        vhtm[k][0][j] =   vhtm[k][nx-1][j];
        vhtm[k][1][j] =   vhtm[k][nx][j];
      }

      for (k=0;k<NZ+1;k++) {
        wd[k][nx+1][j] = wd[k][2][j];
        wd[k][nx+2][j] = wd[k][3][j];
        wd[k][0][j] =   wd[k][nx-1][j];
        wd[k][1][j] =   wd[k][nx][j];
      }
    }

	/* meridional re-entrance            */
#pragma omp for  private(i,k,ii)
    for (i=2;i<=nx;i++) {
     ii = 363 - i;
      for (k=0;k<NZ;k++) {
        uhtm[k][ii][ny+1] = (-1)*uhtm[k][i][ny];
        uhtm[k][ii][ny+2] = (-1)*uhtm[k][i][ny-1];
        vhtm[k][ii][ny+1] = (-1)*vhtm[k][i][ny];
        vhtm[k][ii][ny+2] = (-1)*vhtm[k][i][ny-1];
      }
       for (k=0;k<NZ+1;k++) {
        wd[k][ii][ny+1] = wd[k][i][ny];
        wd[k][ii][ny+2] = wd[k][i][ny-1];
      }
    }


#pragma omp for  private(i,j,k)
  for (k=0;k<NZ;k++)  {
/*  Put the thickness fluxes into uhr and vhr.                  */
    for (j=0;j<=ny+2;j++) {
	for (i=0;i<=nx+2;i++) {

	    uhr[k][i][j] = uhtm[k][i][j]*dt;
	    vhr[k][i][j] = vhtm[k][i][j]*dt;

	    if (h[k][i][j] < EPSILON) {
		h[k][i][j] = 1.0*EPSILON;
	    }

/*   This line calculates the cell volume                       */
        hvol[k][i][j] = DXDYh(i,j)*h[k][i][j];
        hnew[k][i][j] = h[k][i][j];
      }
    }
  }


/* calculate the diapycnal velocities at the interfaces		*/
/*   if we read in the ea, eb and eaml variables                */
/*   Otherwise we read in wd directly                           */

#ifdef ENTRAIN

#pragma omp for  private(i,j)
  for (i=X0;i<=nx+1;i++)                               
      for (j=Y0;j<=ny;j++)
        wd[0][i][j] = nts*eaml[i][j];                        

#pragma omp for  private(i,j,k)
     for (k=1;k<NZ;k++) {
      for (i=X0;i<=nx+1;i++)
	  for (j=Y0;j<=ny;j++)
	      wd[k][i][j] = nts*(ea[k][i][j] - eb[k-1][i][j]); 
      }
#endif

} // omp

#define STANDARD_ADVECTION
//#undef STANDARD_ADVECTION
#ifdef STANDARD_ADVECTION
    /*
    pstage=1;
    print_tr(pstage);
    */

  /* beginning of itt loop */
    for (itt = 0; itt < NUM_ADV_ITER; itt++) {

      /* big loop over k	 */
//ompfail 
#pragma omp parallel 
{
#pragma omp for private(i,j,k,m,minslope,slope,uhh,vhh,fluxtr,hup,hlos,ts2,hlst,hlst1,Ihnew)
      for (k=0;k<NZ;k++)
	{ 
/*    To insure positive definiteness of the thickness at each        */
/*  iteration, the mass fluxes out of each layer are checked each     */
/*  time.  This may not be very efficient, but it should be reliable. */

/* ============================================================ */
/*			first advect zonally			*/
/* ============================================================ */
#ifndef ADV1D
	  for (j=Y1;j<=ny;j++) {

/*   Calculate the i-direction profiles (slopes) of each tracer that  */
/* is being advected.                                                 */
//#pragma omp for  private(i,m,minslope)
	    for (i=X0;i<=nx+1;i++) {
	      for (m=0;m<NTR;m++) {
		minslope = 4.0*((fabs(tr[m][k][i+1][j]-tr[m][k][i][j]) < 
				 fabs(tr[m][k][i][j]-tr[m][k][i-1][j])) ? 
				(tr[m][k][i+1][j]-tr[m][k][i][j]) :
				(tr[m][k][i][j]-tr[m][k][i-1][j]));
		slope[i][m] = umask[i][j]*umask[i-1][j] *
		  (((tr[m][k][i+1][j]-tr[m][k][i][j]) * 
		    (tr[m][k][i][j]-tr[m][k][i-1][j]) < 0.0) ? 0.0 :
		   ((fabs(tr[m][k][i+1][j]-tr[m][k][i-1][j])<fabs(minslope)) ?
		    0.5*(tr[m][k][i+1][j]-tr[m][k][i-1][j]) : 0.5*minslope));
	      }
	    }
            //#pragma omp barrier

/*   Calculate the i-direction fluxes of each tracer, using as much   */
/* the minimum of the remaining mass flux (uhr) and the half the mass */
/* in the cell plus whatever part of its half of the mass flux that   */
/* the flux through the other side does not require.                  */
//#pragma omp for  private(i,m,hup,hlos,ts2)
	    for (i=X0;i<=nx;i++) {
	      if (uhr[k][i][j] == 0.0) {
		uhh[i] = 0.0;
		for (m=0;m<NTR;m++) fluxtr[i][m] = 0.0;
	      }
	      else if (uhr[k][i][j] < 0.0) {

		if (k==0 || k==1 ) {
		  hup = (hvol[k][i+1][j]-DXDYh(i+1,j)*MLMIN);
                } else {
                  hup = (hvol[k][i+1][j]-DXDYh(i+1,j)*EPSILON);
		}

		hlos = D_MAX(0.0,uhr[k][i+1][j]);
		if (((hup + uhr[k][i][j] - hlos) < 0.0) && 
		    ((0.5*hup + uhr[k][i][j]) < 0.0)) {
		  uhh[i] = D_MIN(-0.5*hup,-hup+hlos);
		}
		else uhh[i] = uhr[k][i][j];
		ts2 = 0.5*(1.0 + uhh[i]/hvol[k][i+1][j]);
		for (m=0;m<NTR;m++) {
		  fluxtr[i][m] = uhh[i]*(tr[m][k][i+1][j] - slope[i+1][m]*ts2);
		}
	      }
	      else {

                if (k==0 || k==1 ) {
                  hup = (hvol[k][i][j]-DXDYh(i,j)*MLMIN);
                } else {
                  hup = (hvol[k][i][j]-DXDYh(i,j)*EPSILON);
                }

		hlos = D_MAX(0.0,-uhr[k][i-1][j]);
		if (((hup - uhr[k][i][j] - hlos) < 0.0) && 
		    ((0.5*hup - uhr[k][i][j]) < 0.0)) {
		  uhh[i] = D_MAX(0.5*hup,hup-hlos);
		}
		else uhh[i] = uhr[k][i][j];
		ts2 = 0.5*(1.0 - uhh[i]/hvol[k][i][j]);

		for (m=0;m<NTR;m++) {
		  fluxtr[i][m] = uhh[i]*(tr[m][k][i][j] + slope[i][m]*ts2);
		}
	      }
	    }
            //#pragma omp barrier
/*   Calculate new tracer concentration in each cell after accounting */
/* for the i-direction fluxes.                                        */

	    uhr[k][X0][j] -= uhh[X0];
           // #pragma omp barrier

//#pragma omp for  private(i,m,hlst1,Ihnew)
	    for (i=X1;i<=nx;i++) {

	      if ((uhh[i] != 0.0) || (uhh[i-1] != 0.0)) 
		{
		  uhr[k][i][j] -= uhh[i];
		  hlst1 = hvol[k][i][j];

		  hvol[k][i][j] -= (uhh[i] - uhh[i-1]);
		  Ihnew = 1.0 / hvol[k][i][j];
		  
		  for (m=0;m<NTR;m++) {
		    tr[m][k][i][j] *= hlst1;
		    tr[m][k][i][j] = (tr[m][k][i][j] - 
				      (fluxtr[i][m]-fluxtr[i-1][m])) * Ihnew;
		  }

		}
	    }
          //  #pragma omp barrier
	  } /* j loop */
#endif

/* ============================================================ */
/*			now advect meridionally			*/
/* ============================================================ */
#ifndef ADV1D
	  for (i=X1;i<=nx;i++) {
/*   Calculate the j-direction profiles (slopes) of each tracer that  */
/* is being advected.                                                 */
//#pragma omp for  private(j,m,minslope)
	    for (j=Y0;j<=ny+1;j++) {
	      for (m=0;m<NTR;m++) {
		minslope = 4.0*((fabs(tr[m][k][i][j+1]-tr[m][k][i][j]) <
				 fabs(tr[m][k][i][j]-tr[m][k][i][j-1])) ?
				(tr[m][k][i][j+1]-tr[m][k][i][j]) : 
				(tr[m][k][i][j]-tr[m][k][i][j-1]));
		slope[j][m] = vmask[i][j] * vmask[i][j-1] *
		  (((tr[m][k][i][j+1]-tr[m][k][i][j]) *
		    (tr[m][k][i][j]-tr[m][k][i][j-1]) < 0.0) ? 0.0 :
		   ((fabs(tr[m][k][i][j+1]-tr[m][k][i][j-1])<fabs(minslope)) ?
		    0.5*(tr[m][k][i][j+1]-tr[m][k][i][j-1]) : 0.5*minslope));
	      }
	    }
        //    #pragma omp barrier
  
/*   Calculate the j-direction fluxes of each tracer, using as much   */
/* the minimum of the remaining mass flux (vhr) and the half the mass */
/* in the cell plus whatever part of its half of the mass flux that   */
/* the flux through the other side does not require.                  */
//#pragma omp for  private(j,m,hup,hlos,ts2)
	    for (j=Y0;j<=ny;j++) {
	      if (vhr[k][i][j] == 0.0) { 
		vhh[j] = 0.0;
		for (m=0;m<NTR;m++) fluxtr[j][m] = 0.0;
	      }
	      else if (vhr[k][i][j] < 0.0) {

                if (k==0 || k==1 ) {
                  hup = (hvol[k][i][j+1]-DXDYh(i,j+1)*MLMIN);
                } else {
                  hup = (hvol[k][i][j+1]-DXDYh(i,j+1)*EPSILON);
                }

		hlos = D_MAX(0.0,vhr[k][i][j+1]);
		
		if (((hup + vhr[k][i][j] - hlos) < 0.0) && 
		    ((0.5*hup + vhr[k][i][j]) < 0.0)) {
		  vhh[j] = D_MIN(-0.5*hup,-hup+hlos);
		}
		
		else vhh[j] = vhr[k][i][j];
		ts2 = 0.5*(1.0 + vhh[j]/(hvol[k][i][j+1]));
		
		for (m=0;m<NTR;m++) {
		  fluxtr[j][m] = vhh[j]*(tr[m][k][i][j+1] - slope[j+1][m]*ts2);
		}
	      }
	      else {

                if (k==0 || k==1 ) {
                  hup = (hvol[k][i][j]-DXDYh(i,j)*MLMIN);
                } else {
                  hup = (hvol[k][i][j]-DXDYh(i,j)*EPSILON);
                }

		hlos = D_MAX(0.0,-vhr[k][i][j-1]);
		
		if (((hup - vhr[k][i][j] - hlos) < 0.0) 
		    && ((0.5*hup - vhr[k][i][j]) < 0.0)) {
		  vhh[j] = D_MAX(0.5*hup,hup-hlos);
		}
		
		else vhh[j] = vhr[k][i][j];
		ts2 = 0.5*(1.0 - vhh[j] / (hvol[k][i][j]));
		
		for (m=0;m<NTR;m++) {
		  fluxtr[j][m] = vhh[j]*(tr[m][k][i][j] + slope[j][m]*ts2);
		}
	      }
	    }
       //     #pragma omp barrier

/*   Calculate new tracer concentration in each cell after accounting */
/* for the j-direction fluxes.                                        */

	    vhr[k][i][Y0] -= vhh[Y0];
         //  #pragma omp barrier

//#pragma omp for private(j,m,Ihnew)
	    for (j=Y1;j<=ny;j++) {
	      if ((vhh[j] != 0.0) || (vhh[j-1] != 0.0)) {
		hlst[j] = hvol[k][i][j];
		hvol[k][i][j] -= (vhh[j] - vhh[j-1]);
		Ihnew = 1.0 / hvol[k][i][j];
		vhr[k][i][j] -= vhh[j];
		for (m=0;m<NTR;m++) {
		  tr[m][k][i][j] *= hlst[j];
		  tr[m][k][i][j] = (tr[m][k][i][j] - 
				    fluxtr[j][m] + fluxtr[j-1][m]) * Ihnew;
		}
	      }
	    }
       //     #pragma omp barrier
	  } /* i loop */
#endif

	}			 /* end of big loop over k		*/

/*	calculate new thickness field - to be used for vertical 	*/
/*	tracer advection from updated volumes (vol + fluxes)		*/

#pragma omp for  private(i,j,k)
      for (k=0; k<=NZ-1; k++) {
	for (i=X1; i<=nx; i++) {
          for (j=Y1; j<=ny; j++) {
	    hnew[k][i][j] = hvol[k][i][j]/DXDYh(i,j);
	    
	    if (hnew[k][i][j] < EPSILON) {	
	      hnew[k][i][j] = EPSILON;     
	    }
	  }
	}
      }


/* ============================================================ */
/*			now advect vertically			*/
/* ============================================================ */
#pragma omp for private(i,j,hup,hlos)
      for (j=Y1; j<=ny; j++) {
	 for (i=X1; i<=nx; i++) {
/*      work from top to bottom - by interfaces - interface k is the    */
/*      interface between layer k and layer k-1. net flux at this       */
/*      interface is wd[[k][i][j]= ea[k][i][j] and eb[k-1][i][j]        */
        
/* k=0 */
        
	  if (wd[0][i][j] == 0.0) {
	    wdh[0][i][j] = 0.0;
	  }
	  else if (wd[0][i][j] < 0.0) {
	    hup = hnew[0][i][j] - MLMIN;
	    hlos = D_MAX(0.0, wd[1][i][j]);
	    if (((hup + wd[0][i][j] - hlos) < 0.0) &&
		((0.5*hup + wd[0][i][j]) < 0.0)) {
	      wdh[0][i][j] = D_MIN(-0.5*hup,-hup+hlos);
	    }
	    else wdh[0][i][j] = wd[0][i][j];
	  }
	  else {
	    wdh[0][i][j] = wd[0][i][j];
	  }        
    }
    }

#pragma omp for  private(i,j,hup,hlos)
      for (j=Y1; j<=ny; j++) {
	 for (i=X1; i<=nx; i++) {
/* k=1 */

          if (wd[1][i][j] == 0.0) {
            wdh[1][i][j] = 0.0;
          }
          else if (wd[1][i][j] > 0.0) {
            hup = hnew[0][i][j] - MLMIN;
            hlos = D_MAX(0.0, -wd[0][i][j]);
            if (((hup - wd[1][i][j] - hlos) < 0.0) &&
                ((0.5*hup - wd[1][i][j]) < 0.0)) {
              wdh[1][i][j] = D_MAX(0.5*hup,hup-hlos);
            }
            else wdh[1][i][j] = wd[1][i][j];
          }
          else {
            hup = hnew[1][i][j] - MLMIN;
            hlos = D_MAX(0.0,wd[2][i][j]);
            if (((hup + wd[1][i][j] - hlos) < 0.0) &&
                ((0.5*hup + wd[1][i][j]) < 0.0)) {
              wdh[1][i][j] = D_MIN(-0.5*hup,-hup+hlos);
            }
            else wdh[1][i][j] = wd[1][i][j];
          }
      }
       }

#pragma omp for  private(i,j,hup,hlos)
      for (j=Y1; j<=ny; j++) {
	 for (i=X1; i<=nx; i++) {
/* k=2 */

          if (wd[2][i][j] == 0.0) {
            wdh[2][i][j] = 0.0;
          }
          else if (wd[2][i][j] > 0.0) {
            hup = hnew[1][i][j] - MLMIN;
            hlos = D_MAX(0.0, -wd[1][i][j]);
            if (((hup - wd[2][i][j] - hlos) < 0.0) &&
                ((0.5*hup - wd[2][i][j]) < 0.0)) {
              wdh[2][i][j] = D_MAX(0.5*hup,hup-hlos);
            }
            else wdh[2][i][j] = wd[2][i][j];
          }
          else {
            hup = hnew[2][i][j] - EPSILON;
            hlos = D_MAX(0.0,wd[3][i][j]);
            if (((hup + wd[2][i][j] - hlos) < 0.0) &&
                ((0.5*hup + wd[2][i][j]) < 0.0)) {
              wdh[2][i][j] = D_MIN(-0.5*hup,-hup+hlos);
            }
            else wdh[2][i][j] = wd[2][i][j];
          }
        }
         }

/* k=3 --> NZ-1 */

#pragma omp for  private(i,j,k,hup,hlos)
	  for (k=3; k<=NZ-1; k++) {   	
      for (j=Y1; j<=ny; j++) {
	 for (i=X1; i<=nx; i++) {

	    if (wd[k][i][j] == 0.0) {
	      wdh[k][i][j] = 0.0;
	    }
	    else if (wd[k][i][j] > 0.0) {
	      hup = hnew[k-1][i][j] - EPSILON;
	      hlos = D_MAX(0.0, -wd[k-1][i][j]);
	      if (((hup - wd[k][i][j] - hlos) < 0.0) &&
		  ((0.5*hup - wd[k][i][j]) < 0.0)) {
		wdh[k][i][j] = D_MAX(0.5*hup,hup-hlos);
	      }
	      else wdh[k][i][j] = wd[k][i][j];
	    }
	    else {
	      hup = hnew[k][i][j] - EPSILON;
	      if (k != NZ-1) {
		hlos = D_MAX(0.0,wd[k+1][i][j]);
	      } else {
		hlos = 0.0;
	      }
	      if (((hup + wd[k][i][j] - hlos) < 0.0) &&
		  ((0.5*hup + wd[k][i][j]) < 0.0)) {
		wdh[k][i][j] = D_MIN(-0.5*hup,-hup+hlos);
	      }
	      else wdh[k][i][j] = wd[k][i][j];
	    }
	    
	  } /* k */
	}   /* j */
      }     /* i */

#pragma omp for  private(i,j)
      for (i=X1; i<=nx; i++)
	  for (j=Y1; j<=ny; j++) {
	      ear[0][i][j] = wdh[0][i][j];
	      /* added by Curtis - bottom ebr wasn't set anywhere else */
	      ebr[NZ-1][i][j] = 0;  
	  }

#pragma omp for  private(i,j,k)
      for (k=1;k<=NZ-1;k++) { 
	  for (i=X1; i<=nx; i++) {
	      for (j=Y1; j<=ny; j++) {  
		  ear[k][i][j] =   0.5 * (fabs(wdh[k][i][j]) + wdh[k][i][j]);
		  ebr[k-1][i][j] = 0.5 * (fabs(wdh[k][i][j]) - wdh[k][i][j]);
	      }
	  }
      }     
 
#pragma omp for  private(i,j,k,m,hnew0,bet,gam)
      for (j=Y1; j<=ny; j++) {

	  for (i=X1; i<=nx; i++) {
	      hnew0[i] = hnew[0][i][j];
	      bet[i]=1.0/(hnew[0][i][j] + ebr[0][i][j] + wdh[0][i][j]);

	      for (m=0;m<NTR;m++)
		  tr[m][0][i][j] = bet[i]*(hnew0[i]*tr[m][0][i][j]);
	  }

	  for (k=1;k<=NZ-1;k++) {
	      for (i=X1;i<=nx;i++) {
		  gam[k][i] = ebr[k-1][i][j] * bet[i];

		  bet[i]=1.0/(hnew[k][i][j] + ebr[k][i][j] +
			      (1.0-gam[k][i])*ear[k][i][j]);
		  

		  for (m=0;m<NTR;m++)
		      tr[m][k][i][j] = bet[i] * (hnew[k][i][j]*tr[m][k][i][j] +
						 ear[k][i][j]*(tr[m][k-1][i][j]) );
	      }	      
	  }

	  for (m=0;m<NTR;m++)
	      for (k=NZ-2;k>=0;k--) {
		  for (i=X1;i<=nx;i++) {
		      tr[m][k][i][j] += gam[k+1][i]*tr[m][k+1][i][j];
		  }
       }
      } /*j*/

/* update hvol with diapycnal fluxes */
#pragma omp for  private(i,j,k)
      for (k=0;k<NZ-1;k++) {
	  for (i=X1; i<=nx; i++)
	      for (j=Y1; j<=ny; j++)
		  hnew[k][i][j] += (wdh[k][i][j] - wdh[k+1][i][j]);
      }

#pragma omp for  private(i,j)
      for (i=X1; i<=nx; i++)
	  for (j=Y1; j<=ny; j++)
	      hnew[NZ-1][i][j] += wdh[NZ-1][i][j];

#pragma omp for  private(i,j,k)
      for (k=0;k<=NZ-1;k++)
	  for (i=X1; i<=nx; i++)
	      for (j=Y1; j<=ny; j++) {
		  if (hnew[k][i][j] < EPSILON) hnew[k][i][j] = EPSILON;
		  hvol[k][i][j] = DXDYh(i,j)*hnew[k][i][j];

		  if ( wd[k][i][j] > 0.0 && ( wdh[k][i][j] > wd[k][i][j] ))
		      printf("case 1 wdh[k]\n");

		  else if ( wd[k][i][j] < 0.0 && ( wdh[k][i][j] < wd[k][i][j] ))
		      printf("case 2 wdh[k]\n");
		  wd[k][i][j] -= wdh[k][i][j];	
	      }
 
#else  /* STANDARD_ADVECTION */
      /* big loop over k	 */
//      printf("phos(%d,%d,%d)=%g,uhtm=%g\n",0,190,26,tr[mPHOSPHATE][0][190][26],
//      	     uhtm[0][190][26]);
//      exit(1);

//yanxu: note these are null cycles, so I comment them out
//yanxu      for (k=0;k<NZ;k++) {
	// first advect zonally
//yanxu	for (j=Y1;j<=ny;j++) {
//yanxu	  for (i=X0;i<=nx+1;i++) {
//yanxu	    for (m=0;m<NTR;m++) {
	      //	      fluxtr[i][m] = uhtm[i]*(tr[m][k][i][j]);
//yanxu	    }
//yanxu	  }
//yanxu	}

	// now advect meridionally
// null cycles again, yanxu
//yanxu	for (i=X1;i<=nx;i++) {

//yanxu	}
//yanxu      } /* end of big loop over k */

      /*  calculate new thickness field - to be used for vertical  */
      /*  tracer advection from updated volumes (vol + fluxes)     */
#pragma omp for  private(i,j,k)
      for (k=0; k<=NZ-1; k++) 
	for (i=X1; i<=nx; i++) 
	  for (j=Y1; j<=ny; j++) {
	    hnew[k][i][j] = hvol[k][i][j]/DXDYh(i,j);
	    if (hnew[k][i][j] < EPSILON) hnew[k][i][j] = EPSILON;     
	  }

	  
      // now advect vertically
//yanxu: null cycles
//yanxu      for (j=Y1; j<=ny; j++) {
//yanxu	for (i=X1; i<=nx; i++) {
//yanxu	  
//yanxu	}
//yanxu      }


#endif /* STANDARD_ADVECTION */


      /* zonal re-entrance */
#pragma omp for  private(j,k,m)
    for (j=0;j<NYMEM;j++) {
      for (k=0;k<NZ;k++) {
    
        uhr[k][nx+1][j] = uhr[k][2][j];
        uhr[k][nx+2][j] = uhr[k][3][j];
        uhr[k][0][j]    = uhr[k][nx-1][j];
        uhr[k][1][j]    = uhr[k][nx][j];

        vhr[k][nx+1][j] = vhr[k][2][j];
        vhr[k][nx+2][j] = vhr[k][3][j];
        vhr[k][0][j]    = vhr[k][nx-1][j];
        vhr[k][1][j]    = vhr[k][nx][j];

        hvol[k][nx+1][j] = hvol[k][2][j];
        hvol[k][nx+2][j] = hvol[k][3][j];
        hvol[k][0][j]    = hvol[k][nx-1][j];
        hvol[k][1][j]    = hvol[k][nx][j];

        for (m=0;m<NTR;m++) {
          tr[m][k][nx+1][j] = tr[m][k][2][j];
          tr[m][k][nx+2][j] = tr[m][k][3][j];
          tr[m][k][0][j]    = tr[m][k][nx-1][j];
          tr[m][k][1][j]    = tr[m][k][nx][j];
        }

      }
    }

	/* meridional re-entrance            */

	/* meridional re-entrance            */
#pragma omp for  private(i,k,ii,m)
    for (i=2;i<=nx;i++) {
      ii = 363 - i;
      for (k=0;k<NZ;k++) {
        uhr[k][ii][ny+1] = (-1)*uhr[k][i][ny];
        uhr[k][ii][ny+2] = (-1)*uhr[k][i][ny-1];

        vhr[k][ii][ny+1] = (-1)*vhr[k][i][ny];
        vhr[k][ii][ny+2] = (-1)*vhr[k][i][ny-1];

        hvol[k][ii][ny+1] = hvol[k][i][ny];
        hvol[k][ii][ny+2] = hvol[k][i][ny-1];
  
        for (m=0;m<NTR;m++) {
          tr[m][k][ii][ny+1] = tr[m][k][i][ny];
          tr[m][k][ii][ny+2] = tr[m][k][i][ny-1];
        }
      }
    }

} // omp

#ifdef STANDARD_ADVECTION
# ifdef WRTTS
      printf("itt = %i\n",itt);
      wrts = (double)(itts)+(double)(itt)/11.+0.0001;
      write_ts(wrts);
# endif   // WRTTS

//****************************************************************************************
    }  /* end of temp itt iteration loop */
//****************************************************************************************
#endif
    fract1 = (double)(NTSTEP-itts) / (double)NTSTEP;
    fract2 = 1.0 - fract1;

//#pragma omp parallel for  private(i,j,k) schedule(dynamic)
    for (k=0;k<=NZ-1;k++) {
	for (i=X1; i<=nx; i++) {
	    for (j=Y1; j<=ny; j++) {
# ifdef USE_CALC_H
		h[k][i][j] = hnew[k][i][j];
		if (h[k][i][j] < 0.0)
		  printf("tracadv l 796 - h[%d][%d][%d] = %g\n", k,i,j,
			 h[k][i][j]);	 
# else
		h[k][i][j] = fract1*hstart[k][i][j] + fract2*hend[k][i][j];
//BX 		h[k][i][j] = hend[k][i][j];
# endif
#ifdef HTEST
	  	htest[k][i][j] = hnew[k][i][j]; 
		printf("htest(%d,%d,%d)=%g,hend=%g\n",
		       k,i,j,
//		htest[k][i][j] = h[k][i][j];
#endif
	    }
	}
    }


    //HF
    //	zonal re-entrance
//#pragma omp parallel for  private(j,k) schedule(dynamic)
    for (k=0;k<NZ;k++) {
      for (j=0;j<=NYMEM-1;j++) {
	h[k][nx+1][j] = h[k][2][j];
	h[k][nx+2][j] = h[k][3][j];
	h[k][0][j] =   h[k][nx-1][j];
	h[k][1][j] =   h[k][nx][j];
      }
    }


    //      meridional re-entrance
//#pragma omp parallel for  private(i,k,ii) schedule(dynamic)
    for (i=2;i<=nx;i++) {
      ii = 363 - i;
      for (k=0;k<NZ;k++) {
        h[k][ii][ny+1] = h[k][i][ny];
        h[k][ii][ny+2]   = h[k][i][ny-1];
      }
    }
    //HF-e

# ifdef WRTTS
      printf("End of tracadv\n");
      wrts = (double)(itts)+(double)(itt)/11.+0.0005;
      write_ts(wrts);
# endif   // WRTTS

#ifdef DIFFUSE_TRACER
  if ((KD>0.0) || (KDML>0.0)) {
    diffuse_tracer();        
     }
#endif

     
  pstage=4;
  print_tr(pstage);

  if ((KHTR>0.0)) {
    tracer_hordiff();
  }
 
  pstage=5;
  print_tr(pstage);

  free3d(hvol, NZ);
  free3d(uhr, NZ);
  free3d(vhr, NZ);
  free3d(ear, NZ);
  free3d(ebr, NZ);
  free3d(wdh, NZ);
#if !defined AGE2 && !defined AGE3
  free3d(hnew, NZ);
#endif
}

#ifdef DIFFUSE_TRACER
static void diffuse_tracer()
{
#if (NZ>1)
/*    This subroutine does a fully implicit vertical diffusion        */
/*  of tracer.  Insulating top and bottom b.c.s are used.             */

  double a[NZ+1][NYMEM];   /* a is the coupling coefficient across an */
                           /* interface time integrated over dt, in m.*/
                           /* a times the tracer difference gives the */
                           /* tracer flux across an interface.        */
  double bet[NYMEM];       /* bet and gam are variables used by the   */
  double gam[NZ][NYMEM];   /* tridiagonal solver.                     */

  double z[NYMEM];         /* The distance from the top, normalized   */
                           /* by HMIX, nondimensional.                */
  double I_HMIX = 1.0/HMIX;
  double C2dtKD, C2dtKDML;
  int i, j, k, m;

  C2dtKD = 2.0*KD*dt; C2dtKDML = 2.0*KDML*dt;

	printf("DIFFUSING TRACER - diapycnal\n");

//#pragma omp parallel for private(i,j,k,m,a,z,bet,gam)
  for (i=X1;i<=nx;i++) {
/*    The following loops calculates the diffusive coupling between   */
/*  layers.  Without BULKMIXEDLAYER, elevated diffusion is assumed    */
/*  within HMIX of the surface.                                       */
/*  This is tridag, from Numerical Recipes in C.                      */
    for (j=Y1;j<=ny;j++) {
      if (D[i][j] > MINIMUM_DEPTH) {
#ifdef BULKMIXEDLAYER
      a[1][j] = -C2dtKD / (h[1][i][j] + h[0][i][j]);
#else
      z[j] = h[0][i][j]*I_HMIX;
      a[1][j] = -C2dtKD - C2dtKDML / 
                ((h[1][i][j] + h[0][i][j]) *
                 (1.0 + 0.09*z[j]*z[j]*z[j]*z[j]*z[j]*z[j]) );
#endif
      bet[j]=1.0/(h[0][i][j]-a[1][j]);

      for (m=0;m<NTR;m++)
        tr[m][0][i][j] = bet[j] * h[0][i][j] * tr[m][0][i][j];
    }
    }

    for (k=1;k<=NZ-2;k++) {
      for (j=Y1;j<=ny;j++) {
         if (D[i][j] > MINIMUM_DEPTH) {
#ifdef BULKMIXEDLAYER
        a[k+1][j] = -C2dtKD / (h[k+1][i][j] + h[k][i][j]);
#else
        z[j] += h[k][i][j]*I_HMIX;
        a[k+1][j] = -C2dtKD - C2dtKDML / 
                  ((h[k+1][i][j] + h[k][i][j]) *
                   (1.0 + 0.09*z[j]*z[j]*z[j]*z[j]*z[j]*z[j]) );
/* (           (1.0 + 0.05*z[j]*z[j]*z[j]*z[j]*z[j]*z[j]*z[j]*z[j])); */
#endif
        gam[k][j]=a[k][j]*bet[j];
        bet[j]=1.0/(h[k][i][j] - a[k+1][j] - a[k][j]*(1.0+gam[k][j]));

        for (m=0;m<NTR;m++)
          tr[m][k][i][j]=(h[k][i][j]*tr[m][k][i][j] - 
                          a[k][j]*tr[m][k-1][i][j])*bet[j];
      }
     }
    }
    
    for (j=Y1;j<=ny;j++) {
      if (D[i][j] > MINIMUM_DEPTH) {
      gam[NZ-1][j]=a[NZ-1][j]*bet[j];
      bet[j]=1.0/(h[NZ-1][i][j] - a[NZ-1][j]*(1.0+gam[k][j]));

      for (m=0;m<NTR;m++)
        tr[m][NZ-1][i][j]=(h[NZ-1][i][j]*tr[m][NZ-1][i][j] - 
                           a[NZ-1][j]*tr[m][NZ-2][i][j])*bet[j];
      }
    }

    for (k=(NZ-2);k>=0;k--) {
      for (j=Y1;j<=ny;j++) {
        if (D[i][j] > MINIMUM_DEPTH)
        for (m=0;m<NTR;m++)
          tr[m][k][i][j] -= gam[k+1][j]*tr[m][k+1][i][j];
      }
    }

  }
#endif /* NZ>1 */

}
#endif /* DIFFUSE_TRACER */


static void tracer_hordiff(void)
{

/*   This subroutine does along-coordinate diffusion of all tracers,  */
/*  using the diffusivity KHTR.  Multiple iterations are  used (if    */
/*  necessary) so that there is no limit on the acceptable time       */
/*  increment.                                                        */

  double khdt_rem_x[NXMEM][NYMEM]; /* The value of KHTR*dt which must */
                                /* be accomodated by later iterations */
                                /* in the zonal direction, in m2.     */
  double khdt_rem_y[NXMEM][NYMEM]; /* The value of KHTR*dt which must */
                                /* be accomodated by later iterations */
                                /* in the meridional direction, in m2.*/
  double Coef_x[NXMEM][NYMEM];  /* The coefficient relating zonal     */
                                /* tracer differences to fluxes, m3.  */
  double Coef_y[NXMEM][NYMEM];  /* The coefficient relating meridional*/
                                /* tracer differences to fluxes, m3.  */
  double Coef_x0[NXMEM][NYMEM];  /* The coefficient relating zonal     */
                                /* tracer differences to fluxes, m3.  */
  double Coef_y0[NXMEM][NYMEM];  /* The coefficient relating meridional*/
                                /* tracer differences to fluxes, m3.  */
  double Ihdxdy[NXMEM][NYMEM];  /* The inverse of the volume of fluid */
                                /* in a layer in a grid cell, m-3.    */
  static double Max_khdt_x[NXMEM][NYMEM]; /* The maximum value of     */
                                /* KHTR*dt which can be stably accom- */
                                /* odated within a single iteration in*/
                                /* the zonal-direction, in m2.        */
  static double Max_khdt_y[NXMEM][NYMEM]; /* The maximum value of     */
                                /* KHTR*dt which can be stably accom- */
                                /* odated within a single iteration in*/
                                /* the meridional-direction, in m2.   */
  static double khdt_rem_x0[NXMEM][NYMEM]; /* The value of KHTR*dt which must */
                                /* be accomodated by later iterations */
                                /* in the zonal direction, in m2.     */
  static double khdt_rem_y0[NXMEM][NYMEM]; /* The value of KHTR*dt which must */
                                /* be accomodated by later iterations */
                                /* in the meridional direction, in m2.*/

  extern double umask[NXMEM][NYMEM];   /* _mask are 1 over ocean and 0  */
  extern double vmask[NXMEM][NYMEM];   /* over land on the u & v grids. */

  static int zero_if_first_call = 0;
  int i, j, k, ii, m, domore_k, itt;
  double C0;  /* A work variable with units of m2. */

  printf("DIFFUSING TRACER - isopycnal\n");

  if (zero_if_first_call == 0) {
    for (j=Y1;j<=ny;j++) 
      for (i=X0;i<=nx;i++) 
	      Max_khdt_x[i][j] = 0.125 / (DYu(i,j) * IDXu(i,j) *
           ((IDXDYh(i,j)>IDXDYh(i+1,j)) ? IDXDYh(i,j) : IDXDYh(i+1,j)));

      for (i=X1;i<=nx;i++) 
	    for (j=Y0;j<=ny;j++) 
	      Max_khdt_y[i][j] = 0.125 / (DXv(i,j) * IDYv(i,j) *
	      ((IDXDYh(i,j)>IDXDYh(i,j+1)) ? IDXDYh(i,j) : IDXDYh(i,j+1)));

      for (i=X0;i<=nx;i++) 
	    for (j=Y1;j<=ny;j++) 
	      khdt_rem_x0[i][j] = dt*KHTR*DYu(i,j)*IDXu(i,j)*umask[i][j];

      for (i=X1;i<=nx;i++) 
	    for (j=Y0;j<=ny;j++) 
	      khdt_rem_y0[i][j] = dt*KHTR*DXv(i,j)*IDYv(i,j)*vmask[i][j];
     
      zero_if_first_call = 1;
  }

//#pragma omp parallel for private(k,i,j,khdt_rem_x,khdt_rem_y,C0,Coef_y,Coef_x,Ihdxdy,m,ii,itt,domore_k)
for (k=0;k<NZ;k++) {

      for (i=X0;i<=nx;i++) 
	    for (j=Y1;j<=ny;j++) 
	      khdt_rem_x[i][j] = khdt_rem_x0[i][j];

      for (i=X1;i<=nx;i++) 
	    for (j=Y0;j<=ny;j++) 
	      khdt_rem_y[i][j] = khdt_rem_y0[i][j];

      for (i=X1;i<=nx;i++) 
	    for (j=Y0;j<=ny;j++) 
         Coef_y0[i][j] = 2.0*h[k][i][j]*h[k][i][j+1] / (h[k][i][j]+h[k][i][j+1]);

      for (j=Y1;j<=ny;j++)
        for (i=X0;i<=nx;i++)
          Coef_x0[i][j] = 2.0*h[k][i][j]*h[k][i+1][j] / (h[k][i][j]+h[k][i+1][j]);

      for (j=Y1;j<=ny;j++)
        for (i=X1;i<=nx;i++) 
          Ihdxdy[i][j] = IDXDYh(i,j) / h[k][i][j];

    itt = 0; domore_k = 1;
    while (domore_k) {
      domore_k = 0; itt++;

      for (i=X1;i<=nx;i++) {
	  for (j=Y0;j<=ny;j++) {
          if (khdt_rem_y[i][j] < Max_khdt_y[i][j]) C0 = khdt_rem_y[i][j];
          else {C0 = Max_khdt_y[i][j]; domore_k++;}
          khdt_rem_y[i][j] -= C0;
         Coef_y[i][j] = C0*Coef_y0[i][j];
        }
      }

      for (j=Y1;j<=ny;j++) {
        for (i=X0;i<=nx;i++) {
          if (khdt_rem_x[i][j] < Max_khdt_x[i][j]) C0 = khdt_rem_x[i][j];
          else {C0 = Max_khdt_x[i][j]; domore_k++;}
          khdt_rem_x[i][j] -= C0;
          Coef_x[i][j] = C0*Coef_x0[i][j];
        }
      }

	  for (j=0;j<NYMEM;j++) {
	      Coef_x[nx+1][j] = Coef_x[2][j];
	      Coef_x[nx+2][j] = Coef_x[3][j];
	      Coef_x[0][j]    = Coef_x[nx-1][j];
	      Coef_x[1][j]    = Coef_x[nx][j];
	      Coef_y[nx+1][j] = Coef_y[2][j];
	      Coef_y[nx+2][j] = Coef_y[3][j];
	      Coef_y[0][j]    = Coef_y[nx-1][j];
	      Coef_y[1][j]    = Coef_y[nx][j];
	  }

	  for (i=2;i<=nx;i++) {
	      ii = 363 - i;
	      Coef_x[ii][ny+1] = Coef_x[i][ny];
	      Coef_x[ii][ny+2] = Coef_x[i][ny-1];
	      Coef_y[ii][ny+1] = Coef_y[i][ny];
	      Coef_y[ii][ny+2] = Coef_y[i][ny-1];
	  }

      for (m=0;m<NTR;m++) {

	  for (i=X1;i<=nx;i++) 
	      for (j=Y1;j<=ny;j++) 
		  tr[m][k][i][j] += Ihdxdy[i][j] *
		      (Coef_x[i-1][j] * (tr[m][k][i-1][j] - tr[m][k][i][j]) -
		       Coef_x[i][j] *   (tr[m][k][i][j]   - tr[m][k][i+1][j]) +
		       Coef_y[i][j-1] * (tr[m][k][i][j-1] - tr[m][k][i][j]) -
		       Coef_y[i][j] *   (tr[m][k][i][j]   - tr[m][k][i][j+1]));
	      
	  

//BX-a
//BX   add re-entrace values here
		/* zonal re-entrance		*/
	  for (j=0;j<NYMEM;j++) {
	      tr[m][k][nx+1][j] = tr[m][k][2][j];
	      tr[m][k][nx+2][j] = tr[m][k][3][j];
	      tr[m][k][0][j]    = tr[m][k][nx-1][j];
	      tr[m][k][1][j]    = tr[m][k][nx][j];
	  }
	  /* meridional re-entrance            */
	  for (i=2;i<=nx;i++) {
	      ii = 363 - i;
	      tr[m][k][ii][ny+1] = tr[m][k][i][ny];
	      tr[m][k][ii][ny+2] = tr[m][k][i][ny-1];
	  }
//BX-e
//BX-e
      }   // m

      if (itt > 10) {
	        domore_k = 0;       /* break loop */
      }

    }  /* end while */
} /* end k loop */  
} /* end of subroutine */

void print_tr(int pstage){

  int printwarning=1;
  
  if (printwarning == 1) {
//    printf("tracadv.c: BEGIN PRINT STAGE # %i \n", pstage);
  }
  
}

//BX-a
//BX for debugging
//BX void hvol_integral(void){
void hvol_integral(double ***hvol){
    
  int i, j, k;

  hvolint=0.e0;
//#pragma omp parallel for  private(i,j,k) reduction(+:hvolint) schedule(dynamic)
  for (k=0;k<NZ;k++) {
      for (i=2;i<=NXMEM-1;i++) {
	  for (j=2;j<=NYMEM-1;j++) {
	      if (hvol[k][i][j] > 1e-6) {
		  //BX		  hvolint += h[k][i][j] * dxdyh[i][j];
		  hvolint += hvol[k][i][j];
	      }
	  }
      }
  }

}
void hvol_kintegral(int k, double ***hvol){
    
  int i, j;

  hvolint=0.e0;
//#pragma omp parallel for  private(i,j) reduction(+:hvolint) schedule(dynamic)
  for (i=2;i<=NXMEM-1;i++) {
      for (j=2;j<=NYMEM-1;j++) {
	  //BX if (h[k][i][j] > 1e-6) {
	      //BX	      hvolint += h[k][i][j] * dxdyh[i][j];
	      hvolint += hvol[k][i][j];
	      //BX }
      }
  }

}

//BX tracer_integral is calculated in step.c

void tracer_kintegral(int trnum, int k, double ***hvol){
    
  int i, j;
  double sum;

//  trintegral[trnum]=0.e0;
  sum=0.e0;
//#pragma omp parallel for  private(i,j) reduction(+:sum) schedule(dynamic)
  for (i=2;i<=NXMEM-1;i++) {
      for (j=2;j<=NYMEM-1;j++) {
	  //BX if (h[k][i][j] > 1e-6) {
	      //BX	      trintegral[trnum] += tr[trnum][k][i][j] * h[k][i][j] * dxdyh[i][j];
          //   trintegral[trnum] += tr[trnum][k][i][j] * hvol[k][i][j];
          sum += tr[trnum][k][i][j] * hvol[k][i][j];
	     //BX }
      }
  }

 trintegral[trnum] = sum;

}
//BX-e

