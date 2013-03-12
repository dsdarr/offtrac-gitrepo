#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "init.h"
#include "par_IO.h"

#define NOFLD (NTR+1)

static double dmina1, dmina2;
# define D_MIN(a,b) (dmina1=(a),dmina2=(b),((dmina1) < (dmina2)) ?\
                    (dmina1) : (dmina2))
# define D_MAX(a,b) (dmina1=(a),dmina2=(b),((dmina1) > (dmina2)) ?\
                    (dmina1) : (dmina2))

extern double D[NXMEM][NYMEM];       /* Basin depth, in m.            */
extern double ****tr;

static int num_col;               /* The number of sponge points      */
                                  /* within the computational domain. */
static int fldno = 0;             /* The number of fields which have  */
                                  /* already been registered by calls */
                                  /* to set_up_sponge_field           */
static int nf;                    /* The number of fields for which   */
                                  /* space has been allocated.        */
static int nfNZ;                  /* nfNZ = nf * NZ.                  */
static int *col_addr;             /* An array containing the addresses*/
                                  /* of each of the columns being     */
                                  /* damped in an array of horizontal */
                                  /* extent NXMEM by NYMEM.           */
static double *Iresttime_col;     /* The inverse restoring time of    */
                                  /* each column.                     */
static double (*var[NOFLD])[][NXMEM][NYMEM]; /* An array of pointers  */
                                  /* to the fields which are being    */
                                  /* damped.                          */
static double *Ref_val;           /* The values to which the layers   */
                                  /* are damped.                      */


void initialize_sponge(double Iresttime[NXMEM][NYMEM], int num_fields)
{
/* Arguments: Iresttime - The inverse of the restoring time, in s-1.  */
/*  (in)      num_fields - The number of fields that will later be    */
/*                         registered.                                */

/* This subroutine determines the number of points which are within   */
/* sponges in this computational domain.  Only points that have       */
/* positive values of Iresttime and which hmask (D) indicates are ocean   */
/* points are included in the sponges.  This subroutine then allocates*/
/* the memory for num_fields reference profiles, to be provided later */
/* by calls to set_up_sponge_field.                                   */

  int i, j, col = 0;

  nf = num_fields;
  nfNZ = nf*NZ;
  if (nf > NOFLD) {
    printf("Increse NOFLD to at least %d in sponge.c or decrease the\n"
       "number of fields to be damped in the call to initialize_sponge.\n",
       nf);
    quit(-1000);
  }

  num_col = 0;
  for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++)
    if ((Iresttime[i][j]>0.0) && (D[i][j]>MINIMUM_DEPTH)) num_col++;

  if (num_col > 0) {
    col_addr = calloc(num_col,sizeof(int));
    if (col_addr == NULL) {
      printf("Error allocating memory for col_addr. num_col = %d.\n",
              num_col);
      quit(-1000);
    }
    Iresttime_col = calloc(num_col,sizeof(double));
    if (Iresttime_col == NULL) {
        printf("Error allocating memory for Iresttime_col. num_col = %d.\n",
                num_col); 
        quit(-1000);
    }
    Ref_val = calloc(num_col*nfNZ,sizeof(double));
    if (Ref_val == NULL) {
      printf("Error allocating memory for Ref_val. "
             "num_col = %d, points per column = %d.\n", num_col,nfNZ);
      quit(-1000);
    }

    for (j=Y1;j<=ny;j++) for (i=X1;i<=nx;i++) {
      if ((Iresttime[i][j]>0.0) && (D[i][j]>MINIMUM_DEPTH)) {
        col_addr[col] = i*NYMEM+j;   /* j*NXMEM+i;*/
        Iresttime_col[col] = Iresttime[i][j];
        col++;
      }
    }
  }

}

void set_up_sponge_field(double sp_val[][NXMEM][NYMEM],
                         double (*f_ptr)[][NXMEM][NYMEM], int nlay)
{
/* Arguments: sp_val - The reference profiles of the quantity being   */
/*                     registered.                                    */
/*  (in)      f_ptr - a pointer to the field which will be damped.    */
/*  (in)      nlay - the number of layers in this quantity.           */

/* This subroutine stores the reference profile for the variable      */
/* whose address is given by f_ptr. nlay is the number of layers in   */
/* this variable.  The first call to this subroutine must be to       */
/* register the interface depth profiles, and the second must be to   */
/* register the mixed layer buoyancy (Rml) profiles if BULKMIXEDLAYER */
/* is defined.  Subsequent calls can be made in any order.            */

  int i, j, k, col;

  if (fldno+1 > nf) {
    printf("Only %d fields have been allocated, as per initialize_sponge,\n"
           "but set_up_sponge_field has been called %d times.\n",
           nf, fldno+1);
    quit(-1000);
  }

  for (col=0;col<num_col;col++) {
    i = col_addr[col]/NYMEM; j = col_addr[col]-NYMEM*i;
    for (k=0;k<nlay;k++) Ref_val[col*nfNZ + k*nf + fldno] = sp_val[k][i][j];
    for (k=nlay;k<NZ;k++) Ref_val[col*nfNZ + k*nf + fldno] = 0.0;
  }

  var[fldno] = f_ptr;

  if (nlay!=NZ)
    printf("Danger:  Sponge reference fields require NZ (%d) layers.\n"
           "  A field with %d layers was passed to set_up_sponge_field.\n",NZ,nlay);
  fldno++;

}

void apply_sponge(double h[NZ][NXMEM][NYMEM], double dt,
         double ea[NZ][NXMEM][NYMEM], double eb[NZ][NXMEM][NYMEM])
{
/* Arguments: h -  Layer thickness, in m.                             */
/*            dt - The amount of time covered by this call, in s.     */
/*            ea - an array to which the amount of fluid entrained    */
/*                 from the layer above during this call will be      */
/*                 added, in m.                                       */
/*            eb - an array to which the amount of fluid entrained    */
/*                 from the layer below during this call will be      */
/*                 added, in m.                                       */

/* This subroutine applies damping to the layers thicknesses, mixed   */
/* layer buoyancy, and a variety of tracers for every column where    */
/* there is damping.                                                  */

  double damp;     /* The timestep times the local damping            */
                   /* coefficient.  Nondimensional.                   */
  double e[NZ+1];  /* The interface heights, in m, usually negative.  */
  double w;        /* The thickess of water moving upwards through an */
                   /* interface within 1 timestep, in m.              */
  double wm;       /* wm is w if w is negative and 0 otherwise, in m. */
  double wb;       /* w at the interface below a layer, in m.         */
  double wpb;      /* wpb is wb if wb is positive and 0 otherwise, m. */
  int c, i, j, k, m;

  for (c=0;c<num_col;c++) {
    i = col_addr[c]/NYMEM; j = col_addr[c]-NYMEM*i;
    damp = dt*Iresttime_col[c];

    e[NZ] = -D[i][j];
    for (k=NZ-1;k>=0;k--) e[k] = e[k+1] + h[k][i][j];

    wpb = 0.0; wb = 0.0;
    for (k=NZ-1;k>=0;k--) {
      w = D_MIN(((e[k] - Ref_val[c*nfNZ+k*nf+0]) * damp),
                (wb+h[k][i][j] - EPSILON));
      wm = 0.5*(w-fabs(w));
      for (m=1;m<fldno;m++) {
        (*var[m])[k][i][j] = (h[k][i][j]*(*var[m])[k][i][j] +
                             Ref_val[c*nfNZ+k*nf+m] *
                             (damp*h[k][i][j] + wpb - wm))/
                   (h[k][i][j]*(1.0 + damp) + wpb - wm);
        tr[m-1][k][i][j] = (*var[m])[k][i][j];
      }
      wb = w; wpb = w - wm;
    }

  }
}


