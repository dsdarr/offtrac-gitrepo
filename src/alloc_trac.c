/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *          Memory allocation routine for tracer fields               *
 *                Original code by Hartmut Frenzel                    *
 *              Implemented by Holger Brix - 01NOV07                  *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "init.h"
#include "io.h"


int alloc_trac()
{

  int ix, iz, iv;
  extern double ****tr;

  /* allocation of tr (four levels, all with checks) */
  tr = (double****) malloc(NTR * sizeof(double ***));
  if (tr == NULL) {
      fprintf(stderr, "out of memory in level 1 allocation!\n");
      return(1);
  }
  
  for (iv = 0; iv < NTR; iv++) {
      tr[iv] = (double ***) malloc(NZ * sizeof(double **));
      if (tr[iv] == NULL) {
	  fprintf(stderr, "out of memory in level 2 allocation!\n");
	  return(1);
      }
      
      for (iz = 0; iz < NZ; iz++) {
	  tr[iv][iz] = (double **) malloc(NXMEM * sizeof(double *));
	  if (tr[iv][iz] == NULL) {
	      fprintf(stderr, "out of memory in level 3 allocation!\n");
	      return(1);
	  }
	  
	  for (ix = 0; ix < NXMEM; ix++) {
	      tr[iv][iz][ix] = (double *) malloc(NYMEM * sizeof(double));
	      if (tr[iv][iz][ix] == NULL) {
		  fprintf(stderr, "out of memory in level 4 allocation!\n");
		  return(1);
	      }
	  }
      }
  } /* end of allocation block */

  return 0;
}
