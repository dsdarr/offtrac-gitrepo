#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "init.h"

extern double umask[NXMEM][NYMEM];   /* _mask are 1 over ocean and 0  */
extern double vmask[NXMEM][NYMEM];   /* over land on the u & v grids. */

extern double D[NXMEM][NYMEM];

void initializemasks() {

  int i,j;
  double Dmin;

  Dmin = MINIMUM_DEPTH;  /* value set in init.h */
  Dmin = (Dmin > 2.0*EPSILON) ? Dmin : 2.0*EPSILON;


  for (j=Y0;j<=ny+2;j++) {
    for (i=X0-1;i<=nx+2;i++) {
      if ((D[i][j] <= Dmin) || (D[i+1][j] <= Dmin)) umask[i][j] = 0.0;
      else umask[i][j] = 1.0;
    }
  }

  for (j=Y0-1;j<=ny+2;j++) {
    for (i=X0;i<=nx+2;i++) {
      if ((D[i][j] <= Dmin) || (D[i][j+1] <= Dmin)) vmask[i][j] = 0.0;
      else vmask[i][j] = 1.0;
    }
  }

  /*
  for (j=27;j<=33;j++) {
    printf("D = %g,%g,%g,%g,%g,%g \n", 
	   D[190][j],D[191][j],D[192][j],D[193][j],D[194][j],D[195][j]);
  }
  for (j=27;j<=33;j++) {
    printf("umask = %g,%g,%g,%g,%g,%g \n", 
	   umask[190][j],umask[191][j],umask[192][j],umask[193][j],umask[194][j],umask[195][j]);
  }
  for (j=27;j<=33;j++) {
    printf("vmask = %g,%g,%g,%g,%g,%g \n", 
	   vmask[190][j],vmask[191][j],vmask[192][j],vmask[193][j],vmask[194][j],vmask[195][j]);
  }
  */


}
