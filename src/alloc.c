/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *             Memory allocation routine for 3D fields                *
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

double*** alloc3d(int NZED, int NY, int NX)
{
  int iy, iz;

  double ***arr3d = (double ***) malloc(NZED * sizeof(double **));
  if (arr3d == NULL)
    return NULL;

  for (iz = 0; iz < NZED; iz++) {
    arr3d[iz] = (double **) malloc(NY * sizeof(double *));
    if (arr3d[iz] == NULL)
      return NULL;
  }
  
  arr3d[0][0] = (double *) malloc(NZED * NY * NX * sizeof(double));
  if (arr3d[0][0] == NULL)
    return NULL;
  
  for (iy = 1; iy < NY; iy++)
    arr3d[0][iy] = arr3d[0][0] + iy * NX;
  
  for (iz = 1; iz < NZED; iz++)
    for (iy = 0; iy < NY; iy++)
      arr3d[iz][iy] = arr3d[0][0] + (iz * NY + iy) * NX;

  return arr3d;
}

float*** alloc3d_f(int NZED, int NY, int NX)
{
  int iy, iz;

  float ***arr3d = (float ***) malloc(NZED * sizeof(float **));
  if (arr3d == NULL)
    return NULL;

  for (iz = 0; iz < NZED; iz++) {
    arr3d[iz] = (float **) malloc(NY * sizeof(float *));
    if (arr3d[iz] == NULL)
      return NULL;
  }
  
  arr3d[0][0] = (float *) malloc(NZED * NY * NX * sizeof(float));
  if (arr3d[0][0] == NULL)
    return NULL;
  
  for (iy = 1; iy < NY; iy++)
    arr3d[0][iy] = arr3d[0][0] + iy * NX;
  
  for (iz = 1; iz < NZED; iz++)
    for (iy = 0; iy < NY; iy++)
      arr3d[iz][iy] = arr3d[0][0] + (iz * NY + iy) * NX;

  return arr3d;
}

void free3d(double*** arr3d, int NZED)
{
    int iz;
    free(arr3d[0][0]);
    arr3d[0][0] = NULL;
    for (iz = 0; iz < NZED; iz++) {
	free(arr3d[iz]);
	arr3d[iz] = NULL;
    }
    free(arr3d);
    arr3d = NULL;
}

void free3d_f(float*** arr3d, int NZED)
{
    int iz;
    free(arr3d[0][0]);
    arr3d[0][0] = NULL;
    for (iz = 0; iz < NZED; iz++) {
	free(arr3d[iz]);
	arr3d[iz] = NULL;
    }
    free(arr3d);
    arr3d = NULL;
}


/////////////////////////////////////////////////////////////////////
// T W O - D I M E N S I O N A L   A R R A Y S
/////////////////////////////////////////////////////////////////////

float** alloc2d_f(int NY, int NX)
{
  int iy;

  float **arr2d = (float **) malloc(NY * sizeof(float *));
  if (arr2d == NULL)
    return NULL;

  /* allocate the memory for the array */
  arr2d[0] = (float *) malloc(NY * NX * sizeof(float));
  if (arr2d[0] == NULL)
    return NULL;

  /* assign pointers to rows */
  for (iy = 1; iy < NY; iy++)
    arr2d[iy] = arr2d[0] + iy * NX;

  return arr2d;
}

void free2d_f(float** arr2d, int NY)
{
  free(arr2d[0]);
  arr2d[0] = NULL;

  free(arr2d);
  arr2d = NULL;
}

