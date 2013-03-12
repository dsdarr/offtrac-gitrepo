//HF new utility routines
#include "init.h"

void set_darray2d_zero(double **arr, int NX, int NY) {
	int x, y;
	for (x = 0; x < NX; x++)
		for (y = 0; y < NY; y++)
			arr[x][y] = 0.0;
}

void set_darray3d_zero(double ***arr, int nz, int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr[z][x][y] = 0.0;
}

void set_fix_darray3d_zero(double arr[][NXMEM][NYMEM], int nz) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NXMEM; x++)
			for (y = 0; y < NYMEM; y++)
				arr[z][x][y] = 0.0;
}

void set_fix_darray2d_zero(double arr[NXMEM][NYMEM]) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			arr[x][y] = 0.0;
}

void add_fix_darray2d(double arr1[NXMEM][NYMEM], double arr2[NXMEM][NYMEM]) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			arr1[x][y] += arr2[x][y];
}

void mult_fix_darray2d(double arr[NXMEM][NYMEM], double factor) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			arr[x][y] *= factor;
}

void add_darray3d(double ***arr1, double ***arr2, int nz, int NX, int NY) {
	int x, y, z;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] += arr2[z][x][y];
}

void mult_darray3d(double ***arr, int nz, int NX, int NY, double factor) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr[z][x][y] *= factor;
}
//BX-a
void mult_fix_darray2d_mv(double arr[NXMEM][NYMEM], double factor,
		double D[NXMEM][NYMEM], double mv) {
	int x, y;
	for (x = 0; x < NXMEM; x++)
		for (y = 0; y < NYMEM; y++)
			if (D[x][y] > MINIMUM_DEPTH) {
				arr[x][y] *= factor;
			} else {
				arr[x][y] = mv;
			}
}

void mult_darray3d_mv(double ***arr, int nz, int NX, int NY, double factor,
		double D[NXMEM][NYMEM], double mv) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				if (D[x][y] > MINIMUM_DEPTH) {
					arr[z][x][y] *= factor;
				} else {
					arr[z][x][y] = mv;
				}
}
//BX-e
/*
 void copy_fix_darray2d(double arr1[NXMEM][NYMEM], double arr2[NXMEM][NYMEM])
 {
 int x, y;
 for (x=0; x<NXMEM; x++)
 for (y=0;y<NYMEM; y++)
 arr1[x][y] = arr2[x][y];
 }
 */
void copy_darray3d(double ***arr1, double ***arr2, int nz, int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] = arr2[z][x][y];
}
//BX-a
void copy_fix_darray3d(double ***arr1, double arr2[NZ][NXMEM][NYMEM], int nz,
		int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] = arr2[z][x][y];
}
//BX-e
// begin ashao

void copy_2fix_darray3d(double (*arr1)[NXMEM][NYMEM], double (*arr2)[NXMEM][NYMEM], int nz,
		int NX, int NY) {
	int z, x, y;
	for (z = 0; z < nz; z++)
		for (x = 0; x < NX; x++)
			for (y = 0; y < NY; y++)
				arr1[z][x][y] = arr2[z][x][y];
}
double linear_interp(double x0, double y0, double x1, double y1, double xstar) {
	double ystar;

	// If the two interpolation points are identical, set the output y to input y
	if ( (fabs(y1-y0) < .000000001) || (fabs(x1-x0) < .0000000001) ){
		 ystar=y0;
	}
	// Perform linear interpolation
	else {
		ystar = y0 + (xstar - x0) * (y1 -  y0) / (x1 - x0);
	}

	return ystar;
}






// end ashao
