//HF header file for alloc3d.c

double*** alloc3d(int NZED, int NY, int NX);

float*** alloc3d_f(int NZED, int NY, int NX);

void free3d(double*** arr3d, int NZED);

void free3d_f(float*** arr3d, int NZED);


float** alloc2d_f(int NY, int NX);

void free2d_f(float** arr2d, int NY);
