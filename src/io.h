/* By R. Hallberg 7/99                                                */
/*   This header contains the prototypes and definitions for using    */
/* the I/O subroutines.                                               */

#include <stddef.h>
#define DOUBLE_FILE 0
#define FLOAT_FILE 1
#define NETCDF_FILE 2


struct vardesc {char name[30]; /* The variable name in a NetCDF file. */
                char longname[100]; /* The long name of that variable. */
                char hor_grid; /* The hor. grid:  u, v, h, q, or 1.   */ 
                char z_grid;   /* The vert. grid:  L, i, 2, or 1.     */ 
                char t_grid;   /* The time description: s, a, m, or 1.*/
                char units[50];/* The dimensions of the variable.     */
                char mem_size; /* The size in memory: d or f.         */
                double mval;   /* The missing value */};
/* For example:
struct vardesc vars[2] =
  { {"D","Basin Depth",'h','1','1',"meter", 'd',0.},
    {"taux","Zonal Wind Stress",'u','1','1',"meter2 second-2", 'd',-1.e-6} };
*/

struct varcdfinfo {int id;            /* NetCDF variable ID.          */
                   size_t count[4];   /* The size of the output.      */
                   int x_size;        /* The number of output points  */
                   int y_size;        /* in the x-, y- and z-         */
                   int z_size; };     /* directions.                  */

size_t write_field(int cdfid, FILE *fileptr, struct vardesc vars,
                   struct varcdfinfo varinfo, size_t nrec, double *variable);

void create_file(char filename[], int type, struct vardesc vars[], int novars,
                 FILE **fileptr, int *cdfid, int *timeid, 
                 struct varcdfinfo varinfo[], int output_prec_float);

int name_output_file(char name[], double day, int noexdig,
                     char filepath[], char directory[]);

size_t write_time(int cdfid, FILE *fileptr, int timeid, 
                  size_t nrec, double *day);

size_t read_time(int cdfid, FILE *fileptr, int timeid,
                 size_t nrec, double *day);

void close_file(int *cdfid, FILE **fileptr);

void find_input_file(char name[], char name2[], char fltname[], double day,
                     char directory[], FILE **file, int *cdfid,
                     int *timeid);

int open_input_file(char filename[], FILE **fileptr, int *cdfid, int *timeid);

void read_field(int cdfid, FILE *fileptr, char varname[],
                int size_x, int size_y, int size_z,
                int start_file_x, int start_file_y, int start_file_z,
                int start_array_x, int start_array_y, int start_array_z,
                int recno, double variable[]);                   

double*** alloc3d(int NZED, int NY, int NX);

