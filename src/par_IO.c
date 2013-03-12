/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *  By Robert Hallberg, August 1998                                   *
 *                                                                    *
 *    This program various subroutines which are necessary to make    *
 *  the model run on a parallel computer.  Each of these subroutines  *
 *  except set_up_parallel is also called when not on a parallel      *
 *  computer, but in that case the subroutines are quite simple.      *
 *                                                                    *
 *    collect_double_fields takes whichever field it is given and     *
 *  collects it into io_field_d.                                      *
 *                                                                    *
 *    spread_double_fields distributes one layer's worth of data to   *
 *  the appropriate locations on each processor.                      *
 *                                                                    *
 *    sum_fields does a sum across all processors of the vector       *
 *  indicated by the first argument with a length given by its second *
 *  argument, and returns these sums in that vector.                  *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdio.h>
#include "init.h"
#include <stdlib.h>


#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
#include <mpi.h>
# else
#include <mpp/shmem.h>
# endif
#endif


extern double io_field_d[(NXTOT+1)*(NYTOT+1)]; /* A global field that */
                                     /* is used in all double pre-    */
                                     /* cision I/O.                   */
#ifdef PARALLEL_X
extern int nx;                       /* The number of x-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
extern int X0abs;                    /* The absolute index of X0 on   */
                                     /* the current processor, rel-   */
                                     /* ative to X0 on processor 0.   */
#else
# ifdef PARALLEL_Y
static int X0abs = 0;                /* On 1 processor, X0abs = 0.    */
# endif
#endif
#ifdef PARALLEL_Y
extern int ny;                       /* The number of y-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
extern int Y0abs;                    /* The absolute index of Y0 on   */
                                     /* the current processor, rel-   */
                                     /* ative to Y0 on processor 0.   */
#else
# ifdef PARALLEL_X
static int Y0abs = 0;                /* On 1 processor, Y0abs = 0.    */
# endif
#endif

#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifndef MAXPROC
#  ifndef NXPROC
#   define NXPROC 1
#  endif
#  ifndef NYPROC
#   define NYPROC 1
#  endif
#  define MAXPROC ((NXPROC)*(NYPROC))
# endif
# ifdef CHECKPARALLEL
#  define MAXPROCS ((MAXPROC)+1)
# else
#  define MAXPROCS (MAXPROC)
# endif
extern int pe_here;                  /* The current processor's label.*/
extern int X0abs_proc[MAXPROCS]; /* The absolute index of X0 on each  */
                                 /* processor relative to X0 on       */
                                 /* processor 0.                      */
extern int Y0abs_proc[MAXPROCS]; /* The absolute index of Y0 on each  */
                                 /* processor relative to Y0 on       */
                                 /* processor 0.                      */
extern int nx_proc[MAXPROCS];    /* The number of points in the comp- */
extern int ny_proc[MAXPROCS];    /* utational domain of each processor*/
                                 /* in the x- and y- directions.      */
extern int procs_used;           /* The number of processors actually */
                                 /* used in the present simulation.   */
extern int procs_used_all;       /* procs_used_all equals procs_used  */
                                 /* unless CHECKPARALLEL is used, in  */
                                 /* which case it is 1 greater.       */
#else
#define pe_here 0                    /* Without running in parallel,  */
                                     /* the current processor is 0.   */
#endif


void collect_double_fields(double *field_in, int offx, int offy)
{
/*   Gather a layer's worth of field_in from every processor, and put */
/* it into io_field_d, with offsets offx, offy between the X0,Y0      */
/* corner of the data and the southwest corner of io_field_d.         */
  int i, j;

#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
  {
    if (pe_here == 0) {
      MPI_Request req[MAXPROC];
      MPI_Status stat;
      double pass_field[MAXPROC][NXMEM*NYMEM];
      int p, p2;

      for (p=1;p<procs_used;p++) {
        MPI_Irecv(pass_field[p], NXMEM*NYMEM, MPI_DOUBLE, p, (1000+p), MPI_COMM_WORLD, &req[p]);
      }

      for (i=X0+offx;i<=nx;i++) 
//#pragma _CRI ivdep
        for (j=Y0+offy;j<=ny;j++)
          io_field_d[(i-X0-offx)*(NYTOT+1-offy)+j-Y0-offy] = field_in[NYMEM*i+j];

      for (p2=1;p2<procs_used;p2++) {
        MPI_Waitany(procs_used-1, &req[1], &p, &stat); p++;
        for (i=(((X0abs_proc[p]==0)&&(offx==0))?X0:X1);i<=nx_proc[p];i++) {
          for (j=(((Y0abs==0)&&(offy==0))?Y0:Y1);j<=ny_proc[p];j++) {
            io_field_d[(i+X0abs_proc[p]-X0-offx)*(NYTOT+1-offy) + j+Y0abs_proc[p]-Y0-offy] =
              pass_field[p][NYMEM*i+j];
          }
        }
      }
    }
    else if (pe_here < procs_used) {
      MPI_Send(field_in, NXMEM*NYMEM, MPI_DOUBLE, 0, (1000+pe_here), MPI_COMM_WORLD);
    }
  }
# else
  barrier();
  if (pe_here == 0) {
    for (i=X0+offx;i<=nx;i++) 
//#pragma _CRI ivdep
      for (j=Y0+offy;j<=ny;j++)
        io_field_d[(i-X0-offx)*(NYTOT+1-offy)+j-Y0-offy] = field_in[NYMEM*i+j];
  }
  else {
    for (i=(((X0abs==0)&&(offx==0))?X0:X1);i<=nx;i++) {
      if ((Y0abs == 0) && (offy == 0))
        shmem_double_put(&io_field_d[(X0abs+i-X0-offx)*(NYTOT+1-offy)],
                         field_in+NYMEM*i+Y0,ny+1-Y0,0);
      else
        shmem_double_put(&io_field_d[(X0abs+i-X0-offx)*(NYTOT+1-offy) + Y0abs+1-offy],
                         field_in+NYMEM*i+Y0+offy,ny-Y0,0);
    }
  }
  barrier();

# endif
#else
  for (i=offx;i<=NXTOT;i++) for (j=offy;j<=NYTOT;j++)
    io_field_d[(i-offx)*(NYTOT+1-offy)+j-offy] = field_in[NYMEM*(i+X0)+j+Y0];
#endif
}

void spread_double_fields(double field_out[][NYMEM])
{
/*   Distribute a layer's worth of io_field_d to the appropriate      */
/* processors in field_out.                                           */
  int i, j;

#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
  {
    if (pe_here == 0) {
      MPI_Request req[MAXPROCS];
      MPI_Status stat[MAXPROCS];
      double pass_field[MAXPROCS][NXMEM*NYMEM];
      int p, p2;

      for (p=0;p<procs_used_all;p++) {
        for (i=0;i<=nx_proc[p]-X0;i++) for (j=0;j<=ny_proc[p]-Y0;j++) {
          pass_field[p][NYMEM*i+j] =
              io_field_d[(i+X0abs_proc[p])*(NYTOT+1) + j+Y0abs_proc[p]];
        }
        MPI_Isend(pass_field[p], NXMEM*NYMEM, MPI_DOUBLE, p, (2000+p), MPI_COMM_WORLD, &req[p]);
      }

      for (i=X0;i<=nx;i++) 
//#pragma _CRI ivdep
        for (j=Y0;j<=ny;j++)
          field_out[i][j] = io_field_d[(i-X0)*(NYTOT+1)+j-Y0];

      MPI_Waitall(procs_used_all, req, stat);
    }
    else {
      MPI_Status stat;
      double pass_field[NXMEM*NYMEM];
      MPI_Recv(pass_field, NXMEM*NYMEM, MPI_DOUBLE, 0, (2000+pe_here), MPI_COMM_WORLD, &stat);
      for (i=X0;i<=nx;i++)
//#pragma _CRI ivdep
        for (j=Y0;j<=ny;j++)
          field_out[i][j] = pass_field[(i-X0)*NYMEM+j-Y0];
    }
  }
# else
  barrier();
  if (pe_here == 0) {
    for (i=X0;i<=nx;i++)
//#pragma _CRI ivdep
      for (j=Y0;j<=ny;j++)
        field_out[i][j] = io_field_d[(i-X0)*(NYTOT+1)+j-Y0];

  }
  else {
    for (i=X0;i<=nx;i++) {
      shmem_double_get(&field_out[i][Y0],
                       &io_field_d[(X0abs+i-X0)*(NYTOT+1)+Y0abs],ny-Y0+1,0);
    }
  }
  barrier();
# endif
#else
  for (i=0;i<=NXTOT;i++)
//#pragma _CRI ivdep
   for (j=0;j<=NYTOT;j++) field_out[i+X0][j+Y0]=io_field_d[i*(NYTOT+1)+j];
#endif

}

void spread_double_vector(double vec[], int numpts)
{
/*   This copies the numpts points in vec on processor 0 to every     */
/* vec on every other processor.                                      */
#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
  MPI_Bcast(vec, numpts, MPI_DOUBLE, 0, MPI_COMM_WORLD);
# else
  int i;
  if (pe_here == 0) 
//#pragma _CRI ivdep
    for (i=0;i<numpts;i++) io_field_d[i] = vec[i];
  barrier();
  if (pe_here != 0) {
    shmem_double_get(io_field_d,io_field_d,numpts,0);
//#pragma _CRI ivdep
    for (i=0;i<numpts;i++) vec[i] = io_field_d[i];
  }
  barrier();
# endif
#endif
}

size_t read_layer(double field_out[][NYMEM], FILE *file, int float_vals)
{
/*   Read in one layer's worth data from file and distribute it to    */
/* the appropriate processors in field_out. float_vals is 1 if the    */
/* input file contains floating point numbers and 0 otherwise.        */
  float float_in[NYTOT+1];     /* A column of a layer of a variable.  */
  int floatsize = sizeof(float);
  int i, j, fs = (NXTOT+1)*(NYTOT+1)*sizeof(double);
  size_t err;

/* First PE 0 reads the data from file.                               */

  err = 1;
  if (pe_here == 0) {
    if (float_vals) {
      for (i=0;i<=NXTOT;i++) {
        err *= fread((void *)float_in,floatsize*(NYTOT+1),1,file);
        for (j=0;j<=NYTOT;j++) io_field_d[i*(NYTOT+1)+j] = (double) float_in[j];
      }
    }
    else {
      err = fread((void *)io_field_d,fs,1,file);
    }
  }

/* Now each processor gathers its data from PE 0.                     */

#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
    spread_double_fields(field_out);
# else
  barrier();
  if (pe_here == 0) {
    for (i=X0;i<=nx;i++)
//#pragma _CRI ivdep
      for (j=Y0;j<=ny;j++)
        field_out[i][j] = io_field_d[(i-X0)*(NYTOT+1)+j-Y0];
  }
  else {
    for (i=X0;i<=nx;i++) {
      shmem_double_get(&field_out[i][Y0],
                       &io_field_d[(X0abs+i-X0)*(NYTOT+1)+Y0abs],ny-Y0+1,0);
    }
  }
  barrier();

# endif
#else
  for (i=0;i<=NXTOT;i++) 
//#pragma _CRI ivdep
    for (j=0;j<=NYTOT;j++) field_out[i+X0][j+Y0] = io_field_d[i*(NYTOT+1)+j];
#endif

  return(err);
}

size_t write_layer(double field_in[], FILE *file, int float_vals)
{
/*   Gather a layer's worth of field_in from every processor, and     */
/* write it out to file.  float_vals is 1 if floating point numbers   */
/* are to be written and 0 to write double precision numbers.         */
  float float_out[NYTOT+1];    /* A column of a layer of a variable.  */
  int floatsize = sizeof(float);
  int i, j, fs = (NXTOT+1)*(NYTOT+1)*sizeof(double);
  size_t err;

/* Collect the data to PE 0.                                          */

#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
    collect_double_fields(field_in, 0, 0);
# else
  barrier();
  if (pe_here == 0) {
    for (i=X0;i<=nx;i++)
//#pragma _CRI ivdep
      for (j=Y0;j<=ny;j++)
        io_field_d[(i-X0)*(NYTOT+1)+j-Y0] = field_in[i*NYMEM+j];
  }
  else {
    for (i=((X0abs==0)?X0:X1);i<=nx;i++) {
      if (Y0abs == 0)
        shmem_double_put(&io_field_d[(X0abs+i-X0)*(NYTOT+1)],
                         &field_in[i*NYMEM+Y0],ny-Y0+1,0);
      else
        shmem_double_put(&io_field_d[(X0abs+i-X0)*(NYTOT+1) + Y0abs+1],
                         &field_in[i*NYMEM+Y1],ny-Y0,0);
    }
  }
  barrier();

# endif
#else
  for (i=0;i<=NXTOT;i++) for (j=0;j<=NYTOT;j++)
    io_field_d[i*(NYTOT+1)+j] = field_in[(i+X0)*NYMEM+j+Y0];
#endif

/* Here PE 0 writes out the data.                                     */

  err = 1;
  if (pe_here==0) {
    if (float_vals) {
      for (i=0;i<=NXTOT;i++) {
        for (j=0;j<=NYTOT;j++) float_out[j] = (float) io_field_d[i*(NYTOT+1)+j];
        err *= fwrite((void *)float_out,floatsize*(NYTOT+1),1,file);
      }
    } else {
      err *= fwrite((void *)io_field_d,fs,1,file);
    }
  }

  return(err);
}

void quit(int status) {
/* This subroutine calls the appropriate exit call to terminate the */
/* execution on all processors.                                     */
#if defined(PARALLEL_X) || defined(PARALLEL_Y)
# ifdef MPI_PARALLEL
  MPI_Abort(MPI_COMM_WORLD, status);
# else
  globalexit(status);
# endif
#else
  exit(status);
#endif
}
