/* By R. Hallberg 7/99                                                */
/*   This file contains a subroutine that will create a netcdf file   */
/* or a binary file, a subroutine write a field to that file, a       */
/* subroutine to name a file given a root-name and a time,            */
/* time, binary file, a subroutine to name a file given a root and a   */

#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "init.h"
#ifdef NETCDF_OUTPUT
#include <netcdf.h>
#endif
#include "metrics.h"
#include "io.h"
#include "iocdf.h"
#include "offtrac.h"
#include "par_IO.h"


#ifdef NETCDF_OUTPUT
extern double latq[NYMEM];          /* The latitude of q points.     */
extern double lath[NYMEM];          /* The latitude of h points.     */
extern double lonq[NXMEM];          /* The longitude of q points.    */
extern double lonh[NXMEM];          /* The longitude of h points.    */
# if defined(PARALLEL_Y) && !defined(PARALLEL_IO)
double lath_tot[NYTOT+Y1];           /* The latitude of all points in */
double latq_tot[NYTOT+Y1];           /* the whole domain.             */
# endif
# if defined(PARALLEL_X) && !defined(PARALLEL_IO)
double lonh_tot[NXTOT+X1];           /* The longitude of all points   */
double lonq_tot[NXTOT+X1];           /* in the whole domain.          */
# endif
#endif

double io_field_d[(NXTOT+1)*(NYTOT+1)]; /* A global field that        */
                                     /* is used in all double pre-    */
                                     /* cision I/O.                   */

extern void quit(int status);

#if defined(PARALLEL_X) || defined(PARALLEL_Y)
extern int pe_here;                  /* The current processor's label.*/
extern int procs_used;      /*   The number of processors actually    */
                            /* used in the present simulation.        */
#else
# undef PARALLEL_IO
# define pe_here 0                   /* Without running in parallel,  */
                                     /* the current processor is 0.   */
#endif

#ifdef PARALLEL_X
extern int nx;        //void step_fields(double day);
               /* The number of x-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
extern int X0abs;                    /* The absolute index of X0 on   */
                                     /* the current processor, rel-   */
                                     /* ative to X0 on processor 0.   */
#else
# define X0abs 0                     /* X0abs on processor 0 is 0.    */
#endif

#ifdef PARALLEL_Y
extern int ny;                       /* The number of y-points in the */
                                     /* physical domain calculated by */
                                     /* the current processor.        */
extern int Y0abs;                    /* The absolute index of Y0 on   */
                                     /* the current processor, rel-   */
                                     /* ative to Y0 on processor 0.   */
#else
# define Y0abs 0                     /* Y0abs on processor 0 is 0.    */
#endif

size_t write_field(int cdfid, FILE *fileptr, struct vardesc vars,
                   struct varcdfinfo varinfo, size_t nrec, double *variable)
{

/*   This subrountine will write a full variable, either 1 snapshot of */
/* a timeseries or a single field.  If a netcdf file is used, cdfid    */
/* contains the ID number of the open netcdf file (if non-negative) or */
/* -1 for a non-writing PE if PE 0 is writing a single-threaded netCDF */
/* file.  If a binary file is to be written, cdfid is either -4 (for a */
/* file containing 4-byte floating point numbers) or -8 (for a file    */
/* containing 8-byte floating point numbers).  fileptr points to an    */
/* open binary file.  vars contains the variable descriptions in a     */
/* structure of type vardesc.  varinfo is a structure of type          */
/* varcdfinfo which is obtained from create_file.  nrec is the time    */
/* index number of a netCDF file.  variable is a pointer to the first  */
/* element of the array to be written.                                 */

#ifdef NETCDF_OUTPUT
  if (cdfid > -4) {
    int status;
    char message[144];
    size_t start[4] = {0,0,0,0};

    start[0] = nrec;


    if ((nrec == 0) || (vars.t_grid != '1')) {
      if (vars.hor_grid == '1') {
        if (cdfid >= 0) {
          if (vars.t_grid == '1')  {
            status = nc_put_var_double(cdfid, varinfo.id, variable);
            if (status != NC_NOERR) {
              strcpy(message,"put "); strcat(message,vars.name);
              handle_error(message,status);
            }
          }
          else {
	    printf("writing %s\n", vars.name);
            status = nc_put_vara_double(cdfid, varinfo.id, start,
                                        varinfo.count, variable);
            if (status != NC_NOERR) {
              strcpy(message,"put a "); strcat(message,vars.name);
              handle_error(message,status);
            }
          }
        }
      }
      else {
# if (defined(PARALLEL_Y) || defined(PARALLEL_X)) && !defined(PARALLEL_IO)
        double out_array[(NXTOT+1)*(NYTOT+1)];
        int zindex, i, j, k;

        if (vars.z_grid == '1') {zindex = 3;}
        else if (vars.t_grid == '1') zindex = 0;
        else  zindex = 1;

        varinfo.count[zindex] = 1;
        for (k=0;k<varinfo.z_size;k++) {
          start[zindex] = k;
          collect_double_fields(variable+k*NXMEM*NYMEM,
                      (NXTOT+1-varinfo.x_size),(NYTOT+1-varinfo.y_size));
          if (pe_here == 0) {
            for (i=0;i<varinfo.x_size;i++) for (j=0;j<varinfo.y_size;j++)
              out_array[j*varinfo.x_size+i] = io_field_d[i*varinfo.y_size+j];
            status = nc_put_vara_double(cdfid, varinfo.id, start,
                                        varinfo.count, &out_array[0]);
            if (status != NC_NOERR) {
              sprintf(message,"putting layer %d of ",k); strcat(message,vars.name);
              handle_error(message,status);
            }
          }
        }
# else
/*  With 2-D or 3-D arrays, the data is copied to another field to     */
/*  eliminate the halo regions.                                        */
        double out_array[NXMEM*NYMEM*(NZ+1)];
        int i, j, k, xst, yst, xysize;

        xst = nx-varinfo.x_size+1; yst = ny-varinfo.y_size+1; 
        xysize = varinfo.y_size*varinfo.x_size;
	//printf("writing %s\n", vars.name);
        for (k=0;k<varinfo.z_size;k++) {
          for (i=xst;i<=nx;i++) for (j=yst;j<=ny;j++) {
              out_array[k*xysize + (j-yst)*varinfo.x_size + (i-xst)] = 
                    variable[k*NXMEM*NYMEM + i*NYMEM + j];
	      /* DEBUG if (fabs(variable[k*NXMEM*NYMEM + i*NYMEM + j]) > 1e+35)
		printf("out_array(%d,%d,%d)=%g\n",k,i,j,
		variable[k*NXMEM*NYMEM + i*NYMEM + j]); */
          }
        }
	
        status = nc_put_vara_double(cdfid, varinfo.id, start,
                                    varinfo.count, &out_array[0]);

        if (status != NC_NOERR) {
          strcpy(message,"put m "); strcat(message,vars.name);
          handle_error(message,status);
        }
# endif
      }
    }

    return 1;
  }
  else
#endif
  {
    int k, nlay;
    size_t err = 1;

    if (vars.z_grid == '1') nlay = 1;
    else if (vars.z_grid == 'L') nlay = NZ;
    else if (vars.z_grid == 'i')  nlay = NZ+1;
    else nlay = 0;

    if (vars.hor_grid == '1') {
      if (pe_here == 0) err *= fwrite((void *)variable,8*nlay,1,fileptr);
    }
    else {
      for (k=0;k<nlay;k++) {
        if (cdfid < -4)
          err *= write_layer(&variable[k*NXMEM*NYMEM],fileptr,0);
        else
          err *= write_layer(&variable[k*NXMEM*NYMEM],fileptr,1);
      }
    }

    return err;
    
  }
}

size_t write_time(int cdfid, FILE *fileptr, int timeid, size_t nrec, 
                  double *day)
{
#ifdef NETCDF_OUTPUT
  if (cdfid > -4) {
    if (cdfid >= 0) {
      int status;

      status = nc_put_var1_double(cdfid, timeid, &nrec, day);
      if (status != NC_NOERR)
        handle_error("Writing time ",status);
    }
    return 1;
  }
  else
#endif
  {
    size_t err = 1;
    if (pe_here==0) err *= fwrite((void *)day,8,1,fileptr);
    return err;
  }
}
void create_file(char filename[], int type, struct vardesc vars[], int novars,
                 FILE **fileptr, int *cdfid, int *timeid, 
                 struct varcdfinfo varinfo[], int output_prec_float)
{

#ifdef NETCDF_OUTPUT
  if (type != NETCDF_FILE)
#endif
  {
    if (pe_here == 0) *fileptr = fopen(filename,"w");
    else *fileptr = NULL;

    if (type == FLOAT_FILE) *cdfid = -4;
    else *cdfid = -8;
    return;
  }
#ifdef NETCDF_OUTPUT
  else {
    int lathvid, lonhvid, latqvid, lonqvid;
    int use_lath = 0, use_lonh = 0, use_latq = 0, use_lonq = 0;
    int layervid, layer2vid, intvid;
    int use_layer = 0, use_layer2 = 0, use_int = 0;
    int dayvid, daystvid;

    int k, status, junk, layer[NZ];
# ifdef PARALLEL_IO
    int domain[4];
# endif
    float interface[NZ+1];
    char message[144];
    size_t xhsize, yhsize, xqsize, yqsize;

#define OUTPUT_DEFS
#ifdef OUTPUT_DEFS
    int dummy, kdmlvid, kdvid, khtrvid, tauvid;
    float kdml = KDML; // need to copy from macro to variable
    float kd = KD;     // so pointer can be used in nc_put_var_float call
    float khtr = KHTR;
    extern const double r_bio_tau;
#endif /* OUTPUT_DEFS */

/*  Define the sizes in the file of the various axes.                 */

# ifdef PARALLEL_IO
    xhsize = nx-X0; yhsize = ny-Y0;
# else
    xhsize = NXTOT; yhsize = NYTOT;
# endif
# ifdef REENTRANT
    xqsize = xhsize;
# else
    xqsize = xhsize+1;
# endif
# ifdef REENTRANT_Y
    yqsize = yhsize;
# else
    yqsize = yhsize+1;
# endif


/*   Define the NetCDF file to which output will be written.          */
    if ((strncmp(&filename[strlen(filename)-4],".cdf",4) != 0) &&
        (strncmp(&filename[strlen(filename)-3],".nc",3)))
      strcat(filename,".cdf");
# ifdef PARALLEL_IO
    {
      char suffix[6];
      sprintf(suffix,".%04d",pe_here);
      strcat(filename,suffix);
    }
# elif defined(PARALLEL_Y) || defined(PARALLEL_X)
    if (pe_here != 0) {
      *cdfid = -1;
      for (k=0;k<novars;k++) {
        switch (vars[k].hor_grid) {
          case 'u': varinfo[k].x_size = xqsize; varinfo[k].y_size = yhsize; break;
          case 'v': varinfo[k].x_size = xhsize; varinfo[k].y_size = yqsize; break;
          case 'h': varinfo[k].x_size = xhsize; varinfo[k].y_size = yhsize; break;
          case 'q': varinfo[k].x_size = xqsize; varinfo[k].y_size = yqsize; break;
          case '1': break;
          default: fprintf(stderr,"Unknown horizontal axes[%d] %c.\n",k,vars[k].hor_grid);
        }
        switch (vars[k].z_grid) {
          case 'L': varinfo[k].z_size = NZ; break;
          case 'i': varinfo[k].z_size = NZ+1; break;
          case '2': varinfo[k].z_size = 2; break;
          case '1': varinfo[k].z_size = 1; break;
          default: fprintf(stderr,"Unknown layer axis[%d] %c.\n",k,vars[k].z_grid);
        }
      }

      return;
    }
# endif

    status = nc_create(filename, NC_CLOBBER, cdfid);
     if (status != NC_NOERR) {
      strcpy(message,"Opening "); strcat(message,filename);
      handle_error(message,status);
    }

# ifdef PARALLEL_IO
    {
      status = nc_put_att_int (*cdfid, NC_GLOBAL, "NumFilesInSet",
                                NC_INT, 1, &procs_used);
      if (status != NC_NOERR) handle_error("Att NumFilesInSet",status);
    }
# endif

    status = nc_set_fill(*cdfid, NC_NOFILL, &junk);         /* set nofill */
    if (status != NC_NOERR) handle_error("Set Fill",status);

/*   Define the coordinate variables.                                 */
    {
      int lathid, lonhid, latqid, lonqid;
      int layerid, lay2id, intid;
      int dayid, dayavid = 0, monthid = 0;

      for (k=0;k<novars;k++) 
        if ((vars[k].hor_grid == 'h') || (vars[k].hor_grid == 'u')) {
          use_lath = 1; break;
        }

      if (use_lath == 1) {
        status = nc_def_dim(*cdfid, "lath", yhsize, &lathid);
        if (status != NC_NOERR) handle_error("Dim lath",status);
        status = nc_def_var (*cdfid, "lath", NC_DOUBLE, 1, &lathid, &lathvid);
        if (status != NC_NOERR) handle_error("Var lath",status);
        status = nc_put_att_text(*cdfid, lathvid, "long_name", 9, "Latitude");
        if (status != NC_NOERR) handle_error("Name lath",status);
        status = nc_put_att_text(*cdfid, lathvid, "units", 8, "degrees");
        if (status != NC_NOERR) handle_error("Units lath",status);
# ifdef REENTRANT_Y
        status = nc_put_att_text(*cdfid, lathvid, "modulo", 1, "");
        if (status != NC_NOERR) handle_error("modulo lath",status);
# endif
# if defined(PARALLEL_IO) && defined(PARALLEL_Y)
        domain[0] = 1; domain[1] = NYTOT; 
        domain[2] = Y0abs+1; domain[3] = Y0abs+yhsize;
        status = nc_put_att_int(*cdfid, lathvid, "domain_decomposition",
                                NC_INT, 4, domain);
        if (status != NC_NOERR) handle_error("domain lath",status);
# endif
      }

      for (k=0;k<novars;k++) 
        if ((vars[k].hor_grid == 'q') || (vars[k].hor_grid == 'v')) {
          use_latq = 1; break;
      }

      if (use_latq == 1) {
        status = nc_def_dim(*cdfid, "latq", yqsize, &latqid);
        if (status != NC_NOERR) handle_error("Dim latq",status);
        status = nc_def_var (*cdfid, "latq", NC_DOUBLE, 1, &latqid, &latqvid);
        if (status != NC_NOERR) handle_error("Var latq",status);
        status = nc_put_att_text(*cdfid, latqvid, "long_name", 9, "Latitude");
        if (status != NC_NOERR) handle_error("Name latq",status);
        status = nc_put_att_text(*cdfid, latqvid, "units", 8, "degrees");
        if (status != NC_NOERR) handle_error("Units latq",status);
# ifdef REENTRANT_Y
        status = nc_put_att_text(*cdfid, latqvid, "modulo", 1, "");
        if (status != NC_NOERR) handle_error("modulo latq",status);
# endif
# if defined(PARALLEL_IO) && defined(PARALLEL_Y)
        domain[0] = 1+yhsize-yqsize; domain[1] = NYTOT; 
        domain[2] = Y0abs+1+yhsize-yqsize; domain[3] = Y0abs+yhsize;
        status = nc_put_att_int(*cdfid, latqvid, "domain_decomposition",
                                NC_INT, 4, domain);
        if (status != NC_NOERR) handle_error("domain latq",status);
# endif
      }

      for (k=0;k<novars;k++) 
        if ((vars[k].hor_grid == 'h') || (vars[k].hor_grid == 'v')) {
          use_lonh = 1; break;
        }

      if (use_lonh == 1) {
        status = nc_def_dim(*cdfid, "lonh", xhsize, &lonhid);
        if (status != NC_NOERR) handle_error("Dim lonh",status);
        status = nc_def_var (*cdfid, "lonh", NC_DOUBLE, 1, &lonhid, &lonhvid);
        if (status != NC_NOERR) handle_error("Var lonh",status);
        status = nc_put_att_text(*cdfid, lonhvid, "long_name", 10, "Longitude");
        if (status != NC_NOERR) handle_error("Name lonh",status);
        status = nc_put_att_text(*cdfid, lonhvid, "units", 8, "degrees");
        if (status != NC_NOERR) handle_error("Units lonh",status);
# ifdef REENTRANT
        status = nc_put_att_text(*cdfid, lonhvid, "modulo", 1, "");
        if (status != NC_NOERR) handle_error("modulo lonh",status);
# endif
# if defined(PARALLEL_IO) && defined(PARALLEL_X)
        domain[0] = 1; domain[1] = NXTOT; 
        domain[2] = X0abs+1; domain[3] = X0abs+xhsize;
        status = nc_put_att_int(*cdfid, lonhvid, "domain_decomposition",
                                NC_INT, 4, domain);
        if (status != NC_NOERR) handle_error("domain lonh",status);
# endif
      }

      for (k=0;k<novars;k++) 
        if ((vars[k].hor_grid == 'q') || (vars[k].hor_grid == 'u')) {
          use_lonq = 1; break;
        }

      if (use_lonq == 1) {
        status = nc_def_dim(*cdfid, "lonq", xqsize, &lonqid);
        if (status != NC_NOERR) handle_error("Dim lonq",status);
        status = nc_def_var (*cdfid, "lonq", NC_DOUBLE, 1, &lonqid, &lonqvid);
        if (status != NC_NOERR) handle_error("Var lonq",status);
        status = nc_put_att_text(*cdfid, lonqvid, "long_name", 10, "Longitude");
        if (status != NC_NOERR) handle_error("Name lonq",status);
        status = nc_put_att_text(*cdfid, lonqvid, "units", 8, "degrees");
        if (status != NC_NOERR) handle_error("Units lonq",status);
# ifdef REENTRANT
        status = nc_put_att_text(*cdfid, lonqvid, "modulo", 1, "");
        if (status != NC_NOERR) handle_error("modulo lonq",status);
# endif
# if defined(PARALLEL_IO) && defined(PARALLEL_X)
        domain[0] = 1+xhsize-xqsize; domain[1] = NXTOT; 
        domain[2] = X0abs+1+xhsize-xqsize; domain[3] = X0abs+xhsize;
        status = nc_put_att_int(*cdfid, lonqvid, "domain_decomposition",
                                NC_INT, 4, domain);
        if (status != NC_NOERR) handle_error("domain lonq",status);
# endif
      }

      for (k=0;k<novars;k++) if (vars[k].z_grid == 'L') {
        use_layer = 1; break;
      }
      if (use_layer == 1) {
        status = nc_def_dim(*cdfid, "Layer", (size_t) NZ, &layerid);
        if (status != NC_NOERR) handle_error("Dim Layer",status);
        status = nc_def_var (*cdfid, "Layer", NC_INT, 1, &layerid, &layervid);
        if (status != NC_NOERR) handle_error("Var Layer",status);
        status = nc_put_att_text(*cdfid, layervid, "units", 6, "Layer");
        if (status != NC_NOERR) handle_error("Units layer",status);
      }

      for (k=0;k<novars;k++) if (vars[k].z_grid == 'i') {
          use_int = 1; break;
      }
      if (use_int == 1) {
        status = nc_def_dim(*cdfid, "Interface", (size_t) (NZ+1), &intid);
        if (status != NC_NOERR) handle_error("Dim Interface",status);
        status = nc_def_var (*cdfid, "Interface", NC_FLOAT, 1, &intid, &intvid);
        if (status != NC_NOERR) handle_error("Var Interface",status);
        status = nc_put_att_text(*cdfid, intvid, "units", 10, "Interface");
        if (status != NC_NOERR) handle_error("Units Interface",status);
      }

      for (k=0;k<novars;k++) if (vars[k].z_grid == '2') {
          use_layer2 = 1; break;
      }
      if (use_layer2 == 1) {
        status = nc_def_dim(*cdfid, "lay2", (size_t) 2, &lay2id);
        if (status != NC_NOERR) handle_error("Dim Lay2",status);
        status = nc_def_var (*cdfid, "lay2", NC_INT, 1, &lay2id, &layer2vid);
        if (status != NC_NOERR) handle_error("Var Lay2",status);
        status = nc_put_att_text(*cdfid, layer2vid, "long_name", 6, "Layer");
        if (status != NC_NOERR) handle_error("Name Lay2",status);
        status = nc_put_att_text(*cdfid, layer2vid, "units", 6, "Layer");
        if (status != NC_NOERR) handle_error("Units Lay2",status);
      }


      *timeid = -1;
      for (k=0;k<novars;k++) if (vars[k].t_grid != '1') {
        status = nc_def_dim(*cdfid, "Time", NC_UNLIMITED, &dayid);
        if (status != NC_NOERR) handle_error("Dim Time",status);
        status = nc_def_var (*cdfid, "Time", NC_DOUBLE, 1, &dayid, &dayvid);
        if (status != NC_NOERR) handle_error("Var Time",status);
# ifndef TIMEUNIT
        status = nc_put_att_text(*cdfid, dayvid, "units", 5, "days");
# else
        if (TIMEUNIT < 1800.0)
          status = nc_put_att_text(*cdfid, dayvid, "units", 8, "seconds");
        else if ((TIMEUNIT >= 1800.0) && (TIMEUNIT < 4800.0))
          status = nc_put_att_text(*cdfid, dayvid, "units", 6, "hours");
        else if ((TIMEUNIT >= 4800.0) && (TIMEUNIT < 200000.0))
          status = nc_put_att_text(*cdfid, dayvid, "units", 5, "days");
        else
          status = nc_put_att_text(*cdfid, dayvid, "units", 6, "years");
# endif
        if (status != NC_NOERR) handle_error("Units Time",status);

        dayavid = dayid; /* CORRECT THIS LINE LATER. */
        monthid = dayid; /* CORRECT THIS LINE LATER. */

        *timeid = dayvid;
        break;
      }

      for (k=0;k<novars;k++) if (vars[k].t_grid == 'a') {
        status = nc_def_var (*cdfid, "Time_St", NC_DOUBLE, 1, &dayid, &daystvid);
        if (status != NC_NOERR) handle_error("Var Time",status);
        status = nc_put_att_text(*cdfid, daystvid, "long_name", 31, 
                   "Start time of averaging period");
        if (status != NC_NOERR) handle_error("Name Time_St",status);
# ifndef TIMEUNIT
        status = nc_put_att_text(*cdfid, daystvid, "units", 5, "days");
# else
        if (TIMEUNIT < 1800.0)
          status = nc_put_att_text(*cdfid, daystvid, "units", 8, "seconds");
        else if ((TIMEUNIT >= 1800.0) && (TIMEUNIT < 4800.0))
          status = nc_put_att_text(*cdfid, daystvid, "units", 6, "hours");
        else if ((TIMEUNIT >= 4800.0) && (TIMEUNIT < 200000.0))
          status = nc_put_att_text(*cdfid, daystvid, "units", 5, "days");
        else
          status = nc_put_att_text(*cdfid, daystvid, "units", 6, "years");
# endif
        if (status != NC_NOERR) handle_error("Units Time",status);

        timeid[1] = daystvid;

        break;
      }

#ifdef OUTPUT_DEFS
      status = nc_def_var(*cdfid, "KDML", NC_FLOAT, 0, &dummy, &kdmlvid);
      if (status != NC_NOERR) handle_error("KDML",status);
      status = nc_put_att_text(*cdfid, kdmlvid, "long_name", 30, 
			       "Isopycnal diffusion coeff. ML");
      if (status != NC_NOERR) handle_error("Name KDML",status);

      status = nc_def_var(*cdfid, "KD", NC_FLOAT, 0, &dummy, &kdvid);
      if (status != NC_NOERR) handle_error("KD",status);
      status = nc_put_att_text(*cdfid, kdvid, "long_name", 27, 
			       "Isopycnal diffusion coeff.");
      if (status != NC_NOERR) handle_error("Name KD",status);

      status = nc_def_var(*cdfid, "KHTR", NC_FLOAT, 0, &dummy, &khtrvid);
      if (status != NC_NOERR) handle_error("KHTR",status);
      status = nc_put_att_text(*cdfid, khtrvid, "long_name", 27, 
			       "Diapycnal diffusion coeff.");
      if (status != NC_NOERR) handle_error("Name KHTR",status);

      status = nc_def_var(*cdfid, "r_bio_tau", NC_DOUBLE, 0, &dummy, &tauvid);
      if (status != NC_NOERR) handle_error("r_bio_tau",status);
      status = nc_put_att_text(*cdfid, tauvid, "long_name", 15, 
			       "Restoring time");
      if (status != NC_NOERR) handle_error("Name r_bio_tau",status);
#endif /* OUTPUT_DEFS */


      for (k=0;k<novars;k++) {
        int dims[4], nodims = 0;
	nc_type mem_sz;

        switch (vars[k].t_grid) {
          case 's': dims[0] = dayid; varinfo[k].count[0] = 1; nodims = 1; break;
          case 'a': dims[0] = dayavid; varinfo[k].count[0] = 1; nodims = 1; break;
          case 'm': dims[0] = monthid; varinfo[k].count[0] = 1; nodims = 1; break;
          case '1': nodims = 0; break;
          default: fprintf(stderr,"Unknown time axis[%d] %c.\n",k,vars[k].t_grid);
        }
        switch (vars[k].z_grid) {
          case 'L': dims[nodims] = layerid; varinfo[k].count[nodims] = NZ;
                    varinfo[k].z_size = NZ; nodims++; break;
          case 'i': dims[nodims] = intid; varinfo[k].count[nodims] = NZ+1;
                    varinfo[k].z_size = NZ+1; nodims++; break;
          case '2': dims[nodims] = lay2id; varinfo[k].count[nodims] = 2;
                    varinfo[k].z_size = 2; nodims++; break;
          case '1': varinfo[k].z_size = 1; break;
          default: fprintf(stderr,"Unknown layer axis[%d] %c.\n",k,vars[k].z_grid);
        }
        switch (vars[k].hor_grid) {
          case 'u': dims[nodims] = lathid; varinfo[k].count[nodims] = yhsize;
                    dims[nodims+1] = lonqid; varinfo[k].count[nodims+1] = xqsize;
                    varinfo[k].x_size = xqsize; varinfo[k].y_size = yhsize;
                    nodims+=2; break;
          case 'v': dims[nodims] = latqid; varinfo[k].count[nodims] = yqsize;
                    dims[nodims+1] = lonhid; varinfo[k].count[nodims+1] = xhsize;
                    varinfo[k].x_size = xhsize; varinfo[k].y_size = yqsize;
                    nodims+=2; break;
          case 'h': dims[nodims] = lathid; varinfo[k].count[nodims] = yhsize;
                    dims[nodims+1] = lonhid; varinfo[k].count[nodims+1] = xhsize;
                    varinfo[k].x_size = xhsize; varinfo[k].y_size = yhsize;
                    nodims+=2; break;
          case 'q': dims[nodims] = latqid; varinfo[k].count[nodims] = yqsize;
                    dims[nodims+1] = lonqid; varinfo[k].count[nodims+1] = xqsize;
                    varinfo[k].x_size = xqsize; varinfo[k].y_size = yqsize;
                    nodims+=2; break;
          case '1': break;
          default: fprintf(stderr,"Unknown horizontal axes[%d] %c.\n",k,vars[k].hor_grid);
        }

	if (output_prec_float)
	    mem_sz = NC_FLOAT;
	else if (vars[k].mem_size == 'd') mem_sz = NC_DOUBLE;
        else mem_sz = NC_FLOAT;
        status = nc_def_var(*cdfid, vars[k].name, mem_sz, nodims, 
                            dims, &varinfo[k].id);
        if (status != NC_NOERR) {
          strcpy(message,"Var "); strcat(message,vars[k].name);
          handle_error(message,status);
        }

        status = nc_put_att_text(*cdfid, varinfo[k].id, "long_name", 
                                 strlen(vars[k].longname)+1, vars[k].longname);
        if (status != NC_NOERR) {
          strcpy(message,"Name "); strcat(message,vars[k].name);
          handle_error(message,status);
        }
        status = nc_put_att_text(*cdfid, varinfo[k].id, "units", 
                                 strlen(vars[k].units)+1, vars[k].units);
        if (status != NC_NOERR) {
          strcpy(message,"Units "); strcat(message,vars[k].name);
          handle_error(message,status);
        }
        status = nc_put_att_double(*cdfid, varinfo[k].id, "missing_value", 
                                 NC_DOUBLE, 1, &(vars[k].mval));
        if (status != NC_NOERR) {
          strcpy(message,"missing_value "); strcat(message,vars[k].name);
          handle_error(message,status);
        }
      }

    }

    status = nc_enddef(*cdfid);  /*leave define mode*/
    if (status != NC_NOERR) {
      strcpy(message,"End def "); strcat(message,filename);
      handle_error(message,status);
    }


/*  Write all of the axis variables.                                  */
    {
      double *lathptr, *latqptr, *lonhptr, *lonqptr;
# if (!defined(PARALLEL_Y) || defined(PARALLEL_IO))
      lathptr = &lath[Y1]; latqptr = &latq[Y0];
# else
      lathptr = &lath_tot[Y1]; latqptr = &latq_tot[Y0];
# endif
# if (!defined(PARALLEL_X) || defined(PARALLEL_IO))
      lonhptr = &lonh[X1]; lonqptr = &lonq[X0];
# else
      lonhptr = &lonh_tot[X1]; lonqptr = &lonq_tot[X0];
# endif

      if (use_lath == 1) {
        status = nc_put_var_double(*cdfid, lathvid, lathptr);
        if (status != NC_NOERR) handle_error("Put lath",status);
      }

      if (use_latq == 1) {
# ifdef REENTRANT_Y
        status = nc_put_var_double(*cdfid, latqvid, latqptr+1);
# else
        status = nc_put_var_double(*cdfid, latqvid, latqptr);
# endif
        if (status != NC_NOERR) handle_error("Put latq",status);
      }

      if (use_lonh == 1) {
        status = nc_put_var_double(*cdfid, lonhvid, lonhptr);
        if (status != NC_NOERR) handle_error("Put lonh",status);
      }

      if (use_lonq == 1) {
# ifdef REENTRANT
        status = nc_put_var_double(*cdfid, lonqvid, lonqptr+1);
# else
        status = nc_put_var_double(*cdfid, lonqvid, lonqptr);
# endif
        if (status != NC_NOERR) handle_error("Put lonq",status);
      }
    }

    for (k=0;k<NZ;k++) layer[k] = k+1;
    if (use_layer == 1) {
      status = nc_put_var_int(*cdfid, layervid, layer);
      if (status != NC_NOERR) handle_error("Put layer",status);
    }
  
    if (use_layer2 == 1) {
      status = nc_put_var_int(*cdfid, layer2vid, layer);
      if (status != NC_NOERR) handle_error("Put layer2",status);
    }

    if (use_int == 1) {
      for (k=0;k<=NZ;k++) interface[k] = (float) k + 0.5;
      status = nc_put_var_float(*cdfid, intvid, interface);
      if (status != NC_NOERR) handle_error("Put interface",status);
    }

#ifdef OUTPUT_DEFS
      status = nc_put_var_float(*cdfid, kdmlvid, &kdml);
      if (status != NC_NOERR) handle_error("put KDML",status);

      status = nc_put_var_float(*cdfid, kdvid, &kd);
      if (status != NC_NOERR) handle_error("put KD",status);

      status = nc_put_var_float(*cdfid, khtrvid, &khtr);
      if (status != NC_NOERR) handle_error("put KHTR",status);

      status = nc_put_var_double(*cdfid, tauvid, &r_bio_tau);
      if (status != NC_NOERR) handle_error("put r_bio_tau",status);
#endif /* OUTPUT_DEFS */
  }
#endif

}

void close_file(int *cdfid, FILE **fileptr) {
#ifdef NETCDF_OUTPUT
  if (*cdfid > -4) {
    if (*cdfid >= 0) {
      int status;
      status = nc_close(*cdfid);
      if (status != NC_NOERR) 
        handle_error("Closing file",status);
      *cdfid = -1;
    }
  }
  else
#endif
  {
    if (*fileptr != NULL) fclose(*fileptr);
    *fileptr = NULL;
  }
}

void handle_error (char message[], int status) {
  fprintf(stderr,"PE %d: %s. Error code %d.\n",pe_here,message,status);
#ifdef NETCDF_OUTPUT
  fprintf(stderr,nc_strerror(status)); fprintf(stderr,"\n");
#endif
  quit(status);
}


int name_output_file(char name[], double day, int noexdig,
                     char filepath[], char directory[])
{
/*   This subroutine opens a file with a name that is unique and      */
/* reflects the size and time of the output.  The first two arguments */
/* are the base names that will be used for double and float output   */
/* files.  The time of the output is day.  If noexdig is nonzero,     */
/* noexdig extra digits are added to the time stamp.                  */
/*   This subroutine returns a pointer to the opened file and sets    */
/* the value addressed by floatoutput to 1 for a float file and 0 for */
/* a double file.                                                     */

  char filename[50];        /* The full name of the save file.        */
  char format[43];          /* The format for the file name.          */
  int j, numfig;
  size_t namesize;
#if defined(SAVE_DOUBLE) || defined(NETCDF_OUTPUT)
  const int outsize = 8;
#else
  const int outsize = sizeof(float);
#endif

  if (day > 1e-6) {
    numfig = (int) (floor(log10(1e10*day)) + 1e-3) -
             (int) (floor(log10(1e10*(double)SAVEINT)) + 1e-3) + noexdig;
    numfig = (numfig > 8) ? 8 : ((numfig > 2) ? numfig : 2);
  }
  else numfig = 2;

  if ((outsize != 4) && (outsize != 8)) {
    for (j=0;j<((int) strlen(name));j++) {
      if (islower(name[j])) name[j] = toupper(name[j]);
      else if (isupper(name[j])) name[j] = tolower(name[j]);
    }
  }
  namesize = strlen(name);

/*   The name of the file indicates whether it contains 4-byte or     */
/* 8-byte numbers, by the way the time of the file is written.  It    */
/* also indicates the horizontal size of the fields.  For example, an */
/* 8-byte save file at time 1 will be "save1.00e00.50.80" for 50 x 80 */
/* fields, while the corresponding 4-byte file is "save1.00e0.50.80". */

  sprintf(format,"%s%%%d.%de.%%03d.%%03d",name,numfig+5,numfig);
  sprintf(filename,format,day,NXTOT,NYTOT);
  if (filename[numfig+namesize+3]=='+') {
    for (j=numfig+namesize+3;j<=48;j++) filename[j] = filename[j+1];
    if ((outsize == 4) && (filename[numfig+namesize+3]=='0')) {
      for (j=numfig+namesize+3;j<=47;j++) filename[j] = filename[j+1];
    }
  }
  else if ((outsize == 4) && (filename[numfig+namesize+4]=='0')) {
    for (j=numfig+namesize+4;j<=47;j++) filename[j] = filename[j+1];
  }

  strcpy(filepath, directory);
  strcat(filepath, filename);

#ifdef NETCDF_OUTPUT
  return NETCDF_FILE;
#else
  if (outsize == 8) return 0;
  else return FLOAT_FILE;
#endif

}

void find_input_file(char name[], char name2[], char fltname[], double day,
                     char directory[], FILE **file, int *cdfid,
                     int *timeid)
{
  char infile[150], inpath[300], inname[100];
  char format[30];
  int i, j, err = 1, namesize, try, try_nm, strip_0, flt_try;

/* Read the velocities and thicknesses in from the saved binary files.*/

/* Any non-netcdf file that is opened while stripping the extra 0 out */
/* of the exponent in the file name or while using fltname is assumed */
/* to contain 4-byte binary numbers.  Otherwise 8-byte binary numbers */
/* are in a non-netcdf file.  Any of the 3 name arguments may be      */
/* empty, in which case it is skipped.                                */

  for (try=0;try<=2;try++) {
    for (try_nm=0;try_nm<5;try_nm++) {
      if (try_nm==0) {strcpy(inname, name); strip_0 = 0; flt_try = 0;}
      else if (try_nm==1) {strcpy(inname, name2); strip_0 = 0; flt_try = 0;}
      else if (try_nm==2) {strcpy(inname, name); strip_0 = 1; flt_try = 0;}
      else if (try_nm==3) {strcpy(inname, name2); strip_0 = 1; flt_try = 0;}
      else {strcpy(inname, fltname); strip_0 = 0; flt_try = 1;}

      namesize = strlen(inname); if (namesize == 0) continue;

      for (i=8;i>=2;i--) {
        sprintf(format,"%s%%%d.%de.%%03d.%%03d",inname,i+5,i);
        sprintf(infile,format,day,NXTOT,NYTOT);
        if (infile[i+namesize+3]=='+') {
          for (j=i+namesize+3;j<=148;j++) infile[j] = infile[j+1];
          if (strip_0 && (infile[i+namesize+3]=='0')) {
            for (j=i+namesize+3;j<=47;j++) infile[j] = infile[j+1];
          }
        }
        else if (strip_0 && (infile[i+namesize+4]=='0')) {
          for (j=i+namesize+4;j<=147;j++) infile[j] = infile[j+1];
        }
        strcpy(inpath, directory); strcat(inpath, infile);
        err = open_input_file(inpath,file,cdfid,timeid);
        if (err == 0) {
          if ((flt_try || strip_0) && (*cdfid == -8)) *cdfid = -4;
          break;
        }

        strcat(inpath, ".cdf");
        err = open_input_file(inpath,file,cdfid,timeid);
        if (err == 0) break;
      }

      if (err == 0) break;
    }

    if (err == 0) break;
    else day += day*1.0e-10*(1.0-6.0*(double)try+3.0*(double)(try*try));
  }

  if (err != 0) {
    printf("Unable to find input file.\n");
    printf("Last tried '%s'.\n",inpath);
    exit(-1);
  }
  else printf("Opened file '%s'.\n",inpath);
}


int open_input_file(char filename[], FILE **fileptr, int *cdfid, int *timeid) {

  int status, return_val = 0, i, j = 0;

#ifdef NETCDF_OUTPUT
  *fileptr = NULL;
  status = nc_open(filename, NC_NOWRITE, cdfid);
  if (status != NC_NOERR)
    printf("status after opening %s: %d\n", filename, status);

/*   The indices of all variables starting with "Time" are returned   */
/* in timeid.                                                         */
  for (i=0;;i++) {
    int l;
    char varname[NC_MAX_NAME];

    if (NC_NOERR != nc_inq_varname(*cdfid,i,varname)) break;
    for (l=0;l<(int)strlen(varname);l++) varname[l] = tolower(varname[l]);

    if (strncmp(varname,"time",4) == 0 || strncmp(varname,"month",5) == 0) {
      timeid[0] = i;
      //hf timeid[j] = i; j++; //hf: bug! timeid is a scalar in calls!
    }
  }

  if (status != NC_NOERR)
#endif
  {
    char *e;

    if (pe_here == 0) {
      *fileptr = fopen(filename,"rb"); // hf: this used to be "r"
      if (*fileptr == NULL) return_val = 1;
    }
    else *fileptr = NULL;

/*    If a binary file is opened, cdfid indicates whether the file    */
/*  contains 4-byte or 8-byte numbers.  This is determined from the   */
/*  name of the file, as generated by name output file.  The time     */
/*  of the file is written in the name in a way that contains this    */
/*  information.  For example, a 8-byte save file at time 1 will be   */
/*  "save1.00e00.50.80" for 50 x 80 fields, while the corresponding   */
/*  4-byte file is "save1.00e0.50.80".                                */
    e = strrchr(filename,'e');

/*
if ((*(e+2) == '.') && ((*(e+1) >= '0') && (*(e+1) >= '0'))) *cdfid= -4;
*/

    if (((e != NULL) && (e+3-filename<strlen(filename))) &&
        ( ( (*(e+2) == '.') && (*(e+1) >= '0') && (*(e+1) <= '9') ) ||
          ( (*(e+1) == '-') && (*(e+3) == '.') &&
            (*(e+2) >= '0') && (*(e+2) <= '9') ) ) ) *cdfid = -4; 


    else *cdfid = -8;

    *timeid = 0;
  }

  return return_val;
}

void read_field(int cdfid, FILE *fileptr, char varname[], 
                int size_x, int size_y, int size_z,
                int start_file_x, int start_file_y, int start_file_z,
                int start_array_x, int start_array_y, int start_array_z,
                int recno, double variable[])
{
  char message[500];

#ifdef NETCDF_OUTPUT
  if (cdfid > -4) {
    int i, j, status, varid, ndims, dims_known, dimids[4];
    int xdim, ydim, zdim, tdim, varstart;
    size_t dimlen[4];
    size_t start[8], count[8];
    ptrdiff_t imap[8];
    char dimname[7][20];
    char message[144];

    status = nc_inq_varid (cdfid, varname, &varid);
    if (status != NC_NOERR) {
      strcpy(message,"Unable to get variable ID for ");strcat(message,varname);
      handle_error(message,status);
    }

    status = nc_inq_varndims(cdfid, varid, &ndims);
    if (status != NC_NOERR) handle_error("Number of dimensions", status);
    if (ndims > 4) handle_error("Too many dimensions: ", ndims);

    status = nc_inq_vardimid(cdfid, varid, dimids);
    if (status != NC_NOERR) handle_error("Dimension ids", status);

/*  This is the COARDS order of dimensions.                           */
    if (ndims == 4) {tdim = 0; zdim = 1; ydim = 2; xdim = 3;}
    else if (ndims == 3) {tdim = 4; zdim = 0; ydim = 1; xdim = 2;}
    else if (ndims == 2) {tdim = 4; zdim = 5; ydim = 0; xdim = 1;}
    else if (ndims == 1) {tdim = 4; zdim = 0; ydim = 6; xdim = 7;}
    else {tdim = 0; zdim = 0; ydim = 0; xdim = 0;}

    for (i=0;i<ndims;i++) {
      status = nc_inq_dim(cdfid, dimids[i], dimname[i], &dimlen[i]);
      if (status != NC_NOERR) handle_error("Dimension name", status);
      
      for (j=0;j<(int)strlen(dimname[i]);j++) dimname[i][j] = tolower(dimname[i][j]);

      if ((!strncmp(dimname[i],"lon",3)) || (dimname[i][0] == 'x'))
        {xdim = i; if (tdim == i) tdim = 4; if (ydim == i) ydim = 6; if (zdim == i) zdim = 5;}
      if ((!strncmp(dimname[i],"lat",3)) || (dimname[i][0] == 'y'))
        {ydim = i; if (tdim == i) tdim = 4; if (xdim == i) xdim = 7; if (zdim == i) zdim = 5;}
      if ((!strncmp(dimname[i],"lay",3)) || (dimname[i][0] == 'z'))
        {zdim = i; if (tdim == i) tdim = 4; if (xdim == i) xdim = 7; if (ydim == i) ydim = 6;}
      if ((!strncmp(dimname[i],"tim",3)) || (dimname[i][0] == 't'))
        {tdim = i; if (xdim == i) xdim = 7; if (ydim == i) ydim = 6; if (zdim == i) zdim = 5;}
    }

    dims_known = ((xdim < 4) ? 1 : 0) + ((ydim < 4) ? 1 : 0) +
                 ((zdim < 4) ? 1 : 0) + ((tdim < 4) ? 1 : 0);
    if (dims_known != ndims) {
      sprintf(message,"Unable to interpret dimensions - only %d known:\n",ndims);
      if (tdim < 6) strcat(message,dimname[tdim]); strcat(message,", ");
      if (zdim < 6) strcat(message,dimname[zdim]); strcat(message,", ");
      if (ydim < 6) strcat(message,dimname[ydim]); strcat(message,", ");
      if (xdim < 6) strcat(message,dimname[xdim]); strcat(message,", ");
      strcat(message,"\nOut of:  ");
      for (i=0;i<ndims;i++) {
        strcat(message,dimname[i]); strcat(message,", ");
      }
      strcat(message,"\nFor variable "); strcat(message,varname);
      handle_error(message, varid);
    }

/* At this point the dimensions of the NetCDF field are known.        */


/* Check that the length of the dimensions are sensible.              */
    if ((NXTOT+1) < size_x)
      handle_error("X-dimension to be read too large.",varid);
    if ((NYTOT+1) < size_y)
      handle_error("Y-dimension to be read too large.",varid);
    if ((NZ+1) < size_z)
      handle_error("Z-dimension to be read too large.",varid);

/*   size_x = NXTOT or NXTOT+1, and size_y = NYTOT or NYTOT+1 mean    */
/*  that entire dimesions are read based on the size in the file.     */
    if ((size_x == NXTOT) || (size_x == NXTOT+1)) {
      size_x = dimlen[xdim]; start_array_x = NXTOT+1-size_x;
    }
    if ((size_y == NYTOT) || (size_y == NYTOT+1)) {
      size_y = dimlen[ydim]; start_array_y = NYTOT+1-size_y;
    }

    count[zdim] = (size_z > 0) ? size_z : 1;   count[tdim] = 1;

    start[zdim] = start_file_z;   start[tdim] = recno;

    imap[ydim] = 1;
    imap[xdim] = ((size_y == 0) ? 1 : NYMEM);
    imap[zdim] = ((size_x == 0) ? 1 : NXMEM) * imap[xdim];
    imap[tdim] = NZ*NXMEM*NYMEM;

/*   This block accounts for the possibility that this is a parallel  */
/* simulation.  It makes sure that the current processor only reads   */
/* its own data.                                                      */
    {
      double start_x, end_x, start_y, end_y;

      varstart = start_array_z*imap[zdim];

      if (size_x == 0) {count[xdim] = 1; start[xdim] = start_file_x;}
      else {
        start_x = (start_array_x-X0abs > nx-X0) ? nx+1-X0 :
                   start_array_x-X0abs;
        if (start_x < 0) start_x = 0;
        end_x = (size_x-1+start_array_x-X0abs >= nx-X0) ? nx-X0 : 
                 size_x-1+start_array_x-X0abs;
        if (end_x < 0) end_x = -1;

        count[xdim] = (end_x < start_x) ? 0 : end_x - start_x + 1;
        start[xdim] = start_x + X0abs - start_array_x + start_file_x;
        varstart += (X0+start_x)*imap[xdim];
      }

      if (size_y == 0) {count[ydim] = 1; start[ydim] = start_file_y;}
      else {
        start_y = (start_array_y-Y0abs > ny-Y0) ? ny+1-Y0 :
                   start_array_y-Y0abs;
        if (start_y < 0) start_y = 0;
        end_y = (size_y-1+start_array_y-Y0abs >= ny-Y0) ? ny-Y0 : 
                 size_y-1+start_array_y-Y0abs;
        if (end_y < 0) end_y = -1;

        count[ydim] = (end_y < start_y) ? 0 : end_y - start_y + 1;
        start[ydim] = start_y + Y0abs - start_array_y + start_file_y;
        varstart += Y0+start_y;
      }

    }

    if (count[xdim]*count[ydim] > 0) {
      
      status = nc_get_varm_double(cdfid, varid, start, count, NULL, 
                                  imap, &variable[varstart]);
      if (status != NC_NOERR) {
        sprintf(message,"Reading field %s",varname);
        handle_error(message,status);
      } 
   }

  } else
#endif
  {
/*  Here a field is read from a binary file.                          */
    int k, flt;
    size_t err = 1;
    
    flt = (cdfid == -4) ? 1 : 0;

    if (start_array_z+size_z > NZ) {
      strcpy(message,"Vertical offset + size too large to read ");strcat(message,varname);
      handle_error(message,-1);
    }

    if (((size_x == NXTOT) || (size_x == NXTOT+1)) &&
        ((size_y == NYTOT) || (size_y == NYTOT+1)))
    {
      for (k=start_array_z;k<size_z+start_array_z;k++) 
        err *= read_layer(variable+k*NXMEM*NYMEM,fileptr,flt);
    }
    else if (((size_x == 0) && (size_y == 0))) {
      if (pe_here == 0) err *= fread((void *)variable,8*size_z,1,fileptr);
      spread_double_vector(variable,size_z);
    }
    else {
      strcpy(message,"Only full layers or single columns are supported with binary files.  Variable ");
      strcat(message,varname);
      handle_error(message,-1);
    }
  }
}

size_t read_time(int cdfid, FILE *fileptr, int timeid, size_t nrec, 
                 double *day)
{
#ifdef NETCDF_OUTPUT
  if (cdfid > -4) {
    if (cdfid >= 0) {
      int status;
      char message[144];

      status = nc_get_var1_double(cdfid, timeid, &nrec, day);
      if (status != NC_NOERR) {
        sprintf(message,"Reading time %d",timeid);
        handle_error(message,status);
      }
    }
    return 1;
  }
  else
#endif
  {
    size_t err = 1;
    if (pe_here==0) err *= fread((void *)day,8,1,fileptr);
    return err;
  }
}
