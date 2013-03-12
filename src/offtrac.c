/********+*********+*********+*********+*********+*********+*********+*
 *                                                                    *
 *                OFFTRAC - off-line tracer advection code            *
 *                                                                    *
 *                    David Darr - April/May 2002                     *
 *                                                                    *
 *                    OCMIP biogeochemistry added                     *
 *                     cdeutsch 2004                                  *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <netcdf.h>

#include "init.h"
#define MAINFILE
#include "metrics.h"
#include "io.h"
#include "util.h"
#include "alloc_arrays.h"
#include "alloc_trac.h"
#include "initialize.h"
#include "iocdf.h"
#include "read.h"

/*-------------------------------------------------------------------*
 *                                                                   *
 *     define variables and subroutines
 *                                                                   *
 *-------------------------------------------------------------------*/

/*  To add a variable, increase NOVARS (in init.h)  */
/*  then add information to struct vardesc (below) */
//BX  31AUG07  added missing value - last fields
//BX  31OCT07  Changed some variables to double precision (needed for restart)
//BX           All variables were float ('f') before
//BX           This is a bug fix as memory restrictions do not allow to create
//BX           a separate vardesc field for restarts.

//BX-a
//BX  missing values are set here and to the same value in "vars" declaration
//BX  they will be written on their respective mn_ fields with mult_.._mv calls
//HF the parameter and the macro need to be the same!
const double misval = -1.e+6;
#define MISVAL -1.e+6

const double r_bio_tau = 1.0 / (30.0 * 86400.0); // restoring time step


//BX-e
struct vardesc vars[NOVARS] =
	{
		{
			"D", "Basin Depth", 'h', '1', '1', "meter", 'f', 0
		}, // 0
		{
			"mn_u",
			"Zonal Velocity",
			'u',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 1
		{
			"mn_v",
			"Meridional Velocity",
			'v',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 2
		{
			"mn_h", "Layer Thickness", 'h', 'L', 's', "meter", 'd', 9.e-10
		}, // 3
		{
			"mn_uhtm",
			"zonal vel - cumul",
			'u',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 4
		{
			"mn_vhtm",
			"merid vel - cumul",
			'v',
			'L',
			's',
			"meter second-1",
			'f',
			0
		}, // 5
		{
			"mn_ea",
			"downward entrainment",
			'h',
			'L',
			's',
			"meters",
			'd',
			0
		}, // 6
		{
			"mn_eb", "upward entrainment", 'h', 'L', 's', "meters", 'd', 0
		}, // 7
		{
			"mn_eaml",
			"downward ML detrain",
			'h',
			'1',
			's',
			"meters",
			'd',
			0
		}, // 8
		{
			"mn_age",
			"age tracer",
			'h',
			'L',
			's',
			"years",
			'd',
			-3.17097930758233E-14
		}// 9

	};

void alloc_fields(void);
void set_metrics(void);

/* Begin edit DT */
void step_fields(int iyear, int itts, int imon, int iterno); // modified ashao
/*  End edit DT  */
void initialize_sponge(int from_file, char filename[]);

// mn_ refers to mean value for output interval, WRINT (see init.h)

double ***u;
double ***v;
double h[NZ][NXMEM][NYMEM];
double hend[NZ][NXMEM][NYMEM];
#ifdef HTEST
double htest[NZ][NXMEM][NYMEM];
#endif
double hstart[NZ][NXMEM][NYMEM];
#if defined AGE2 || defined AGE3
//  for debugging
//double hnew[NZ][NXMEM][NYMEM];
double ***hnew;
#endif
double depth[NZ][NXMEM][NYMEM];
double rml[2][NXMEM][NYMEM];
double ***mn_u;
double ***mn_v;
double ***mn_h;
double ***mn_rml;

double ***uhtm;
double ***vhtm;
double ***ea;
double ***eb;
double eaml[NXMEM][NYMEM];
double ***mn_uhtm;
double ***mn_vhtm;
double ***mn_ea;
double ***mn_eb;
double mn_eaml[NXMEM][NYMEM];
double wd[NZ + 1][NXMEM][NYMEM];

double Temptm[NZ][NXMEM][NYMEM];
double Salttm[NZ][NXMEM][NYMEM];
double salt_woa[NXMEM][NYMEM];
double windspeed[NXMEM][NYMEM];
double xkw[NXMEM][NYMEM];
double fice[NXMEM][NYMEM];
double atmpres[NXMEM][NYMEM];

double areagr[NXMEM][NYMEM];
double D[NXMEM][NYMEM];
double ****tr;
double trintegral[NTR];
double trprofile[NZ];
double trwarn[NTR][2];

double umask[NXMEM][NYMEM];
double vmask[NXMEM][NYMEM];

int ksub[NXMEM][NYMEM];
int ksub_clim[NXMEM][NYMEM];

int month;
int mlen[12];
int icnt;

double qlat[NYMEM];
double hlat[NYMEM];

#ifdef ENTRAIN
double ea_init[NZ][NXMEM][NYMEM];
double eaml_init[NXMEM][NYMEM];
#endif

/* Begin added DT */
int beginyear;
int theyear;
int lastsave;
/*  End added DT  */


#ifdef AGE
//HF  double age[NZ][NXMEM][NYMEM];
double ***age;
double age_init[NZ][NXMEM][NYMEM];
double ***mn_age;
int mAGE;
#endif

#ifdef SMOOTH_REST
//double ***h_init;
double h_init[NZ][NXMEM][NYMEM];
#endif



// begin ashao
double currtime;
int iterno;
#ifdef HINDCAST
	int usehindcast = 1;
#else
	int usehindcast = 0;
#endif

int hindindex = 0;
// Calculate what iteration to start/stop reading hindcast
int starthindindex = (BEGHIND-BEGYEAR)*12; //cmon added later
int numhindmonths = (ENDHIND-BEGHIND+1)*12;
// end ashao

double k0, k1, k2, kb, kw, ks, kf, k1p, k2p, k3p, ksi;
double ff, htotal, htotal2, xco2, xacc, x1, x2;
double bt, st, ft, sit, pt, dic1, ta, co2starair;

#ifdef SPONGE
double Idamp[NXMEM][NYMEM];
double num_fields;
double sp_e[NZ][NXMEM][NYMEM];
# if defined(AGE)
double sp_age[NZ][NXMEM][NYMEM];
# endif
#endif

double dt;

double junk[(NZ + 1) * (NXMEM) * (NYMEM)];

double *var[NOVARS];
long varsize[NOVARS];
int flags[NOVARS];
int rflags[NOVARS];
char directory[75];
char fname[75];
char inittrac[200];
char kw_varname[100]; //ashao: Specify the name of the kw variable in the gasx input file
#ifdef WRTTS
struct varcdfinfo varinfo[NOVARS];
int varmap[NOVARS];
int itts; /* tracer time step counter */
FILE *fn;
char output_filename[200];
double wrts;
int nvar = 0, cdfid, timeid[2];
size_t nrec = 0;
#endif

/*-------------------------------------------------------------------*
 *
 *     begin main executable
 *
 *-------------------------------------------------------------------*/

int main(void)
    {
    double inmon, tmon;
    double mon, nxt, lst;
    double *iyr;
    double *nyr;
    double dyr, day;
    double ndyr;
    double *dy;
    double dmon[12];
    double dbmon;
    int imon, inxt, ilst;
# ifndef WRTTS
    int itts; /* tracer time step counter */
# endif
    int nmnfirst;
#ifdef SEPFILES
    double smon, snxt;
    int ismon, isnxt, ihnxt;
#endif

#ifndef WRTTS
    size_t nrec = 0;
#endif
    int err, i, j, k;
    int cmon;
    int nmn;
    double frac;
    static int m;
#ifndef WRTTS
    int varmap[NOVARS];

    FILE *fn;
    char output_filename[200];
#endif
    char run_name[200];
    char restart_filename[200];
    struct vardesc var_out[NOVARS];
#ifndef WRTTS
    struct varcdfinfo varinfo[NOVARS];
    int nvar = 0, cdfid, timeid[2];
#endif

    extern int flags[NOVARS];
    extern int rflags[NOVARS];

    //BX-a  for testing only
    int status;
    char message[144];
    //BX-e

    //BX  allocate tracer fields
    err = alloc_arrays();
    if (err)
	printf("Error allocating arrays.\n");

    err = alloc_trac();
    if (err)
	printf("Error allocating tracer field.\n");

    iyr = malloc(sizeof(double));
    nyr = malloc(sizeof(double));
    dy = malloc(sizeof(double));

    mlen[0] = 31; /* January      */
    mlen[1] = 28; /* February     */
    mlen[2] = 31; /* March        */
    mlen[3] = 30; /* April        */
    mlen[4] = 31; /* May          */
    mlen[5] = 30; /* June         */
    mlen[6] = 31; /* July         */
    mlen[7] = 31; /* August       */
    mlen[8] = 30; /* September    */
    mlen[9] = 31; /* October      */
    mlen[10] = 30; /* November     */
    mlen[11] = 31; /* December     */

    dmon[0] = 0.0;

    for (i = 1; i <= 11; i++)
	{
	dmon[i] = dmon[i - 1] + mlen[i - 1];
	}

    /*----------------------------------*
     *
     *     get user input
     *
     *----------------------------------*/

    {
    char line[100];
    int scan_count, done = 1;

    printf("Enter directory to use to read.\n");
    fgets(directory, sizeof(directory), stdin);

    directory[strlen(directory) - 1] = '\0';
    k = strlen(directory);
    if (k > 0)
	if (directory[k - 1] != '/')
	    {
	    directory[k] = '/';
	    directory[k + 1] = '\0';
	    printf("Using directory %s first.\n", directory);
	    }

    strcat(directory, fname);
    strcpy(fname, directory);
    printf("file path = %s\n", fname);

    while (done)
	{

	printf(
		"\nEnter the starting month and the total months to integrate.\n");

	fgets(line, sizeof(line), stdin);
	line[strlen(line) - 1] = '\0';
	scan_count = sscanf(line, "%lg, %lg,", &inmon, &tmon);
	if (scan_count == 2)
	    {
	    if (inmon < 0 || tmon < 0)
		printf("Negative values not allowed\n");
	    else
		done = 0;
	    }
	else
	    printf("Incorrect number of values, %d, read.\n", scan_count);
	}

    printf("\ninitial month = %g \n", inmon);
    printf("final month = %g \n", inmon + tmon - 1);
    printf("total months = %g \n\n", tmon);

    /*-----------------------------------
     *
     *     set output print flags to 0
     *
     *     added restart flags as a restart bug fix until
     *     memory restriction problem is solved 31OCT07 BX
     *
     *----------------------------------*/

    for (i = 0; i <= NOVARS - 1; i++)
	flags[i] = 0;
    for (i = 0; i <= NOVARS - 1; i++)
	rflags[i] = 0;

    flags[1] = 0;
    flags[2] = 0; /* u,v */
    rflags[1] = 0;
    rflags[2] = 0; /* u,v */
    flags[0] = 0;
    flags[3] = 0; /* D, h */
    rflags[0] = 0;
    rflags[3] = 0; /* D, h */
    flags[4] = 0;
    flags[5] = 0; /* uhtm, vhtm */
    rflags[4] = 0;
    rflags[5] = 0; /* uhtm, vhtm */
    flags[6] = 0;
    flags[7] = 0;
    flags[8] = 0; /* ea, eb, eaml */
#ifdef ENTRAIN
    rflags[6]=1; rflags[7]=1; rflags[8]=1; /* ea, eb, eaml */
#endif
    flags[18] = 0; /* ML potential density */
    rflags[18] = 0; /* ML potential density */
#ifdef AGE
    flags[9] = 1;
    rflags[9] = 1; /* ideal age tracer*/
#endif

    printf("Enter base name for output\n");

    fgets(run_name, sizeof(run_name), stdin);
    run_name[strlen(run_name) - 1] = '\0';
    sprintf(output_filename, "%s.%04g.nc", run_name, inmon + tmon - 1);
    printf("Create NetCDF output file '%s'.\n", output_filename);



    } // end of block "get user input"

    //DT
    lastsave = -1;
    //DT

    /*-----------------------------------
     *
     *     allocate and initialize fields
     *
     *----------------------------------*/

    read_grid();
    printf("Done reading grid or metric file.\n");

    set_metrics();
    printf("Done setting metrics.\n");

    /* Copy the variable descriptions to a list of the actual output variables. */
    for (i = 0; i < NOVARS; i++)
	if (flags[i] > 0)
	    {
	    var_out[nvar] = vars[i];
	    varmap[i] = nvar;
	    nvar++;
	    }
    // force float precision output with last argument
    printf("Making NETCDF %s file\n", output_filename);
    create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
	    timeid, varinfo, 1);
    // don't force
    // create_file(output_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid, timeid, varinfo, 0);
    printf("Closing file \n");
    close_file(&cdfid, &fn);

    /* Allocate the memory for the fields to be calculated.		*/
    alloc_fields();

    /* initialize tracer pointers					*/

    for (m = 0; m < NOVARS; m++)
	{
	if (flags[m])
	    for (k = 0; k < varsize[m]; k++)
		var[m][k] = 0.0;
	}

    /********** set means to zero                                   */
    /*  2d variables */
    if (flags[8])
	set_fix_darray2d_zero(mn_eaml);


    /*  3d variables */
    if (flags[1])
	set_darray3d_zero(mn_u, NZ, NXMEM, NYMEM);
    if (flags[2])
	set_darray3d_zero(mn_v, NZ, NXMEM, NYMEM);
    if (flags[3])
	set_darray3d_zero(mn_h, NZ, NXMEM, NYMEM);
    if (flags[4])
	set_darray3d_zero(mn_uhtm, NZ, NXMEM, NYMEM);
    if (flags[5])
	set_darray3d_zero(mn_vhtm, NZ, NXMEM, NYMEM);
    if (flags[6])
	set_darray3d_zero(mn_ea, NZ, NXMEM, NYMEM);
    if (flags[7])
	set_darray3d_zero(mn_eb, NZ, NXMEM, NYMEM);
#ifdef AGE
    if (flags[9])
	set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);

#endif
    printf("Reading bathymetry, D.\n");

    // initialize D to be all ocean first

    for (i = 0; i <= NXMEM - 1; i++)
	{
	for (j = 0; j <= NYMEM - 1; j++)
	    {
	    D[i][j] = MINIMUM_DEPTH;
	    }
	}
#ifndef USE_CALC_H
    printf("calling read_D\n");
    read_D();

    for (i = 0; i <= NXMEM - 1; i++)
	for (j = 0; j <= NYMEM - 1; j++)
	    if (D[i][j] < 10.0)
		D[i][j] = MINIMUM_DEPTH;
#endif

    read_grid();
    set_metrics();

    dyr = inmon / 12;
    smon = (double) ((int) inmon % NMONTHS);
    mon = 12 * modf(dyr, iyr);

    printf("\n initial mon = %g - %g - %g \n\n", inmon, smon, mon);


    imon = (int) (mon + 0.00001);
    /* Begin edit DT */
    int iyear = floor((inmon + imon) / 12);
    theyear = iyear;
#ifdef RESTART
    theyear++;
    iyear++;
#endif
    /*  End edit DT  */
    inxt = imon + 1;

    ismon = (int) (smon + 0.00001) % NMONTHS;

    isnxt = (ismon+1) % NMONTHS;
    ihnxt = (ismon+2) % NMONTHS;

    dbmon = (double) inmon;
    lst = 12 * modf((dbmon - 1 + 12) / 12, iyr);
    ilst = (int) (lst + 0.00001);

    // ashao: Read in next month's isopycnal thickness fields
    // (will be copied at the beginning of the integration)
    // Done this way so that if hindcasts are used the physical fields switch smoothly
    // to/from climatological fields
    //BX  files are not in regular order (uvw are midmonth and h starts with last month)

        currtime = BEGYEAR;
    if (usehindcast) {
        // Check to see if simulation started within the hindcast years
        if ( (currtime >= BEGHIND) || (currtime < (ENDHIND+1) ) ) {
        	hindindex=inmon;
        	read_h(ismon,hend,"hind");
        }
    }
    else {
    	read_h(ismon, hend,"clim");
    }
    // for files in regular order (h before uvw) use code here
    //BX
    // read_uvw(imon,1);
    //BX
    // read_h(imon,inxt);



#ifdef USE_CALC_H
    z_sum(h, D);
#endif
    printf("\nSetting up and initializing\n");

#ifdef RESTART
    initialize(inmon,run_name);
#else
    initialize(imon);
#endif
    nmn = 0;

    //HF the next line should be in init.h (and be optional!)
#undef OUTPUT_IC
#ifdef OUTPUT_IC
    /*-----------------------------------
     *
     *     write tracer initial conditions
     *
     *----------------------------------*/

    printf("Writing initial condition variables out to netCDF\n\n",cmon);

    //  open netcdf file for each writing
    status = nc_open(output_filename, NC_WRITE, &cdfid);
    if (status != NC_NOERR)
	{
	strcpy(message,"Opening file"); strcat(message,output_filename);
	handle_error(message,status);
	}

    err = write_time(cdfid,fn,timeid[0],nrec, dy);
    if (err == -1) printf("Error writing day.\n");

    if (flags[3])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_h[k][i][j] += h[k][i][j];
	}

#ifdef AGE
    if (flags[9]) add_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_age[k][i][j] += age[k][i][j];
     }*/
#endif /* AGE */


    for (m=0;m<NOVARS;m++) if (flags[m])
	{
	err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]],
		nrec,var[m]);
	if (err == -1) printf("Error writing %s.\n",vars[m].name);
	}
    //  close file after each writing
    close_file(&cdfid, &fn);
    printf("netcdf record = %d\n",nrec++);
#endif /* OUTPUT_IC */




    /*-------------------------------------------------------------------*
     *
     *     begin integration loop
     *
     *-------------------------------------------------------------------*/

    mon = 0.0; /* reiniti */
    nmnfirst = 1;


    for (cmon = inmon; cmon < inmon + tmon; cmon++)
	{
	nmn++;

	dyr = cmon / 12.0;
	ndyr = (cmon + 1) / 12.0;

	dbmon = (double) cmon;
	lst = 12 * modf((dbmon - 1 + 12) / 12, iyr);
	ilst = (int) (lst + 0.00001);
	printf("double-mon=%g,lastmon=%g,ilst=%i\n", dbmon, lst, ilst);


	smon = (double) ((int) cmon % NMONTHS);
	snxt = (double) ((int) (cmon + 1) % NMONTHS);

	mon = 12.0 * modf(dyr, iyr);
	nxt = 12.0 * modf(ndyr, nyr);

	printf("the current month is %i-%g-%g \n", cmon, smon, mon);
	printf("the current year is %g \n", *iyr);

	imon = (int) (mon + 0.00001);
	inxt = (int) (nxt + 0.00001);

	ismon = (int) (smon + 0.00001);
	isnxt = (int) (snxt + 0.00001);

	day = (*iyr) * 365.0 + dmon[imon];
	*dy = day;

	printf("the current day is -%g- \n", *dy);
	printf("the current ismon/smon is -%i-%g- \n", ismon, smon);
	printf("the next smonth/mean index: -%i- \n", isnxt);

	dt = ((double) mlen[imon]);
	dt = dt * 86400 / (double) NTSTEP;


	for (itts = 1; itts <= NTSTEP; itts++)
	    {
	    printf("\nSub time step number= %i \n", itts);
	    /*-----------------------------------
	     *
	     *     get physical fields and bc's
	     *
	     *----------------------------------*/

	    printf("Month %d timestamp Start: %4.2f End: %4.2f\n",cmon,currtime,currtime+dt/31536000);
	    currtime += dt / 31536000; //ashao: advance currenttime
	    copy_2fix_darray3d(hstart,hend,NZ,NXMEM,NYMEM);
	    copy_2fix_darray3d(h,hstart,NZ,NXMEM,NYMEM);
	    if (usehindcast) {
	        if ( ( cmon >= starthindindex) && ( hindindex <= (numhindmonths-1) )) {
			
			printf("Reading in UVW from hindcast\n");
	        	read_uvw(hindindex,"hind");
		        read_h(hindindex+1,hend,"hind");
	        	hindindex++;
	        }
		else {
		
		printf("Reading in UVW from climatology\n");
	    	read_uvw(isnxt,"clim");
	    	read_h(isnxt, hend,"clim");
		}
	    }
	    else {
		printf("Reading in UVW from climatology\n");
	    	read_uvw(isnxt,"clim");
	    	read_h(isnxt, hend,"clim");
	    }
	    printf("Month- hstart:%d hend:%d\n",ismon,isnxt);

	    /*-----------------------------------
	     *
	     *     integrate 1 time step
	     *
	     *----------------------------------*/

	    printf("step fields - day = %g\n\n", day);
	    step_fields(iyear, itts, imon, iterno); // modified ashao


	    /*-------------------------------------------------------------------*
	     *
	     *     end integration loop
	     *
	     *-------------------------------------------------------------------*/

	    /*-----------------------------------
	     *
	     *     calculate tracer averages
	     *
	     *----------------------------------*/

	    printf("calculate tracer averages\n");

	    if (flags[8])
		add_fix_darray2d(mn_eaml, eaml);
	    if (flags[1])
		add_darray3d(mn_u, u, NZ, NXMEM, NYMEM);
	    if (flags[2])
		add_darray3d(mn_v, v, NZ, NXMEM, NYMEM);
	    if (flags[3])
		{
		for (k = 0; k < NZ; k++)
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    {
			    mn_h[k][i][j] += h[k][i][j];
			    }
		}
	    if (flags[4])
		add_darray3d(mn_uhtm, uhtm, NZ, NXMEM, NYMEM);
	    if (flags[5])
		add_darray3d(mn_vhtm, vhtm, NZ, NXMEM, NYMEM);
	    if (flags[6])
		add_darray3d(mn_ea, ea, NZ, NXMEM, NYMEM);
	    if (flags[7])
		add_darray3d(mn_eb, eb, NZ, NXMEM, NYMEM);
#ifdef AGE
	    //DT
	    //    if (flags[9]) add_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
	    if (flags[9])
	    {
	    	for (k = 0; k < NZ; k++)
	    	{
	    		for (i = 0; i < NXMEM; i++)
	    		{
	    			for (j = 0; j < NYMEM; j++)
	    			{
	    				if (age[k][i][j] > 0.0)
	    				{
	    					{
	    						mn_age[k][i][j] += age[k][i][j];
	    					}
	    				}}}}}

	    				//DT
#endif
	    if (flags[18])
		{
		for (k = 0; k < 2; k++)
		    for (i = 0; i < NXMEM; i++)
			for (j = 0; j < NYMEM; j++)
			    mn_rml[k][i][j] += rml[k][i][j];
		}


	    /* calculate the mean */

	    if (nmn == WRINT && itts == 1)
		{
		printf("***nmn= %i, itts= %i, nmnfirst= %i\n", nmn, itts,
			nmnfirst);
		frac = 1.0 / ((double) nmn * (double) NTSTEP);

		if (flags[8])
		    mult_fix_darray2d(mn_eaml, frac);
		if (flags[1])
		    mult_darray3d(mn_u, NZ, NXMEM, NYMEM, frac);
		if (flags[2])
		    mult_darray3d(mn_v, NZ, NXMEM, NYMEM, frac);
		if (flags[3])
		    mult_darray3d(mn_h, NZ, NXMEM, NYMEM, frac);
		if (flags[4])
		    mult_darray3d(mn_uhtm, NZ, NXMEM, NYMEM, frac);
		if (flags[5])
		    mult_darray3d(mn_vhtm, NZ, NXMEM, NYMEM, frac);
		if (flags[6])
		    mult_darray3d(mn_ea, NZ, NXMEM, NYMEM, frac);
		if (flags[7])
		    mult_darray3d(mn_eb, NZ, NXMEM, NYMEM, frac);
#ifdef AGE
		if (flags[9])
		    mult_darray3d(mn_age, NZ, NXMEM, NYMEM, frac);
#endif

		if (flags[18])
		    mult_darray3d_mv(mn_rml, 2, NXMEM, NYMEM, frac, D, misval);

		/*-----------------------------------
		 *
		 *     write tracer averages and reset to 0
		 *
		 *----------------------------------*/
		printf("Writing month %i variables out to netCDF\n\n", cmon);

		status = nc_open(output_filename, NC_WRITE, &cdfid);
		if (status != NC_NOERR)
		    {
		    strcpy(message, "Opening file");
		    strcat(message, output_filename);
		    handle_error(message, status);
		    }

		err = write_time(cdfid, fn, timeid[0], nrec, dy);
		if (err == -1)
		    printf("Error writing day.\n");

		for (m = 0; m < NOVARS; m++)
		    if (flags[m]==1)
			{
			printf("m = %d\n",m);
			err = write_field(cdfid, fn, vars[m],
				varinfo[varmap[m]], nrec, var[m]);
			if (err == -1)
			    printf("Error writing %s.\n", vars[m].name);
			}
		close_file(&cdfid, &fn);

		/*    reset all means to zero */
		nmn = 0;
		if (flags[8])
		    set_fix_darray2d_zero(mn_eaml);
		if (flags[1])
		    set_darray3d_zero(mn_u, NZ, NXMEM, NYMEM);
		if (flags[2])
		    set_darray3d_zero(mn_v, NZ, NXMEM, NYMEM);
		if (flags[3])
		    set_darray3d_zero(mn_h, NZ, NXMEM, NYMEM);
		if (flags[4])
		    set_darray3d_zero(mn_uhtm, NZ, NXMEM, NYMEM);
		if (flags[5])
		    set_darray3d_zero(mn_vhtm, NZ, NXMEM, NYMEM);
		if (flags[6])
		    set_darray3d_zero(mn_ea, NZ, NXMEM, NYMEM);
		if (flags[7])
		    set_darray3d_zero(mn_eb, NZ, NXMEM, NYMEM);
		if (flags[18])
		    set_darray3d_zero(mn_rml, 2, NXMEM, NYMEM);
#ifdef AGE
		if (flags[9])
		    set_darray3d_zero(mn_age, NZ, NXMEM, NYMEM);

#endif
		// begin ashao
		// end ashao

		printf("netcdf record = %d\n", nrec + 1);
		nrec++;

		} /*  end if nmn==WRITEINT */

	    /*-------------------------------------------------------------------*
	     *
	     *     end integration loop
	     *
	     *-------------------------------------------------------------------*/

	    } /* end itts loop over NTSTEP */
	} /* end while */

    //BX-a
    /*-----------------------------------
     *
     *     write restart file
     *
     *----------------------------------*/

    //   First write the up-to-date tracer values to the mean fields.
    //   In order to save memory the instantaneous values for the
    //   restart will be written on mean fields in the restart file.

    printf("start of array copying\n");

#ifdef ENTRAIN
    copy_fix_darray2d(mn_eaml, eaml);
#endif
#ifdef ENTRAIN
    copy_darray3d(mn_ea, ea, NZ, NXMEM, NYMEM);
    //    copy_darray3d(mn_eb, Salttm, NZ, NXMEM, NYMEM);
#endif
    //HF #ifdef SMOOTH_REST
    copy_fix_darray3d(mn_h, h, NZ, NXMEM, NYMEM);
    //HF #endif
#ifdef AGE
    copy_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
#endif
    for (k = 0; k < NZ; k++)
	{
	for (i = 0; i < NXMEM; i++)
	    {
	    for (j = 0; j < NYMEM; j++)
		{
		} /* for j loop */
	    } /* for i loop */
	} /* for k loop */

    //  Second, create restart file name and open file

    sprintf(restart_filename, "restart.%s.%04d.nc", run_name, cmon);
    printf("Writing NetCDF restart file '%s'.\n\n", restart_filename);

    /* Copy the variable descriptions to a list of the actual restart variables. */
    nvar = 0;
    for (i = 0; i < NOVARS; i++)
	if (rflags[i] > 0)
	    {
	    var_out[nvar] = vars[i];
	    varmap[i] = nvar;
	    nvar++;
	    }

    // do NOT force float precision output with last argument
    create_file(restart_filename, NETCDF_FILE, var_out, nvar, &fn, &cdfid,
	    timeid, varinfo, 0);

    for (m = 0; m < NOVARS; m++)
	if (rflags[m])
	    {
	    err = write_field(cdfid, &fn, vars[m], varinfo[varmap[m]], 0,
		    var[m]);
	    if (err == -1)
		printf("Error writing %s.\n", vars[m].name);
	    }

    close_file(&cdfid, &fn);

    printf("\n programme termine normallement. \n");
    return (0);
    }

/*-------------------------------------------------------------------*
 *
 *     begin subroutines
 *                                                                   *
 *-------------------------------------------------------------------*/

#ifdef WRTTS
void write_ts(double wrts)
    {
    int i,j,k,m;
    int status;
    char message[144];
    int err;
    double *wrdy;
#ifdef WRTTS
    wrdy = malloc(sizeof(double));
#endif
    /********** set means to zero                                   */
    /*  2d variables */
    if (flags[8]) set_fix_darray2d_zero(mn_eaml);
    /*---------------------------------------
     *
     *     write tracer values on output field
     *
     *--------------------------------------*/

    printf("write output field\n");

    if (flags[8]) add_fix_darray2d(mn_eaml, eaml);
    //HF outer if, inner for loops, instead of vice versa!
    //HF think about writing appropriate subroutine(s)!!!
    if (flags[1]) copy_darray3d(mn_u, u, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_u[k][i][j] = u[k][i][j];
     } */
    if (flags[2]) copy_darray3d(mn_v, v, NZ, NXMEM, NYMEM);
    /* {
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_v[k][i][j] = v[k][i][j];
     }*/
    if (flags[3])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_h[k][i][j] = h[k][i][j];
	}
    if (flags[4]) copy_darray3d(mn_uhtm, uhtm, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_uhtm[k][i][j] = uhtm[k][i][j];
     }*/
    if (flags[5]) copy_darray3d(mn_vhtm, vhtm, NZ, NXMEM, NYMEM);
    /*{
     for (k=0;k<NZ;k++)
     for (i=0;i<NXMEM;i++)
     for (j=0;j<NYMEM;j++)
     mn_vhtm[k][i][j] = vhtm[k][i][j];
     }*/
    if (flags[6])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_ea[k][i][j] = ea[k][i][j];
	}
    if (flags[7])
	{
	for (k=0;k<NZ;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_eb[k][i][j] = Salttm[k][i][j];
	}
#ifdef AGE
    if (flags[9]) copy_darray3d(mn_age, age, NZ, NXMEM, NYMEM);
#endif
    if (flags[18])
	{
	for (k=0;k<2;k++)
	    for (i=0;i<NXMEM;i++)
		for (j=0;j<NYMEM;j++)
		    mn_rml[k][i][j] = rml[k][i][j];
	}
	    printf("Writing variables for sub timestep %i out to netCDF\n\n",itts);
	    //  open netcdf file for each writing
	    status = nc_open(output_filename, NC_WRITE, &cdfid);
	    if (status != NC_NOERR)
		{
		strcpy(message,"Opening file"); strcat(message,output_filename);
		handle_error(message,status);
		}
	    *wrdy = wrts;
	    err = write_time(cdfid,fn,timeid[0],nrec, wrdy);
	    if (err == -1) printf("Error writing day.\n");

	    for (m=0;m<NOVARS;m++) if (flags[m])
		{
		err = write_field(cdfid, fn, vars[m], varinfo[varmap[m]],
			nrec,var[m]);
		if (err == -1) printf("Error writing %s.\n",vars[m].name);
		}
	    close_file(&cdfid, &fn);

	    printf("netcdf record = %d\n",nrec++);
	    }
#endif //WRTTS
	void alloc_fields(void)
	    {

	    int m;

	    extern double *var[NOVARS];
	    extern double junk[(NZ + 1) * (NXMEM) * (NYMEM)];
	    extern long varsize[NOVARS];
	    extern int flags[NOVARS];

	    for (m = 0; m < NOVARS; m++)
		{

		//HF if ( m>3 )    {
		//HF: added mn_h
		if (m >= 3)
		    {
		    switch (vars[m].z_grid)
			{
		    case 'L':
			varsize[m] = NZ * (NXMEM) * (NYMEM);
			break;
		    case 'i':
			varsize[m] = (NZ + 1) * (NXMEM) * (NYMEM);
			break;
		    case '2':
			varsize[m] = 2 * (NXMEM) * (NYMEM);
			break;
		    case '1':
			varsize[m] = (NXMEM) * (NYMEM);
			break;
		    default:
			printf("Unknown layer axis[%d] %c.\n", m, vars[m].z_grid);
			}

		    var[m] = calloc((size_t) varsize[m], sizeof(double));
		    if (var[m] == NULL)
			{
			printf("Unable to allocate memory for var[%d].\n", m);
			exit(-21);
			}

		    //	temporarily disabling this output
		    //      printf("Allocated memory for var[%d],%s.\n\n",m,vars[m].name);

		    }
		else
		    {
		    var[m] = junk;
		    varsize[m] = 0;
		    }
		}

	    var[0] = &D[0][0];
	    var[1] = &mn_u[0][0][0];
	    var[2] = &mn_v[0][0][0];
	    var[3] = &mn_h[0][0][0];
	    var[4] = &mn_uhtm[0][0][0];
	    var[5] = &mn_vhtm[0][0][0];
	    var[6] = &mn_ea[0][0][0];
	    var[7] = &mn_eb[0][0][0];
	    var[8] = &mn_eaml[0][0];
	    var[18] = &mn_rml[0][0][0];

#ifdef AGE
	    var[9] = &mn_age[0][0][0];

#endif



	    // end ashao

	    }
