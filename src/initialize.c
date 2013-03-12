/********+*********+*********+*********+*********+*********+*********+*
 *         Initialize                                                 *
 *                                                                    *
 ********+*********+*********+*********+*********+*********+*********+*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "init.h"
#include "offtrac.h"
#include "alloc.h"
#include "read.h"
#ifdef SMOOTH_REST
#include "step.h"
#endif

extern const double misval;


extern void initializemasks(void);

extern double ****tr;
extern double umask[NXMEM][NYMEM];   /* _mask are 1 over ocean and 0  */
extern double vmask[NXMEM][NYMEM];   /* over land on the u & v grids. */

extern double D[NXMEM][NYMEM];
extern double lonh[NXMEM],lath[NYMEM];
extern double dxh[NYMEM],dy;

extern double ***u;
extern double ***v;
extern double h[NZ][NXMEM][NYMEM];

extern double ***uhtm;
extern double ***vhtm;

#ifdef AGE
extern double age_init[NZ][NXMEM][NYMEM];
//extern double age[NZ][NXMEM][NYMEM];
extern double ***age;
extern int mAGE;
#endif

#ifdef SMOOTH_REST
//BXextern double ***h_init;
extern double h_init[NZ][NXMEM][NYMEM];
extern double hstart[NZ][NXMEM][NYMEM];
double tr_in[NZ];
double z_in[NZ];
double depth[NZ][NXMEM][NYMEM];
double hdepth[NZ][NXMEM][NYMEM];
double tmp_in;
#endif


/*   Begin added DT    */
extern int beginyear;
extern int theyear;
extern int lastsave;
extern int estTTD;
/*   End added DT      */


// end ashao




//BX-a for restart reading
#ifdef ENTRAIN
extern double ea_init[NZ][NXMEM][NYMEM];
extern double eaml_init[NXMEM][NYMEM];
#endif
//BX-e

#ifdef SPONGE
extern double Idamp[NXMEM][NYMEM];
extern int num_fields;
extern double sp_e[NZ][NXMEM][NYMEM];
# if defined (AGE) || defined(CFC) || defined(SF6) // ashao
extern double sp_age[NZ][NXMEM][NYMEM];
# endif
#endif

#ifdef RESTART
void initialize(int imon, char *run_name)
#else
void initialize(int imon)
#endif
{
/* Arguments: mon - Time of the start of the run in months.             */
/* Arguments: inmon - Time of the restart of the run in months.         */

    int i, j, k, m;

#ifdef SPONGE

    num_fields=1;  /* that e thing */
# ifdef AGE
    num_fields++;
# endif

// end ashao

    printf("num_fields = %d\n\n",num_fields);
    read_sponge();

    initialize_sponge(Idamp,num_fields);
    set_up_sponge_field(sp_e,NULL,NZ);
 
#endif  /* SPONGE */

    initializemasks();   /* create land/ocean masks based on topography array (D) */

    m=0; 

    //  printf("Initialize 173 for month %i.\n",imon);
    //  printf("Initialize 174 for run %s.\n",run_name);


#ifdef RESTART
    read_tracer_init(imon,run_name);  /* reading all tracer initial values from restart file*/
#else
    /*  for (i=0;i<=NXMEM-1;i++)
    for (j=0;j<=NYMEM-1;j++)
    printf("l202 ini D[%d][%d]=%g\n",i,j,D[i][j]); */

    read_tracer_init(imon);  /* reading all tracer initial values */
#endif
    // printf("Initialize 173 for month %i.\n",imon);

#ifdef SMOOTH_REST
    //BX interpolate restart files to new depth field associated 
    //BX with time step wrap around at the end of a climatological forcing file
    printf("Restart using smoothed transition of depth fields. \n");
    z_depth(h_init,depth);
    z_depth(hstart,hdepth);
    //printf("Initialize 199,   h_init= %g\n", h_init[20][90][146]);
    //printf("Initialize 199,   hstart= %g\n", hstart[20][90][146]);
    //printf("Initialize 199,    depth= %g\n", depth[20][90][146]);
    //printf("Initialize 199,   hdepth= %g\n", hdepth[20][90][146]);
    for (i=X1;i<=nx;i++) {
      for (j=Y1;j<=ny;j++) {
	  if (D[i][j]>MINIMUM_DEPTH) {
	      for (k=0;k<NZ;k++)
		  z_in[k]  = depth[k][i][j];
# ifdef ENTRAIN
	      for (k=0;k<NZ;k++) tr_in[k] = ea_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  ea_init[k][i][j] = tmp_in;
	      }
	      //BX not action neccessary for eaml_init - is only 2D
# endif
# ifdef AGE
	      for (k=0;k<NZ;k++) tr_in[k] = age_init[k][i][j];
	      for (k=0;k<NZ;k++) {
		  tmp_in = lin_interpp(hdepth[k][i][j], tr_in, 
		  				z_in, 0, NZ);
		  if (tmp_in < 0.e0) tmp_in = 0.;
		  age_init[k][i][j] = tmp_in;
	      }
# endif


//end ashao
	  } else {   //BX if(D>MINIMUM_DEPTH)
#ifdef WE_DONT_USE_THIS
	      for (k=0;k<NZ;k++) {
# ifdef ENTRAIN
		  ea_init[k][i][j] = misval;
# endif
# ifdef AGE
		  age_init[k][i][j] = misval;
# endif
}
#endif /* NOT USED */
	  }
      }
  }
#endif /* SMOOTH_REST */

    //BX printf("Initialize 239, dic_init= %g\n", dic_init[0][90][146]);


#ifdef AGE

    printf("initialize tracer internally\n\n");

    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<NZ;k++) {
		age_init[k][i][j] = 0.0;
		uhtm[k][i][j] = 0.0;
		vhtm[k][i][j] = 0.0;
		u[k][i][j] = 0.0;
		v[k][i][j] = 0.0;
	    }
	}
    }


    for (i=0;i<=NXMEM-1;i++) {
      for (j=0;j<=NYMEM-1;j++) {
	if (D[i][j] > MINIMUM_DEPTH) {
	  for (k=0; k<NZ; k++) {
	    age[k][i][j] = age_init[k][i][j];  /* in years */
	    tr[m][k][i][j] = age_init[k][i][j]*365.0*86400.0; /* in seconds */
	  }
	} else {
	  for (k=0; k<NZ; k++) {
	    age[k][i][j] = misval;
	    tr[m][k][i][j] = misval;
	  }                            
	}
      }
    }
    

# ifdef SPONGE
    set_up_sponge_field(sp_age,age,NZ);
# endif

    mAGE=m;
    printf("index of AGE tracer m = %d\n\n",mAGE);
    m++;
#endif   // AGE



// begin ashao
// end ashao



#ifdef MERGED_ML
    //BX merge 1st-2nd and 3rd-4th layer for auxiliary fields
    for (i=0;i<=NXMEM-1;i++) {
	for (j=0;j<=NYMEM-1;j++) {
	    for (k=0;k<=2;k=k+2) {
# ifdef AGE
	    age[k][i][j]       = (age[k][i][j]*h[k][i][j]+
				  age[k+1][i][j]*h[k+1][i][j])/
	      (h[k][i][j] + h[k+1][i][j]);
	    age[k+1][i][j]     = age[k][i][j];
# endif
// end ashao

	    }
	}
    }
# ifdef SMOOTH_REST
    //BX Needed here as fields have been altered
    merge_ml_tr();
# endif
#endif
    /* zonal, meridional re-entrance    */
    for (m=0;m<NTR;m++) {
	for (k=0;k<NZ;k++) {
	    for (j=0;j<=NYMEM-1;j++) {
		tr[m][k][0][j] = tr[m][k][nx-1][j];
		tr[m][k][1][j] = tr[m][k][nx][j];	
		tr[m][k][nx+1][j] = tr[m][k][2][j];
		tr[m][k][nx+2][j] = tr[m][k][3][j];
	    }
	}
	for (k=0;k<NZ;k++) {
	    for (i=2;i<=nx;i++) {
		tr[m][k][363-i][ny+1] = tr[m][k][i][ny];
		tr[m][k][363-i][ny+2] = tr[m][k][i][ny-1];
	    }
	}
    }
}
