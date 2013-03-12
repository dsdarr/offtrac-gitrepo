#include <stddef.h>
#include <stdio.h>

/* The following structure is used in various subroutines to contain  */
/* a number of run-time changeable parameters.                        */

struct params {double dt;      /* The time step, in s.                */
               double dtbt;    /* The barotropic time step, in s.     */
               double be;   /* A nondimensional number from 0.5 to 1. */
               double bebt; /* A nondimensional number from 0 to 1.   */
               double daymax;  /* The final day of the simulation.    */
               double daysave; /* The first day to save output.       */
               double saveint; /* The number of days between saves.   */
               int saves_per_file; /* The number of time leves to     */
                               /* save to each file.                  */
               double energysavedays; /* The interval between writing */
                               /* the energies of the run.            */
               char energyfile[50];  /* The name of the energy file.  */
               char restartfile[50]; /* The name of the restart file. */
               double restint; /* The interval between saves of the   */
                               /* the restart file.                   */
               char inputdir[50]; /* The directory in which NetCDF    */
                               /* input files might be found.         */
               double Kh;      /* The Lapacian horizontal viscosity   */
                               /* in m2 s-1.                          */
               double Ah;      /* The biharmonic horizontal viscosity */
                               /* in m4 s-1.                          */
               double Kh_vel_scale; /* This velocity times the grid   */
                               /* spacing gives the Laplacian         */
                               /* viscosity, in m s-1.                */
               double Ah_vel_scale; /* This velocity times the grid   */
                               /* spacing cubed gives the biharmonic  */
                               /* viscosity, in m s-1.                */
               double Smag_Lap_const; /* The nondimensional Lapacian  */
                               /* Smagorinsky constant.               */
               double Smag_bi_const; /* The nondimensional biharmonic */
                               /* Smagorinsky constant.               */
               double Khth;    /* The interface depth diffusivity in  */
                               /* m2 s-1.                             */
               int iord;       /* The number of iterations of MPDATA. */
               int ntstep;     /* The number of time steps between    */
                               /* tracer updates or diabatic forcing. */
               double Hmix;    /* The mixed layer thickness in m.     */
               double Kvml;    /* The mixed layer vertical viscosity  */
                               /* in m2 s-1.                          */
               double Kdml;    /* The diapycnal diffusivity in the    */
                               /* mixed layer in m2 s-1.              */
               double Kv;      /* The vertical viscosity in m2 s-1.   */
               double Kd;      /* The diapycnal diffusivity in m2 s-1.*/
               double Hbbl;    /* The static bottom boundary layer    */
                               /* thickness, in m.                    */
               double Kvbbl;   /* The linear viscosity in the bottom  */
                               /* boundary layer, in m2 s-1.          */
               double cdrag;   /* The quadratic drag coefficient.     */
               double drag_bg_vel; /* An assumed unresolved background*/
                               /* velocity for calculating the bottom */
                               /* drag, in m s-1.                     */
               double BBL_thick_min; /* The minimum bottom boundary   */
                               /* layer thickness, in m.  This might  */
                               /* be Kv / (cdrag * drag_bg_vel) to    */
                               /* give Kv as the minimum near-bottom  */
                               /* viscosity.                          */
               double KhTr;    /* The along-isopycnal tracer          */
                               /* diffusivity in m2 s-1.              */
               double maxvel;  /* Velocity components greater than    */
                               /* maxvel, in m s-1, are truncated.    */
               int maxtrunc;   /* The number of truncations per time  */
                               /* step at which the run is stopped.   */
               int save_aves;  /* If defined (i.e. nonzero), time     */
                               /* averaging is enabled.               */
               double maxcpu;  /* The maximum amount of cpu time per  */
                               /* processor for which HIM should run  */
                               /* before saving a restart file and    */
                               /* quiting with a return value that    */
                               /* indicates that a further execution  */
                               /* is required to complete the simu-   */
                               /* lation, in wall-clock seconds.      */
              };


/* Subroutines in offtrac.c */
void alloc_fields(void);

/* Subroutine in tracadv.c */
void tracer(int itts);
//BX-TS void tracer(void);

/* Subroutine in initialize.c */
#ifdef RESTART
void initialize(int imon, char *run_name);
#else
void initialize(int imon);
#endif

/* Subroutines in par_IO.c */
void collect_double_fields(double *field_in, int offx, int offy);
void spread_double_fields(double field_out[NXMEM][NYMEM]);
void spread_double_vector(double vec[], int numpts);
size_t write_layer(double field_in[], FILE *file, int float_vals);
size_t read_layer(double field_out[NXMEM][NYMEM], FILE *file, int float_vals);

/* Subroutine in set_metrics.c */
void set_metrics(void);

/* Subroutine in masks.c */
void initializemasks(void);

/* Subroutines in off_par.c */
void set_up_parallel(void);
void pass_vars(double var[][NXMEM][NYMEM], int nz, int pts_out, int nrows,
               double bounds, int dosync, int safemem, int dirflag);
void pass_var_2d(double var[NXMEM][NYMEM], int nz, int pts_out, int nrows,
                 double bounds, int dosync, int safemem, int dirflag);  
void spread_string(char name[], int length);
void sum_fields(double vec[], int numsum);
void sum_int_fields(int vec[], int numsum);
void check_field(double field_in[NXMEM][NYMEM], char write_string[],
                 int i0, int i1, int j0, int j1);
void quit(int status);

#define TO_NORTH 1
#define TO_SOUTH 2
#define TO_EAST 4
#define TO_WEST 8
