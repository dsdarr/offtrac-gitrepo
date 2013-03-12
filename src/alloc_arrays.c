#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "init.h"
#include "io.h"
#include "alloc.h"
#include "par_IO.h"


int alloc_arrays()
{
    extern double ***u;
    extern double ***v;
    extern double ***mn_u;
    extern double ***mn_v;
    extern double ***mn_h;
    extern double ***uhtm;
    extern double ***vhtm;
    extern double ***ea;
    extern double ***eb;
    extern double ***mn_rml;
    extern double ***mn_uhtm;
    extern double ***mn_vhtm;
    extern double ***mn_ea;
    extern double ***mn_eb;
#ifdef AGE
    extern double ***age;
    extern double ***mn_age;
#endif
    /*#ifdef SMOOTH_REST
    extern double ***h_init;
    #endif*/


    if (! (u = alloc3d(NZ,NXMEM,NYMEM))) handle_error("u");
    if (! (v = alloc3d(NZ,NXMEM,NYMEM))) handle_error("v");
    if (! (mn_u = alloc3d(NZ,NXMEM,NYMEM))) handle_error("mn_u");
    if (! (mn_v = alloc3d(NZ,NXMEM,NYMEM))) handle_error("mn_v");
    if (! (mn_h = alloc3d(NZ,NXMEM,NYMEM))) handle_error("mn_h");
    uhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(uhtm == NULL) {
	fprintf(stderr,"not enough memory for uhtm!\n");
	return 1;
    }
    vhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(vhtm == NULL) {
	fprintf(stderr,"not enough memory for vhtm!\n");
	return 1;
    }
    ea = alloc3d(NZ,NXMEM,NYMEM);
    if(ea == NULL) {
	fprintf(stderr,"not enough memory for ea!\n");
	return 1;
    }    eb = alloc3d(NZ,NXMEM,NYMEM);
    if(eb == NULL) {
	fprintf(stderr,"not enough memory for eb!\n");
	return 1;
    }
    mn_rml = alloc3d(2,NXMEM,NYMEM);
    if(mn_rml == NULL) {
	fprintf(stderr,"not enough memory for mn_rml!\n");
	return 1;
    }
    mn_uhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_uhtm == NULL) {
	fprintf(stderr,"not enough memory for mn_uhtm!\n");
	return 1;
    }
    mn_vhtm = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_vhtm == NULL) {
	fprintf(stderr,"not enough memory for mn_vhtm!\n");
	return 1;
    }
    mn_ea = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_ea == NULL) {
	fprintf(stderr,"not enough memory for mn_ea!\n");
	return 1;
    }
    mn_eb = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_eb == NULL) {
	fprintf(stderr,"not enough memory for mn_eb!\n");
	return 1;
    }
    
#ifdef AGE
    age = alloc3d(NZ,NXMEM,NYMEM);
    if(age == NULL) {
	fprintf(stderr,"not enough memory for age!\n");
	return 1;
    }
    mn_age = alloc3d(NZ,NXMEM,NYMEM);
    if(mn_age == NULL) {
	fprintf(stderr,"not enough memory for mn_age!\n");
	return 1;
    }
#endif /* AGE */

    /*#ifdef SMOOTH_REST
    h_init = alloc3d(NZ,NXMEM,NYMEM);
    if(h_init == NULL) {
	fprintf(stderr,"not enough memory for h_init!\n");
	return 1;
    }
    #endif*/


  return 0;
}


int alloc_tracadv(int NZED, int NY, int NX, double ***vartrac)
{

    vartrac = alloc3d(NZED,NY,NY);
    if(vartrac == NULL) {
	fprintf(stderr,"not enough memory for tracer variables!\n");
	return 1;
    }

  return 0;
}

int alloc_error(char* arr_name)
{
  fprintf(stderr,"not enough memory for %s!\n", arr_name);
  quit(2);
  return 2;
}
