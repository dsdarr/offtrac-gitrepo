void read_var3d( char inpath[200], char varname[200], int imon, double*** data);
void read_fieldave(int imon);
void read_D();
void read_uvw(int imon, char *fieldtype);
void read_h(int imon,double (*hread)[NXMEM][NYMEM], char *fieldtype);
void read_fields(int imon,int itts);
void read_biotic_bc(int imon,int itts);
void read_ts(int imon,int itts);
void read_clim(int imon,int inxt,int ilst);
#ifdef RESTART
void read_tracer_init(int imon, char *run_name);
#else
void read_tracer_init(int imon);
#endif
void read_buoy(int imon);
void read_sponge(void);
void read_grid();
#ifdef LEV_OXY
void read_oxy_ic(void);
#endif
