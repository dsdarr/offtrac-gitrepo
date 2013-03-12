
void spread_double_vector(double vec[], int numpts);

void quit(int status);

size_t read_layer(double field_out[][NYMEM], FILE *file, int float_vals);

size_t write_layer(double field_in[], FILE *file, int float_vals);
