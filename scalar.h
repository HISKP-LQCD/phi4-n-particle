#include "types.h"

int lattice_point(int t, int x, int y, int z);
int initialize_field(scalar_field *p_aux);
int copy_field_point(scalar_field *p_old, scalar_field *p_new, int t, int x, int y, int z);
int copy_field(scalar_field *p_old, scalar_field *p_new);
void update_field_point(scalar_field *p_aux, int t, int x, int y, int z);
int fprint_field(scalar_field phiaux, long long n_conf);
int fread_field(const char *filename, scalar_field *p_phi);

  
