
int inittype4_readdpca(char *inprefix, double *h_densities, double *h_potentials, double *eF);
double npart_TF(int nxyz, double *U, double mu, double dxyz);

// ZBRENT functions
#define ITMAX 100
typedef struct
{
    double mu[ITMAX+2]; // values of chemical potentials for each iteration
    double N[ITMAX+2]; // particle number for corresponding potential
    int status[ITMAX+2]; // status=0 - not computed, status=1 - request for computation, status=2 - computed
} zbrent_data_t;

void reset_zbrent_data(zbrent_data_t *zbrent_data);
double get_point_from_zbrent_data(zbrent_data_t *zbrent_data);
void set_point_in_zbrent_data(zbrent_data_t *zbrent_data, double value);
double zbrent2(zbrent_data_t *zbrent_data, double *x1, double *x2, double tol, int *ierr);
