#include "gpe_engine_api.h"
#include "gpe_engine.h"

void gpe_get_lattice_api(int *_nx, int *_ny, int *_nz)
{
    gpe_get_lattice(_nx, _ny, _nz);
}
int gpe_create_engine_api(double alpha, double beta, double dt, double npart)
{
    int ierr;
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    return ierr;
}
int gpe_set_user_params_api(int size, double *params)
{
    int ierr;
    gpe_exec( gpe_set_user_params(size, params), ierr );
    return ierr;
}
int gpe_set_psi_api(double t, double complex * psi)
{
    int ierr;
    gpe_exec( gpe_set_psi(t, (Complex *)psi), ierr );
    return ierr;
}
int gpe_normalize_psi_api()
{
    int ierr;
    gpe_exec( gpe_normalize_psi(), ierr );
    return ierr;
}
int gpe_energy_api(double *t, double *ekin, double *eint, double *eext)
{
    int ierr;
    gpe_exec( gpe_energy(t, ekin, eint, eext), ierr );
    return ierr;
}
int gpe_evolve_api(int nt)
{
    int ierr;
    gpe_exec( gpe_evolve(nt), ierr );
    return ierr;
}
int gpe_get_psi_api(double *t, double complex * psi)
{
    int ierr;
    gpe_exec( gpe_get_psi(t, (Complex *)psi), ierr );
    return ierr;
}
int gpe_destroy_engine_api()
{
    int ierr;
    gpe_exec( gpe_destroy_engine(), ierr );
    return ierr;
}