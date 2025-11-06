/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 31-07-2020
 * */

#include "wslda_potdens.h"

#ifndef __WSLDA_WRITEVARS__
#define __WSLDA_WRITEVARS__

int create_wdata_metadata(metadata_t *input, int datadim, double t0, double dt, int spinsymmetry, wdata_metadata *wdmd);
int clear_files(metadata_t *input, wdata_metadata *wdmd);
int write_wdata_metadata_file(metadata_t *input, wdata_metadata *wdmd, char *codename);
int write_measurments(wdata_metadata *wdmd, MPI_Comm mpi_comm, char *codetype, int it, wslda_density h_densities, wslda_potential h_potentials);
int write_measurments_subset(wdata_metadata *wdmd, MPI_Comm mpi_comm, char *codetype, int it, wslda_density h_densities);
int create_wtxt_file_for_wf(const char * prefix, int iogroup,
                           int nwf, int nx, int ny, int nz, double dx, double dy, double dz,
                           double kF, double *mu, double ec, double beta, int codedim
                           );
#ifdef TDWSLDA
void set_ptr_d_delta(void * ptr);
#endif

#endif
