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
#ifdef TDWSLDA
void set_ptr_d_delta(void * ptr);
#endif

#endif
