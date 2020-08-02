/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 31-07-2020
 * */

#ifndef __WSLDA_WRITEVARS__
#define __WSLDA_WRITEVARS__

int create_wdata_metadata(metadata_t *input, int datadim, double t0, double dt, int spinsymmetry, wdata_metadata *wdmd);
int clear_files(metadata_t *input, wdata_metadata *wdmd);
int write_wdata_metadata_file(metadata_t *input, wdata_metadata *wdmd, char *codename);
int write_measurments(wdata_metadata *wdmd, MPI_Comm mpi_comm, char *codetype, double *h_densities, double *h_potentials, double *h_potentials_ext);
#endif
