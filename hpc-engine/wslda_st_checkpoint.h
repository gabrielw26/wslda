/**
 * W-SLDA Toolkit
 * 
 * Warsaw University of Technology, Faculty of Physics (2021)
 * 
 * @author Gabriel Wlazlowski 
 * */ 

#ifndef __WSLDA_ST_CHECKPOINT__
#define __WSLDA_ST_CHECKPOINT__
/**
 * This file defines functions for reading and writing checkpoint
 * */

#define WSLDA_ST_CHECKPOINT_DAT 1
#define WSLDA_ST_CHECKPOINT_OLD 2
#define WSLDA_ST_CHECKPOINT_NOFILE 2

#define ST_CHECKPOINT_2D_TO_3D 33
#define ST_CHECKPOINT_1D_TO_3D 34
#define ST_CHECKPOINT_1D_TO_2D 35
#define ST_CHECKPOINT_RESIZE 36


int wslda_stcheckpoint_format(int codedim);

int wslda_st_required_operations(int codedim, int *intepolation, int *resize);

int wslda_st_checkpoint_convert(int operation, int fileidx, int codedim, int *it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              );

int wslda_st_read_checkpoint(int fileidx, int codedim, int *it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              );

int wslda_st_write_checkpoint(int codedim, int it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              );
#endif
