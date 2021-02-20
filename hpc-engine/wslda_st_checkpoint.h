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

int wslda_st_read_checkpoint();


int wslda_st_write_checkpoint(int codedim, int it, 
                              int nconsts, double *consts, 
                              int npot, double *h_potentials, 
                              int ndens, double *h_densities, 
                              int nenergy, double *energy,
                              int nbroy, double **dens_in, double **dens_out
                              );
#endif
