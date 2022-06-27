#ifndef __GET_EXT_POTENTIALS__
#define __GET_EXT_POTENTIALS__

// TO MAKE CONFORM WITH WSLDA TOOLKIT
#ifndef TDWSLDA
#define TDWSLDA
int get_v_ext(int datadim, int spin, int it, double *data); // implemented on cuda
int get_delta_ext(int datadim, int it, void *deltain, void *data); // implemented on cuda
int get_velocity_ext(int datadim, int spin, int it, double *data); // implemented on cuda
#endif

#ifdef TDWSLDA
void * ptr_d_delta; //pointer to delta on device
void set_ptr_d_delta(void * ptr) {ptr_d_delta=ptr;}
#endif

#endif