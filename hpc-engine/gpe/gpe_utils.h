#ifndef __GPE_UTILS__
#define __GPE_UTILS__

int gpe_compute(int type);

// TO MAKE COMPATIBILE WITH WSLDA TOOLKIT
#ifndef TDWSLDA
#define TDWSLDA
int get_v_ext(int datadim, int spin, int it, double *data); // implemented on cuda
void * ptr_d_delta; //pointer to delta on device 
int get_delta_ext(int datadim, int it, void *deltain, void *data); // implemented on cuda
void set_ptr_d_delta(void * ptr) {ptr_d_delta=ptr;}
int get_velocity_ext(int datadim, int spin, int it, double *data); // implemented on cuda
#endif

#endif