#ifndef __GPE_UTILS__
#define __GPE_UTILS__

#include <complex.h>
#define Complex double complex

void set_gpu_device(int device);
void alloc_host_memory(uint nxyz, Complex **psi);
void free_host_memory(void *psi);
void set_initial_wave_function(uint nxyz, Complex *psi);
void read_initial_wave_function(uint nxyz, Complex *psi);
void write_to_binary_file(uint nxyz, Complex *psi);
void write_to_txt_file(uint nx, uint ny, uint nz, Complex *psi);
void print_header();
void print_intial_results(double time, double npart, double etot, double ekin, double eint, double eext);
void print_results(double time, double npart, double etot, double ekin, double eint, double eext, double diff, double rt);

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