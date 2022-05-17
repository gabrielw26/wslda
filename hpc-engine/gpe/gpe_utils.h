#ifndef __GPE_UTILS__
#define __GPE_UTILS__

#include <cufft.h>
typedef cufftDoubleComplex Complex;

#include "wslda_errors.h"
#define host_exec( cmd )                                                         \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "CPU ERROR: cannot execute: %s\n" , #cmd) ;\
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        report_error(ierr, stderr);                                             \
        return( EXIT_FAILURE ) ;                                                \
    } }

void set_initial_wave_function(uint nxyz, Complex *psi);
void read_initial_wave_function(uint nxyz, Complex *psi);
void write_to_binary_file(uint nxyz, Complex *psi);
void write_to_txt_file(uint nx, uint ny, uint nz, Complex *psi);
void print_header();
void print_intial_results(double time, double npart, double etot, double ekin, double eint, double eext);
void print_results(double time, double npart, double etot, double ekin, double eint, double eext, double diff, double rt);
void read_of_input_parameters(char* execcmd, int argc , char ** argv);
int parse_command_line_and_get_idx_of_input_file(int argc , char ** argv);
void read_input_file(int idx, char ** argv);
void save_extradata_to_file_with_outprefix(size_t size, void *extra_data, char* outprefix);

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