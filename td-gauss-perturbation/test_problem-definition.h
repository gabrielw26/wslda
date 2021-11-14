/**
 * W-SLDA Toolkit
 * Engine version: 2021.11.12
 * */

/**
 * Switch function - performs switch in time interval [0-T]
 * */
__device__ __host__ inline double switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}

__device__ __host__ inline double smooth_step(double t, double step_start, double step_stop, double T, double alpha)
{
    if(t<=step_start || t>=step_stop) return 0.0; 
    if(t>=step_start+T && t<=step_stop-T) return 1.0;
    if(t>step_start && t<step_start+T) return switch_function(t-step_start, T, alpha);
    else return 1.0-switch_function(t-step_stop+T, T, alpha);
}

/** 
 * EXTERNAL POTENTIAL V_ext
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2),
 *           NOTE: in case of 1d code iy=0
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 *           NOTE: in case of 1d and 2d codes iz=0
 * @param it iteration number
 * @param spin spin indicator, value from set {SPINA,SPINB}
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the external potential V_spin(x,y,z)
 * */
__device__ double v_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data)
{
    double x = DX*(ix-NX/2);
    double y = DY*(iy-NY/2);     // for 1d code iy will be always 0
    double z = DZ*(iz-NZ/2);     // for 1d and 2d codes iz will be always 0
    double t = dc_t0 + dc_dt*it; // time

    // ADD HERE FORMULA FOR V_ext(r)
    double _w = params[10]; // width, in lattice unit
    double _a = params[11]; // amplitude, in eF unit
    double s = smooth_step(t, params[12], params[13], params[14], 1.0);
    
    double V_ext;
    V_ext = -1.0*_a*s*exp(-0.5*pow(x/_w,2)); // attractive
    
    if(spin==SPINB && params[15]>0.5) // change sign
        V_ext*=-1.0;
        
    return V_ext; 
}
 
/** 
 * EXTERNAL PAIRING POTENTIAL Delta_ext
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2)
 *           NOTE: in case of 1d code iy=0
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 *           NOTE: in case of 1d and 2d codes iz=0
 * @param it iteration number
 * @param delta - value of delta computed self-consistently for given iteration it. 
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of external pairing potential Delta_{ext}(x,y,z)
 * Complex type is equivalent to thrust::complex<double>
 * for more info see: https://thrust.github.io/doc/structthrust_1_1complex.html
 * */
__device__ Complex delta_ext(int ix, int iy, int iz, int it, Complex delta, double *params, size_t extra_data_size, void *extra_data)
{
//     double x = DX*(ix-NX/2);
//     double y = DY*(iy-NY/2);     // for 1d code iy will be always 0
//     double z = DZ*(iz-NZ/2);     // for 1d and 2d codes iz will be always 0
//     double t = dc_t0 + dc_dt*it; // time

    // ADD HERE FORMULA FOR Delta_ext(r)
    Complex D_ext = Complex(0.0,0.0);

    return D_ext; 
}

/** 
 * EXTERNAL VELOCITY FIELD vec[v]_ext
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2)
 *           NOTE: in case of 1d code iy=0
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 *           NOTE: in case of 1d and 2d codes iz=0
 * @param it iteration number
 * @param spin spin indicator, value from set {SPINA,SPINB}
 * @param coordinate - Cartesian coordinate of the external velocity vector that should be computed, value from set {XAXIS, YAXIS, ZAXIS}
 *                     NOTE: for 1d code only XAXIS is requested, for 2d code XAXIS and YAXIS are requested.
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the external velocity vector v_ext(x,y,z)
 * */
__device__ double velocity_ext(int ix, int iy, int iz, int it, int spin, int coordinate, double *params, size_t extra_data_size, void *extra_data)
{
//     double x = DX*(ix-NX/2);
//     double y = DY*(iy-NY/2);     // for 1d code iy will be always 0
//     double z = DZ*(iz-NZ/2);     // for 1d and 2d codes iz will be always 0
//     double t = dc_t0 + dc_dt*it; // time

    // ADD HERE FORMULAS FOR vec{v}_ext=(vx, vy, vz)
    double v_ext;
    if(coordinate==XAXIS) v_ext=0.0;
    if(coordinate==YAXIS) v_ext=0.0;
    if(coordinate==ZAXIS) v_ext=0.0;

    return v_ext; 
}

/** 
 * THIS FUNCTION IS CALLED AT THE BEGINNING OF SIMULATION.
 * After loading params array from input file, the parameters are processed by this routine.
 * @param params array of size MAX_USER_PARAMS with parameters from input file. 
 * @param kF typical Fermi momentum scale of the problem. 
 *           kF=referencekF if the referencekF tag is indicated in the input file, 
 *           otherwise to kF value is assigned according formula kF=(3*pi^2*n)^{1/3}, where n corresponds to density in the box center.
 *           Note that it is passed via a pointer, to access/modify it use kF[0] or (*kF).
 * @param mu array with chemical potentials: mu[SPINA] and mu[SPINB]. 
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
extern "C" void process_params(double *params, double *kF, double *mu, size_t extra_data_size, void *extra_data)
{
    // PROCESS INPUT FILE PARAMETERS 
    double eF = 0.5*kF[0]*kF[0];
    params[10] /= kF[0]; // sigma
    params[11] *= eF; // amplitude
    params[12]/=eF; 
    params[13]/=eF;
    params[14]/=eF;

}

/**
 * THIS FUNCTION IS CALLED AFTER EACH COMPUTATION OF POTENTIALS
 * IT IS CALLED ONLY IF ENABLE_MODIFY_POTENTIALS IS DEFINED (in predefines.h)
 * Before each diagonalization process, user can modify arbitrarily potentials
 * @param it iteration number
 * @param h_densities structure with densities, see (wiki) documentation for list of fields
 *                    NOTE: densities structure is processed by modify_densities(...) function before call ot this function.
 * @param h_potentials struture with potentials, see (wiki) documentation for list of fields
 * Global variabls deliver by tdwslda_functionals_framework_enable.h are:
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
#include "tdwslda_functionals_framework_enable.h" // DO NOT REMOVE!
__global__ void modify_potentials(int it, wslda_density h_densities, wslda_potential h_potentials)
{    
    size_t ixyz= threadIdx.x + blockIdx.x * blockDim.x; // compute for this point
    int ix, iy, iz, i;
    
    if(ixyz<BLOCKLENGTH)
    {
        // decode_ixyz2ixiyiz(ixyz,ix,iy,iz, i);
        // // Now ix, iy, iz keeps lattice coordinate. 
        
        // h_potentials.V_a[ixyz] stores value of spin-up particles mean-field potential for coordinate (x,y,z)
        // and similarly for other potentials
        // ... below you can modify them your wish ...
        // 
    }
}
#include "tdwslda_functionals_framework_disable.h" // DO NOT REMOVE!

/**
 * This function provides size of extra_data array, in bytes.
 * The extra_data of specified size will be allocated by the main process.
 * This function is thread-safe.
 * @param params with input file parameters. 
 *              NOTE: the array contains bare input file values, not processed by process_params()!
 * @return size of the extra_data array that needs to be allocated, if 0 then extra_data will not be allocated.
 * */
extern "C" size_t get_extra_data_size(double *params)
{
    return 0;
}

/**
 * This function loads data into extra_data array.
 * This function is thread-safe.
 * @param size size of array computed using function get_extra_data_size()
 * @param extra_data pointer to array that should be filled with data
 * @param params with input file parameters. 
 *               NOTE: the array contains bare input file values, not processed by process_params()!
 * @return 0 if load is successful, otherwise return error code. If nonzero value is returned the main code will terminate.
 * */
extern "C" int load_extra_data(size_t size, void *extra_data, double *params)
{
    return 0;
}



