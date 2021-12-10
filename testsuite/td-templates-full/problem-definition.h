
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
//     double x = DX*(ix-NX/2);
//     double y = DY*(iy-NY/2);     // for 1d code iy will be always 0
//     double z = DZ*(iz-NZ/2);     // for 1d and 2d codes iz will be always 0
//     double t = dc_t0 + dc_dt*it; // time

    // ADD HERE FORMULA FOR V_ext(r)
    double V_ext = 0.0;

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

}

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

/**
 * Scattering length, in code units.
 * This function is meaningful only in the case of BDG or SLDAE functionals.
 * For SLDA and ASLDA the scattering length is assumed to be infinite, and the function is ignored.
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2),
 *           NOTE: in case of 1d code iy=0
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 *           NOTE: in case of 1d and 2d codes iz=0
 * @param it time iteration, use dc_t0 + dc_dt*it to compute corresponding time 
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the scattering length a(x,y,z,t).
 * */
__device__ double scattering_length(int ix, int iy, int iz, int it, double *params, size_t extra_data_size, void *extra_data)
{
    return dc_sclgth; // by default return value from input file.
    
    // however, here you can define your own prescription
    // it can be time and position dependent
    // ...
}


