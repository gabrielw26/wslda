//
// This file provides function for computation of external potential
//

#ifndef __PCA_UEXT__
#define __PCA_UEXT__

__constant__ double dc_params[MAX_USER_PARAMS]; // array with params from input file

// =================================================================================
// ============================== EXTERNAL POTENTIAL ===============================
// =================================================================================

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
 * For testing
 * */
__device__  double u_ext_test(int ix, int iy, int iz, int it, int spin)
{

    double _x = (double)(ix) - 1.0*(NX/2);
    double _y = (double)(iy) - 1.0*(NY/2);
    double _z = (double)(iz) - 1.0*(NZ/2);
    double time=dc_t0 + dc_dt*it;
    
    #define eF_a 0.5
    #define T_START 10.0
    #define T_STOP  100.0
    #define T_SWITCH 25.0
    #define AMPLITUDE 2.0
    #define SIGX 2.7
    #define SIGY 2.8
    #define SIGZ 2.9

    double gauss = AMPLITUDE * eF_a *
                   smooth_step(time, T_START/eF_a, T_STOP/eF_a, T_SWITCH/eF_a, 1.0) *
                   exp(-1.*_x*_x/(2.*SIGX*SIGX) -1.*_y*_y/(2.*SIGY*SIGY) -1.*_z*_z/(2.*SIGZ*SIGZ) );
    
//     return gauss;

    // you can also play with other options
    if(spin==SPINA) 
        return gauss;
    else
        return -1.0*gauss;
    
    #undef eF_a
    #undef T_START
    #undef T_STOP
    #undef T_SWITCH
    #undef AMPLITUDE
    #undef SIGX
    #undef SIGY
    #undef SIGZ
}

/**                                                                                                                                                                                                               
 * 0 - r1                                                                                                                                                                                                         
 * 1 - r2                                                                                                                                                                                                         
 * 2 - V0                                                                                                                                                                                                         
 * */
__device__ double u_ext_tube(int ix, int iy, int iz, int it, int spin)
{
    double _x = ((double)(ix) - 1.0*(NX/2));
    double _y = ((double)(iy) - 1.0*(NY/2));
    double _r = hypot(_x,_y);

    if(_r<dc_params[0])
        return 0.0;
    else if(_r<dc_params[1])
        return dc_params[2]*switch_function(_r-dc_params[0], dc_params[1]-dc_params[0], 1.0);
    else
        return dc_params[2];

}
void process_u_ext_tube(double *params, double kF)
{
    
    double eF = kF*kF/2.0;
    params[2] = params[2] * eF;
}

/**                                                                                                                                                                                                               
 * 0 - amplitude, in eF on entry                                                                                                                                                                                                        
 * 1 - time when the amplitude is maximal, in eF^-1 on entry                                                                                                                                                                                                         
 * 2 - width of perturbation in time, in eF^-1 on entry    
 * 3 - wave-length                                                                                                                                                                                                     
 * */
__device__ double u_sin(int ix, int iy, int iz, int it, int spin)
{
    double _x = (double)(ix)*DX - 1.0*(LX/2);
    double _t = dc_t0 + dc_dt*it;
    
    double A = dc_params[0]*exp(-0.5*pow((_t-dc_params[1])/dc_params[2], 2)); 
    double u= A * sin(2.*M_PI*_x/dc_params[3]);
    return u;
}
void process_u_sin(double *params, double kF)
{
    
    double eF = kF*kF/2.0;
    params[0] *= eF;
    params[1] /= eF;
    params[2] /= eF;
}


// *******************************************************************************************
// **************************** FUNCTIONS CALLED BY THE CODE *********************************
// *******************************************************************************************
/** 
 * THIS FUNCTION IS CALLED BY KERLNELS FROM '(c)(c)pca_kernels.cu'
 * @param ix - coordinate x,  in range [0,NX)
 * @param iy - coordinate y,  in range [0,NY), in case of ccpca code it is always set as iy=0 
 * @param iz - coordinate z,  in range [0,NZ), in case of ccpca and cpca code it is always set as iz=0
 * @param it - index if time step, time is computed as time = dc_t0 + dc_dt*it
 * @param spin - spin selector, spin in {SPINA,SPINB}
 * @return value of external potential V_{spin}(x,y,z,t)
 * */
__device__  double u_ext(int ix, int iy, int iz, int it, int spin)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0; // no external potential
#endif

    return u_sin(ix, iy, iz, it, spin);
}

/** 
 * THIS FUNCTION IS CALLED BY KERLNELS FROM '(c)(c)pca_kernels.cu' ONLY IF ENABLE_DELTA_EXT IS ACTIVE!
 * @param ix - coordinate x,  in range [0,NX)
 * @param iy - coordinate y,  in range [0,NY), in case of ccpca code it is always set as iy=0 
 * @param iz - coordinate z,  in range [0,NZ), in case of ccpca and cpca code it is always set as iz=0
 * @param it - index if time step, time is computed as time = dc_t0 + dc_dt*it
 * @param delta - value of delta computed self-consitently for given point in time. Note that simulation will conserve particle number only if arg[delta] = arg[Delta_{ext}], otherwise the particle number conervation will be violated. 
 * @return value of exterrnal pairing potential Delta_{ext}(x,y,z,t)
 * */
__device__  Complex delta_ext(int ix, int iy, int iz, int it, Complex delta)
{
#ifdef UNIFORM_TEST_MODE
    return Complex(0.0,0.0); // no external potential
#endif

    return Complex(0.0,0.0);
}

/** 
 * THIS FUNCTION IS CALLED FROM 'pca.c'
 * AFTER LOADING params ARRAY FROM INPUT FILES.
 * AFTER PROCESSING THE params ARE LOADED TO CONST MEMORY ON GPU.
 * THE PARAMS ARE VISIBLE IN dc_params
 * @param params array of size MAX_USER_PARAMS with parameters
 * @param kF typical Fermi momentum of the problem
 * @param mu chemical potentials, mu[SPINA] and mu[SPINB]
 * */
extern "C" void process_params(double *params, double kF, double *mu)
{
#ifndef UNIFORM_TEST_MODE
    process_u_sin(params,kF);
#endif
}

#ifdef WORK_IN_ROTATING_FRAME
/** 
 * THIS FUNCTION IS CALLED FROM 'pca.c'
 * IT SETS CURRENT VALUES OF dc_Omega_a and dc_Omega_a;
 * @param params array of size MAX_USER_PARAMS with parameters
 * @param kF typical Fermi momentum of the problem
 * @param time value of time in code units
 * @param Omega_a returns to main code value of Omega_a (OUTPUT)
 * @param Omega_b returns to main code value of Omega_b (OUTPUT)
 * */
extern "C" int set_Omega(double *params, double kF, double time, double *Omega_a, double *Omega_b)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0; // no external potential
#endif

//     if(time<params[20])
//     {
//         *Omega_a = 0.0;
//         *Omega_b = 0.0;
//     }
//     else if(time<=params[21])
//     {
//         double s = switch_function(time-params[20], params[21]-params[20], 1.0);
//         *Omega_a = s*params[22];
//         *Omega_b = s*params[22];        
//     
//         // and copy to register if needed
//         if( cudaMemcpyToSymbol(dc_Omega_a, Omega_a, sizeof(double))!= cudaSuccess ) return 1;
//         if( cudaMemcpyToSymbol(dc_Omega_b, Omega_b, sizeof(double))!= cudaSuccess ) return 2;
//     }
//     else
//     {
//         *Omega_a = params[22];
//         *Omega_b = params[22];        
//     }
    
    return 0;
}
#endif

#endif
