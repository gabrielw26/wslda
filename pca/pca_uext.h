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

// -------------------------------------------------------------------------------------------
// ---------------------------VORTEX RECONNECTION STUDIES-------------------------------------
// -------------------------------------------------------------------------------------------
__device__ double borderfun(double v0, double alpha, double width, double xmax, double x) {
    if(fabs(x) <= xmax - (width + 2.)) return 0.;
    else if(fabs(x) >= xmax - 2.) return v0;
    else {
        if(x >= 0) return ( v0 * ( .5 + .5 * tanh(alpha * tan(M_PI * (x - xmax + width + 2) / width - M_PI / 2.)) ) );
        else return ( v0 * ( .5 + .5 * tanh(alpha * tan(M_PI * (-x - xmax + width + 2) / width - M_PI / 2.)) ) );
    }
}

__device__ double u_ext_VR_study(int ix, int iy, int iz, int it, int spin)
{
    double _ix = ((double)(ix) - 1.0*(NX/2))*DX;
    double _iy = ((double)(iy) - 1.0*(NY/2))*DY;
    double _iz = ((double)(iz) - 1.0*(NZ/2))*DZ;
    
    double borderpot = 0.;
    
    borderpot += borderfun(2., 2., 7., LX / 2., _ix);
    borderpot += borderfun(2., 2., 7., LY / 2., _iy);
    borderpot += borderfun(2., 2., 7., LZ / 2., _iz);
    return borderpot;

}
void process_params_u_ext_VR_study(double *params, double kF)
{
    // no process
}

// -------------------------------------------------------------------------------------------
// ----------------------------- HIGGS MODE STUDIES ------------------------------------------
// -------------------------------------------------------------------------------------------
/**
 * Smooth at the boundries HO potential
 * 0 - omega_x^2 / 2, in input file value at the boundry in units of eF_a
 * 1 - omega_y^2 / 2, in input file value at the boundry in units of eF_a
 * 2 - Umaxx - calculated automatically
 * 3 - Umaxy - calculated automatically
 * 4 - swich time, in number of iterations - IGNORE
 * ------------- stirer ---------------
 * 5 - speed of dragging, in input file in units of vF
 * 6 - distane from corner to start turning on potential
 * 7 - A - aplitude of stirer, in input file in units of eF
 * 8 - start time, in input in units of eF
 * 9 - switch time, in input in units of eF
 * 10 - none
 * 11 - gaussian width, in input file sigma
 * */
#define xs (3.0)
#define ys (3.0)
#define zs (3.0)
#define HO(x, omega2p2) (omega2p2*x*x)
__device__ double u_ext_smooth_HO(int ix, int iy, int iz, int it, int spin)
{
    double _x = (double)(ix)*DX - 1.0*(LX/2);
    double _y = (double)(iy)*DY - 1.0*(LY/2);
    double _z = (double)(iz)*DZ - 1.0*(LZ/2);
    
    // ---------------------------- HO ----------------------------
    double s;
    
    double hox;
    if(_x>=(double)(-LX/2)+xs && _x<=(double)(LX/2)-xs) 
    {
        hox = HO(_x, dc_params[0]);
    }
    else if(_x<(double)(-LX/2)+xs) 
    {   
        s = switch_function(_x+(double)(LX/2), xs, 1.0);
        hox = s*HO(_x, dc_params[0]) + (1.0-s)*dc_params[3];
    }
    else
    {
        s = switch_function(_x-((double)(LX/2)-xs), xs, 1.0);
        hox=(1.0-s)*HO(_x, dc_params[0]) + s*dc_params[3];
    }
    
    double hoy;
    if(_y>=(double)(-LY/2)+ys && _y<=(double)(LY/2)-ys) 
    {
        hoy = HO(_y, dc_params[1]);
    }
    else if(_y<(double)(-LY/2)+ys) 
    {   
        s = switch_function(_y+(double)(LY/2), ys, 1.0);
        hoy = s*HO(_y, dc_params[1]) + (1.0-s)*dc_params[4];
    }
    else
    {
        s = switch_function(_y-((double)(LY/2)-ys), ys, 1.0);
        hoy=(1.0-s)*HO(_y, dc_params[1]) + s*dc_params[4];
    }
    
    double hoz;
    if(_z>=(double)(-LZ/2)+zs && _z<=(double)(LZ/2)-zs) 
    {
        hoz = HO(_z, dc_params[2]);
    }
    else if(_z<(double)(-LZ/2)+zs) 
    {   
        s = switch_function(_z+(double)(LZ/2), zs, 1.0);
        hoz = s*HO(_z, dc_params[2]) + (1.0-s)*dc_params[5];
    }
    else
    {
        s = switch_function(_z-((double)(LZ/2)-zs), zs, 1.0);
        hoz=(1.0-s)*HO(_z, dc_params[2]) + s*dc_params[5];
    }
        
    return hox + hoy + hoz;

}
void process_params_u_ext_smooth_HO(double *params, double kF)
{
//     kF=  1.0;
    double eF = 0.5*kF*kF;
    
    // double Rx = 0.87*LX/2;
    double omega_x = 0.0212856534696; //  = sqrt(2.*eF / (Rx*Rx) );
    double omega_y = omega_x * (151./91);
    double omega_z = omega_x * (235./91);

    params[0] = omega_x * omega_x / 2.0;
    params[1] = omega_y * omega_y / 2.0;
    params[2] = omega_z * omega_z / 2.0;
    
//     Umax=HO(xs) + 1.0*(HO(0)-HO(xs))
    params[3] = HO((double)(-LX/2)+xs, params[0]) + 1.0*(HO((double)(-LX/2), params[0])-HO((double)(-LX/2)+xs, params[0]));
    params[4] = HO((double)(-LY/2)+ys, params[1]) + 1.0*(HO((double)(-LY/2), params[1])-HO((double)(-LY/2)+ys, params[1]));
    params[5] = HO((double)(-LZ/2)+zs, params[2]) + 1.0*(HO((double)(-LZ/2), params[2])-HO((double)(-LZ/2)+zs, params[2]));
    
    // internal parameters for coefficient GAMMA switch
    params[21]/=eF;
    params[22]/=eF;
    params[23]/=eF;
    params[24]/=eF;
    
    // noise timing
    params[10]*=eF;
    params[11]/=eF;
    params[12]/=eF;
}
#undef xs
#undef ys
#undef zs
#undef HO

/**
 * dc_params[10] - amplitude
 * dc_params[11] - t_kick
 * dc_params[12] - sigma
 * dc_params[13] - lambda_x
 * dc_params[14] - lambda_y
 * dc_params[15] - lambda_z
 * */
__device__ double u_ext_higgs_noise(int ix, int iy, int iz, int it, int spin)
{
    double _x = (double)(ix)*DX - 1.0*(LX/2);
    double _y = (double)(iy)*DY - 1.0*(LY/2);
    double _z = (double)(iz)*DZ - 1.0*(LZ/2);
    double _t = dc_t0 + dc_dt*it;
    
    double A = dc_params[10]*exp(-0.5*pow((_t-dc_params[11])/dc_params[12], 2)); 
    double u= A * sin(2.*M_PI*_x/dc_params[13])*sin(2.*M_PI*_y/dc_params[14])*sin(2.*M_PI*_z/dc_params[15]);
    
    return u;
}

// *******************************************************************************************
// **************************** FUNCTIONS CALLED BY THE CODE *********************************
// *******************************************************************************************
/** 
 * THIS FUNCTION IS CALLED BY KERLNELS FROM 'pca_kernels.cu'
 * */
__device__  double u_ext(int ix, int iy, int iz, int it, int spin)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0; // no external potential
#endif

//     return 0.0; // Higgs mode in uniform system
    return u_ext_higgs_noise(ix, iy, iz, it, spin);
//     return u_ext_smooth_HO(ix, iy, iz, it, spin);
}

/** 
 * THIS FUNCTION IS CALLED FROM 'pca.c'
 * AFTER LOADING params ARRAY FROM INPUT FILES.
 * AFTER PROCESSING THE params ARE LOADED TO CONST MEMORY ON GPU.
 * THE PARAMS ARE VISIBLE IN dc_params
 * @param params array of size MAX_USER_PARAMS with parameters
 * @param kF typical Fermi momentum of the problem
 * */
extern "C" void process_params(double *params, double kF)
{
#ifndef UNIFORM_TEST_MODE
    process_params_u_ext_smooth_HO(params, kF);
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
