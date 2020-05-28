// Author: Gabriel Wlazlowski

// Code implements external potential - in similar way as in pca_uext

#ifndef __S3DPCA_UEXT__
#define __S3DPCA_UEXT__

/**
 * Switch function - performs switch in time interval [0-T]
 * */
inline double switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}

inline double smooth_step(double t, double step_start, double step_stop, double T, double alpha)
{
    if(t<=step_start || t>=step_stop) return 0.0; 
    if(t>=step_start+T && t<=step_stop-T) return 1.0;
    if(t>step_start && t<step_start+T) return switch_function(t-step_start, T, alpha);
    else return 1.0-switch_function(t-step_stop+T, T, alpha);
}

// =================================================================================
// ============================== EXTERNAL POTENTIAL ===============================
// =================================================================================

double u_ext_HO(int ix, int iy, int iz, int it, int spin)
{
    double _x = ((double)(ix) - 1.0*(NX/2))*DX;
    double _y = ((double)(iy) - 1.0*(NY/2))*DY;
    double _z = ((double)(iz) - 1.0*(NZ/2))*DZ;
    
    double swtch=1.0;
//     if(1.0*it<dc_params[2]) swtch=switch_function(1.0*it, dc_params[2], 1.0);
    
    return swtch*(dc_params[0]*_x*_x + dc_params[1]*_y*_y + dc_params[2]*_z*_z);
}
void process_params_u_ext_HO(double *params, double kF)
{
    double eF = 0.5*kF*kF;
    params[0] = 4.0*params[0]*eF/(DX*DX*NX*NX);
    params[1] = 4.0*params[1]*eF/(DY*DY*NY*NY);
    params[2] = 4.0*params[2]*eF/(DZ*DZ*NZ*NZ);
}
void imprint_vortex_along_z(int it, double *h_densities, double *h_potentials, double *extra_data)
{
    // pontentials - decode
    double complex *delta = (double complex *)(h_potentials +  2*NXYZ);
    
    int ix, iy, iz, ixyz=-1;
    
    double _x, _y, _r;
    double arg, abs_delta;
    
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++) // iterate over all lattice points
    {
        ixyz++;
        
        double _x = ((double)(ix) - 1.0*(NX/2))*DX;
        double _y = ((double)(iy) - 1.0*(NY/2))*DY;
        _r = sqrt(_x*_x + _y*_y); // distance from the center
    
        if(_r<1.0*(LX/2.)) 
        {
            abs_delta = cabs(delta[ixyz]);
            arg = atan2(_y,_x);             
            delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);
        }            
    }
}

// -------------------------------------------------------------------------------------------
// ---------------------------VORTEX RECONNECTION STUDIES-------------------------------------
// -------------------------------------------------------------------------------------------
double borderfun(double v0, double alpha, double width, double xmax, double x) {
    if(fabs(x) <= xmax - (width + 2.)) return 0.;
    else if(fabs(x) >= xmax - 2.) return v0;
    else {
        if(x >= 0) return ( v0 * ( .5 + .5 * tanh(alpha * tan(M_PI * (x - xmax + width + 2) / width - M_PI / 2.)) ) );
        else return ( v0 * ( .5 + .5 * tanh(alpha * tan(M_PI * (-x - xmax + width + 2) / width - M_PI / 2.)) ) );
    }
}

double u_ext_VR_study(int ix, int iy, int iz, int it, int spin)
{
    double _ix = ((double)(ix) - 1.0*(NX/2))*DX;
    double _iy = ((double)(iy) - 1.0*(NY/2))*DY;
    double _iz = ((double)(iz) - 1.0*(NZ/2))*DZ;
    
    double borderpot = 0.;
    
    borderpot += borderfun(2., 2., 7., LX / 2., _ix);
    borderpot += borderfun(2., 2., 7., LY / 2., _iy);
    borderpot += borderfun(2., 2., 7., LZ / 2., _iz);
    
    if(dc_params[12]<0.1)  return borderpot; // dc_params[12] codes if activate this mode

    // addition for studies of vortex polarized-unpolarized reconnection
    // add part chemical potantial dependent
    double mu_avg = 0.5*(dc_params[10]+dc_params[11]); // average chemical potential
    double mu_add;
    if(spin==SPINA)
      mu_add=dc_params[10]-mu_avg;
    else
      mu_add=dc_params[11]-mu_avg;
    
    #define SWTCH 5.0
    if(_ix<-SWTCH) // no modify at left
      return borderpot;
    else if(_ix<SWTCH) // switch region
      return borderpot + mu_add*switch_function(_ix+SWTCH, 2.*SWTCH, 1.0);
    else
      return borderpot + mu_add;
    

}
void process_params_u_ext_VR_study(double *params, double kF)
{
    // no process
}

void modify_potentials_u_ext_VR_study(int it, double *h_densities, double *h_potentials, double *extra_data)
{
    // here extra_data contains phase of the order parameter to be imprinted
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NXYZ);
    double *V_b = (double *)(h_potentials +  1*NXYZ);
    double complex *delta = (double complex *)(h_potentials +  2*NXYZ);
    
    int ixyz;
    double abs_delta, arg;
    for(ixyz=0; ixyz<NXYZ; ixyz++)
    {
        abs_delta = cabs(delta[ixyz]);
        arg = extra_data[ixyz]; // take phase from external source
        
        // phase imprint
        delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg); 
    }
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
double u_ext_smooth_HO(int ix, int iy, int iz, int it, int spin)
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
//     double kF=  1.0;
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
    
}
#undef xs
#undef ys
#undef zs
#undef HO

void modify_potentials_phase_random(int it, double *h_densities, double *h_potentials)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=0;
    double _x, _y, _r;
    double arg, abs_delta;
    srand(123);
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        abs_delta = cabs(delta[ixyz]);
        arg = (double)rand() / (double)RAND_MAX;
        arg = (2.0*arg-1.0)*M_PI;

        // phase imprint
        delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);         
        
        ixyz++;
    }
}


// *******************************************************************************************
// **************************** FUNCTIONS CALLED BY THE CODE *********************************
// *******************************************************************************************

/** 
 * THIS FUNCTION IS CALLED BY MAIN ENGINE.
 * After loading params array from input files, the parameters are processed by this routine.
 * The routine is executed at beginning of each iteration.
 * @param params array of size MAX_USER_PARAMS with parameters from input file. 
 * @param kF typical Fermi momentum scale of the problem. 
 *           kF=referencekF if the referencekF tag is indicated in the input file, 
 *           otherwise to kF value (3*pi^2*n)^{1/3} is assigned, where n corresponds to density in the box center
 * @param mu array with chemical potentials: mu[SPINA] and mu[SPINB]. 
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void process_params(double *params, double kF, double *mu, size_t extra_data_size, void *extra_data)
{
#ifndef UNIFORM_TEST_MODE
    // no processing
    
    // test - HO
//     process_params_u_ext_HO(params, kF);
    
    // studies of VR
//     process_params_u_ext_VR_study(params, kF);
    
    // Higgs mode
    return process_params_u_ext_smooth_HO(params, kF);
#endif
}

/**
 * THIS FUNCTION IS CALLED AFTER EACH EXTRATION OF densities AND potentials
 * Use this function is you want to impose additional constraints on densities and/or potentials.
 * @param it iteration number
 * @param h_densities array with densities, see (wiki) documentation for decoding prescription
 * @param h_potentials array with potentials, see (wiki) documentation for decoding prescription
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void modify_potentials(int it, double *h_densities, double *h_potentials, double *params, size_t extra_data_size, void *extra_data)
{
#ifdef UNIFORM_TEST_MODE
    if(dc_params[31]>0.5 && it<=1) modify_potentials_phase_random(it, h_densities, h_potentials); 
#endif 
    
#ifndef UNIFORM_TEST_MODE
    // no modify
    
    // VR studies
//     modify_potentials_u_ext_VR_study(it, h_densities, h_potentials, extra_data);
    
    // imprint vortex line 
//     imprint_vortex_along_z(it, h_densities, h_potentials, extra_data);
#endif
}

/** 
 * THIS FUNCTION IS CALLED BY MAIN ENGINE.
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2)
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 * @param it iteration number
 * @param spin spin indicator, value from set {SPINA,SPINB}
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the external potential V_spin(x,y,z)
 * */
double u_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0; // no external potential
#endif
    
//     return u_ext_HO(ix, iy, iz, it, spin);
    
    // studies of VR
//     return u_ext_VR_study(ix, iy, iz, it, spin);
    
    // Higgs mode
    return u_ext_smooth_HO(ix, iy, iz, it, spin);
}

/** 
 * THIS FUNCTION IS CALLED BY MAIN ENGINE.
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2)
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 * @param it iteration number
 * @param delta - value of delta computed self-consitently for given iteration it. 
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of external pairing potential Delta_{ext}(x,y,z)
 * */
double complex delta_ext(int ix, int iy, int iz, int it, double complex delta, double *params, size_t extra_data_size, void *extra_data)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0 + I*0.0; // no external potential
#endif

    // ADD HERE YOUR CODE
    
    return 0.0 + I*0.0;
}

/** 
 * THIS FUNCTION IS CALLED BY MAIN ENGINE.
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2)
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 * @param it iteration number
 * @param spin spin indicator, value from set {SPINA,SPINB}
 * @param coordinate - Cartesian coordinate of the external velocity vector that should be computed, value from set {XAXIS, YAXIS, ZAXIS}
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the external velocity vector v_ext(x,y,z)
 * */
double vector_vext(int ix, int iy, int iz, int it, int spin, int coordinate, double *params, size_t extra_data_size, void *extra_data)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0; // no external velocity field
#endif 

    // ADD HERE YOUR CODE
    
    return 0.0;
}

/**
 * This function provides size of extra_data array, in bytes.
 * The extra_data of specified size will be allocated by the main process.
 * This function is thread-safe.
 * @param array with input file parameters. NOTE: the array contains bare input file values, not processed by process_params()!
 * @return size of the extra_data array than needs to be allocated, if 0 then extra_data then extra_data will not be allocated.
 * */
size_t get_extra_data_size(double *params)
{
    return 0;
}

/**
 * This function loads data into extra_data array.
 * This function is thread-safe.
 * @param size size of array compute using function get_extra_data_size()
 * @param extra_data pointer to array that should be filled with data
 * @param array with input file parameters. NOTE: the array contains bare input file values, not processed by process_params()!
 * @return 0 if load is successful, otherwise return error code. If nonzero value is returned the main code will terminate.
 * */
int load_extra_data(size_t size, void *extra_data, double *params)
{
    return 0;
}
#endif

