
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
double v_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data)
{
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
 * */
double complex delta_ext(int ix, int iy, int iz, int it, double complex delta, double *params, size_t extra_data_size, void *extra_data)
{
    // ADD HERE FORMULA FOR Delta_ext(r)
    double complex D_ext = 0.0 + I*0.0;

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
double velocity_ext(int ix, int iy, int iz, int it, int spin, int coordinate, double *params, size_t extra_data_size, void *extra_data)
{
    // ADD HERE FORMULAS FOR vec{v}_ext=(vx, vy, vz)
    double v_ext;
    if(coordinate==XAXIS) v_ext=0.0;
    if(coordinate==YAXIS) v_ext=0.0;
    if(coordinate==ZAXIS) v_ext=0.0;

    return v_ext; 
}

/** 
 * THIS FUNCTION IS CALLED DURING THE SELF-CONSISTENT PROCESS.
 * After loading params array from input file, the parameters are processed by this routine.
 * The routine is executed at beginning of each iteration.
 * @param params array of size MAX_USER_PARAMS with parameters from input file. 
 * @param kF typical Fermi momentum scale of the problem. 
 *           kF=referencekF if the referencekF tag is indicated in the input file, 
 *           otherwise to kF value is assigned according formula kF=(3*pi^2*n)^{1/3}, where n corresponds to density in the box center
 * @param mu array with chemical potentials: mu[SPINA] and mu[SPINB]. 
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void process_params(double *params, double kF, double *mu, size_t extra_data_size, void *extra_data)
{
    // PROCESS INPUT FILE PARAMETERS 

}

/**
 * THIS FUNCTION IS CALLED DURING THE SELF-CONSISTENT PROCESS.
 * Before each diagonalization process, user can modify arbitrarily densities
 * @param it iteration number
 * @param h_densities array with densities, see (wiki) documentation for decoding prescription
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void modify_densities(int it, double *h_densities, double *params, size_t extra_data_size, void *extra_data)
{
    // DENSITIES DECODING
    double *rho_a = (double *)(h_densities +  0*BLOCKSIZE);
    double *rho_b = (double *)(h_densities +  1*BLOCKSIZE);
    double *tau_a = (double *)(h_densities +  2*BLOCKSIZE);
    double *tau_b = (double *)(h_densities +  3*BLOCKSIZE);
    double complex *nu = (double complex *)(h_densities +  4*BLOCKSIZE);
    double *j_a_x = (double *)(h_densities +  6*BLOCKSIZE);
    double *j_a_y = (double *)(h_densities +  7*BLOCKSIZE);
    double *j_a_z = (double *)(h_densities +  8*BLOCKSIZE);
    double *j_b_x = (double *)(h_densities +  9*BLOCKSIZE);
    double *j_b_y = (double *)(h_densities + 10*BLOCKSIZE);
    double *j_b_z = (double *)(h_densities + 11*BLOCKSIZE);    
    
    // DETERMINE VARIANT OF THE CODE
    int lNX, lNY, lNZ; // local sizes
    int ix, iy, iz, ixyz;
    if(BLOCKSIZE==NX      ) {lNX=NX; lNY=1 ; lNZ=1 ;} // 1D code
    if(BLOCKSIZE==NX*NY   ) {lNX=NX; lNY=NY; lNZ=1 ;} // 2D code
    if(BLOCKSIZE==NX*NY*NZ) {lNX=NX; lNY=NY; lNZ=NZ;} // 3D code
    
    if(params[31]>0.5 && it<=1) // add noise
    {
        srand(123);
        double arg, abs_nu;
        
        // ITERATE OVER ALL POINTS
        ixyz=0;
        for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
        {
            double x = DX*(ix-lNX/2);
            double y = DY*(iy-lNY/2); // for 1d code y will be 0
            double z = DZ*(iz-lNZ/2); // for 1d and 2d codes z will be 0
            
            // rho_a[ixyz] stores value of spin-up particles densities for coordinate (x,y,z)
            // and similarly for other densities
            // ... below you can modify at your wish ...
            abs_nu = cabs(nu[ixyz]);
            arg = (double)rand() / (double)RAND_MAX;
            arg = (2.0*arg-1.0)*M_PI;

            // phase imprint
            nu[ixyz]=abs_nu*cos(arg) + I*abs_nu*sin(arg);
          
            ixyz++; // go to next point,  it should be last line of the triple loop
        }
    }
}

/**
 * THIS FUNCTION IS CALLED DURING THE SELF-CONSISTENT PROCESS.
 * Before each diagonalization process, user can modify arbitrarily potentials
 * @param it iteration number
 * @param h_densities array with densities, see (wiki) documentation for decoding prescription
 *                    NOTE: densities array is processed by modify_densities(...) function.
 * @param h_potentials array with potentials, see (wiki) documentation for decoding prescription
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void modify_potentials(int it, double *h_densities, double *h_potentials, double *params, size_t extra_data_size, void *extra_data)
{
    // DENSITIES DECODING
    double *rho_a = (double *)(h_densities +  0*BLOCKSIZE);
    double *rho_b = (double *)(h_densities +  1*BLOCKSIZE);
    double *tau_a = (double *)(h_densities +  2*BLOCKSIZE);
    double *tau_b = (double *)(h_densities +  3*BLOCKSIZE);
    double complex *nu = (double complex *)(h_densities +  4*BLOCKSIZE);
    double *j_a_x = (double *)(h_densities +  6*BLOCKSIZE);
    double *j_a_y = (double *)(h_densities +  7*BLOCKSIZE);
    double *j_a_z = (double *)(h_densities +  8*BLOCKSIZE);
    double *j_b_x = (double *)(h_densities +  9*BLOCKSIZE);
    double *j_b_y = (double *)(h_densities + 10*BLOCKSIZE);
    double *j_b_z = (double *)(h_densities + 11*BLOCKSIZE); 
    
    // POTENTIALS DECODING
    double *V_a = (double *)(h_potentials +  0*BLOCKSIZE);
    double *V_b = (double *)(h_potentials +  1*BLOCKSIZE);
    double complex *delta = (double complex *)(h_potentials +  2*BLOCKSIZE); 
    
    // DETERMINE VARIANT OF THE CODE
    int lNX, lNY, lNZ; // local sizes
    int ix, iy, iz, ixyz;
    if(BLOCKSIZE==NX      ) {lNX=NX; lNY=1 ; lNZ=1 ;} // 1D code
    if(BLOCKSIZE==NX*NY   ) {lNX=NX; lNY=NY; lNZ=1 ;} // 2D code
    if(BLOCKSIZE==NX*NY*NZ) {lNX=NX; lNY=NY; lNZ=NZ;} // 3D code
    
    // ITERATE OVER ALL POINTS
    ixyz=0;
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {
        double x = DX*(ix-lNX/2);
        double y = DY*(iy-lNY/2); // for 1d code y will be 0
        double z = DZ*(iz-lNZ/2); // for 1d and 2d codes z will be 0
        
        // V_a[ixyz] stores value of spin-up particles mean-field potential for coordinate (x,y,z)
        // and similarly for other potentials
        // ... below you can modify at your wish ...
        
        
        ixyz++; // go to next point,  it should be last line of the triple loop
    }
}


/**
 * This function provides size of extra_data array, in bytes.
 * The extra_data of specified size will be allocated by the main process.
 * This function is thread-safe.
 * @param params with input file parameters. 
 *              NOTE: the array contains bare input file values, not processed by process_params()!
 * @return size of the extra_data array that needs to be allocated, if 0 then extra_data will not be allocated.
 * */
size_t get_extra_data_size(double *params)
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
int load_extra_data(size_t size, void *extra_data, double *params)
{
    return 0;
}



