
/** 
 * EXTERNAL POTENTIAL V_ext
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2),
 *           NOTE: in case of 1d code iy=0
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 *           NOTE: in case of 1d and 2d codes iz=0
 * @param it iteration number
 * @param spin spin indicator, - IGNORE FOR GPE
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the external potential V_spin(x,y,z)
 * */
// analog of inline __device__  double gpe_external_potential(uint ix, uint iy, uint iz, uint it)
inline __device__ double v_ext(int ix, int iy, int iz, int it, int spin, double *params, size_t extra_data_size, void *extra_data)
{
    // harmonic trap:
    // V(x,y,z) = 0.5*(omega_x*x)^2 + 0.5*(omega_y*y)^2 + 0.5*(omega_z*z)^2
    
    // frequencies are passed through d_user_param array
    double omega_x = params[0];
    double omega_y = params[1];
    double omega_z = params[2];
    
    // coordinate with respect to center of the box
    double _ix = (double)(ix) - 1.0*(NX/2);
    double _iy = (double)(iy) - 1.0*(NY/2);
    double _iz = (double)(iz) - 1.0*(NZ/2);


    double trap =   0.5*_ix*_ix*omega_x*omega_x
                  + 0.5*_iy*_iy*omega_y*omega_y 
                  + 0.5*_iz*_iz*omega_z*omega_z;
    
    return trap;
}

/** 
 * THIS FUNCTION IS CALLED DURING THE SELF-CONSISTENT PROCESS.
 * After loading params array from input file, the parameters are processed by this routine.
 * The routine is executed at beginning of each iteration.
 * @param params array of size MAX_USER_PARAMS with parameters from input file. 
 * @param kF typical Fermi momentum scale of the problem. 
 *           kF=referencekF if the referencekF tag is indicated in the input file, 
 *           otherwise to kF value is assigned according formula kF=(3*pi^2*n)^{1/3}, where n corresponds to maximal density.
 *           You can also set kF at request in this function using (*kF)=myvalue;
 * @param mu array with chemical potentials: mu[SPINA] and mu[SPINB]. 
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * */
void process_params(double *params, double *kF, double *mu, size_t extra_data_size, void *extra_data)
{
    // PROCESS INPUT FILE PARAMETERS 
}

/**
 * This function computes Fermi momentum, which is used as the reference value. 
 * Other reference scales are set automatically to: eF=kF^2/2, Effg=(3/5)*N*eF (N-total number of particles)
 * For more details see: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Reference%20scales
 * NOTE units are: hbar=m=k_b=1
 * @param it iteration number
 * @param h_densities structure with densities, see (wiki) documentation for list of fields
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of Fermi momentum for your problem
 * */
double referencekF(int it, wslda_density h_densities, double *params, size_t extra_data_size, void *extra_data)
{
    return 0; // TODO
    // if(input->referencekF>0.0) return input->referencekF; // take it from input file
    
    // // define here your prescription for computing kF
    // // ...
    // // default: extract max density and use it for definition of kF
    // double max_dens=0.0, kF;
    // int ixyz;
    // for(ixyz=0; ixyz<h_densities.nx*h_densities.ny*h_densities.nz; ixyz++) 
    //     if(h_densities.rho_a[ixyz]+h_densities.rho_b[ixyz]>max_dens) max_dens=h_densities.rho_a[ixyz]+h_densities.rho_b[ixyz];
    
    // // depending on dimensionality of the problem
    // if(NY==1 && NZ==1) kF = 0.5*M_PI*max_dens;                // 1D
    // else if(NZ==1)     kF = pow(2.0*M_PI*max_dens,1./2.);     // 2D
    // else               kF = pow(3.*M_PI*M_PI*max_dens,1./3.); // 3D
    
    // return kF;
}

/**
 * Function changes wave function.
 * This function is called before each integration step.
 * NOTE: This function assumes that norm is not changed after modification. 
 * @param ix - x coordinate, ix=0,1,...,d_nx-1, where d_nx is global variable 
 * @param iy - y coordinate, iy=0,1,...,d_ny-1, where d_ny is global variable 
 * @param iz - z coordinate, iy=0,1,...,d_nz-1, where d_ny is global variable 
 * @param it - time value, ie. time = d_t0 + it*d_dt, d_t0 and d_dt are global variables
 * @param psi - psi(ix, iy, iz, it) - psi is normalized, i.e. int n(r) d^3r = npart, where n(r) computed according gpe_density(psi)
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if modification is successful, otherwise return error code. If nonzero value is returned the main code will terminate.
 * */
inline __device__  int gpe_modify_psi(int ix, int iy, int iz, int it, Complex *psi, double *params, size_t extra_data_size, void *extra_data)
{
    // no change
    return 0;
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

/**
 * Scattering length, in code units.
 * This function is meaningful only in the case of BDG or SLDAE functionals.
 * For SLDA and ASLDA the scattering length is assumed to be infinite, and the function is ignored.
 * @param ix x-coordinate from range [0,NX), to convert to Cartesian use: x = DX*(ix-NX/2)
 * @param iy y-coordinate from range [0,NY), to convert to Cartesian use: y = DY*(iy-NY/2),
 *           NOTE: in case of 1d code iy=0
 * @param iz z-coordinate from range [0,NZ), to convert to Cartesian use: z = DZ*(iz-NZ/2)
 *           NOTE: in case of 1d and 2d codes iz=0
 * @param it iteration number
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return value of the scattering length a(x,y,z,t).
 * */
double scattering_length(int ix, int iy, int iz, int it, double *params, size_t extra_data_size, void *extra_data)
{
    return 0; // TODO
    // return input->sclgth; // by default return value from input file.
    
    // however, here you can define your own prescription
    // it can be time and position dependent
    // ...
}

/**
 * ------------------------ FOR FUNCTIONAL == CUSTOMEDF ------------------------
 * */

/**
 * This function computes internal energy in case is CUSTOMEDF functional is selected.
 * Otherwise the function is ignored.
 * For more info see wiki pages. 
 * @param it iteration number
 * @param rho - density, computed according gpe_density(psi)
 * @param energy  (OUTPUT)
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if computation is successful, otherwise return error code. If nonzero value is returned the main code will terminate.
 * */
inline __device__ int compute_energy_gpe(int it, double rho, double *energy, double *params, size_t extra_data_size, void *extra_data)
{
    // TODO call scattering_length()
    // Density energy functional for unitary Fermi gas
    // see: Phys. Rev. A 90, 043638 (2014)
    (*energy) = 0.37*0.6*rho*pow(3.0*M_PI*M_PI*rho, 2.0/3.0)/2.; // unitary limit
    return 0;
}

/**
 * This function computes potentials defining Hamiltonian in case is CUSTOMEDF functional is selected.
 * Otherwise the function is ignored.
 * For more info see wiki pages. 
 * @param it iteration number
 * @param rho - density, computed according gpe_density(psi)
 * @param dEDFdn (OUTPUT) updated values as output (INPUT/OUTPUT) 
 * @param params array of input parameters, before call of this routine the params array is processed by process_params() routine
 * @param extra_data_size size of extra_data in bytes, if extra_data size=0 the optional data is not uploaded
 * @param extra_data optional set of data uploaded by load_extra_data()
 * @return 0 if computation is successful, otherwise return error code. If nonzero value is returned the main code will terminate.
 * */
inline __device__ int compute_potentials_gpe(int it, double rho, double *dEDFdn, double *params, size_t extra_data_size, void *extra_data)
{
    // TODO call scattering_length()
    // see: Phys. Rev. A 90, 043638 (2014)
    (*dEDFdn) = 0.37*pow(3.0*M_PI*M_PI*rho, 2.0/3.0)/2.0; // unitary limit
    return 0;
}

