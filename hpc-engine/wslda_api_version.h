/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * @author Gabriel Wlazlowski
 * @date 21.02.2022
 * */  

// -------------------------------- 20220221 --------------------------------
#if API_VERSION<20220221

#ifdef API_PROBLEM_DEFINITION
double referencekF(int it, wslda_density h_densities, double *params, size_t extra_data_size, void *extra_data)
{
    if(input->referencekF>0.0) return input->referencekF; // take it from input file
    
    // define here your prescription for computing kF
    // ...
    // default: extract max density and use it for definition of kF
    double max_dens=0.0, kF;
    int ixyz;
    for(ixyz=0; ixyz<h_densities.nx*h_densities.ny*h_densities.nz; ixyz++) 
        if(h_densities.rho_a[ixyz]+h_densities.rho_b[ixyz]>max_dens) max_dens=h_densities.rho_a[ixyz]+h_densities.rho_b[ixyz];
    
    // depending on dimensionality of the problem
    if(NY==1 && NZ==1) kF = 0.5*M_PI*max_dens;                // 1D
    else if(NZ==1)     kF = pow(2.0*M_PI*max_dens,1./2.);     // 2D
    else               kF = pow(3.*M_PI*M_PI*max_dens,1./3.); // 3D
    
    return kF;
}
#endif

#ifdef API_LOGGER
// empty
#endif

#endif

// -------------------------------- 20221120 --------------------------------
#if API_VERSION<20221120

#ifdef API_PROBLEM_DEFINITION
void modify_energies(int it, wslda_density h_densities, wslda_potential h_potentials, double *energy, double *params, size_t extra_data_size, void *extra_data)
{
    // DETERMINE LOCAL SIZES OF ARRAYS (CODE DIMENSIONALITY DEPENDENT)
    int lNX=h_densities.nx, lNY=h_densities.ny, lNZ=h_densities.nz; // local sizes
    int ix, iy, iz, ixyz;

    // extract volume element
    double volume_element=0.0;
    if(h_densities.datadim==3) volume_element=DX*DY*DZ;
    if(h_densities.datadim==2) volume_element=DX*DY*LZ;
    if(h_densities.datadim==1) volume_element=DX*LY*LZ;

    // ITERATE OVER ALL POINTS TO INTEGRATE OVER ALL POINTS
    double myE_contrib=0.0;
    ixyz=0;
    for(ix=0; ix<lNX; ix++) for(iy=0; iy<lNY; iy++) for(iz=0; iz<lNZ; iz++)
    {
        double x = DX*(ix-lNX/2);
        double y = DY*(iy-lNY/2); // for 1d code y will be always 0
        double z = DZ*(iz-lNZ/2); // for 1d and 2d codes z will be always 0

        // compute your contribution to the energy
        // myE_contrib += (...)*volume_element;

        ixyz++; // go to the next point, it should be the last line of the triple loop
    }

    // Add contributio to desired tag, for example
    // energy[EPOT]+=myE_contrib;
}
#endif


#ifdef API_LOGGER
double energy_unit(double kF, double *mu, double *npart,
           double *params, size_t extra_data_size, void *extra_data)
{
    double Effg;
    double eF = kF*kF/2.0; // Fermi energy
    double N = npart[SPINA]+npart[SPINB]; // total number of particles

    // depending on dimensionality of the problem
    if(NY==1 && NZ==1) Effg=(1./3.)*N*eF;   // 1D
    else if(NZ==1)     Effg=(1./2.)*N*eF;   // 2D
    else               Effg=(3./5.)*N*eF;   // 3D

    return Effg;
}
#endif

#endif

// -------------------------------- 20230426 --------------------------------
#if API_VERSION<20230426

#ifdef API_PROBLEM_DEFINITION
// empty
#endif

#ifdef API_LOGGER
// empty
#endif

#endif

// -------------------------------- 20231218 --------------------------------
#if API_VERSION<20231218

#ifdef API_PROBLEM_DEFINITION
// empty
#endif

#ifdef API_LOGGER

int add_custom_variable_to_wdata_metadata(wdata_metadata *wdmd,
           double *params, size_t extra_data_size, void *extra_data)
{
    return 0;
}

int write_custom_variable_to_wdata_set(wdata_metadata *wdmd,
           int it,
           wslda_density h_densities, wslda_potential h_potentials,
           double kF, double *mu,
           double *params, size_t extra_data_size, void *extra_data)
{
    return 0;
}

#endif

#endif
