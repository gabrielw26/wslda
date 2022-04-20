/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * @author Gabriel Wlazlowski
 * @date 21.02.2022
 * */  

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
#endif
