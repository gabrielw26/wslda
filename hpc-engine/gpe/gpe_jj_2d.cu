/***************************************************************************
 *   Copyright (C) 2015 by                                                 *
 *   WARSAW UNIVERSITY OF TECHNOLOGY                                       *
 *   FACULTY OF PHYSICS                                                    *
 *   NUCLEAR THEORY GROUP                                                  *
 *   See also AUTHORS file                                                 *
 *                                                                         *
 *   This file is a part of GPE for GPU project.                           *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/ 
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <math.h>
#include <complex>      // std::complex

/***************************************************************************/ 
/**************************** GPE HEADERS **********************************/
/***************************************************************************/
#include "gpe_engine.h"
#include "gpe_timing.h"

// For writing measurements
#include "gpe_io.h"

typedef std::complex<double> cplx;

/***************************************************************************/ 
/**************************** USER DEFINED *********************************/
/***************************************************************************/
// See: gpe_user_defined.h
#include "gpe_user_defined.h"

void compute_population_imbalance(Complex *psi, double *rho, double *outNL, double *outNR, double *outphaseL, double *outphaseR)
{
    int ix, iy, iz, ixyz=0;
    double NL=0.0;
    double NR=0.0;
    double phaseL=0.0;
    double phaseR=0.0;
    
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        if(1.0*(ix-NX/2)+0.5<0.0)  
        {
            NL+=rho[ixyz];
            cplx aaa;
            phaseL+=rho[ixyz]*atan2(psi[ixyz].y,psi[ixyz].x);
        }
        else                       
        {
            NR+=rho[ixyz];
            phaseR+=rho[ixyz]*atan2(psi[ixyz].y,psi[ixyz].x);
        }
        ixyz++;
    }
    
    double z=(NL-NR)/(NL+NR);
    phaseL/=(NL*M_PI);
    phaseR/=(NR*M_PI);
    
    // fix
    ixyz= NZ/2 + NZ*NY/2 + NZ*NY*(NX/2-5);
    phaseL=atan2(psi[ixyz].y,psi[ixyz].x)/M_PI;
    ixyz= NZ/2 + NZ*NY/2 + NZ*NY*(NX/2+5);
    phaseR=atan2(psi[ixyz].y,psi[ixyz].x)/M_PI;   
    
    // save results 
    *outNL = NL;
    *outNR = NR;
    *outphaseL = phaseL;
    *outphaseR = phaseR;
}

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
int main( int argc , char ** argv ) 
{
    // This code does simulation of Josephson juction experiment
    // see: Phys. Rev. Lett. 124, 045301 (2020)
    
    if(argc != 3) 
    {
        printf("%s deviceNo prefix\n", argv[0]);
        return 1;
    }
    
    // SETTINGS
    double alpha;
    double beta;
    double dt=0.025;
//     double npart=60000.0*2;
    double npart=348685.847482;
    int steps=100;
    int device = atoi(argv[1]);
    
    char *prefix = argv[2];
    char file_name[128];
//     sprintf(prefix, "jj.V%2.1f.w%2.1f.p%2.1f.psi", atof(argv[2]), atof(argv[3]), atof(argv[4]));
//     sprintf(prefix, "torm");
    printf("# OUTPUT FILES PREFIX: `%s`\n", prefix);
    
    int ierr;
    cudaError err;
    
    printf("# TAKING DEVICE: %d\n", device);
    err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        return 1;
    }
    
    int nx, ny, nz;
    gpe_get_lattice(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    
    uint nxyz=nx*ny*nz;
    uint ixyz;
    
    double ekin, eint, eext, etot, etot_prev;
    double time;
    double rt;
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , 2*sizeof(double)*nxyz*2, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    double *rho; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &rho , sizeof(double)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
       
    double *jx = (double *)(psi); // reuse memory
    double *jy = jx+nx*ny*nz; // reuse memory
    double *jz = jx+2*nx*ny*nz; // reuse memory
    
    // ***************** Imaginary time projection **************************
    printf("# IMAGINARY TIME PROJECTION - WITH BARRIER\n");
    alpha=0.0;
    beta=1.0;
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
            
    // Set initial wave function
    for(ixyz=0; ixyz<nxyz; ixyz++) { psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; } // Fill with initial values
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;
    
    // Normalize state
    gpe_exec( gpe_normalize_psi(), ierr );
    
    // user defined params
    #define MY_USER_PARAMS 12
    double params[MY_USER_PARAMS];
    
    // estimate frequency of trap to get desired eF 
        
    double lx = 1.*NX/48.; // according: Phys. Rev. Lett. 124, 045301 (2020)
    double omegax = 1./(lx*lx);
    params[0] = 0.5 * omegax * omegax;
    
//     double Nsigma = 0.5*npart;
    double Nsigma = 60000.0;
    double omega_ho = omegax*pow(OMEGA_YX*OMEGA_ZX, 1./3.);
    double EFho=pow(6.*Nsigma, 1./3.)*omega_ho;
    double KFho = sqrt(2.0*EFho);
    params[10] = 1./(4.6*KFho); // according: Phys. Rev. Lett. 124, 045301 (2020)
    printf("# EFho=%f, KFho=%f, a=1/(4.6*kF)=%f\n", EFho, KFho, params[10]);
    double eF=EFho; 
    double pF=sqrt(2.*eF);
    
//     double mu=114.0*omegax;
    double mu = 0.4*EFho;
    double V0_per_mu=1.0; // according: Phys. Rev. Lett. 124, 045301 (2020)
    double cohlen = 2.0 / (M_PI * 0.5 * KFho);
    printf("# cohlen=%f\n", cohlen);
//     double cohlen = 0.067*lx; // according: Phys. Rev. Lett. 124, 045301 (2020)
    double w=4.0*cohlen; // according: Phys. Rev. Lett. 124, 045301 (2020)
    params[1] = V0_per_mu*mu; // barrier height
    params[2] = 2.0/(w*w); 
    params[3] = 0.0005*mu;

    params[4] = 0.0; // 2D case
    
    printf("# BARRIER PARAMETERS: V0_per_mu=%f, w=%f, w/cohlen=%f, EFho=%f\n", V0_per_mu, w, w/cohlen, EFho);
    
    // Copy parameters to engine
    // Copy table to d_user_param on device side.
    gpe_exec( gpe_set_user_params(MY_USER_PARAMS, params), ierr );
            
    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %14s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "diff", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    fflush(stdout);
    
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 
        gpe_exec( gpe_evolve(steps), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(); // get time
        
        etot_prev=etot;
        etot = ekin + eint + eext;
        double diff=(etot_prev-etot)/npart;
        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);  
        fflush(stdout);
        
        if(fabs(diff)<1.0e-9) break;
    }
   
    // Get psi and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;
        
    int ix, iy, iz;
    char prefix_cross[128];
    sprintf(prefix_cross, "%s.cross.nb", prefix);
    printf("# Writing sections to text file : `%s`\n", prefix_cross);
    FILE * fout = fopen(prefix_cross, "w");
    iy=ny/2;
    iz=nz/2;
    for(ix=0; ix<nx; ix++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, 2.0*(psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, 2.0*(psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, 2.0*(psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fclose(fout);
    
    // Find max density
    double rho0=0.0;
    for(ixyz=0; ixyz<nxyz; ixyz++) if(rho0<2.0*(psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y)) rho0=2.0*(psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y);
    double pf0 = pow(3.*M_PI*M_PI*rho0,1./3.);
    double ef0 = 0.5*pf0*pf0;    

    printf("# MAX INFO: rho=%f, eF=%f, pF=%f, V0/eF=%f\n", rho0, ef0, pf0, params[1]/ef0);
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;           
    gpe_exec( gpe_get_density(&time, rho), ierr ) ;

    
    double NL, NR, phaseL, phaseR;
    double z;
    compute_population_imbalance(psi, rho, &NL, &NR, &phaseL, &phaseR); z=(NL-NR)/(NL+NR);
        
//     gpe_exec( gpe_get_currents(&time, jx, jy, jz), ierr ) ;
// //         for(ixyz=0; ixyz<3*nxyz; ixyz++) jx[ixyz]*=0.5;
//     sprintf(file_name, "%s_current_a.dpca", prefix);
//     gpe_exec( create_measurement_file_with_header(file_name, nx, ny, nz, 1.0, 1.0, 1.0, ef0, 0.0, dt*steps), ierr) ;
//     gpe_exec( add_measurement_entry(file_name, jx, 3*sizeof(double)*nx*ny*nz), ierr) ; 

            
    // ***************** Real time evolution **************************
    printf("# REAL TIME EVOLUTION\n");
    double T_perid = 2.0*M_PI/omegax;
    double steps_perid = T_perid / dt; 
    int points_per_perid=100;
//     steps = steps_perid / points_per_perid;
    int max_periods=5;
    printf("# OSCILLATION PERIOD: %f (STEPS: %.2f)\n", T_perid*eF, steps_perid);
    printf("# PRINTING EVERY %d STEPS\n", steps);
    
    // Before saving - reduce data amount
    iz=NZ/2;
    double testnpart=0.0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        psi[iy+ny*ix] = psi[ixyz];
        rho[iy+ny*ix] = rho[ixyz];
        testnpart+=rho[ixyz];
    }
//     printf("# TESTNPART=%f, targe-npart=%f\n", testnpart, testnpart*NZ);
    
    sprintf(file_name, "%s_delta.dpca", prefix);
    gpe_exec( create_measurement_file_with_header(file_name, nx, ny, 1, 1.0, 1.0, 1.0, eF, 0.0, dt*steps), ierr) ;
    gpe_exec( add_measurement_entry(file_name, psi, sizeof(Complex)*nx*ny), ierr) ;
    
    sprintf(file_name, "%s_density_a.dpca", prefix);
    gpe_exec( create_measurement_file_with_header(file_name, nx, ny, 1, 1.0, 1.0, 1.0, eF, 0.0, dt*steps), ierr) ;
    gpe_exec( add_measurement_entry(file_name, rho, sizeof(double)*nx*ny), ierr) ;
    
    // Reset time
    gpe_exec( gpe_set_time(0.0), ierr );
    
    // Update parameters
    params[3] = 0.0; // no linear potential
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(MY_USER_PARAMS, params), ierr );
    
    // switch from imaginary time to real time evolution
    alpha=1.0;
    beta=0.0;
    gpe_exec( gpe_change_alpha_beta(alpha, beta), ierr);
    
    // evolution
    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", "t*eF", "etot", "NL", "NR", "(NL+NR)", "z", "phaseL","phaseR","Dphase","comp.time");
    time=0.0; // reset time
    rt = 0.0; // reset 
    printf("%8.2f %12.8f %12.6f %12.6f %12.6f %12.8f %12.4f %12.6f %12.6f %12.8f\n",time*eF, etot/npart, NL, NR, NL+NR, z, phaseL, phaseR, phaseL-phaseR, rt);  
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 
        gpe_exec( gpe_evolve(steps), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        // Add mesurement to file
        gpe_exec( gpe_get_psi(&time, psi), ierr ) ;        
        gpe_exec( gpe_get_density(&time, rho), ierr ) ;
        
        compute_population_imbalance(psi, rho, &NL, &NR, &phaseL, &phaseR); z=(NL-NR)/(NL+NR);
        
        // Before saving - reduce data amount
        iz=NZ/2;
        for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++)
        {
            ixyz = iz + nz*iy + nz*ny*ix;
            psi[iy+ny*ix] = psi[ixyz];
            rho[iy+ny*ix] = rho[ixyz];
        }
        sprintf(file_name, "%s_delta.dpca", prefix);
        gpe_exec( add_measurement_entry(file_name, psi, sizeof(Complex)*nx*ny), ierr) ;
        
        sprintf(file_name, "%s_density_a.dpca", prefix);
        gpe_exec( add_measurement_entry(file_name, rho, sizeof(double)*nx*ny), ierr) ; 
        
        rt = e_t(); // get time

        etot = ekin + eint + eext;
        printf("%8.2f %12.8f %12.6f %12.6f %12.6f %12.8f %12.4f %12.6f %12.6f %12.8f\n",time*eF, etot/npart, NL, NR, NL+NR, z, phaseL, phaseR, phaseL-phaseR, rt);  
        fflush(stdout);

//         if(time>T_perid*max_periods) break;
        if(time*eF>2000) break;
    }
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
   
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}
