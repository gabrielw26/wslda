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
#include "predefines.h"
#include "gpe_timing.h"
//#include "pca_utils.h"
//#include "wslda_toolkit.h"
// -> include "pca_utils.h", "gpe_user_defined.h"
//extern int wsldapid;

/***************************************************************************/ 
/************************* MAIN FUNCTION  **********************************/
/***************************************************************************/
extern "C" int gpe_real()
{
    // SETTINGS
    double alpha=1.0;
    double beta=0.0;
    double dt=0.025;
    double npart=1000.0;
    const int device=0;  

    int ierr;
    cudaError err;
    
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
    uint ix, iy, iz, ixyz;
    
    double ekin, eint, eext, etot;
    double time;
    double rt;
  
    // TODO
    // ***************** Real time evolution **************************
    printf("# REAL TIME EVOLUTION\n");
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    
    // Set initial wave function - read it from file psi.dat (see gpe_imag.cu)
    printf("# Reading psi from file\n");
    FILE * psiFile;
    psiFile = fopen ("psi.dat", "rb");
    size_t readok = fread (psi , sizeof(Complex)*nxyz, 1, psiFile);
    if (readok != 1)
    {
        printf("Reading error\n");
        return 1;
    }
    fclose (psiFile);  
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    
    // Prepare user defined parameters
    double params[3];
    params[0] = 0.01; // omega_x
    params[1] = 0.10; // omega_y
    params[2] = 0.11; // omega_z
    
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(3, params), ierr );
    
    // Copy wave function to GPU
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;

    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    
    // evolve in real time
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 100 steps forward
        gpe_exec( gpe_evolve(100), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(0); // get time
        
        etot = ekin + eint + eext;

        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, rt);  
        
        if(time>100.) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;    
            
    // write to binary file
    printf("# Writing psi to file\n");
    psiFile = fopen ("psi2.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);   
    
    // Write to txt file |Psi(x,y,x)|^2 - sections along x, y and z axis
    FILE * fout = fopen("psi2.txt", "w");
    iy=ny/2;
    iz=nz/2;
    for(ix=0; ix<nx; ix++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fclose(fout);
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}

extern "C" int gpe_imag()
{
    // SETTINGS
    double alpha=0.0;
    double beta=1.0;
    double dt=0.025;
    double npart=1000.0;
    int steps=100;
    int device=0;
    
    int ierr;
    cudaError err;
    
    err=cudaSetDevice( device );
    if(err != cudaSuccess) 
    {
        printf("this is error: %d\n", err);
        printf("Error: Cannot cudaSetDevice(%d)!\n", device);
        return 1;
    }
    
    int nx, ny, nz;
    gpe_get_lattice(&nx, &ny, &nz);
    printf("# GPE engine compiled for lattice: %d x %d x %d\n", nx, ny, nz);
    
    uint nxyz=nx*ny*nz;
    uint ix, iy, iz, ixyz;
    
    double ekin, eint, eext, etot, etot_prev;
    double time;
    double rt;
  
    // ***************** Imaginary time projection **************************
    printf("# IMAGINARY TIME PROJECTION\n");
    
    // CPU memory for wave function - pinned for fast transfers
    Complex *psi; // Complex type defined in gpe_engine.h - structure with two doubles x and y for real and imaginary parts
    err=cudaHostAlloc( &psi , sizeof(Complex)*nxyz, cudaHostAllocDefault );
    if(err != cudaSuccess) 
    {
        printf("Error: Cannot allocate memory!\n");
        return 1;
    }
    
    // Set initial wave function - constant
    for(ixyz=0; ixyz<nxyz; ixyz++) { psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; } // Fill with initial values
    
    // Create engine
    gpe_exec( gpe_create_engine(alpha, beta, dt, npart), ierr );
    
    // Prepare user defined parameters
    double params[3];
    params[0] = 0.01; // omega_x
    params[1] = 0.10; // omega_y
    params[2] = 0.11; // omega_z
    
    // Copy parameters to engine
    gpe_exec( gpe_set_user_params(3, params), ierr );
    
    // Copy wave function to GPU and normalize
    gpe_exec( gpe_set_psi(0.0, psi), ierr ) ;
    gpe_exec( gpe_normalize_psi(), ierr );

    // For nice printing
    printf("#%7s %12s %12s %12s %12s %12s %14s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "diff", "comp.time");
    
    // initial energy
    gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
    etot = ekin + eint + eext;
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
    
    // evolve in imaginary time until convergence is achieved
    while(1)
    {
        b_t(); // reset timer
        
        // Evolve 100 steps forward
        gpe_exec( gpe_evolve(100), ierr ); 
        
        // Compute energy 
        gpe_exec( gpe_energy(&time, &ekin, &eint, &eext), ierr );
        
        rt = e_t(0); // get time
        
        etot_prev=etot;
        etot = ekin + eint + eext;
        double diff=(etot_prev-etot)/npart; // diference in energy per particle
        printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);  
        
        if(fabs(diff)<1.0e-9) break;
    }    
    
    // Get wave function and save it to file
    gpe_exec( gpe_get_psi(&time, psi), ierr ) ;    
            
    // write to binary file
    FILE * psiFile;
    psiFile = fopen ("psi.dat", "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);   
    
    // Write to txt file |Psi(x,y,x)|^2 - sections along x, y and z axis
    FILE * fout = fopen("psi.txt", "w");
    iy=ny/2;
    iz=nz/2;
    for(ix=0; ix<nx; ix++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", ix-nx/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n");
    ix=nx/2;
    iz=nz/2;
    for(iy=0; iy<ny; iy++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iy-ny/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fprintf(fout, "\n\n"); 
    ix=nx/2;
    iy=ny/2;
    for(iz=0; iz<nz; iz++)
    {
        ixyz = iz + nz*iy + nz*ny*ix;
        fprintf(fout, "%6d %12.6g\n", iz-nz/2, (psi[ixyz].x*psi[ixyz].x + psi[ixyz].y*psi[ixyz].y) );
    }
    
    fclose(fout);
    
    // Destroy engine
    gpe_exec( gpe_destroy_engine(), ierr) ;
    
    // Clear memory
    cudaFreeHost(psi);
    
    return 0;
}