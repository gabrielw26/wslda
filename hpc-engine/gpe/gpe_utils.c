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
#include <string.h>
#include <tgmath.h>
#include <complex.h>

#include "predefines.h"
#include "wdata.h"
#include "pca_utils.h"
#include "gpe_utils.h"

void set_initial_wave_function(uint nxyz, Complex *psi)
{
    uint ixyz;
    for(ixyz=0; ixyz<nxyz; ixyz++) 
    { 
        psi[ixyz].x = 1.0; psi[ixyz].y = 0.0; 
    }
}

int read_initial_wave_function(uint nxyz, Complex *psi, double* t0)
{
    char* psiFilename = (char*)malloc(strlen(input->inprefix) * sizeof(char));
    strcpy(psiFilename, input->inprefix);
    strcat(psiFilename, ".wtxt");
    printf("# Reading psi from wdata set: `%s`\n", psiFilename);
    wdata_metadata mdin;

    int ierr = wdata_parse_metadata_file(psiFilename, &mdin);
    if (ierr != 0)
    {
        printf("Cannot read metadata file! ERROR: #%d\n", ierr);
        return 1;
    }
    *t0 = mdin.t0 + (mdin.cycles - 1) * mdin.dt; 
    ierr = wdata_read_cycle(&mdin, "psi", mdin.cycles-1, psi);
    if (ierr != 0)
    {
        printf("ERROR: Cannot read psi!\n");
        return 1;
    }
    return 0;
}

void write_to_binary_file(uint nxyz, Complex *psi)
{
    FILE * psiFile;
    printf("# Writing psi to file\n");

    char* psiFilename = (char*)malloc(strlen(input->outprefix) * sizeof(char));
    strcpy(psiFilename, input->outprefix);
    switch (input->gpe_mode)
    {
    // output file for imaginary
    case 0:
        strcat(psiFilename, "_psi.dat");
        break;
    // output file for real
    case 1:
        strcat(psiFilename, "_psi2.dat");
        break;
    default:
        break;
    }

    psiFile = fopen (psiFilename, "wb");
    fwrite (psi , sizeof(Complex)*nxyz, 1, psiFile);
    fclose (psiFile);    
}

void write_to_txt_file(uint nx, uint ny, uint nz, Complex *psi)
{
    FILE * fout;
    uint ix, iy, iz, ixyz;

    char* psiFilename = (char*)malloc(strlen(input->outprefix) * sizeof(char));
    strcpy(psiFilename, input->outprefix);
    switch (input->gpe_mode)
    {
    // output file for imaginary
    case 0:
        strcat(psiFilename, "_psi.txt");
        break;
    // output file for real
    case 1:
        strcat(psiFilename, "_psi2.txt");
        break;
    default:
        break;
    }

    fout = fopen(psiFilename, "w");

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
}


void print_header_image()
{
    printf("#%7s %12s %12s %12s %12s %12s %12s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "diff", "comp.time");
}

void print_header_real()
{
    printf("#%7s %12s %12s %12s %12s %12s %12s\n", "time", "etot", "ekin", "eint", "eext", "(eint+eext)", "comp.time");
}

void print_results_image(double time, double npart, double etot, double ekin, double eint, double eext, double diff, double rt)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.6g %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, diff, rt);
}

void print_results_real(double time, double npart, double etot, double ekin, double eint, double eext, double rt)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f %12.4f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart, rt);  
}

void print_intial_results(double time, double npart, double etot, double ekin, double eint, double eext)
{
    printf("%8.2f %12.8f %12.8f %12.8f %12.8f %12.8f\n",time, etot/npart, ekin/npart, eint/npart, eext/npart, (eint+eext)/npart);  
}

void read_of_input_parameters(char* execcmd, int argc , char ** argv)
{
    int i;
    strcpy( execcmd , argv[ 0 ] ) ;
    for( i = 1 ; i < argc ; i++ ) 
    {
        strcat( execcmd , " " ) ; 
        strcat( execcmd , argv[ i ] ) ;
    }
}

int parse_command_line_and_get_idx_of_input_file(int argc , char ** argv)
{
    int i = readcmd( argc , argv ) ;
    if( i == -1 )
    {
        printf( "TERMINATING! NO INPUT FILE.\n" ) ;
        return(EXIT_FAILURE);
    }
    return i;
}

void read_input_file(int idx, char ** argv)
{
    // Read input file
    // Info from file is loaded into metadata structure
    int j = parse_input_file(argv[idx]);
    if ( j == 0 )
    {
        printf("PROBLEM WITH INPUT FILE: `%s`.\n" , argv[ idx ] ) ;
        exit(EXIT_FAILURE);      
    }
        
    // Input file tags are accessible through pointer `input`
}

void save_extradata_to_file_with_outprefix(size_t size, void *extra_data, char* outprefix)
{
    char fname[1024];
    sprintf(fname, "%s_extra_data.dat", outprefix);
    wprintf("# SAVING EXTRA_DATA TO FILE: %s\n", fname);
    FILE * f = fopen(fname, "wb");
    fwrite(extra_data, size, 1, f);
    fclose(f);
}

int malloc_extra_data(size_t extra_data_size, void *extra_data)
{
    wprintf("# EXTRA_DATA IS ACTIVE.\n");
    wprintf("# ALLOCATING EXTRA_DATA OF SIZE %ld B.\n", extra_data_size); fflush(stdout);
    if ( ( extra_data = (void *) malloc( extra_data_size ) ) == NULL  )
    {                                                             
        wfprintf( stderr , "error: cannot malloc()! Exiting!\n") ; 
        wfprintf( stderr , "error: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; 
        return( EXIT_FAILURE ) ; 
    }
    wprintf("# EXECUTING: load_extra_data(%zu, extra_data, input->params)\n", extra_data_size);
    return 0;
}