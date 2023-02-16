/**
 * W-SLDA Toolkit
 *
 * Simple code showing how to compute observables like the energy of the flow and the condensation energy for 3d data
 *
 * Copy this file to your project folder and compile using:
 *    gcc -std=gnu99 compute-observables-3d.c -I. -I$WSLDA/hpc-engine -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdata -lm -o compute-observables-3d
 *
 * */

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include "wdata.h"

void something_to_cheer_you_up_pid0(FILE *stream) {} // otherwise, pca_macro.h does not compile
#include "pca_macro.h"

int main( int argc , char ** argv )
{
    char *inprefix;
    if(argc != 2)
    {
        printf("%s prefix\n", argv[0]);
        return 1;
    }

    char file_name[256], report[256];
    sprintf(file_name, "%s.wtxt", argv[1]);         // <-- input
    sprintf(report, "%s.observables.txt", argv[1]); // <-- output

    // create metadata handler
    wdata_metadata md;
    int ierr;

    // read metadata from file
    printf("# Reading file `%s`\n", file_name);
    ierr = wdata_parse_metadata_file(file_name, &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;}

    // get basic info
    if(md.datadim!=3)
    {
        printf("# Error: this code does computation for 3D data [datadim=%d]!\n", md.datadim);
        return 1;
    }
    int nx=md.nx, ny=md.ny, nz=md.nz;
    int nxyz = nx*ny*nz;
    double dx=md.dx, dy=md.dy, dz=md.dz;
    double dxyz = dx*dy*dz;
    int ix, iy, iz, ixyz;
    int nom=md.cycles;
    double dt=md.dt, t0=md.t0, eF, kF;
    eF = wdata_getconst_value(&md, "eF");
    kF = wdata_getconst_value(&md, "kF");
    printf("# Lattice: %d x %d x %d with spacing %f x %f x %f\n", nx, ny, nz, dx, dy, dz);

    // arrays

    // order parameter
    double complex *delta;
    cppmallocl(delta, nxyz, double complex);

    // density, I assume here that the system spin symmetric
    double *rho;
    cppmallocl(rho, nxyz, double);

    double *j;
    cppmallocl(j, nxyz*3, double);

    // prepare file for the report
    printf("# Creating file: %s\n", report);
    FILE *fout = fopen(report, "w");
    fprintf(fout, "# COLUMNS\n");
    fprintf(fout, "# 1: cycle\n");
    fprintf(fout, "# 2: time*eF\n");
    fprintf(fout, "# 3: N\n");
    fprintf(fout, "# 4: Eflow=\\int j(r)^2/2n(r) d^3r\n");
    fprintf(fout, "# 5: Econd=\\int (3/8) \\Delta(r)^2*n(r)/eF(r) d^3r\n");
    fprintf(fout, "# CONSTS\n");
    fprintf(fout, "# kF=%f\n", kF);
    fprintf(fout, "# eF=%f\n", eF);
    int inom;
    for(inom=0; inom<md.cycles; inom++) // for each record
    {
        double time = md.t0+md.dt*inom; // time
        wdata_read_cycle(&md, "delta", inom, delta);
        wdata_read_cycle(&md, "rho_a", inom, rho);
        wdata_read_cycle(&md, "j_a", inom, j);
        for(ixyz=0; ixyz<nxyz; ixyz++) rho[ixyz]*=2.0; // I assume the system is spin symmetric
        for(ixyz=0; ixyz<nxyz*3; ixyz++) j[ixyz]*=2.0; // I assume the system is spin symmetric

        // observables
        double Ntot = 0.0; // \int n d^3r, particle number
        for(ixyz=0; ixyz<nxyz; ixyz++)
            Ntot+=rho[ixyz] * dxyz;

        double Eflow = 0.0; // \int j^2/2n d^3r, flow energy
        for(ixyz=0; ixyz<nxyz; ixyz++)
            Eflow+=(pow(j[ixyz+0*nxyz],2)+pow(j[ixyz+1*nxyz],2)+pow(j[ixyz+2*nxyz],2)) / (2.0*rho[ixyz]+1.0e-9) * dxyz;

        double Econd = 0.0; // \int (3/8) \Delta^2*n/eF d^3r, BCS condensation energy
        for(ixyz=0; ixyz<nxyz; ixyz++)
        {
            double leF = 0.5*pow(3.0*M_PI*M_PI*rho[ixyz], 2./3.); // local Fermi energy
            double ldelta=cabs(delta[ixyz]);
            Econd+=(3./8.)*pow(ldelta,2)*rho[ixyz]/(leF+1.0e-9) * dxyz;
        }

        fprintf(fout, "%6d %12.6f %12.6f %12.6f %12.6f\n", inom, time*eF, Ntot, Eflow, Econd);
    }

    // close file
    fclose(fout);
    printf("# Done\n");

    return 0;
}

