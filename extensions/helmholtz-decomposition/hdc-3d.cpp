// W-SLDA Toolkit
//
// Documentation: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Helmholtz-decomposition-code
//
// Author: Andrea Barresi (Warsaw University of Technology, 2020)
//
// Compilation 			     =>  mpic++ hd-3d.cpp -o wslda-hd-3d -L../lib/wdata -I../lib/wdata/c -lwdata -lfftw3 -lm
// Execution (Unitarity)	=> mpirun -n 50 ./mpihh /home2/scratch/td-qt-48-study/run2.wtxt 0
//

// very small number
#define EPS 1.0e-14

#include <fftw3.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sstream>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include <math.h>
#include <numeric>
#include <complex>
#include <fstream>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <iomanip>
#include <mpi.h>

// Custom libraries
#include "wdata.h"

#define cppmallocl(pointer, size, type)                                     \
    if ((pointer = (type *)malloc((size) * sizeof(type))) == NULL)          \
    {                                                                       \
        fprintf(stderr, "error: cannot malloc()! Exiting!\n");              \
        fprintf(stderr, "error: file=`%s`, line=%d\n", __FILE__, __LINE__); \
        return -1;                                                          \
    }
    
#define S 70000
#define Master 0
static double t_gettimeofday ;
static struct timeval s;
using namespace std;						
typedef complex<double> comp;		

double Omega(double *k1, double k, double dk);
void b_t();
double e_t();

int main(int argc, char *argv[])
{
      comp I = comp(0.,1.);	double t=2.*M_PI;     //cout<<endl;

	// Create metadata handler
    wdata_metadata md;
 	md.datadim=3; // 3D data

	// Initilize MPI
    int ip = 0, np = 0;                 	// Basic MPI indicators
    MPI_Init(&argc, &argv);              
    MPI_Comm_size(MPI_COMM_WORLD, &np);  	// Number of processes 
    MPI_Comm_rank(MPI_COMM_WORLD, &ip); 	// id of process st 0 <= pid < np 
	MPI_Status status, check;
	int count;
	b_t();
  
	// Compile by giving a in-filename prefix
	if(ip==0)	cout<<endl;
	if(ip==0)
		if(argc!=3)
		{
			printf("Usage:\n\t %s [Path to .wtxt file] [Spin component]\n\n", argv[0]);
			printf("Spin components:\n\t 0 - a (spin up) \n\t 1 - b (spin down) \n\t 2 - a+b (both components) \n%s");
			exit(EXIT_FAILURE);
		}

	// Component naming
	char si[8], sc[8];
	char *pipo = argv[2];
	char u;
	int c=atoi(pipo);
	if(c==0)		{ u='a'; sprintf(si,"wi_a"); sprintf(sc,"wc_a"); }
	else if(c==1)	{ u='b'; sprintf(si,"wi_b"); sprintf(sc,"wc_b"); }
	else if(c==2)	{ u='t'; sprintf(si,"wi_t"); sprintf(sc,"wc_t"); }

	// Outputs file to Dwarf repository; change path for different machines
	char file_name[256], file_name2[256];
	char fwi[256];
	char fwc[256];

	// Lattice
	int NX, NY, NZ;
	double *kx, *ky, *kz, *kp;
	double k; 
	int ierr, nom, inom=1, track=0;
	double t0, dt, eF, cycles;
    
	if(ip==0)
	{
		// Read metadata file
		ierr = wdata_parse_metadata_file(argv[1], &md);
		if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;}

		// Clear wi & wc if already present
		wdata_variable test;
		ierr =  wdata_get_variable(&md, si, &test);
		if(ierr==0) 	
        {
            printf("\t Variables already exists! Cleaning data sets...\n");
            wdata_clear_file(&md, si);
            wdata_clear_file(&md, sc);
        }
		else		// Add Helmholtz variables to data set ???
		{	
            printf("\t Adding new variables [%s, %s] to w-data set.\n", sc, si);
            wdata_add_comment_to_metadata_file(argv[1], "Variables added by wslda-hdc-3d");
			wdata_variable vwc = {"wc", "vector", "none", "wdat"};
            sprintf(vwc.name,sc);
            wdata_add_variable(&md, &vwc);
            wdata_add_var_to_metadata_file(argv[1], &vwc); 
			wdata_variable vwi = {"wi", "vector", "none", "wdat"};
            sprintf(vwi.name,si);
            wdata_add_variable(&md, &vwi);
            wdata_add_var_to_metadata_file(argv[1], &vwi);
		}
	
		// Some checks just in case
		printf("\t Number of measurements = %d\n", md.cycles);printf("\t Starting time 		= %.0lf \n", md.t0);
		printf("\t Lattice size 	= %d x %d x %d \n", md.nx, md.ny, md.nz);
		printf("\t Ending time 		= %.3lf \n", md.t0+md.cycles*md.dt);
		printf("\t Component = %c \n", u);
		printf("\t Variables being produced = %s & %s \n\n", si, sc);
		if(c>2)
		{
			cout<<"\n Please select a suitable component for inspection and run again.\n"<<endl;
			printf("Spin components:\n\t 0 - a (spin up) \n\t 1 - b (spin down) \n\t 2 - a+b (both components) \n\n");
			return 0;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast (&md, sizeof(wdata_metadata), MPI_BYTE, 0, MPI_COMM_WORLD);
    if(md.datadim!=3)
    {
        if(ip==0) 
        {   
            printf("The code supports only 3D data!\n");
            printf("For more info see wiki pages:\n");
            printf("https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Helmholtz-decomposition-code!\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();	
        return 0;
    }
	int nxyz = md.nx * md.ny * md.nz;
	double dk=2.*M_PI/(md.dx*md.nx);
	cppmallocl(kx, nxyz, double);
	cppmallocl(ky, nxyz, double);
	cppmallocl(kz, nxyz, double);
	cppmallocl(kp, nxyz, double);

	for(int iz=0; iz<md.nz; iz++)
	{
		for(int iy=0; iy<md.ny; iy++)				
		{
			for(int ix=0; ix<md.nx; ix++)
			{
				if (ix<md.nx/2)
					kx[iz + iy*md.nz + ix*md.ny*md.nz] = (t * double(ix) / (md.nx*md.dx));
				else
					kx[iz + iy*md.nz + ix*md.ny*md.nz] = (t * double(ix-md.nx) / (md.nx*md.dx));

				if (iy<md.ny/2)
					ky[iz + iy*md.nz + ix*md.ny*md.nz] = (t * double(iy) / (md.ny*md.dy));
				else
					ky[iz + iy*md.nz + ix*md.ny*md.nz] = (t * double(iy-md.ny) / (md.ny*md.dy));

				if (iz<md.nz/2)
					kz[iz + iy*md.nz + ix*md.ny*md.nz] = (t * double(iz) / (md.nz*md.dz));
				else
					kz[iz + iy*md.nz + ix*md.ny*md.nz] = (t * double(iz-md.nz) / (md.nz*md.dz));
				
				kp[iz + iy*md.nz + ix*md.ny*md.nz] = sqrt(kx[iz + iy*md.nz + ix*md.ny*md.nz]*kx[iz + iy*md.nz + ix*md.ny*md.nz] + 
											ky[iz + iy*md.nz + ix*md.ny*md.nz]*ky[iz + iy*md.nz + ix*md.ny*md.nz] + 
											kz[iz + iy*md.nz + ix*md.ny*md.nz]*kz[iz + iy*md.nz + ix*md.ny*md.nz]);
			}
		}
	}
	
	// DENSITY
	double *rho_a;
	cppmallocl(rho_a, nxyz, double);	
	double *rho_b;
	cppmallocl(rho_b, nxyz, double);

	// CURRENTS
	double *j_a;
	cppmallocl(j_a, 3 * nxyz, double);
	double *j_b;
	cppmallocl(j_b, 3 * nxyz, double);
	
	// Currents' components
	double *j_a_x = j_a + 0*nxyz;
	double *j_a_y = j_a + 1*nxyz;
	double *j_a_z = j_a + 2*nxyz;
	
	double *j_b_x = j_b + 0*nxyz;
	double *j_b_y = j_b + 1*nxyz;
	double *j_b_z = j_b + 2*nxyz;

	double *Q, q=0, *D, d=0;
	cppmallocl(Q,  nxyz, double);
	cppmallocl(D,  nxyz, double);

	// Weighted velocities and components
	double *waxr, *wayr, *wazr, *waxi, *wayi, *wazi;
	cppmallocl(waxr,  nxyz, double);
	cppmallocl(wayr,  nxyz, double);
	cppmallocl(wazr,  nxyz, double);
	cppmallocl(waxi,  nxyz, double);
	cppmallocl(wayi,  nxyz, double);
	cppmallocl(wazi,  nxyz, double);
	double *lxr = (double*) waxr;
	double *lyr = (double*) wayr;
	double *lzr = (double*) wazr;
	double *lxi = (double*) waxi;
	double *lyi = (double*) wayi;
	double *lzi = (double*) wazi;
	
	fftw_complex *wax, *way, *waz;
	cppmallocl(wax, 2* nxyz, fftw_complex);
	cppmallocl(way, 2* nxyz, fftw_complex);
	cppmallocl(waz, 2* nxyz, fftw_complex);
	comp *lx = (comp*) wax;		
	comp *ly = (comp*) way;
	comp *lz = (comp*) waz;						
	
	fftw_complex *temppx, *temppy, *temppz;
	cppmallocl(temppx, 2* nxyz, fftw_complex);
	cppmallocl(temppy, 2* nxyz, fftw_complex);
	cppmallocl(temppz, 2* nxyz, fftw_complex);
	comp *px = (comp*) temppx;
	comp *py = (comp*) temppy;
	comp *pz = (comp*) temppz;

	fftw_complex *tempfx, *tempfy, *tempfz;					
	cppmallocl(tempfx, 2* nxyz, fftw_complex);
	cppmallocl(tempfy, 2* nxyz, fftw_complex);
	cppmallocl(tempfz, 2* nxyz, fftw_complex);
	comp *fx = (comp*) tempfx;	
	comp *fy = (comp*) tempfy;
	comp *fz = (comp*) tempfz;			

	fftw_complex *tgx, *tgy, *tgz, *tp;
	cppmallocl(tgx, 2* nxyz, fftw_complex);							
	cppmallocl(tgy, 2* nxyz, fftw_complex);	
	cppmallocl(tgz, 2* nxyz, fftw_complex);
	comp *gx = (comp*) tgx;	
	comp *gy = (comp*) tgy;
	comp *gz = (comp*) tgz;

	fftw_complex *tempcx, *tempcy, *tempcz;							
	cppmallocl(tempcx, 2* nxyz, fftw_complex);							
	cppmallocl(tempcy, 2* nxyz, fftw_complex);	
	cppmallocl(tempcz, 2* nxyz, fftw_complex);
	comp *wcx = (comp*) tempcx;	
	comp *wcy = (comp*) tempcy;
	comp *wcz = (comp*) tempcz;					

	fftw_complex *tempix, *tempiy, *tempiz;					
	cppmallocl(tempix, 2* nxyz, fftw_complex);	
	cppmallocl(tempiy, 2* nxyz, fftw_complex);	
	cppmallocl(tempiz, 2* nxyz, fftw_complex);
	comp *hx = (comp*) tempix;
	comp *hy = (comp*) tempiy;
	comp *hz = (comp*) tempiz;
    
	fftw_complex *tempvx, *tempvy, *tempvz;							
	cppmallocl(tempvx, 2* nxyz, fftw_complex);	
	cppmallocl(tempvy, 2* nxyz, fftw_complex);	
	cppmallocl(tempvz, 2* nxyz, fftw_complex);
	comp *wix = (comp*) tempvx;
	comp *wiy = (comp*) tempvy;
	comp *wiz = (comp*) tempvz;
	double Ec=0, Ei=0, Eck=0, Eik=0, Ek=0, U=0, ecv=0, eiv=0, EcSpec=0, EiSpec=0, OmegaVolume=0;
	double *ecvprint, *eivprint;
	int print_lght=(int)(M_PI/(md.dx*dk));
	cppmallocl(ecvprint, print_lght, double);
	cppmallocl(eivprint, print_lght, double);
	
	double *wixr, *wiyr, *wizr, *wcxr, *wcyr, *wczr;
	cppmallocl(wixr, nxyz, double);
	cppmallocl(wiyr, nxyz, double);
	cppmallocl(wizr, nxyz, double);
	cppmallocl(wcxr, nxyz, double);
	cppmallocl(wcyr, nxyz, double);
	cppmallocl(wczr, nxyz, double);

	fftw_complex *tempox, *tempoy, *tempoz;							
	cppmallocl(tempox, 2* nxyz, fftw_complex);	
	cppmallocl(tempoy, 2* nxyz, fftw_complex);	
	cppmallocl(tempoz, 2* nxyz, fftw_complex);
	comp *wkix = (comp*) tempvx;
	comp *wkiy = (comp*) tempoy;
	comp *wkiz = (comp*) tempoz;

	double *wi, *wc, *wia, *wca, *wib, *wcb;
	comp *wki, *wkc;
	cppmallocl(wi, 3*nxyz, double);
	cppmallocl(wc, 3*nxyz, double);
	cppmallocl(wki, 3*nxyz, comp);
	cppmallocl(wkc, 3*nxyz, comp);
	cppmallocl(wia, 3*nxyz, double);
	cppmallocl(wca, 3*nxyz, double);
	cppmallocl(wib, 3*nxyz, double);
	cppmallocl(wcb, 3*nxyz, double);

	// Create files for Energies and Spectra and print headers
	double time, Y=M_PI/(md.dx*dk), runtime, g;  
	clock_t IT, ET;
	//MPI_File es, spec;
	FILE *es, *spec;
	char buf[512];	
	int tint, min;
    
	sprintf(file_name, "%s_hh_%c.txt", md.prefix, u);
	sprintf(file_name2, "%s_hhspec_%c.txt", md.prefix, u);
// 	snprintf(bif, 42, "#   Time*eF \t Ec \t\t Ei \t EcSpec \t EiSpec \n\n");
	if(ip==0)
	{
        printf("\t Creating file: %s\n", file_name);
		es = fopen(file_name, "w");
        fprintf(es, "# %10s %12s %12s %12s %12s\n", "time", "Ec", "Ei", "EcSpec", "EiSpec");
        fclose(es);
        
        printf("\t Creating file: %s\n", file_name2);
		es = fopen(file_name2, "w");
        fclose(es);
	}

	for(inom=0; inom<md.cycles; inom++)												
	{
        if(inom%np != ip) continue; 

		// printf("HERE: ip=%d, np=%d, line=%d, frame=%d\n", ip, np, __LINE__, inom);
      
		Ec=0, Ei=0, Eck=0, Eik=0, Ek=0, U=0, q=0;
		wdata_get_time(&md, inom, &time);

		// Read binary files
		if(ip==0)
// 			if(inom%p==0)
			{
				printf("\t Reading frame [%d], currently at ", inom);
				cout<<setprecision(2)<<double(inom)*100./double(md.cycles)<<"%"<<endl;
			}

		ierr = wdata_read_cycle(&md, "rho_a", inom, rho_a);
		ierr = wdata_read_cycle(&md, "rho_b", inom, rho_b);
		ierr = wdata_read_cycle(&md, "j_a", inom, j_a);
		ierr = wdata_read_cycle(&md, "j_b", inom, j_b);
		// Good so far

		// Choose spin component(s) to investigate
		if(c==0)									// Only for "a" component
		{
			for(int ixyz=0; ixyz<nxyz; ixyz++)
				waxr[ixyz] = j_a_x[ixyz]/sqrt(rho_a[ixyz]+EPS);

			for(int ixyz=0; ixyz<nxyz; ixyz++)
				wayr[ixyz] = j_a_y[ixyz]/sqrt(rho_a[ixyz]+EPS);

			for(int ixyz=0; ixyz<nxyz; ixyz++)
				wazr[ixyz] = j_a_z[ixyz]/sqrt(rho_a[ixyz]+EPS);
		}

		if(c==1)									// Only for "b" component
		{
			for(int ixyz=0; ixyz<nxyz; ixyz++)
				waxr[ixyz] = j_b_x[ixyz]/sqrt(rho_b[ixyz]+EPS);

			for(int ixyz=0; ixyz<nxyz; ixyz++)
				wayr[ixyz] = j_b_y[ixyz]/sqrt(rho_b[ixyz]+EPS);

			for(int ixyz=0; ixyz<nxyz; ixyz++)
				wazr[ixyz] = j_b_z[ixyz]/sqrt(rho_b[ixyz]+EPS);
		}

		if(c==2)									// "a+b" components
		{
			for(int ixyz=0; ixyz<nxyz; ixyz++)
				waxr[ixyz] = (j_a_x[ixyz]+ j_b_x[ixyz]) / (sqrt(rho_a[ixyz] + rho_b[ixyz] + EPS)); // + sqrt(rho_b[ixyz]));

			for(int ixyz=0; ixyz<nxyz; ixyz++)
				wayr[ixyz] = (j_a_y[ixyz]+ j_b_y[ixyz]) / (sqrt(rho_a[ixyz] + rho_b[ixyz] + EPS)); // (sqrt(rho_a[ixyz]) + sqrt(rho_b[ixyz]));

			for(int ixyz=0; ixyz<nxyz; ixyz++)
				wazr[ixyz] = (j_a_z[ixyz]+ j_b_z[ixyz]) / (sqrt(rho_a[ixyz] + rho_b[ixyz] + EPS)); // (sqrt(rho_a[ixyz]) + sqrt(rho_b[ixyz]));
		}

		// Initialize w(r) in coordinate space
		for(int ixyz=0; ixyz<nxyz; ixyz++) 							// Check w0
		{
			lxr[ixyz] = waxr[ixyz];
			lxi[ixyz] = 0.;
			lx[ixyz] = comp(lxr[ixyz], lxi[ixyz]);

			lyr[ixyz] = wayr[ixyz];
			lyi[ixyz] = 0.;
			ly[ixyz] = comp(lyr[ixyz], lyi[ixyz]);

			lzr[ixyz] = wazr[ixyz];
			lzi[ixyz] = 0.;
			lz[ixyz] = comp(lzr[ixyz], lzi[ixyz]);
		}

		// Obtain w(k) in momentum space
		fftw_plan Q1x = fftw_plan_dft_3d(md.nx, md.ny, md.nz, wax, temppx, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(Q1x);
		fftw_destroy_plan(Q1x);
		fftw_plan Q1y = fftw_plan_dft_3d(md.nx, md.ny, md.nz, way, temppy, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(Q1y);
		fftw_destroy_plan(Q1y);
		fftw_plan Q1z = fftw_plan_dft_3d(md.nx, md.ny, md.nz, waz, temppz, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(Q1z);
		fftw_destroy_plan(Q1z);

		// Make f(k) = (k*tA(k)/k^2)*k = w_c(k)
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			if(kp[ixyz]==0)
			{
				fx[ixyz] = 0.;
				fy[ixyz] = 0.;
				fz[ixyz] = 0.;
			}
			else
			{
				fx[ixyz] = (kx[ixyz] * px[ixyz] + ky[ixyz] * py[ixyz] + kz[ixyz] * pz[ixyz]) / (abs(kp[ixyz])*abs(kp[ixyz])*md.ny*md.nx*md.nz) * kx[ixyz];		// component fx and fy
				fy[ixyz] = (kx[ixyz] * px[ixyz] + ky[ixyz] * py[ixyz] + kz[ixyz] * pz[ixyz]) / (abs(kp[ixyz])*abs(kp[ixyz])*md.ny*md.nx*md.nz) * ky[ixyz];
				fz[ixyz] = (kx[ixyz] * px[ixyz] + ky[ixyz] * py[ixyz] + kz[ixyz] * pz[ixyz]) / (abs(kp[ixyz])*abs(kp[ixyz])*md.ny*md.nx*md.nz) * kz[ixyz];
			}
		}

		// Copy f(k)
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			gx[ixyz] = fx[ixyz];
			gy[ixyz] = fy[ixyz];
			gz[ixyz] = fz[ixyz];
		}

		// Anti-FT => Return w_c(r) in coordinate space
		fftw_plan L1x = fftw_plan_dft_3d(md.nx, md.ny, md.nz, tempfx, tempcx, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(L1x);
		fftw_destroy_plan(L1x);
		fftw_plan L1y = fftw_plan_dft_3d(md.nx, md.ny, md.nz, tempfy, tempcy, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(L1y);
		fftw_destroy_plan(L1y);
		fftw_plan L1z = fftw_plan_dft_3d(md.nx, md.ny, md.nz, tempfz, tempcz, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(L1z);
		fftw_destroy_plan(L1z);

		// w(k) - f(k) = w_i(k)
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			hx[ixyz] = (px[ixyz]/double(md.nx*md.ny*md.nz) - gx[ixyz]);
			hy[ixyz] = (py[ixyz]/double(md.nx*md.ny*md.nz) - gy[ixyz]);
			hz[ixyz] = (pz[ixyz]/double(md.nx*md.ny*md.nz) - gz[ixyz]);
		}

		// Anti-FT => Return w_i(r) in coordinate space
		fftw_plan L2x = fftw_plan_dft_3d(md.nx, md.ny, md.nz, tempix, tempvx, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(L2x);
		fftw_destroy_plan(L2x);
		fftw_plan L2y = fftw_plan_dft_3d(md.nx, md.ny, md.nz, tempiy, tempvy, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(L2y);
		fftw_destroy_plan(L2y);
		fftw_plan L2z = fftw_plan_dft_3d(md.nx, md.ny, md.nz, tempiz, tempvz, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(L2z);
		fftw_destroy_plan(L2z);

		// Real part of w_i & w_c components
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			wixr[ixyz] = real(wix[ixyz]);		wcxr[ixyz] = real(wcx[ixyz]);
			wiyr[ixyz] = real(wiy[ixyz]);		wcyr[ixyz] = real(wcy[ixyz]);
			wizr[ixyz] = real(wiz[ixyz]);		wczr[ixyz] = real(wcz[ixyz]);
		}

		// Integrals to get energies (coordinate space) => Ei & Ec
		// Integrals to get energies (Fourier space) => Eik & Eck
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			// Ek = Ek + (abs(px[ixyz])*abs(px[ixyz]) + abs(py[ixyz])*abs(py[ixyz]) + abs(pz[ixyz])*abs(pz[ixyz])) /(2.*md.nx*md.ny*md.nz);
			Ek = Ek + (abs(waxr[ixyz])*abs(waxr[ixyz]) + abs(wayr[ixyz])*abs(wayr[ixyz]) + abs(wazr[ixyz])*abs(wazr[ixyz])) /2.;

			Eik = Eik + (abs(hx[ixyz])*abs(hx[ixyz]) + abs(hy[ixyz])*abs(hy[ixyz]) + abs(hz[ixyz])*abs(hz[ixyz]))*md.nx*md.ny*md.nz /2.;
			Ei = Ei + (abs(wix[ixyz])*abs(wix[ixyz]) + abs(wiy[ixyz])*abs(wiy[ixyz]) + abs(wiz[ixyz])*abs(wiz[ixyz])) /2.;

			Eck = Eck + (abs(gx[ixyz])*abs(gx[ixyz]) + abs(gy[ixyz])*abs(gy[ixyz]) + abs(gz[ixyz])*abs(gz[ixyz]))*md.nx*md.ny*md.nz /2.;
			Ec = Ec + (abs(wcx[ixyz])*abs(wcx[ixyz]) + abs(wcy[ixyz])*abs(wcy[ixyz]) + abs(wcz[ixyz])*abs(wcz[ixyz])) /2.;
		}

		// Vectors for irrotational and compressive components
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			wi[ixyz+0*nxyz] = wixr[ixyz];
            	wi[ixyz+1*nxyz] = wiyr[ixyz];
            	wi[ixyz+2*nxyz] = wizr[ixyz];

			wc[ixyz+0*nxyz] = wcxr[ixyz];
            	wc[ixyz+1*nxyz] = wcyr[ixyz];
            	wc[ixyz+2*nxyz] = wczr[ixyz];
		}

		// Copy w_i & w_c in Fourier space for spectral analysis
		for(int ixyz=0; ixyz<nxyz; ixyz++)
		{
			wki[ixyz+0*nxyz] = hx[ixyz];
            	wki[ixyz+1*nxyz] = hy[ixyz];
            	wki[ixyz+2*nxyz] = hz[ixyz];

			wkc[ixyz+0*nxyz] = gx[ixyz];
            	wkc[ixyz+1*nxyz] = gy[ixyz];
            	wkc[ixyz+2*nxyz] = gz[ixyz];
		}		// Good so far

		// Set up MPI ordering

		// Spectral analysis => eiv & ecv
		k=0, eiv=0, ecv=0, EiSpec=0, EcSpec=0;
		for(int l=1; l<=M_PI/(md.dx*dk); l++)	// From 0 to pi/dx for each k
		{
			k = dk*l;	eiv=0.;	ecv=0.; 	OmegaVolume=0.;
			for(int ixyz=0; ixyz<nxyz; ixyz++)
			{
				double k1[3] = {kx[ixyz], ky[ixyz], kz[ixyz]};
				eiv += (abs(hx[ixyz])*abs(hx[ixyz]) + abs(hy[ixyz])*abs(hy[ixyz]) + abs(hz[ixyz])*abs(hz[ixyz])) * Omega(k1, k, dk) * (1./2.) *md.nx*md.ny*md.nz;		// Spectra
				ecv += (abs(gx[ixyz])*abs(gx[ixyz]) + abs(gy[ixyz])*abs(gy[ixyz]) + abs(gz[ixyz])*abs(gz[ixyz])) * Omega(k1, k, dk) * (1./2.) *md.nx*md.ny*md.nz;
				// Omega_volume stands for integral over angles
				OmegaVolume += Omega(k1, k, dk);
			}

			// Rescale integral to match 4pi
			eiv *= 4.*M_PI/OmegaVolume;
			ecv *= 4.*M_PI/OmegaVolume;
			// Linearization of the integral
			eiv *= k*k;
			ecv *= k*k;
			// Additional missing coefficient
			eiv /= 2.*M_PI/(md.nx*md.dx) * 2.*M_PI/(md.ny*md.dy) * 2.*M_PI/(md.nz*md.dz);
			ecv /= 2.*M_PI/(md.nx*md.dx) * 2.*M_PI/(md.ny*md.dy) * 2.*M_PI/(md.nz*md.dz);
			// Integration over dk => Energies from spectra EiSpec & EcSpec
			EiSpec += eiv*2.*M_PI/md.nx;											// Energies from spectra
			EcSpec += ecv*2.*M_PI/md.nx;

			// Print spectra data to .txt
			ecvprint[l-1] = ecv;
			eivprint[l-1] = eiv;
		}
		// printf("HERE: ip=%d, np=%d, line=%d, frame=%d\n", ip, np, __LINE__,inom);

		// Print energies data to .txt
		int mpiflag=1, j;
		MPI_Status MPIStat;
		if(ip!=0) 	// Wait until ip-1 process finishes his job
		{
			// printf("Recv: ip=%d, np=%d, line=%d, frame\n", ip, np, __LINE__, inom);
			j = MPI_Recv(&mpiflag, 1, MPI_INT, ip-1, 99, MPI_COMM_WORLD, &MPIStat);
		}

		FILE *f = fopen(file_name, "a");
		fprintf(f, "%12.3f %12.4f %12.4f %12.4f %12.4f \n", time, Ec, Ei, EcSpec, EiSpec);
		fclose(f);
		FILE *g = fopen(file_name2, "a");
		fprintf(g, "# t*eF = %12.6f \n", time);
        fprintf(g, "# %10s %16s \t %16s \n", "k", "ecv(k)", "eiv(k)");
		for(int l=1; l<=M_PI/(md.dx*dk); l++)
			fprintf(g, "%12.4f %16.10g \t %16.10g \n", l*dk, ecvprint[l-1], eivprint[l-1]);
		fprintf(g, "\n\n");
		fclose(g);

		int ierr2, ierr3; 		// Something weird happens here, didn't find a fix other than removing the security check
		ierr2 = wdata_write_cycle(&md, sc, wc);
		if(ierr2!=0) { printf("ierr=%d, ERROR: Cannot add w_c!\n", ierr2); return 1;}
		ierr3 = wdata_write_cycle(&md, si, wi);
		if(ierr3!=0) { printf("ierr=%d, ERROR: Cannot add w_i!\n", ierr3); return 1;}

		if(inom==md.cycles-1)
            	break; 	//MPI_Abort(MPI_COMM_WORLD,0);

		if(ip!=np-1)
		{
			// printf("Send: ip=%d, np=%d, line=%d, frame\n", ip, np, __LINE__, inom);
			j = MPI_Send(&mpiflag, 1, MPI_INT, ip+1, 99, MPI_COMM_WORLD);
		}

		//if(inom>50)	break;
		// printf("HERE: ip=%d, np=%d, line=%d, frame\n", ip, np, __LINE__, inom);
	}

	double t1 = e_t();
	if(ip==0)	
	{
		tint = (int)t1;
		min = tint/60;
		printf("\n \tCompleted in %4.0fs (%2dm%2ds)! \n", t1, min, tint-min*60);
	}

	MPI_Finalize();	
	return 0;
}

double Omega(double *k1, double k, double dk)
{
	double L = sqrt(k1[0]*k1[0] + k1[1]*k1[1] + k1[2]*k1[2]);

    	if(L>=k-dk/2. &&  L<k+dk/2.) 
		return 1.0;
   	else return 0.0;
}

void b_t() 
{ 
  	gettimeofday( &s , NULL ) ;
  	t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec ;
}

double e_t() 
{
    	gettimeofday( &s , NULL ) ;
    	return s.tv_sec + 1e-6 * s.tv_usec - t_gettimeofday ;
}
