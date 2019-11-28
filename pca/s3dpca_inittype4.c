#include <stdlib.h>
#include <stddef.h>
#include <complex.h>
#include <fftw3.h>
#include "pca_settings.h"
#include "pca_macro.h"
#include "s3dpca_me.h"
#include "s3dpca_inittype4.h"

// extern functions from pca.io file
int read_measurement_file_header(const char * file_name, 
                                        int *nx, int *ny, int *nz, 
                                        double *dx, double *dy, double *dz,
                                        double *eF, double *t0, double *dt,
                                        int *number_of_entries
                                       );
int read_measurement_entry(const char * file_name, int entry_number, void * array, size_t size); 

/**
 * Function returns interpolated value of function for coordinate (x,y,z).
 * Fourier transform coefficients are reqired
 * */
double complex interpolate(double x, double y, double z, double complex *fcoeffs, double *kkx, double *kky, double *kkz, int nx, int ny, int nz)
{
    int ix, iy, iz, ixyz;
    double complex r = 0.0 + I*0.0;
    ixyz=0;
    for(ix=0; ix<nx; ix++) for(iy=0; iy<ny; iy++) for(iz=0; iz<nz; iz++) 
    {
        r += fcoeffs[ixyz]*cexp(I*( kkx[ix]*x + kky[iy]*y + kkz[iz]*z));
        ixyz++;
    }
    
    return r;
}

/**
 * Function initializes data based on dpca files.
 * Function automatically do interpolations if lattice if destination lattice is different that originial one
 * @param eF reference value of eF (OUTPUT)
 * @return 0 -OK, otherwise ERROR
 * */
int inittype4_readdpca(char *inprefix, double *h_densities, double *h_potentials, double *eF)
{
    // Open files input data
    char file_name[256];
    int nx, ny, nz, nxyz;
    double dx, dy, dz;
    double t0, dt;
    int nom;
    int ierr;
    int i, j, ixyz, ix, iy, iz;
    
    // read density file
    sprintf(file_name, "%s_density_a.dpca", inprefix);
    printf("# INITYPE[4]: Reading file `%s`\n", file_name);
    ierr = read_measurement_file_header(file_name,  &nx, &ny, &nz, &dx, &dy, &dz, eF, &t0, &dt, &nom);
    if(ierr!=0) return ierr;
    printf("# INITYPE[4]: Lattice for input data: %d x %d x %d with lattice spacing %.2f x %.2f x %.2f\n", nx, ny, nz, dx, dy, dz);
    printf("# INITYPE[4]: Target lattice        : %d x %d x %d with lattice spacing %.2f x %.2f x %.2f\n", NX, NY, NZ, DX, DY, DZ);
    
    // check is the target lattice is the same as intial
    if(nx!=NX || ny!=NY || nz!=NZ || dx!=DX || dy!=DY || dz!=DZ) return -99;
    nxyz=nx*ny*nz;
    
    // densities - decode - target arrays
    double *rho_a = (double *)(h_densities +  0*NXYZ); // <- from dpca
    double *rho_b = (double *)(h_densities +  1*NXYZ); // <- from dpca
    double *tau_a = (double *)(h_densities +  2*NXYZ); // <- fill using TF formula
    double *tau_b = (double *)(h_densities +  3*NXYZ); // <- fill using TF formula
    double complex *nu = (double complex *)(h_densities +  4*NXYZ);
    double *j_a_x = (double *)(h_densities +  6*NXYZ); // <- from dpca
    double *j_a_y = (double *)(h_densities +  7*NXYZ); // <- from dpca
    double *j_a_z = (double *)(h_densities +  8*NXYZ); // <- from dpca
    double *j_b_x = (double *)(h_densities +  9*NXYZ); // <- from dpca
    double *j_b_y = (double *)(h_densities + 10*NXYZ); // <- from dpca
    double *j_b_z = (double *)(h_densities + 11*NXYZ); // <- from dpca
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NXYZ);
    double *V_b = (double *)(h_potentials +  1*NXYZ);
    double complex *delta = (double complex *)(h_potentials +  2*NXYZ); // <- from dpca
        
    // read from files
    ierr = read_measurement_entry(file_name, nom-1, rho_a, sizeof(double)*nx*ny*nz); 
    if(ierr!=0) return ierr;
    
    sprintf(file_name, "%s_density_b.dpca", inprefix);
    printf("# INITYPE[4]: Reading file `%s`\n", file_name);
    ierr = read_measurement_entry(file_name, nom-1, rho_b, sizeof(double)*nx*ny*nz); 
    if(ierr!=0) return ierr;  
    
    sprintf(file_name, "%s_delta.dpca", inprefix);
    printf("# INITYPE[4]: Reading file `%s`\n", file_name);
    ierr = read_measurement_entry(file_name, nom-1, delta, sizeof(double complex)*nx*ny*nz); 
    if(ierr!=0) return ierr;  
    
    sprintf(file_name, "%s_current_a.dpca", inprefix);
    printf("# INITYPE[4]: Reading file `%s`\n", file_name);
    ierr = read_measurement_entry(file_name, nom-1, j_a_x, sizeof(double)*nx*ny*nz*3); 
    if(ierr!=0) return ierr; 
    
    sprintf(file_name, "%s_current_b.dpca", inprefix);
    printf("# INITYPE[4]: Reading file `%s`\n", file_name);
    ierr = read_measurement_entry(file_name, nom-1, j_b_x, sizeof(double)*nx*ny*nz*3); 
    if(ierr!=0) return ierr; 
    
    // fill tau using TF formula
    double rho, eFt;
    for(ixyz=0; ixyz<nxyz; ixyz++)
    {
        rho=rho_a[ixyz];
        eFt = 0.5 * pow(6.0*M_PI*M_PI*rho, 2./3.);
        tau_a[ixyz] = 0.6 * rho * eFt;
        
        rho=rho_b[ixyz];
        eFt = 0.5 * pow(6.0*M_PI*M_PI*rho, 2./3.);
        tau_b[ixyz] = 0.6 * rho * eFt;
        
        // rest arrays to zeros
        nu[ixyz] = 0.0 + I*0.0;
        V_a[ixyz] = 0.0;
        V_b[ixyz] = 0.0;
    }
    
    return 0;
}

// Function computes particle number in T-F approximation
double npart_TF(int nxyz, double *U, double mu, double dxyz)
{
    double kF, rho=0.0, M=1.0;
    int i;
    for(i=0; i<nxyz; i++) if((mu-U[i])>0.0)
    {
        kF=sqrt(2.*M*(mu-U[i]));
        rho+=kF*kF*kF/(6.*M_PI*M_PI);
    }
    
    return rho*dxyz;
}

// ZBRENT
void reset_zbrent_data(zbrent_data_t *zbrent_data)
{
    int i;
    for(i=0; i<ITMAX+2; i++) zbrent_data->status[i]=0;
}

double get_point_from_zbrent_data(zbrent_data_t *zbrent_data)
{
    int i;
    for(i=0; i<ITMAX+2; i++) if(zbrent_data->status[i]==1) return zbrent_data->mu[i];    
    
    return 0.0;
}

void set_point_in_zbrent_data(zbrent_data_t *zbrent_data, double value)
{
    int i;
    for(i=0; i<ITMAX+2; i++) 
    {
        if(zbrent_data->status[i]==1) 
        {
            zbrent_data->status[i]=2;
            zbrent_data->N[i]=value;
            break;
        }
    }    
}

int zbrent_data_iters(zbrent_data_t *zbrent_data)
{
    int i;
    for(i=0; i<ITMAX+2; i++) if(zbrent_data->status[i]==1) return i;    
    
    return 0;
}

// Machine floating-point precision.
#define EPS 1.0e-16
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
// Using Brent's method, find the root of a function func known to lie between x1 and x2. The
// root, returned as zbrent, will be refined until its accuracy is tol .
// Taken from c recepies
double zbrent2(zbrent_data_t *zbrent_data, double *x1, double *x2, double tol, int *ierr)
{
    int iter;
    double a,b,c,d,e=0.0,min1,min2;
    double fa,fb,fc,p,q,r,s,tol1,xm;
     a=(*x1);
     b=(*x2);
     c=(*x2);
    
    // function for f(a)
    if(zbrent_data->status[0]==2)
    {
        if(fabs(a-zbrent_data->mu[0])>1.0e-8) {*ierr=-1; return 0.0; }
        fa=zbrent_data->N[0];
    }
    else
    {
        zbrent_data->mu[0]=a; // for this point
        zbrent_data->status[0]=1; // request computation
        *ierr=1; 
        return 0.0;
    }
        
    // function for f(b)
    if(zbrent_data->status[1]==2)
    {
        if(fabs(b-zbrent_data->mu[1])>1.0e-8) {*ierr=-2; return 0.0; }
        fb=zbrent_data->N[1];
    }
    else
    {
        zbrent_data->mu[1]=b; // for this point
        zbrent_data->status[1]=1; // request computation
        *ierr=2; 
        return 0.0;
    }
        
    if (fa > 0.0 && fb > 0.0) // lower window
    {
        zbrent_data->mu[1]=zbrent_data->mu[0]; 
        zbrent_data->N[1]=zbrent_data->N[0];
        zbrent_data->status[1]=2; // it is computed
        
        d=b-a; // delta
        b=a;
        a=a-d;
        *x1=a;
        *x2=b;
        
        zbrent_data->mu[0]=a; 
        zbrent_data->N[0]=-1.0;
        zbrent_data->status[0]=1; // request computation  
        
        *ierr=1; 
        return 0.0;
    }
    
    if (fa < 0.0 && fb < 0.0) // upper window
    {
        zbrent_data->mu[0]=zbrent_data->mu[1]; 
        zbrent_data->N[0]=zbrent_data->N[1];
        zbrent_data->status[0]=2; // it is computed
        
        d=b-a; // delta
        a=b;
        b=b+d;
        *x1=a;
        *x2=b;
        
        zbrent_data->mu[1]=b; 
        zbrent_data->N[1]=-1.0;
        zbrent_data->status[1]=1; // request computation  
        
        *ierr=2; 
        return 0.0;
    }
    
    fc=fb;
    for (iter=1;iter<=ITMAX;iter++) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c=a; //Rename a, b, c and adjust bounding interval d.
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        tol1=2.0*EPS*fabs(b)+0.5*tol; //Convergence check.
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0) { *ierr=0; return b; }
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            s=fb/fa; //Attempt inverse quadratic interpolation.
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            if (p > 0.0) q = -q; // Check whether in bounds.
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2)) {
                e=d; //Accept interpolation.
                d=p/q;
            } else {
                d=xm; //Interpolation failed, use bisection.
                e=d;
            }
        } else { //Bounds decreasing too slowly, use bisection.
            d=xm;
            e=d;
        }
        a=b; //Move last best guess to a.
        fa=fb;
        if (fabs(d) > tol1) //Evaluate new trial root.
            b += d;
        else
            b += SIGN(tol1,xm);
        
        // function for f(b)
        if(zbrent_data->status[iter+1]==2)
        {
            if(fabs(b-zbrent_data->mu[iter+1])>1.0e-8) {*ierr=-1*(iter+1); return b; }
            fb=zbrent_data->N[iter+1];
        }
        else
        {
            zbrent_data->mu[iter+1]=b; // for this point
            zbrent_data->status[iter+1]=1; // request computation
            *ierr=iter+2; 
            return b;
        }
    }
    *ierr=-998;
    return b; //Never get here.
} 

