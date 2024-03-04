// Author: Gabriel Wlazlowski
// This is simple code for vortex position extraction
//
// NOTE: the code requires GSL - GNU Scientific Library
//
// To compile:
//    g++ simple_vortex_track_2d.cpp -o simple_vortex_track_2d -I$WSLDA/hpc-engine -I$WSLDA/lib/wdata/c -L$WSLDA/lib/wdata -lwdata -O3 -std=c++98 -lfftw3 -lgsl -lgslcblas -lm


// standard libraries
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "wdata.h"

#include "gsl/gsl_multimin.h"

#define cppmallocl(pointer,size,type)                                           \
    if ( ( pointer = (type *) malloc( (size) * sizeof( type ) ) ) == NULL )     \
    {                                                                           \
        fprintf( stderr , "ERROR: cannot malloc()! Exiting!\n") ;               \
        fprintf( stderr , "ERROR: file=`%s`, line=%d\n", __FILE__, __LINE__ ) ; \
        return -1 ;                                                             \
    }

#define file_operationl( cmd )                                                  \
    { ierr=cmd;                                                                 \
    if(ierr)                                                                    \
    {                                                                           \
        fprintf( stderr , "FILE ERROR:: cannot execute: %s\n" , #cmd);          \
        fprintf( stderr , "file=`%s`, line=%d\n" ,__FILE__,__LINE__) ;          \
        fprintf( stderr , "Error=%d\nExiting!\n" ,ierr) ;                       \
        return( EXIT_FAILURE ) ;                                                \
    } }

// ------------------------
class Vector3D
{
public:

	static const Vector3D X_UNIT;// = Vector3D(1.0f, 0.0f, 0.0f);
	static const Vector3D Y_UNIT;// = Vector3D(0.0f, 1.0f, 0.0f);
	static const Vector3D Z_UNIT;// = Vector3D(0.0f, 0.0f, 1.0f);
	static const Vector3D ZERO;// = Vector3D(0.0f, 0.0f, 0.0f);

	static const float PI = 3.14159265358979323846264338327950288f;


    float X;
    float Y;
    float Z;

    float degToRad(float deg) {
            return deg * 3.14159265358979323846264338327950288f / 180;
    }


    Vector3D() {
        this->X = 0;
        this->Y = 0;
        this->Z = 0;
    }

    Vector3D(float x, float y, float z) {
        this->X = x;
        this->Y = y;
        this->Z = z;
    }

    Vector3D(Vector3D* other) {
        this->X = other->X;
        this->Y = other->Y;
        this->Z = other->Z;
    }

    float getX(){return this->X;}
    float getY(){return this->Y;}
    float getZ(){return this->Z;}

    void setX(float x){this->X = x;}
    void setY(float y){this->Y = y;}
    void setZ(float z){this->Z = z;}

    float getAngleRad(Vector3D from) {
    	float len = this->length()*from.length();
    	if(len == 0)
    		return 0;
    	return acosf(this->dot(from)/(this->length()*from.length()));
    }

    float getAngleXRad(Vector3D from) {
    	return atan2f(this->Z, this->Y) - atan2f(from.Z, from.Y);
    }

    float getAngleYRad(Vector3D from) {
    	return atan2f(this->X, this->Z) - atan2f(from.X, from.Z);
    }

    float getAngleZRad(Vector3D from) {
    	return atan2f(this->Y, this->X) - atan2f(from.Y, from.X);
    }

    /////////////////////////////////////////////// +++ //////////////////////////////////////////////////////////
    Vector3D operator++() {
        this->X++; this->Y++; this->Z++;
        return *this;
    }

    void operator+=(Vector3D v) {
        this->X += v.X; this->Y += v.Y; this->Z += v.Z;
    }

    Vector3D operator+(Vector3D v) {
        Vector3D v2;
        v2.X = this->X + v.X; v2.Y = this->Y + v.Y; v2.Z = this->Z + v.Z;
        return v2;
    }

    /////////////////////////////////////////////// --- //////////////////////////////////////////////////////////
    Vector3D operator--() {
        this->X--; this->Y--; this->Z--;
        return *this;
    }

    void operator-=(Vector3D v) {
        this->X -= v.X; this->Y -= v.Y; this->Z -= v.Z;
    }

    Vector3D operator-(Vector3D v) {
        Vector3D v2;
        v2.X = this->X - v.X; v2.Y = this->Y - v.Y; v2.Z = this->Z - v.Z;
        return v2;
    }

    /////////////////////////////////////////////// *** //////////////////////////////////////////////////////////

    bool operator*=(float scalar) {
        this->X *= scalar; this->Y *= scalar; this->Z *= scalar;
    }

    Vector3D operator*(float scalar){
        Vector3D v2;
        v2.X = this->X * scalar; v2.Y = this->Y * scalar; v2.Z = this->Z * scalar;
        return v2;
    }

    bool operator*=(Vector3D v) {
        this->X *= v.X; this->Y *= v.Y; this->Z *= v.Z;
    }

    Vector3D operator*(Vector3D v){
        Vector3D v2;
        v2.X = this->X * v.X; v2.Y = this->Y * v.Y; v2.Z = this->Z * v.Z;
        return v2;
    }

    /////////////////////////////////////////////// /// //////////////////////////////////////////////////////////

     bool operator/=(float scalar) {
        this->X /= scalar; this->Y /= scalar; this->Z /= scalar;
    }

     Vector3D operator/(float scalar){
         Vector3D v2;
        v2.X = this->X / scalar; v2.Y = this->Y / scalar; v2.Z = this->Z / scalar;
        return v2;
    }

    bool operator/=( Vector3D v) {
        this->X /= v.X; this->Y /= v.Y; this->Z /= v.Z;
    }

     Vector3D operator/( Vector3D v){
        Vector3D v2;
        v2.X = this->X / v.X; v2.Y = this->Y / v.Y; v2.Z = this->Z / v.Z;
        return v2;
    }

    /////////////////////////////////////////////// ^^^ //////////////////////////////////////////////////////////

    bool operator^=(float power) {
        this->X = pow(this->X, power); this->Y = pow(this->Y, power); this->Z = pow(this->Z, power);
    }

    Vector3D operator^(float power){
        Vector3D v2;
        v2.X = pow(this->X, power); v2.Y = pow(this->Y, power); v2.Z = pow(this->Z, power);
        return v2;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    float length(){
        return sqrt(this->X * this->X + this->Y * this->Y +  this->Z * this->Z);
    }

    float sqLength(){
    	return this->X * this->X + this->Y * this->Y +  this->Z * this->Z;
    }

    Vector3D unit() {
        return *(this) / this->length();
    }
    
    void normalize(){
        float l=this->length();
        this->X/=l;
        this->Y/=l;
        this->Z/=l;
    }

    Vector3D cross(Vector3D other) {
    	Vector3D* a = this;
    	Vector3D* b = &other;
        return Vector3D(
                a->Y*b->Z - a->Z*b->Y,
                a->Z*b->X - a->X*b->Z,
                a->X*b->Y - a->Y*b->X
        		);
    }

    float dot(Vector3D other) {
    	return this->X*other.X + this->Y*other.Y + this->Z*other.Z;
    }

    void rotateXrad(float theta)
    {
        Vector3D newv;
        newv.X = this->X;
        newv.Y = this->Y * cos(theta) + this->Z * -sin(theta);
        newv.Z = this->Y * sin(theta) + this->Z * cos(theta);

        (*this) = newv;
    }

    void rotateXdeg(float theta)
    {
        this->rotateXrad(this->degToRad(theta));
    }


    void rotateYrad(float theta)
    {
        Vector3D newv;
        newv.X = this->X * cos(theta) + this->Z * sin(theta);
        newv.Y = this->Y;
        newv.Z = this->Z * -sin(theta) + this->Z * cos(theta);

        (*this) = newv;
    }
    void rotateYdeg(float theta)
    {
        this->rotateYrad(this->degToRad(theta));
    }

    void rotateZrad(float theta)
    {
        Vector3D newv;
        newv.X = this->X * cos(theta) + this->Y * -sin(theta);
        newv.Y = this->X * sin(theta) + this->Y * cos(theta);
        newv.Z = this->Z;

        (*this) = newv;
    }
    void rotateZdeg(float theta)
    {
        this->rotateZrad(this->degToRad(theta));
    }


    void rotateAroundVectorRad(Vector3D vect2, float theta)
    {
        Vector3D newv;
        Vector3D unit = vect2.unit();
        //theta = Math.toRadians(theta);
        float q0 = cos(theta/2);
        float q1 = sin(theta/2)*unit.X;
        float q2 = sin(theta/2)*unit.Y;
        float q3 = sin(theta/2)*unit.Z;

        // column vect
        newv.X = (q0*q0 + q1*q1 - q2*q2 - q3*q3)*this->X +               2*(q2*q1 - q0*q3) * this->Y +               2*(q3*q1 + q0*q2) * this->Z;
        newv.Y =               2*(q1*q2 + q0*q3)*this->X + (q0*q0 - q1*q1 + q2*q2 - q3*q3) * this->Y +               2*(q3*q2 - q0*q1) * this->Z;
        newv.Z =               2*(q1*q3 - q0*q2)*this->X +               2*(q2*q3 + q0*q1) * this->Y + (q0*q0 - q1*q1 - q2*q2 + q3*q3) * this->Z;

        (*this) = newv;
    }
    void rotateAroundVectorDeg(Vector3D vect2, float theta)
    {
        this->rotateAroundVectorRad(vect2, this->degToRad(theta));
    }


        void rotateGLrad(float i,float j,float k)
    {

        Vector3D vy(0,1,0);
        Vector3D vz(0,0,1);

        this->rotateXrad(i);
        vy.rotateXrad(i);
        vz.rotateXrad(i);

        this->rotateAroundVectorRad(vy,j);
        vz.rotateAroundVectorRad(vy,j);

        this->rotateAroundVectorRad(vz,k);
    }
    void rotateGLdeg(float i,float j, float k)
    {
        /*
            glRotatef(i,1,0,0);
            glRotatef(j,0,1,0);
            glRotatef(k,0,0,1);
        */
        this->rotateGLrad(this->degToRad(i),this->degToRad(j),this->degToRad(k));
    }


    std::string toString() {
        std::stringstream strs;
        strs << "[" << this->X << "," << this->Y << "," << this->Z << "]";
        std::string str = strs.str();
        return str;
    }

    Vector3D copy() {
    	return Vector3D(this);
    }

    static Vector3D* trilaterate(
    		Vector3D P1,Vector3D P2, Vector3D P3,
    		float r1, float r2, float r3){
    	Vector3D* ret = new Vector3D[2];

    	Vector3D temp1 = P2-P1;
    	Vector3D e_x = temp1.unit();
    	Vector3D temp2 = P3-P1;
    	float i = e_x.dot(temp2);
        Vector3D temp3 = temp2 - e_x*i;
        Vector3D e_y = temp3.unit();
        Vector3D e_z = e_x.cross(e_y);
        float d = (P2-P1).length();
        float j = e_y.dot(temp2);
        float x = (r1*r1 - r2*r2 + d*d) / (2.0*d);
        float y = (r1*r1 - r3*r3 -2*i*x + i*i + j*j) / (2.0*j);
        float temp4 = r1*r1 - x*x - y*y;
        if(temp4<0)
            throw 42;
        float z = sqrt(temp4);

        ret[0] = P1 + e_x*x + e_y*y + e_z*z;
        ret[1] = P1 + e_x*x + e_y*y - e_z*z;

    	return ret;
    }
};

/**
 * Funtion returns inerpolated value of function for coordinate (x,y,z).
 * Fourier transform coefficients are reqired
 * */
double complex interpolate(double complex *fcoeffs, double x, double y, int nx, int ny)
{
    int ix, iy, ixyz;
    double kx, ky;
    double complex r = 0.0 + I*0.0;
    
    // multiply by momentum
    ixyz=0;
    for(ix=0; ix<nx; ix++)
    {
        for(iy=0; iy<ny; iy++)
        {
            // extract momentum
            if(ix<nx/2) kx=2.*M_PI/(( double )nx) * ( double )(ix   );
            else        kx=2.*M_PI/(( double )nx) * ( double )(ix-nx);

            
            if(iy<ny/2) ky=2.*M_PI/(( double )ny) * ( double )(iy   );
            else        ky=2.*M_PI/(( double )ny) * ( double )(iy-ny); 
            
            r += fcoeffs[ixyz]*cexp(I*( kx*x + ky*y ));
            
            ixyz++;
        }
    }
    
    return r;
}

/** 
 * Function computes phase difference between two argumenents
 * To be consistent with computation of velocity field by formula v=j/rho and v=grad Phi
 * return value is from range (-pi,+pi)
 * */
double pdiff(double arg1, double arg2)
{
    
    double p;
    if(arg1>=arg2) p=arg1-arg2;
    else           p=2.0*M_PI - (arg2-arg1);
    
    if(p>M_PI) p-=2.*M_PI;
    return p;
    
}

/**
 * @return 1 - pahase decreas monotnicaly 
 * */
int check_if_monotonic_down(const int size, double * args)
{
    int i;
    double argdiffs[size];
    for(i=0; i<size; i++) argdiffs[i]=pdiff(args[i],args[(i+1)%size]);
//     for(i=0; i<size; i++) printf("check_if_monotonic_down: %d %f\n", i, argdiffs[i]);
    for(i=0; i<size; i++) if(argdiffs[i]<0.0) return 0;
        
    return 1; 
}


// ----------------------- GSL MINIMIZATION ------------------------
int _gslnx;
int _gslny;

double my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double complex *p = (double complex *)params;

  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);

  double complex r = interpolate(p, x, y, _gslnx, _gslny);
  #define cnorm(a) (creal(a)*creal(a) + cimag(a)*cimag(a))
  return cnorm(r);
}


// ======================================================================================
// ============================= MAIN PART ==============================================
// ======================================================================================
int main( int argc , char* argv[] )
{
    if(argc!=5)
    {
        printf("Usage: %s wtxt x y tsign\n", argv[0]);
        printf("\twtxt: name of wtxt file\n");
        printf("\t(x,y): initial coordinate of vortex to be tracked\n");
        printf("\tsign: vortex sign, +1 or -1\n");
        return 1;
    }
    
    char file_name[256];
    int nx, ny, nz, nxyz;
    int ix, iy, iz, ixyz;
    int i, ierr, status;
    int inom, nom;
    double dx, dy, dz, dt, t0, eF;
    
    // expected position of the vortex
    Vector3D rg=Vector3D(atof(argv[2]),atof(argv[3]),0.0);
    const int Nphi = 8; // number of points used for circle discretization
    const float Rcircle_max=5.0; // Radius of circle used for checking if delta rotates
    const float Rcircle_start=2.5;
    const float Rcircle_min=0.05;
    const float dRcircle=0.05;
    double vortex_sign =atof(argv[4]);
    
    
    float phis[Nphi];
    float dphi = 2.*M_PI / Nphi;
    for(i=0; i<Nphi; i++) phis[i]=dphi*i;
    double argvalues[Nphi];
    
    float Rcircle = Rcircle_start;
    
    // create metadata handler
    wdata_metadata md;
    
    // read metadata from file
    printf("# VORTEX-TRACKER: Reading file `%s`\n", argv[1]);
    ierr = wdata_parse_metadata_file(argv[1], &md);
    if(ierr!=0) {printf("Cannot read metadata file!\n"); return 1;}
    
    nx = md.nx; ny = md.ny; nz = md.nz;
    dx = md.dx; dy = md.dy; dz = md.dz;
    eF = wdata_getconst_value(&md, "eF");
    t0 = md.t0;
    dt = md.dt;
    nom = md.cycles;
    
    printf("# VORTEX-TRACKER: Lattice for input data: %d x %d x %d with lattice spacing %.2f x %.2f x %.2f\n", nx, ny, nz, dx, dy, dz);
    nz=1; // add hoc
    nxyz = nx*ny*nz;
    
    // fftw
    #define USE_FFTW_PLANNER FFTW_ESTIMATE
    double complex *delta_k;
    cppmallocl(delta_k,nx*ny,double complex);
    fftw_plan plan_f = fftw_plan_dft_2d(nx, ny, delta_k, delta_k, FFTW_FORWARD, USE_FFTW_PLANNER);
    fftw_plan plan_b = fftw_plan_dft_2d(nx, ny, delta_k, delta_k, FFTW_BACKWARD, USE_FFTW_PLANNER);
        
    double complex *rho_a_k, *rho_b_k;
    cppmallocl(rho_a_k,nx*ny,double complex);
    cppmallocl(rho_b_k,nx*ny,double complex);
    
    // order parameter
    double complex *delta;
    cppmallocl(delta, nxyz, double complex);
    
    double *rho_a, *rho_b;
    cppmallocl(rho_a, nxyz, double);
    cppmallocl(rho_b, nxyz, double);
    
    // --------------------------- gsl -------------------------
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    double size;
    
    _gslnx=nx;
    _gslny=ny;
        
    // File header
    printf("# STARTING POINT: (%f,%f)\n", atof(argv[2]),atof(argv[3])); 
    printf("# VORTEX SIGN: %f\n", vortex_sign); 
    printf("# COLUMNS:\n");
    printf("# 1: measurement id\n");
    printf("# 2: time*eF\n");
    printf("# 3: n_a(x,y)\n");
    printf("# 4: n_b(x,y)\n");
    printf("# 5: zero delta (x)\n");
    printf("# 6: zero delta (y)\n");
//     printf("# 11: zero delta (r)\n");
//     printf("# 12: zero delta (theta)\n");
    for(inom=0; inom<nom; inom++)
    {
//         printf("# inom=%d...\n", inom);
        wdata_read_cycle(&md, "delta", inom, delta);
        wdata_read_cycle(&md, "rho_a", inom, rho_a);
        wdata_read_cycle(&md, "rho_b", inom, rho_b);
        
        double time=t0+dt*inom;

        // prepare data for interpolations
        for(ixyz=0; ixyz<nxyz; ixyz++) delta_k[ixyz]=rho_a[ixyz];
        fftw_execute(plan_f);  
        for(ixyz=0; ixyz<nxyz; ixyz++) rho_a_k[ixyz]=delta_k[ixyz]/nxyz; 
        
        for(ixyz=0; ixyz<nxyz; ixyz++) delta_k[ixyz]=rho_b[ixyz];
        fftw_execute(plan_f);  
        for(ixyz=0; ixyz<nxyz; ixyz++) rho_b_k[ixyz]=delta_k[ixyz]/nxyz; 
        

        for(ixyz=0; ixyz<nxyz; ixyz++) delta_k[ixyz]=delta[ixyz];
        fftw_execute(plan_f);  
        for(ixyz=0; ixyz<nxyz; ixyz++) delta_k[ixyz]=delta_k[ixyz]/nxyz;     
        

//         printf("# inom2=%d...\n", inom);
        // these parts are taken from vdetect code
        
        // Find basis plane perpendicular to psv (up to arbitrary rotation)
        Vector3D bX = Vector3D(-1.0, 0.0, 0.0);
        Vector3D bY = Vector3D(0.0, 1.0, 0.0);
        
        // phases over circle
        double xc=0.0; // initial value
        double yc=0.0; // initial value
        Vector3D rcircle;
        
        // Increase radius of circle to have sigularity inside the circle
        while(1) // Increase until singularity is inside
        {
            // compute arguments
            for(i=0; i<Nphi; i++)
            {
                rcircle = rg + bX*(xc + Rcircle*cos(vortex_sign*phis[i])) + bY*(yc + Rcircle*sin(vortex_sign*phis[i]));
                argvalues[i] = carg( interpolate(delta_k, rcircle.getX(), rcircle.getY(), nx, ny) );
            }
                            
            // check if monotonic decreas
            status=check_if_monotonic_down(Nphi, argvalues);
            if(status==1 || Rcircle>Rcircle_max) break;
            Rcircle+=dRcircle;
        } /*while(1) // Increase until singularity is inside*/
                
//         printf("Rcircle=%f\n", Rcircle);
        
        
        // If cannot find singularity in vicinity of the point guess point
        if(Rcircle>Rcircle_max)
        {
            printf("# ===> INTIAL POINT rg=%s IS OF BAD QUALITY!\n", rg.toString().c_str());
            return 1; 
        }
        
        // decreas value of Radius
        int RcircleIter=0;
        while(1) // until convergence
        {
            
            
            while(1)
            {
                // decreas until singularity is inside
                // compute arguments
                for(i=0; i<Nphi; i++)
                {
                    rcircle = rg + bX*(xc + Rcircle*cos(vortex_sign*phis[i])) + bY*(yc + Rcircle*sin(vortex_sign*phis[i]));
                    argvalues[i] = carg( interpolate(delta_k, rcircle.getX(), rcircle.getY(), nx, ny) );
                }
                
//                 for(i=0; i<Nphi; i++) 
//                 {
//                     rcircle = rg + bX*(xc + Rcircle*cos(phis[i])) + bY*(yc + Rcircle*sin(phis[i]));
//                     printf("%3d %24s %12.6f %12.6f %12.6f\n", i, rcircle.toString().c_str(), Rcircle, argvalues[i], pdiff(argvalues[i],argvalues[(i+1)%Nphi]));
//                 }
                            
                // check if monotonic decreas
                status=check_if_monotonic_down(Nphi, argvalues);
                if(status==0 || Rcircle<=Rcircle_min) break;
                Rcircle-=dRcircle;
                
                
            }
//             printf("Rcircle=%f\n", Rcircle);
//             return 1;
            
            if(Rcircle<=Rcircle_min) break; // we have solution with given accuracy
            
            // shift circle
            double xcshift=0.0, ycshift=0.0;
            
            for(i=0; i<Nphi; i++) if(pdiff(argvalues[i],argvalues[(i+1)%Nphi])<0.0) 
            {
                xcshift += dRcircle*cos(vortex_sign*phis[i]) + dRcircle*cos(vortex_sign*phis[(i+1)%Nphi]);
                ycshift += dRcircle*sin(vortex_sign*phis[i]) + dRcircle*sin(vortex_sign*phis[(i+1)%Nphi]);
            }
            // normalize vector
            xcshift/=sqrt(xcshift*xcshift + ycshift*ycshift);
            ycshift/=sqrt(xcshift*xcshift + ycshift*ycshift);
            xc+=dRcircle*xcshift;
            yc+=dRcircle*ycshift;
//             printf("xc=%f (xcshift=%f), yc=%f (ycshift=%f)\n", xc, xcshift, yc, ycshift);
            
            RcircleIter++;
            if(RcircleIter>(int)(Rcircle_max/dRcircle)) { printf("# ERROR2...\n"); return 2;};
        }

        Vector3D vortex_core = rg + bX*(xc) + bY*(yc);
        
        // compute density of particles at vortex position
        double na[3], nb[3];
        double _x, _y, _xs, _ys;

        // ------------------ gsl minization ----------------------
        /* Starting point */
        x = gsl_vector_alloc (2);
        gsl_vector_set (x, 0, vortex_core.getX());
        gsl_vector_set (x, 1, vortex_core.getY());

        /* Set initial step sizes to 1 */
        ss = gsl_vector_alloc (2);
        gsl_vector_set_all (ss, Rcircle_min);

        /* Initialize method and iterate */
        minex_func.n = 2;
        minex_func.f = my_f;
        minex_func.params = delta_k;

        s = gsl_multimin_fminimizer_alloc (T, 2);
        gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
        
        iter=0;
        do
        {
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);

            if (status) break;

            size = gsl_multimin_fminimizer_size (s);
            status = gsl_multimin_test_size (size, 1e-3);

        }
        while (status == GSL_CONTINUE && iter < 500);

        _xs = gsl_vector_get (s->x, 0); // VORTEX POSITION
        _ys = gsl_vector_get (s->x, 1); // VORTEX POSITION
        gsl_vector_free(x);
        gsl_vector_free(ss);
        gsl_multimin_fminimizer_free (s);
        
        // no shift
        na[0]=creal(interpolate(rho_a_k, _xs, _ys, nx, ny));
        nb[0]=creal(interpolate(rho_b_k, _xs, _ys, nx, ny));
        
        
        printf("%6d %9.4f %12.9f %12.9f %9.4f %9.4f\n", 
            inom, 
            time*eF,
            na[0], nb[0],
             _xs, _ys
              );
        
        vortex_core.setX(_xs);
        vortex_core.setY(_ys);
        rg = vortex_core; // use it as new value for starting point
                       
    }
    
    printf("# VORTEX-TRACKER: Done.\n");

    return( EXIT_SUCCESS ) ;
}
