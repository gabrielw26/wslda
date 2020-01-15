// Author: Gabriel Wlazlowski

// Code implements external potential - in similar way as in pca_uext

#ifndef __KZPCA_UEXT__
#define __KZPCA_UEXT__


extern double *dc_params; // array with params from input file 

/**
 * Switch function - performs switch in time interval [0-T]
 * */
inline double switch_function(double t, double T, double alpha)
{
    return 0.5*( 1.0+tanh( alpha*tan( M_PI_2*( 2.0*t/T-1.0 ) ) ) );
}

inline double smooth_step(double t, double step_start, double step_stop, double T, double alpha)
{
    if(t<=step_start || t>=step_stop) return 0.0; 
    if(t>=step_start+T && t<=step_stop-T) return 1.0;
    if(t>step_start && t<step_start+T) return switch_function(t-step_start, T, alpha);
    else return 1.0-switch_function(t-step_stop+T, T, alpha);
}

// =================================================================================
// ============================== EXTERNAL POTENTIAL ===============================
// =================================================================================
/**
 * 0 - r1
 * 1 - r2
 * 2 - V0
 * */
double u_ext_tube(int ix, int iy, int it, int spin)
{
    double _x = ((double)(ix) - 1.0*(NX/2))*DX;
    double _y = ((double)(iy) - 1.0*(NY/2))*DY;
    double _r = sqrt(_x*_x + _y*_y);
    
    if(_r<dc_params[0]) 
        return 0.0;
    else if(_r<dc_params[1]) 
        return dc_params[2]*switch_function(_r-dc_params[0], dc_params[1]-dc_params[0], 1.0);
    else 
        return dc_params[2]; 
    
}
void process_u_ext_tube(double *params, double kF)
{
    double eF = kF*kF/2.0;
    params[2] = params[2] * eF;
}
void modify_potentials_phase_imprint2(int it, double *h_densities, double *h_potentials, double *extra_data)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=0;
    double arg, abs_delta;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        double _x = ((double)(ix) - 1.0*(NX/2))*DX;
        double _y = ((double)(iy) - 1.0*(NY/2))*DY;
    
        abs_delta = cabs(delta[ixyz]);
        arg = atan2(_y,_x);
        
        delta[ixyz] = cexp(I*arg) * abs_delta;
        ixyz++;
    }
}

void modify_potentials_imprint_ferron(int it, double *h_densities, double *h_potentials, double *extra_data)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=0;
    double abs_delta;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        double _x = ((double)(ix) - 1.0*(NX/2))*DX;
        double _y = ((double)(iy) - 1.0*(NY/2))*DY;
        double _r = sqrt(_x*_x + _y*_y);
    
        abs_delta = cabs(delta[ixyz]);

        if(_r<dc_params[0])      delta[ixyz] = -1.0 * abs_delta + I*0.0;
        else if(_r>dc_params[1]) delta[ixyz] =  1.0 * abs_delta + I*0.0;
        // else pass
            
        ixyz++;
    }
}

/**
 * 0 - omega_x^2 / 2, in input file value at the boundry in units of eF_a
 * 1 - omega_y^2 / 2, in input file value at the boundry in units of eF_a
 * 2 - swich time, in number of iterations
 * */
double u_ext_HO(int ix, int iy, int it, int spin)
{
    double _x = ((double)(ix) - 1.0*(NX/2))*DX;
    double _y = ((double)(iy) - 1.0*(NY/2))*DY;
    
    double swtch=1.0;
    // if(1.0*it<dc_params[2]) swtch=switch_function(1.0*it, dc_params[2], 1.0);
    
    return swtch*(dc_params[0]*_x*_x + dc_params[1]*_y*_y);
}
void process_params_u_ext_HO(double *params, double kF)
{
    double eF = 0.5*kF*kF;
    params[0] = 4.0*params[0]*eF/(DX*DX*NX*NX);
    params[1] = 4.0*params[1]*eF/(DY*DY*NY*NY);
    
    // params[0] = omega_x^2 / 2;
    double _omega = sqrt(2.0*params[0]);
    // rotation as fraction of the trap frequency
    dc_Omega_a = params[6] * _omega;
    dc_Omega_b = params[6] * _omega;
}

/**
 * Smooth at the boundries HO potential
 * 0 - omega_x^2 / 2, in input file value at the boundry in units of eF_a
 * 1 - omega_y^2 / 2, in input file value at the boundry in units of eF_a
 * 2 - Umaxx - calculated automatically
 * 3 - Umaxy - calculated automatically
 * 4 - swich time, in number of iterations
 * */
#define xs (0.08*(double)(LX))
#define ys (0.08*(double)(LY))
#define HO(x, omega2p2) (omega2p2*x*x)
double u_ext_smooth_HO(int ix, int iy, int it, int spin)
{
    double _x = ((double)(ix) - 1.0*(NX/2))*DX;
    double _y = ((double)(iy) - 1.0*(NY/2))*DY;
        
    double swtch=1.0;
//     if(1.0*it<dc_params[4]) swtch=switch_function(1.0*it, dc_params[4], 1.0);
    
    // ---------------------------- HO ----------------------------
    double s;
    
    double hox;
    if(_x>=(double)(-LX/2)+xs && _x<=(double)(LX/2)-xs) 
    {
        hox = HO(_x, dc_params[0]);
    }
    else if(_x<(double)(-LX/2)+xs) 
    {   
        s = switch_function(_x+(double)(LX/2), xs, 1.0);
        hox = s*HO(_x, dc_params[0]) + (1.0-s)*dc_params[2];
    }
    else
    {
        s = switch_function(_x-((double)(LX/2)-xs), xs, 1.0);
        hox=(1.0-s)*HO(_x, dc_params[0]) + s*dc_params[2];
    }
    
    double hoy;
    if(_y>=(double)(-LY/2)+ys && _y<=(double)(LY/2)-ys) 
    {
        hoy = HO(_y, dc_params[1]);
    }
    else if(_y<(double)(-LY/2)+ys) 
    {   
        s = switch_function(_y+(double)(LY/2), ys, 1.0);
        hoy = s*HO(_y, dc_params[1]) + (1.0-s)*dc_params[3];
    }
    else
    {
        s = switch_function(_y-((double)(LY/2)-ys), ys, 1.0);
        hoy=(1.0-s)*HO(_y, dc_params[1]) + s*dc_params[3];
    }
        
    return swtch*(hox + hoy);
    
//     double swtch2=0.0;
//     swtch2 = smooth_step((double)it, dc_params[10], dc_params[11], dc_params[12], 1.0);
//     
//     return swtch*(hox + hoy + swtch2*dc_params[5]*_y*_y*_y + swtch2*dc_params[5]*_x*_x*_x);

}
void process_params_u_ext_smooth_HO(double *params, double kF)
{
    double eF = 0.5*kF*kF;
    params[0] = 4.0*params[0]*eF/(DX*DX*NX*NX);
    params[1] = 4.0*params[1]*eF/(DY*DY*NY*NY);
    
    // params[0] = omega_x^2 / 2;
    double _omega = sqrt(2.0*params[0]);
    // rotation as fraction of the trap frequency
    dc_Omega_a = params[6] * _omega;
    dc_Omega_b = params[6] * _omega;
 
//     Umax=HO(xs) + 1.0*(HO(0)-HO(xs))
    params[2] = HO((double)(-LX/2)+xs, params[0]) + 1.0*(HO((double)(-LX/2), params[0])-HO((double)(-LX/2)+xs, params[0]));
    params[3] = HO((double)(-LY/2)+ys, params[1]) + 1.0*(HO((double)(-LY/2), params[1])-HO((double)(-LY/2)+ys, params[1]));
    
    
// // //     kF = 1.417086; // use fixed value
// //     kF = 1.0000;
// //     
// //     double eF = 0.5*kF*kF;
// //     params[0] = 4.0*params[0]*eF/(NX*NX); // m*omega_x^2 / 2
// // //     params[1] = 4.0*params[1]*eF/(NY*NY);
// //     double asym=0.015;
// //     params[1] = params[0]*pow((1.-asym)/(1.+asym),2); // m*omega_y^2 / 2
// //     
// // 
// //     
// // // //     // rotating frame    
// //     dc_Omega_a=0.0;
// //     dc_Omega_b=0.0;
// // //     dc_Omega_a = params[6]*sqrt(2.*params[0]); // take this fraction of angular momentum
// // //     dc_Omega_b = params[6]*sqrt(2.*params[0]); // take this fraction of angular momentum
// // //     
// // //     // See supplement of https://arxiv.org/abs/1404.1038 for this formulas
// // //     double mu = 0.37 * eF; // use local density approximation to estimate mu
// // //     double R = sqrt(mu / params[1]); // Thomas-Fermi radius
// // //     double C0 = params[1] / R;
// // //     params[5] = params[5] * C0; 
}

void modify_potentials_u_ext_smooth_HO(int it, double *h_densities, double *h_potentials)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=0;
    double _x, _y, _r;
    double arg, abs_delta;
    if(dc_params[21] < 1.0) for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        _x = ((double)(ix) - 1.0*(NX/2) + 1.0)*DX;
        _y = ((double)(iy) - 1.0*(NY/2) + 2.0)*DY;
        _r = sqrt(_x*_x + _y*_y);
        
        if(_r<0.01)
        {
            delta[ixyz]=0.0 + I*0.0;
        }
        else if(_r<14.1)
        {
            abs_delta = cabs(delta[ixyz]);
            arg = atan2(_y, _x);
    
            // phase imprint
            delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);         
        }

        
        ixyz++;
    }
}
#undef xs
#undef ys
#undef HO

void modify_potentials_phase_random(int it, double *h_densities, double *h_potentials)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=0;
    double _x, _y, _r;
    double arg, abs_delta;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        abs_delta = cabs(delta[ixyz]);
        arg = (double)rand() / (double)RAND_MAX;
        arg = (2.0*arg-1.0)*M_PI;

        // phase imprint
        delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);         
        
        ixyz++;
    }
}

void modify_potentials_phase_imprint(int it, double *h_densities, double *h_potentials, double *extra_data)
{
    // densities - decode 
    double *rho_a = (double *)(h_densities +  0*NX*NY);
    double *rho_b = (double *)(h_densities +  1*NX*NY);
    double *tau_a = (double *)(h_densities +  2*NX*NY);
    double *tau_b = (double *)(h_densities +  3*NX*NY);
    double complex *nu = (double complex *)(h_densities +  4*NX*NY);
    double *j_a_x = (double *)(h_densities +  6*NX*NY);
    double *j_a_y = (double *)(h_densities +  7*NX*NY);
    double *j_a_z = (double *)(h_densities +  8*NX*NY);
    double *j_b_x = (double *)(h_densities +  9*NX*NY);
    double *j_b_y = (double *)(h_densities + 10*NX*NY);
    double *j_b_z = (double *)(h_densities + 11*NX*NY);
    
    // pontentials - decode
    double *V_a = (double *)(h_potentials +  0*NX*NY);
    double *V_b = (double *)(h_potentials +  1*NX*NY);
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=0;
    double arg, abs_delta;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        abs_delta = cabs(delta[ixyz]);
        arg = extra_data[ixyz];

        // phase imprint
        delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);         
        
        ixyz++;
    }
}

/**
 * Flat well for QT studies
 * 0 - R_min - inside this radius potential is 0.0
 * 1 - R_max - outside this radius potential is V_max 
 * 2 - V_max - value of the potential ouside R_max tube, in units of eF
 * 3 - 
 * 4 - swich time, in number of iterations
 * 5 - number of iterations with imprint
 * */

double u_ext_QT_tube(int ix, int iy, int it, int spin)
{
    double _x = (double)(ix) - 1.0*(NX/2);
    double _y = (double)(iy) - 1.0*(NY/2);
    double _r = sqrt(_x*_x + _y*_y); // distance from the center
        
    double swtch=1.0;
    if(1.0*it<dc_params[4]) swtch=switch_function(1.0*it, dc_params[4], 1.0);
    
    // ---------------------------- potential ----------------------------
    double _u = 0.0;
    if(_r>=dc_params[1])      _u = dc_params[2];
    else if(_r>=dc_params[0]) _u = dc_params[2]*switch_function(_r-dc_params[0], dc_params[1]-dc_params[0], 1.0);

    
    return swtch*_u;
}
void process_params_u_ext_QT_tube(double *params, double kF)
{
    kF=1.0;
    double eF = 0.5*kF*kF;
    
    params[2]*=eF;
    
    // rotating frame
    dc_Omega_a = 0.0; // 0.025;
    dc_Omega_b = 0.0; // 0.025;
}
void modify_potentials_u_ext_QT_tube_lda(int it, double *h_densities, double *h_potentials)
{
  // pontentials - decode
  double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
  double *rho_a = (double *)(h_densities +  0*NX*NY);
  double *rho_b = (double *)(h_densities +  1*NX*NY);
    
  int ixyz;
  double eF;
  for(ixyz=0; ixyz<NXYZ; ixyz++) 
    {
      eF = 0.5 * pow(3.0*M_PI*M_PI*(rho_a[ixyz]+rho_b[ixyz]),2./3.);
      delta[ixyz] = 0.5* eF + I*0.0; // local density approximation
    }
}

void modify_potentials_u_ext_vortex_imprint(int it, double *h_densities, double *h_potentials)
{
    // pontentials - decode
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ixyz=0, ix, iy;
    double _x, _y, _r;
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++)
    {
        _x = (double)(ix) - 1.0*(NX/2);
        _y = (double)(iy) - 1.0*(NY/2);
        _r = sqrt(_x*_x + _y*_y);
            
        if(_r<0.1)
        {
            delta[ixyz]=0.0 + I*0.0;
        }
        else if(_r<dc_params[0])
        {
            double abs_delta = cabs(delta[ixyz]);
            double arg = atan2(_y, _x);
    
            // phase imprint
            delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);         
        }
    
        ixyz++;
    }
}

/**
 * Flat well for QT studies, formula follows experimental box-shape potential
 * Potential in radial direction has form: V(r) = alphap*(r/betap)**16.0
 * for more see gnuplot script "pot_flat.gp"
 * PARAMTERS:
 * 0 - alphap, in input file value of potential at the boundry of the trap, in units of eF
 * 1 - betap, in lattice units
 * 2 - Umax - computed atomaticaly
 * 3 - 
 * 4 - swich time, in number of iterations
 * 5 - number of iterations with imprint
 * MAKRO-parameters
 * rs - distance where potential changes into smooting form
 * */
#define rs (0.925*NX/2)
#define HO(r,alphap,betap) (alphap*pow(r / betap, 16.0))
double u_ext_QT_tube2(int ix, int iy, int it, int spin)
{
    double _x = (double)(ix) - 1.0*(NX/2);
    double _y = (double)(iy) - 1.0*(NY/2);
    double _r = sqrt(_x*_x + _y*_y); // distance from the center
        
    double swtch=1.0;
    if(1.0*it<dc_params[4]) swtch=switch_function(1.0*it, dc_params[4], 1.0);
    
    // ---------------------------- potential ----------------------------
    double _u, _t;
    if(_r<=rs) 
    {
        _u = HO(_r, dc_params[0], dc_params[1]);
    }
    else       
    {   
        _t = switch_function(_r-rs,NX/2-rs,1.0);
        _u = (1.0-_t)*HO(_r, dc_params[0], dc_params[1]) + _t*dc_params[2];
    }

    
    return swtch*_u;

}
void process_params_u_ext_QT_tube2(double *params, double kF)
{
    // in converget state I should get eF in the cloud center about 0.5 (kF=1.0)
    double eF = 0.5;
    // alphap = params[0]
    // betap = params[1]
    params[0]*=eF;
    params[0] = params[0] * pow( params[1] / rs, 16.0);
    // params[2] = Umax=HO(rs) + 1.0*(HO(R)-HO(rs))
    params[2] = HO(rs, params[0], params[1]) + 1.0*(HO(NX/2, params[0], params[1])-HO(rs, params[0], params[1]));
    
//     // rotating frame
//     dc_Omega_a = 0.03;
//     dc_Omega_b = 0.03;
    
    // rotating frame
    dc_Omega_a = 0.0006;
    dc_Omega_b = 0.0006;
}


void modify_potentials_u_ext_QT_tube(int it, double *h_densities, double *h_potentials)
{
    // pontentials - decode
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=-1;
    double ixtr, iytr;
    double tr_a=8.; // side of triangle    
    double imprintR = 0.35*tr_a; // sqrt(3.)*tr_a/4.;
    int oddeven;
    
    double _x, _y, _r;
    double arg, abs_delta;
    
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) // iterate over all lattice points
    {
        ixyz++;
        
        
        // check proximity to traingular lattice point
        for(ixtr=0.0; ixtr<1.0*NX; ixtr+=tr_a) for(iytr=0.0; iytr<1.0*NY; iytr+=sqrt(3.0)*tr_a/2.0)
        {
            // relative coordinate
            oddeven = (int)( round(iytr/( sqrt(3.0)*tr_a/2.0 )) );
            if(oddeven%2 == 0) _x = (double)(ix) - ixtr;
            else               _x = (double)(ix) - (ixtr+0.5*tr_a); 
            _y = (double)(iy) - iytr;
            
            // distance
            _r = sqrt(_x*_x + _y*_y);
            
            if(_r<0.01)
            {
                delta[ixyz]=0.0 + I*0.0;
            }
            else if(_r<imprintR)
            {
                abs_delta = cabs(delta[ixyz]);
                arg = atan2(_y,_x);
        
                // phase imprint
                delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);  
            }            
        }
        
    }
}

void imprint_vortex_u_ext_QT_tube(int it, double *h_densities, double *h_potentials)
{
    // pontentials - decode
    double complex *delta = (double complex *)(h_potentials +  2*NX*NY);
    
    int ix, iy, ixyz=-1;
    
    double _x, _y, _r;
    double arg, abs_delta;
    
    for(ix=0; ix<NX; ix++) for(iy=0; iy<NY; iy++) // iterate over all lattice points
    {
        ixyz++;
        
        _x = (double)(ix) - 1.0*(NX/2);
        _y = (double)(iy) - 1.0*(NY/2);
        _r = sqrt(_x*_x + _y*_y); // distance from the center
    
        if(_r<1.0*(NX/2)) 
        {
            abs_delta = cabs(delta[ixyz]);
            arg = atan2(_y,_x);             
            delta[ixyz]=abs_delta*cos(arg) + I*abs_delta*sin(arg);
        }            
    }
}
#undef rs
#undef HO

// *******************************************************************************************
// **************************** FUNCTIONS CALLED BY THE CODE *********************************
// *******************************************************************************************
/** 
 * THIS FUNCTION IS CALLED BY KERLNELS FROM 'pca_kernels.cu'
 * */
double u_ext(int ix, int iy, int it, int spin)
{
#ifdef UNIFORM_TEST_MODE
    return 0.0; // no external potential
#endif

//     return 0.0; // no external potential
    
    return u_ext_tube(ix, iy, it, spin);

//       return u_ext_HO(ix, iy, it, spin); 
      
//     return u_ext_smooth_HO(ix, iy, it, spin); // glitch studies
//      return u_ext_QT_tube(ix, iy, it, spin);
//     return u_ext_QT_tube2(ix, iy, it, spin);
}

/** 
 * THIS FUNCTION IS CALLED FROM 'pca.c'
 * AFTER LOADING params ARRAY FROM INPUT FILES.
 * AFTER PROCESSING THE params ARE LOADED TO CONST MEMORY ON GPU.
 * THE PARAMS ARE VISIBLE IN dc_params
 * @param params array of size MAX_USER_PARAMS with parameters
 * @param kF typical Fermi momentum of the problem
 * */
void process_params(double *params, double kF)
{
#ifndef UNIFORM_TEST_MODE
    // no processing
    process_u_ext_tube(params, kF);
//      process_params_u_ext_HO(params, kF);
//      process_params_u_ext_smooth_HO(params, kF); // glitch studies
//       process_params_u_ext_QT_tube(params, kF);
//     process_params_u_ext_QT_tube2(params, kF);
#endif
}

/**
 * THIS FUNCTIONS IS CALLED AFTER EACH EXECUTION `recompute_potentials`
 * */
void modify_potentials(int it, double *h_densities, double *h_potentials, double *extra_data)
{
#ifndef UNIFORM_TEST_MODE
//     if(it<dc_params[5]) modify_potentials_imprint_ferron(it, h_densities, h_potentials, extra_data);
    
    modify_potentials_phase_imprint2(it, h_densities, h_potentials, extra_data);
    
    //no modify
//     if(it<dc_params[5]) modify_potentials_phase_random(it, h_densities, h_potentials);
//     if(it<dc_params[5]) modify_potentials_phase_imprint(it, h_densities, h_potentials, extra_data);
    
//     if(it<dc_params[5])
//         modify_potentials_u_ext_smooth_HO(it, h_densities, h_potentials);
    
//      if(it<dc_params[5])
//      {
//        //         imprint_vortex_u_ext_QT_tube(it, h_densities, h_potentials);
// // //         modify_potentials_u_ext_QT_tube(it, h_densities, h_potentials);
// //         modify_potentials_u_ext_QT_tube_lda(it, h_densities, h_potentials);
//      
//         modify_potentials_u_ext_vortex_imprint(it, h_densities, h_potentials);
//      }
#endif
}
#endif

