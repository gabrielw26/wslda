/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 15-10-2020
 * 
 * This file constains funtions for generation of ky and kz plane waves for 2D and 1D codes
 * */ 

#ifndef __WSLDA_WAVEVECTORS__
#define __WSLDA_WAVEVECTORS__


// ------------------------------- IMPROVED INTERFACE FOR 1D CODES -------------------------------
typedef struct
{
    int cnt; // number of unique modes for 1d case
    double *ky; 
    double *kz;
    int *weight; // degeneracy of the mode
} wslda_kmodes_1d;

static wslda_kmodes_1d _kmodes_1d = {0, NULL, NULL, NULL}; // local structure

int fill_wslda_kmodes_1d(double *kky, double *kkz)
{
    double *lky, *lkz, *lk2;
    int *lweight;
    int lcnt=0;
    
    // take memory
    cppmallocl(lky,NY*NZ,double);
    cppmallocl(lkz,NY*NZ,double);
    cppmallocl(lk2,NY*NZ,double);
    cppmallocl(lweight,NY*NZ,int);
    
    // take only unique kx^2+ky^2
    int iz, iy, i;
    int hasit;
    double k2;
    for(iy=0; iy<NY; iy++) for(iz=0; iz<NZ; iz++)
    {
        k2=kky[iy]*kky[iy] + kkz[iz]*kkz[iz];
        
        if(iy==NY/2) continue; // state without (+,-) pair, SKIP for 1D mode 
        if(iz==NZ/2) continue; // state without (+,-) pair, SKIP for 1D mode 
        
        hasit=0;
        for(i=0; i<lcnt; i++) if(fabs(lk2[i]-k2)<1.0e-9) {hasit=1; break;}
        
        if(hasit==1)
        {
            lweight[i]++;
//             if(wsldapid==0) printf("->>>> HAS %12.6f %12.6f %12.6f %6d\n", kky[iy], kkz[iz], k2, i);
        }
        else
        {
            lweight[lcnt]=1;
            lky[lcnt]=kky[iy];
            lkz[lcnt]=kkz[iz];
            lk2[lcnt]=k2;
            lcnt++;
//             if(wsldapid==0) printf("->>>> NEW %12.6f %12.6f %12.6f %6d\n", kky[iy], kkz[iz], k2, lcnt-1);
        }
    }
    
    // copy to destination buffer
    cppmallocl(_kmodes_1d.ky,lcnt,double);
    cppmallocl(_kmodes_1d.kz,lcnt,double);
    cppmallocl(_kmodes_1d.weight,lcnt,int);
    _kmodes_1d.cnt=lcnt;
        
    int isum=0;
    for(i=0; i<lcnt; i++)
    {
        _kmodes_1d.ky[i]=lky[i];
        _kmodes_1d.kz[i]=lkz[i];
        _kmodes_1d.weight[i]=lweight[i];
        
        // check internal sum
        isum+=lweight[i];
    }
    
    // free memory
    free(lky);
    free(lkz);
    free(lk2);
    free(lweight);
    
//     printf("->>>> lcnt=%d isum=%d NY*NZ=%d\n", lcnt, isum, (NY-1)*(NZ-1));
    if(isum!=(NY-1)*(NZ-1)) return WSLDA_ERR_INTRISTIC_ERROR;
    
    return WSLDA_OK;
}


// ------------------------------- OLD INTERFACE FOR 1D AND 2D CODES -------------------------------

/**
 * Structure containg info about k mode for 2D and 1D codes only
 * */
typedef struct
{
    double ky; // if codedim==2 then ky is set to zero
    double kz;
    int weight; // degeneracy of the mode
} wslda_kmode; 

/**
 * Count number of plane waves to be considerd in the calculatons
 * */
int count_number_of_k_modes(double *kkx, double *kky, double *kkz, int codedim, int *k_modes)
{
    int kvecs_to_consder = 0;
    int iy=0; 
    int iz=0;
    int lNZ=NZ;
    int lNY=NY;
    
    if(codedim==1) // apply improved interface
    {
        int ierr=WSLDA_OK;
        if(_kmodes_1d.cnt==0) ierr=fill_wslda_kmodes_1d(kky, kkz); 
        kvecs_to_consder = _kmodes_1d.cnt;
        
        *k_modes = kvecs_to_consder; // save result
        return ierr;
    }
    
    // otherwise, use old method
    if(codedim==2) lNY=1;
    
    // take only positive energy states
    for(iz=0; iz<lNZ; iz++) for(iy=0; iy<lNY; iy++)
        {
            if     (kky[iy]>1.0e-12       && kkz[iz]>1.0e-12      ) kvecs_to_consder += 1; // 2*2; // one solution x (-ky, +ky) x (-kz, +kz)
            else if(kky[iy]>1.0e-12       && fabs(kkz[iz])<1.0e-12) kvecs_to_consder += 1; // 2*1; // one solution x (-ky, +ky) x (  kz=0  )
            else if(fabs(kky[iy])<1.0e-12 && kkz[iz]>1.0e-12      ) kvecs_to_consder += 1; // 1*2; // one solution x (  ky=0  ) x (-kz, +kz)
            else if(fabs(kky[iy])<1.0e-12 && fabs(kkz[iz])<1.0e-12) kvecs_to_consder += 1; // 1*1; // one solution x (  ky=0  ) x (  kz=0  )
        }
        
    *k_modes = kvecs_to_consder; // save result
    return WSLDA_OK;
}

/**
 * Fill structure of k_modes with data
 * */
int create_k_modes(double *kkx, double *kky, double *kkz, int codedim, wslda_kmode *k_modes)
{
    int kvecs_to_consder = -1;
    int iy=0; 
    int iz=0;
    int lNZ=NZ;
    int lNY=NY;
    
    if(codedim==1) // apply improved interface
    {
        int ierr=WSLDA_OK;
        if(_kmodes_1d.cnt==0) ierr=fill_wslda_kmodes_1d(kky, kkz); 

        int i;
        for(i=0; i<_kmodes_1d.cnt; i++)
        {
            k_modes[i].ky=_kmodes_1d.ky[i];
            k_modes[i].kz=_kmodes_1d.kz[i];
            k_modes[i].weight=_kmodes_1d.weight[i];
        }
            
        return ierr;
    }
    
    // otherwise, use old method
    if(codedim==2) lNY=1;
    
    // take only positive energy states
    for(iz=0; iz<lNZ; iz++) for(iy=0; iy<lNY; iy++)
        {
            if     (kky[iy]>1.0e-12       && kkz[iz]>1.0e-12      ) 
            {
                kvecs_to_consder += 1; 
                k_modes[kvecs_to_consder].ky=kky[iy];
                k_modes[kvecs_to_consder].kz=kkz[iz];
                k_modes[kvecs_to_consder].weight=2*2; // (-ky, +ky) x (-kz, +kz)
            }
            else if(kky[iy]>1.0e-12       && fabs(kkz[iz])<1.0e-12)
            {
                kvecs_to_consder += 1; 
                k_modes[kvecs_to_consder].ky=kky[iy];
                k_modes[kvecs_to_consder].kz=kkz[iz];
                k_modes[kvecs_to_consder].weight=2*1; // (-ky, +ky) x (  kz=0  )
            }
            else if(fabs(kky[iy])<1.0e-12 && kkz[iz]>1.0e-12      )
            {
                kvecs_to_consder += 1; 
                k_modes[kvecs_to_consder].ky=kky[iy];
                k_modes[kvecs_to_consder].kz=kkz[iz];
                k_modes[kvecs_to_consder].weight=1*2; //  (  ky=0  ) x (-kz, +kz)
            }
            else if(fabs(kky[iy])<1.0e-12 && fabs(kkz[iz])<1.0e-12)
            {
                kvecs_to_consder += 1; 
                k_modes[kvecs_to_consder].ky=kky[iy];
                k_modes[kvecs_to_consder].kz=kkz[iz];
                k_modes[kvecs_to_consder].weight=1*1; //  (  ky=0  ) x (  kz=0  )
            }
        }
        
    return WSLDA_OK;
}

/**
 * Create wave vectors, x-direction
 * */
int create_kkx(double *kkx)
{
    int i,j;
    
    for ( i = 0 ; i <= NX / 2 - 1 ; i++ ) {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) i ; }
    j = - i ;
    for ( i = NX / 2 ; i < NX ; i++ ) 
    {
        kkx[ i ] = 2. * ( double ) M_PI / LX * ( double ) j ;
        j++ ;
    }
    
    return WSLDA_OK;
}

/**
 * Create wave vectors, y-direction
 * */
int create_kky(double *kky)
{
    int i,j;
    
    for ( i = 0 ; i <= NY / 2 - 1 ; i++ ) {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) i ; }
    j = - i ;
    for ( i = NY / 2 ; i < NY ; i++ ) 
    {
        kky[ i ] = 2. * ( double ) M_PI / LY * ( double ) j ;
        j++ ;
    }
    
    return WSLDA_OK;
}

/**
 * Create wave vectors, z-direction
 * */
int create_kkz(double *kkz)
{
    int i,j;
    
    for ( i = 0 ; i <= NZ / 2 - 1 ; i++ ) {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) i ; }
    j = - i ;
    for ( i = NZ / 2 ; i < NZ ; i++ ) 
    {
        kkz[ i ] = 2. * ( double ) M_PI / LZ * ( double ) j ; 
        j++ ;
    }
    
    return WSLDA_OK;
}

/**
 * Get weight for density computation in 1d
 * */
int get_weight_1d(double ky, double kz)
{
    int wcnt = 1; 
    if(fabs(ky)>1.0e-12) wcnt*=2; // account for -ky and +ky 
    if(fabs(kz)>1.0e-12) wcnt*=2; // account for -kz and +kz 
    return wcnt;
}

/**
 * Get weight for density computation in 1d
 * */
int get_weight_2d(double kz)
{
    int wcnt = 1;  
    if(fabs(kz)>1.0e-12) wcnt*=2; // account for -kz and +kz 
    return wcnt;
}


#endif

