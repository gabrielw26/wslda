
#define params dc_params
#define extra_data dc_extra_data
#define extra_data_size dc_extra_data_size

#if CODEDIM==1
#define NUMBER_ELEMENT NX
#define VOLUME_ELEMENT (DXYZ*NY*NZ)
#elif CODEDIM==2
#define NUMBER_ELEMENT NXY
#define VOLUME_ELEMENT (DXYZ*NZ)
#elif CODEDIM==3
#define NUMBER_ELEMENT NXYZ
#define VOLUME_ELEMENT (DXYZ)
#else
#error CODEDIM NOT SPECIFIED
#endif

// -------------------- generic decoding functions ------------------------------------
#if CODEDIM==3
    
#define  decode_ixyz2ixiyiz(ixyz,_ix,_iy,_iz,i) \
    ixyz2ixiyiz(ixyz,_ix,_iy,_iz,i)

#elif CODEDIM==2
    
#define  decode_ixyz2ixiyiz(ixyz,_ix,_iy,_iz,i) \
    ixy2ixiy2d(ixyz,_ix,_iy)                    \
    _iz=0;
    
#elif CODEDIM==1
    
#define  decode_ixyz2ixiyiz(ixyz,_ix,_iy,_iz,i) \
    _ix=ixyz;                                   \
    _iy=0;                                      \
    _iz=0;
    
#else
    
#error CODEDIM NOT SPECIFIED
    
#endif
    
