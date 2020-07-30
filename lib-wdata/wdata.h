/**
 * W-SLDA Toolkit
 * */


#ifndef __W_DATA_LIB__
#define __W_DATA_LIB__

#define MD_CHAR_LGTH 256
#define MD_VARNAME_LGTH 32
#define WDATA_MAX_NVARS 32

typedef struct
{
    char name[MD_VARNAME_LGTH];
    char type[MD_VARNAME_LGTH];
    char unit[MD_VARNAME_LGTH];
} wdata_variable;

typedef struct
{
    char name[MD_VARNAME_LGTH];
    char linkto[MD_VARNAME_LGTH];
} wdata_link;

typedef struct
{
    int NX;
    int NY;
    int NZ;
    double DX;
    double DY;
    double DZ;
    int datadim;               /// dimension of block size 1=NX, 2=NX*NY, 3=NX*NY*NZ
    char prefix[MD_CHAR_LGTH]; /// prefix for files belonging to this data set, binary files have names prefix_variable.wdat
    int cycles;                /// number of cycles (measurments)
    double t0;                 /// time value for the first mesurments
    double dt;                 /// time interval between cycles    
    int nvars;                 /// number of variables
    int nlinks;                /// number of links
    
    // variables
    wdata_variable vars[WDATA_MAX_NVARS];
    
    // links
    wdata_link links[WDATA_MAX_NVARS];
    
} wdata_metadata; 


int wdata_parse_metadata_file(const char * file_name, wdata_metadata *md);
void wdata_print_metadata(wdata_metadata *md, FILE *out);
void wdata_print_variable(wdata_variable *md, FILE *out);
void wdata_print_link(wdata_link *md, FILE *out);
void wdata_add_variable(wdata_metadata *md, wdata_variable *var);
void wdata_add_link(wdata_metadata *md, wdata_link *link);
int wdata_get_blocksize(wdata_metadata *md);
size_t wdata_get_blocksize_bytes(wdata_metadata *md, wdata_variable *var);
int wdata_add_datablock(wdata_metadata *md, wdata_variable *var, void *data);
int wdata_write_cycle(wdata_metadata *md, const char *varname, void *data);
int wdata_read_cycle(wdata_metadata *md, const char *varname, int cycle, void *data);
void wdata_get_filename(wdata_metadata *md, wdata_variable *var, char *file_name);
int wdata_get_variable(wdata_metadata *md, const char *varname, wdata_variable *var);

#endif



