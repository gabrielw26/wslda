/**
 * W-SLDA Toolkit
 * */


#ifndef __W_DATA_LIB__
#define __W_DATA_LIB__

#include <stdio.h>

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
    char name[MD_VARNAME_LGTH];
    double value;
} wdata_const;

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
    int nconsts;               /// number of constants
    
    // variables
    wdata_variable vars[WDATA_MAX_NVARS];
    
    // links
    wdata_link links[WDATA_MAX_NVARS];

    // constants
    wdata_const consts[WDATA_MAX_NVARS];
    
    // auxliary vars
    int issetwrkdir; 
    char wrkdir[MD_CHAR_LGTH]; // working directory
    
} wdata_metadata; 

/** 
 * Function reads metadata
 * and puts values into struct wdata_metadata
 * @return 0: ok, 1: Cannot open metadata file,
 * */
int wdata_parse_metadata_file(const char * file_name, wdata_metadata *md);

void wdata_print_metadata(wdata_metadata *md, FILE *out);

void wdata_print_variable(wdata_variable *md, FILE *out);

void wdata_print_link(wdata_link *md, FILE *out);

void wdata_print_const(wdata_const *md, FILE *out);

void wdata_add_variable(wdata_metadata *md, wdata_variable *var);

void wdata_add_link(wdata_metadata *md, wdata_link *link);

void wdata_add_const(wdata_metadata *md, wdata_const *_const);

int wdata_get_blocksize(wdata_metadata *md);

size_t wdata_get_blocksize_bytes(wdata_metadata *md, wdata_variable *var);

/**
 * Low level function. It adds block of data to binary file.
 * @param md metadata for data set
 * @param var variable to be added to file with name `prefix`_`varname`.wdat
 * @param data pointer to binary data (INPUT)
 * @return 0: ok; 1: cannot open binary file; 2: cannot add datablock to file
 * */
int wdata_add_datablock(wdata_metadata *md, wdata_variable *var, void *data);

/**
 * Function adds new block to data file.
 * @param md metadata for data set
 * @param varname name of variable, can be from list of vars or links
 * @param data pointer to binary data (INPUT)
 * @return 0: ok; 1: cannot open binary file; 2: cannot add datablock to file, 11: variable is not defined
 * */
int wdata_write_cycle(wdata_metadata *md, const char *varname, void *data);

/**
 * Function reads block of data from file
 * @param md metadata for data set
 * @param varname name of variable, can be from list of vars or links
 * @param data pointer to binary data (OUTPUT)
 * @return 0: ok; 1: cannot open binary file; 2: cannot read data block from file; 3: cannot shift pointer;  11: variable is not defined
 * */
int wdata_read_cycle(wdata_metadata *md, const char *varname, int cycle, void *data);

void wdata_get_filename(wdata_metadata *md, wdata_variable *var, char *file_name);

/**
 * Functions extracts variable corresponding to given name
 * @param md metadata for data set
 * @param varname name of variable, can be from list of vars or links
 * @param var pointer to variable from md structure (OUTPUT)
 * @return 0: ok; 1: cannot find variable
 * */
int wdata_get_variable(wdata_metadata *md, const char *varname, wdata_variable *var);

int wdata_get_const(wdata_metadata *md, const char *constname, wdata_const *_const);

/**
 * @return value of constant from metadata structure, if const is not found zero is returned
 * */
double wdata_getconst_value(wdata_metadata *md, const char *constname);

/**
 * Functions sets value of constant. If constant was not added before it adds it and sets its value.
 * */
void wdata_setconst(wdata_metadata *md, const char *constname, double constvalue);

/**
 * Function checks if binary file exists for variable.
 * @return 1 if binary file for this variable exists, otherwise 0
 */
int wdata_file_exists(wdata_metadata *md, const char *varname);

/**
 * Function removes binary file associated with variable varname.
 * */
void wdata_clear_file(wdata_metadata *md, const char *varname);

/**
 * Removes all binary data files, except metadata file, and sets cycles to zero.
 * */
void wdata_clear_database(wdata_metadata *md);

/**
 * Functions increases number of cycles.
 * It is equivalent to: md->cycles++;
 * */
int wdata_add_cycle(wdata_metadata *md);

/**
 * Function writes metadata to file.
 * @param md metadata to be written
 * @param filename name of file to write md, if filename is empty string then default name of file will be used `prefix`.wtxt
 * @return 0: ok, 1: cannot create file
 * */
int wdata_write_metadata_to_file(wdata_metadata *md, const char * filename);

/**
 * Function sets working dir for given metadata.
 * Working dir defines directory where all binary files will be stored.
 * Default is working dir of code that executes wdata functions.
 * */
void wdata_set_working_dir(wdata_metadata *md, const char * wkrdir);


#endif



