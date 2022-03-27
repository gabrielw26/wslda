
#ifndef __PREDEFINES__
#define __PREDEFINES__
/**
 * Define lattice size and lattice spacing
 * */
#define NX 256
#define NY 32
#define NZ 32

#define DX 1.0
#define DY 1.0
#define DZ 1.0

/**
 * Maximal number of parameters in params array
 * */
#define MAX_USER_PARAMS 32 


/**
 * Machine file. 
 * This file contains info about machine that will be used in the computation process.
 * You can specify the file in the following ways: 
 * - copy `machine.h` file to the current directory, see templates folder for various examples,
 * - specify the folder with `machine.h` file via -I option in Makefile  
 * - use system variable WSLDA_MACHINE to specify the folder with `machine.h` file, for example:
 *   export WSLDA_MACHINE=...
 **/ 
#include "machine.h"

/**
 * Files: predefines.h, problem-definition.h, logger.h are assumed to be compatible with this API version
 * For list of API versions see: https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/API-version
 * */
#define API_VERSION 20220221
#endif
