/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 15-10-2020
 * 
 * Declaration file for errors reporting
 * */ 
#ifndef __WSLDA_ERRORS__
#define __WSLDA_ERRORS__

#include <stdio.h>
void report_error(int errcode, FILE *stream);

// no error
#define WSLDA_OK  0

// errors
#define WSLDA_ERR_CANNOT_CREATE_DIR 10001
#define WSLDA_ERR_CANNOT_OPEN_FILE 10002
#define WSLDA_ERR_CANNOT_OVERWRITE 10003
#define WSLDA_ERR_BINARY_FILE_CORRUPTED 10004

#endif
