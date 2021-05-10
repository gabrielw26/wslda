/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 08-03-2021
 * 
 * Static variables of td codes
 * */ 

#ifndef __TDWSLDA_STATIC_VARS__
#define __TDWSLDA_STATIC_VARS__

static double dc_ec;
static double dc_t0;
static int dc_np;
static int dc_nwfip;
static int dc_nwf;

#define TDWSLDA_SET_STATIC_VARS \
  dc_ec=ec; dc_t0=t0; dc_np=np; dc_nwfip=nwfip; dc_nwf=nwf; 

#endif


