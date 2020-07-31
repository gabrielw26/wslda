/**
 * W-SLDA Toolkit
 * Author: Gabriel Wlazlowski
 * Creation date: 31-07-2020
 * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "wdata.h"

#include "pca_settings.h"
#include "pca_utils.h"

/**
 * @param input input structure (INPUT)
 * @param datadim data dimensonality (INPUT)
 * @param t0 initial  value of time
 * @param dt time increment between measurments
 * @param spinsymmetry flag inicating if system is spin-symmetric
 * @param wdmd metadata for wdata format (OUTPUT)
 * @return exit status 0: ok, 
 * */
int create_wdata_metadata(metadata_t *input, int datadim, double t0, double dt, int spinsymmetry, wdata_metadata *wdmd)
{
    // create artificial data for visulisation in visit
    wdata_metadata tmd = {NX, NY, NZ, DX, DY, DZ, 0, "none", 0, 0.0, 0.0, 0, 0};
    tmd.datadim=datadim;
    sprintf(tmd.prefix, "%s", md.outprefix);
    tmd.t0=t0;
    tmd.dt=dt;
    
    // add variables
    int lnvars;
    char lvars[MAX_WRITEVARS][MAX_VARNAME_LGTH]; // and their names
    int i;
    
    // TODO
    // simple copy of variables to list
    lnvars = md.nwritevars;
    for(i=0; i<lnvars; i++) strcpy(lvars[i],md.writevars[i]);
    
    

    for(i=0; i<lnvars; i++)
    {        
        if(strcmp (lvars[i],"density") == 0)
        {
            wdata_variable va = {"density_a", "real", "none"};
            wdata_variable vb = {"density_b", "real", "none"};
            wdata_link l = {"density_b", "density_a"};
            wdata_add_variable(&tmd, &va);
            if(spinsymmetry==0) wdata_add_variable(&tmd, &vb);
            else                wdata_add_link(&tmd, &l);
        }
    }
    
    *wdmd = tmd; // copy to output buffers
    return 0;
}


