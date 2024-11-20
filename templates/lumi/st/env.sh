# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  

# for LUMI system
module load PrgEnv-gnu
module load LUMI  partition/G
module load rocm
module load cray-fftw
module load craype-accel-amd-gfx90a
module load lumi-CrayPath

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/project_465000995/share/libs/elpa-IX2024/elpa_git/lib
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/project_465000995/share/cold-atoms/lib/wdata:/project/project_465000995/share/cold-atoms/lib/wderiv:/project/project_465000995/share/cold-atoms/lib/winterp
export PATH=$PATH:/project/project_465000995/share/cold-atoms/lib/wdata/bin:/project/project_465000995/share/cold-atoms/tools/bin
export WSLDA_MACHINE=$WSLDA/templates/lumi

cat $WSLDA/VERSION.h
