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

export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH:/project/project_465002810/share/libs/lib
export ELPALIB=/project/project_465002810/share/libs/lib
export ELPAINCLUDE=/project/project_465002810/share/libs/include/elpa-2026.02.001

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/project_465002810/share/cold-atoms/lib/wdata:/project/project_465002810/share/cold-atoms/lib/wderiv:/project/project_465002810/share/cold-atoms/lib/winterp
export PATH=$PATH:/project/project_465002810/share/cold-atoms/lib/wdata/bin:/project/project_465002810/share/cold-atoms/tools/bin
export WSLDA_MACHINE=$WSLDA/platform/lumi

cat $WSLDA/VERSION.h
