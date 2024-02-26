# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  

# for LUMI system
module load PrgEnv-gnu/8.4.0
module load LUMI/23.09  partition/G
module load rocm/5.6.1
module load cray-fftw
module load craype-accel-amd-gfx90a

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/project_465000708/share/libs/lib

cat $WSLDA/VERSION.h
