# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  

# for LUMI system
module load PrgEnv-gnu/8.3.3
module load LUMI/22.06
module load rocm/5.1.4
module load cray-fftw
module load craype-accel-amd-gfx90a

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/project_465000708/share/elpa-libsci/PrgEnv-gnu-8.3.3/LUMI-22.06/rocm-5.1.4/lib

cat $WSLDA/VERSION.h
