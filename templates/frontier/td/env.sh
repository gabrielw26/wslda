# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  

export WSLDA=/lustre/orion/proj-shared/nfi125/cold-atoms/

export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:$LD_LIBRARY_PATH

export MPICH_GPU_SUPPORT_ENABLED=1

# for FRONTIER system
module purge
module load PrgEnv-cray
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module load cray-fftw
module load cray-pmi/6.1.8
module load amd-mixed/5.3.0
module load xpmem/2.6.2-2.5_2.22__gd067c3f.shasta
module load perftools-base

echo $WSLDA
cat $WSLDA/VERSION.h
