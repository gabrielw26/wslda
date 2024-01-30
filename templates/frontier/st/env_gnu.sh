# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  

# for FRONTIER system
module purge
module load rocm
module load ums/default
module load ums002/default
module load hip/5.3.0-rocm-pwsxi
module load craype-accel-amd-gfx90a
module load cray-python
module load cray-libsci
module load cpe
module load DefApps/default
module load autoconf/2.69
module load automake
module load libtool
module load cray-mpich
module load craype-x86-trento
module load cray-fftw

export WSLDA=/lustre/orion/proj-shared/nfi125/cold-atoms/

export LD_LIBRARY_PATH=/lustre/orion/proj-shared/nfi125/elpa/elpaloc/lib:$LD_LIBRARY_PATH

cat $WSLDA/VERSION.h
