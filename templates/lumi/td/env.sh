# source env.sh

module load LUMI
module load craype-accel-amd-gfx90a
module load rocm
module load rocmlibs
module load rocThrust
module load cray-fftw
module list

cat $WSLDA/VERSION.h
