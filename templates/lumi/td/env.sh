# source env.sh

module load LUMI
module load rocm
module load craype-accel-amd-gfx90a
module load cray-fftw
module list

cat $WSLDA/VERSION.h
