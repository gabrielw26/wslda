# source env.sh

# module load LUMI
# module load rocm
# module load rocmlibs
# module load rocThrust

module load LUMI/22.06
module load rocm/5.1.4
module load craype-accel-amd-gfx90a
module load cray-fftw
module list

cat $WSLDA/VERSION.h
