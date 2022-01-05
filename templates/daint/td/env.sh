# source env.sh

module load daint-gpu
module load cudatoolkit
module load craype-accel-nvidia60
module load cray-fftw

module list
cat $WSLDA/VERSION.h
