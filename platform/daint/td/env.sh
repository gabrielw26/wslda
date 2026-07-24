# source env.sh

module load daint-gpu
module load cudatoolkit/10.2.89_3.28-2.1__g52c0314
module load craype-accel-nvidia60
module load cray-fftw

module list
cat $WSLDA/VERSION.h
