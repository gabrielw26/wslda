# source env.sh

module load cuda/10.2.89
module load intel/19.1.0.166
module load intel-mpi/19.6.166
module load fftw/3.3.6
module list

export WSLDADIR=/gs/hs1/hp190063/share/wslda/
cat $WSLDA/VERSION.h
