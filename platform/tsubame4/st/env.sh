# source env.sh

module load cuda/12.3.2
module load intel/2024.0.2
module load intel-mpi/2021.11
module load openmpi/5.0.2-intel
module load fftw/3.3.10-intel
module list

export MKLROOT=/apps/t4/rhel9/isv/intel/2024.0
export LD_LIBRARY_PATH=$MKLROOT/lib:/gs/bs/hp190063/share/libs/lib:$LD_LIBRARY_PATH
export WSLDA=/gs/bs/hp190063/share/cold-atoms/
export WSLDA_MACHINE=$WSLDA/templates/tsubame4
cat $WSLDA/VERSION.h
