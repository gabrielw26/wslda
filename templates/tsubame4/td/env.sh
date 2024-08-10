# source env.sh

module load openmpi/5.0.2-gcc
module load fftw/3.3.10-gcc
module list

export WSLDA=/gs/bs/hp190063/share/cold-atoms/
export WSLDA_MACHINE=$WSLDA/templates/tsubame4
cat $WSLDA/VERSION.h
