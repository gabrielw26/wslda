# source env.sh

source scl_source enable devtoolset-7 python27
module load cuda/9.0
module load openmpi-gcc721-Cuda90/3.1.1
module list

export WSLDA_MACHINE=$WSLDA/templates/dwarf
cat $WSLDA/VERSION.h
