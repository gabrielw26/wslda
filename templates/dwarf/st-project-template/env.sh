# source env.sh

source scl_source enable devtoolset-7 python27
module load cuda/9.0
module load openmpi-gcc721-Cuda90/3.1.1
module load lapack/391
module load scalapack-gcc721-cuda90-openmpi311/210
module load elpa-gcc721/202005
module list

cat $WSLDA/VERSION.h
