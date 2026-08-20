
source scl_source enable devtoolset-8 rh-python36
module load cuda/9.0
module load openmpi-gcc721-Cuda90/3.1.1
module load lapack/391
module load scalapack-gcc721-cuda90-openmpi311/210
module load elpa-gcc721/202005
# show modules
module list

# path to W-SLDA engine
export WSLDA=/home/gabrielw/tmp/GW-branch/cold-atoms/
