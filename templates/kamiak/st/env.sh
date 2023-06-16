# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  
# export WSLDA=...
module load openmpi/4.1.5
module load liblas/1.8.0
module load openblas/0.3.0
module load fftw/3.3.10
module load cuda/11.0.3
module list 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/apps/blas/3.8.0/gcc/7.3.0/lib:/opt/apps/lapack/3.8.0/gcc/7.3.0/lib/:/opt/apps/scalapack/2.1.0/gcc/7.3.0/openmpi/4.0.3/lib/:/opt/apps/elpa/2018.11.001/gcc/7.3.0/openmpi/4.0.3/lib

cat $WSLDA/VERSION.h

