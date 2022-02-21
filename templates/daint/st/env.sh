module load daint-gpu
module swap PrgEnv-cray PrgEnv-gnu/6.0.10
# module load cudatoolkit/11.2.0_3.39-2.1__gf93aa1c
module load cudatoolkit/10.2.89_3.28-2.1__g52c0314
module load craype-accel-nvidia60
module load cray-fftw
export LD_LIBRARY_PATH=/project/pr125/share/elpa-2020.11.001/lib:$LD_LIBRARY_PATH

module list
cat $WSLDA/VERSION.h

