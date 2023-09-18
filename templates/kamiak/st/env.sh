# 
# W-SLDA Toolkit
# Warsaw University of Technology
#  
# export WSLDA=...
module load intel/20.2
module load fftw/3.3.10
module load cuda/11.0.3
module list 
export ELPA_HOME=/data/lab/forbes/apps/elpa
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ELPA_HOME/lib:$MKL_HOME/lib/intel64

cat $WSLDA/VERSION.h

