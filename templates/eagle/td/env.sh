# source env.sh

module load openmpi/3.1.4_gcc620
module unload cudatoolkit
module load cuda/11.2.1_460.32.03
module list
cat $WSLDA/VERSION.h
