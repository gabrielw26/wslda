# source env.sh

module load impi/2020.4.912
module load mkl/2020.0.4
module load cuda/11.2.1_460.32.03
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/mnt/storage_2/project_data/grant_518/share/elpa/lib/ 
module list
cat $WSLDA/VERSION.h
