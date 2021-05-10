#!bin/bash

source scl_source enable devtoolset-8 rh-python36
module load openmpi-gcc721-Cuda90
module load cuda/9.0
source env.sh
python3 main.py
