#!/bin/bash

# This script converts the code to HIP environment.
# To run it, you need to have installed:
#   hipify-perl
# Make sure you have exported the path to the tool to PATH variable.
# For example:
#   export PATH=$PATH:/opt/HIPIFY/bin

echo "---> Converting WSLDA to HIP environment"
set -x

# Create & clear location for hip files
cd hpc-engine
mkdir -p hip
rm -rf hip/*

# Convert CUDA files
for cufile in *.cu; do
    hipify-perl -hip-kernel-execution-syntax -o hip/${cufile%.cu}.cpp $cufile
done

set +x
echo "---> Conversion of WSLDA to HIP is DONE."
