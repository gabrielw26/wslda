#!/bin/bash

# This script converts the code to HIP environment.
# To run it, you need to have installed:
#   hipify-perl
# Make sure you have exported the path to the tool to PATH variable.
# For example:
#   export PATH=$PATH:/opt/HIPIFY/bin

# SETTINGS
HIPIFY_CMD="hipify-perl -hip-kernel-execution-syntax"

echo "---> Converting WSLDA to HIP environment."
set -x

# Create & clear location for hip files
cd hpc-engine
mkdir -p hip
rm -rf hip/*

# Convert CUDA files
for cufile in *.cu; do
    $HIPIFY_CMD -o hip/${cufile%.cu}.cpp $cufile
done

# Other files
$HIPIFY_CMD -o hip/wslda_cuda_utils.hpp wslda_cuda_utils.h
$HIPIFY_CMD -o hip/tdwslda_energy_generic.hpp tdwslda_energy_generic.h
$HIPIFY_CMD -o hip/tdwslda_ode_integrator.hpp tdwslda_ode_integrator.h

set +x
echo "---> Conversion of WSLDA to HIP is DONE."
