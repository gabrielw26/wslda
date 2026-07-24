#!/bin/bash

# Set the C compiler to use GNU99 standard
export CC="gcc -std=gnu99"
# export CC="cc -std=gnu99"


# Test if WSLDA is already set
if [ -z "$WSLDA" ]; then
    echo "# WSLDA is not set. Setting it to the current directory."
    export WSLDA=$(pwd)
else
    echo "# WSLDA is already set to $WSLDA. Using the existing value."
fi

# go to WSLDA root directory
echo "# Changing to WSLDA root directory: $WSLDA"
cd $WSLDA

# libs
# wderiv
echo "# Making lib/wderiv.."
make -C lib/wderiv

# winterp
echo "# Making lib/winterp.."
make -C lib/winterp

# wdata
echo "# Making lib/wdata.."
make -C lib/wdata

# wbox
echo "# Making lib/wbox.."
make -C lib/wbox

# tools
echo "# Making tools.."
make -C tools

echo "# =========================================================================="
echo "# =========================== SET ENVIRONMENT =============================="
echo "# =========================================================================="
echo "#  Add to your .bashrc"
echo " "
cmd=export\ WSLDA=$(pwd)
echo $cmd
cmd=export\ LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$(pwd)/lib/wdata:$(pwd)/lib/wderiv:$(pwd)/lib/winterp:$(pwd)/lib/wbox
echo $cmd
cmd=export\ PATH=\$PATH:$(pwd)/lib/wdata/bin:$(pwd)/tools/bin
echo $cmd
echo " "
