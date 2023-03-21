#!/bin/bash
export CC="cc -std=gnu99"
export MPICPP="cc -lstdc++"

echo "Installing libs..."

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
make -C lib/wdata

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
cmd=export\ LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$(pwd)/lib/wdata:$(pwd)/lib/wderiv:$(pwd)/lib/winterp
echo $cmd
cmd=export\ PATH=\$PATH:$(pwd)/lib/wdata/bin:$(pwd)/tools/bin
echo $cmd
echo " "
