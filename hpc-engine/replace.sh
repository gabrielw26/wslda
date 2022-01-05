# This script can be used Find&Replace opration applied to many files

# find -name "what" -exec sed -i 's+find+replace+g' {} +

# find -name "Makefile*" -exec sed -i 's+$(WSLDADIR)hpc-engine+$(WSLDADIR)/hpc-engine+g' {} +
# find -name "Makefile*" -exec sed -i 's+(must end with /)++g' {} +
# find -name "Makefile*" -exec sed -i 's+$(WSLDA)/+$(WSLDA)+g' {} +

find -name "Makefile" -exec sed -i 's+instalation+installation+g' {} +
find -name "Makefile" -exec sed -i 's+(for read only)+(for read-only)+g' {} +
find -name "Makefile" -exec sed -i 's+temporarry+temporary+g' {} +
