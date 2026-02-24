# This script can be used Find&Replace opration applied to many files

# find -name "what" -exec sed -i 's+find+replace+g' {} +

# find -name "Makefile*" -exec sed -i 's+$(WSLDADIR)hpc-engine+$(WSLDADIR)/hpc-engine+g' {} +
# find -name "Makefile*" -exec sed -i 's+(must end with /)++g' {} +
# find -name "Makefile*" -exec sed -i 's+$(WSLDA)/+$(WSLDA)+g' {} +

# find -name "*.cu" -exec sed -i 's+amb+abm+g' {} +
# find -name "*.c" -exec sed -i 's+amb+abm+g' {} +
# find -name "*.h" -exec sed -i 's+amb+abm+g' {} +

find -name "*.sh" -exec sed -i 's+465001656+465002810+g' {} +
find -name "Makefile*" -exec sed -i 's+465001656+465002810+g' {} +
find -name "README*" -exec sed -i 's+465001656+465002810+g' {} +

