# This script can be used Find&Replace opration applied to many files

# find -name "what" -exec sed -i 's+find+replace+g' {} +

# find -name "Makefile*" -exec sed -i 's+$(WSLDADIR)hpc-engine+$(WSLDADIR)/hpc-engine+g' {} +
# find -name "Makefile*" -exec sed -i 's+(must end with /)++g' {} +
# find -name "Makefile*" -exec sed -i 's+$(WSLDA)/+$(WSLDA)+g' {} +

find -name "*.cu" -exec sed -i 's+hff_mu+hkf_mu+g' {} +
find -name "*.c" -exec sed -i 's+hff_mu+hkf_mu+g' {} +
find -name "*.h" -exec sed -i 's+hff_mu+hkf_mu+g' {} +
find -name "input.txt" -exec sed -i 's+hff_mu+hkf_mu+g' {} +

find -name "*.cu" -exec sed -i 's+hff_T+hkf_T+g' {} +
find -name "*.c" -exec sed -i 's+hff_T+hkf_T+g' {} +
find -name "*.h" -exec sed -i 's+hff_T+hkf_T+g' {} +
find -name "input.txt" -exec sed -i 's+hff_T+hkf_T+g' {} +
