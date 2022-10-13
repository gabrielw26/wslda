# This script can be used Find&Replace opration applied to many files

# find -name "what" -exec sed -i 's+find+replace+g' {} +

# find -name "Makefile*" -exec sed -i 's+$(WSLDADIR)hpc-engine+$(WSLDADIR)/hpc-engine+g' {} +
# find -name "Makefile*" -exec sed -i 's+(must end with /)++g' {} +
# find -name "Makefile*" -exec sed -i 's+$(WSLDA)/+$(WSLDA)+g' {} +

find -name "*.cu" -exec sed -i 's+SLDAE_FORCE_A1+SLDA_FORCE_A1+g' {} +
find -name "*.c" -exec sed -i 's+SLDAE_FORCE_A1+SLDA_FORCE_A1+g' {} +
find -name "*.h" -exec sed -i 's+SLDAE_FORCE_A1+SLDA_FORCE_A1+g' {} +
# find -name "input.txt" -exec sed -i 's+hff_mu+hkf_mu+g' {} +

