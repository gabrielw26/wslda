# 2022.xx.xx
* Optimization of the HIP version of the code.
* The td codes can boost the performance by exploiting GPU-aware MPI.
* Performance improvement of SLDAE functional.
* New API_VERSION=20221120.
* Minor bug fixes related to the stability, improved database of test and templates.
* User can select integration scheme for time-dependent problems via predefines.h
* Improved parser in the input file

# 2022.08.27
* Porting code to [HIP](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Converting%20the%20toolkit%20to%20HIP).
* New library [W-Box](https://gitlab.fizyka.pw.edu.pl/wtools/wbox) has been incorporated into the W-SLDA.
* Minor bug fixes related to the stability of the code.

# 2022.06.27
* Speed-up of diagonalization process with ELPA for spin-symmetric mode.
* Minor bug fixes.
* Added Extended Thomas Fermi Model [experimental].

# 2022.04.05
* Update of SLDAE functional (redefinition of the effective mass, now A=alpha; added a new flag that forces the effective mass to be 1.0).
* Redefinition of delta_ext, now it is consistent with https://arxiv.org/abs/2201.07626.
* Automatic adjustment of energy cut-off has been improved.
* Improvement: code can now compute in HF mode (delta=0) and fully polarized gas (Nb=0).
* Added tool for generating [custom SLDAE functional](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Custom-SLDAE-functional).
* Added [Monitoring of conservation laws](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Monitoring-of-conservation-laws).

# 2022.01.05
* A new functional SLDAE has been implemented.
* Reorganization of templates.
* New user-defined file has been added: machine.h.
* New functionality: Time and position-dependent scattering length has been implemented.
* New functionality: Periodic checkpoints in td codes.
* New functionality: [Tracking of selected states](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Tracking%20of%20selected%20states)
* Framework for defining CUSTOMEDF in td codes has been added. 
* Improved performance of td codes (see *Implementing parallel reduction for multiple arrays* in Snippets). 
* Various minor updates and bug fixes. 

# 2021.09.05 
* [w-deriv](https://gitlab.fizyka.pw.edu.pl/wtools/wderiv) lib was included into W-SLDA. 
* [w-interp](https://gitlab.fizyka.pw.edu.pl/wtools/winterp) lib was included into W-SLDA. 
* Update of [w-data](https://gitlab.fizyka.pw.edu.pl/wtools/wdata) format to version 0.2.0.
* New tools were added: wdata2checkpoint, [Helmholtz decomposition code](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Helmholtz%20decomposition%20code)
* Reading of wave-functions via `td` codes was improved.
* Documentation of input files was improved. 

# 2021.05.25 
* Improved [Results reproducibility](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Results%20reproducibility).
* Improved error reporting system. 
* Better integration with [ELPA Library](https://elpa.mpcdf.mpg.de/).
* Simplified compilation process. 
* New feature: Automated testing system ([testsuite](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Testsuite))
* Added computation of entropy in st codes. 
* Improved stability of Broyden algorithm in st codes.

# 2021.03.15
* Improved performance of 1D codes
* Support for [regularization schemes of the pairing field](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Regularization%20schemes%20of%20the%20pairing%20field)
* New feature: [automatic interpolations](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Automatic%20interpolations) 
* Minor improvements and bug fixes 

# 2021.02.09 
* Added Extensions: small codes supporting computation process and data analysis.
* Improvements of data plugin for VisIt: it provides now full information about computation settings in File Information window.
* Better system of error reporting.
* Better control of uniform solver.

# 2020.12.28
* BdG functional has been extended to Fermi-Fermi mixtures
* Better support for mass-imbalanced systems
* Minor improvements: change of default values for some input parameters

# 2020.12.14  
* Initial release of W-SLDA Toolkit
