# 2022.xx.xx
* Update of SLDAE functional (redefinition of the effective mass, now A=alpha; added a new flag that forces the effective mass to be 1.0)
* Redefinition of delta_ext, now it is consistent with https://arxiv.org/abs/2201.07626
* Automatic adjustment of energy cut-off has been improved

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
