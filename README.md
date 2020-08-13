# W-SLDA Toolkit

### Implemented functionalities
* Solvers: static equations and time dependent equations
* Dimensionality: 3D and 2D and 1D
* Functionals: BdG, (A)SLDA
* Zero temperature and finite temperature
* Spin balanced and spin imbalanced systems
* User defined external potential, external pairing potential, external velocity field
* Broyden mixing
* Quantum friction
* Integration with VisIt application

### Codes
* **td-wslda-3d** - code for solving time-dependent density functional equations in 3D Cartesian mesh. (status: **up to date**)
* **td-wslda-2d** - constrained-pca code for solving time-dependent density functional equations in 3D Cartesian mesh with imposed constraint that system is uniform in _z_ direction (status: **up to date**)
* **td-wslda-1d** - constrained-constrained-pca code for solving time-dependent density functional equations in 3D Cartesian mesh with imposed constraint that system is uniform in _z_  and _y_ directions. (status: **up to date**)
* **st-wslda-3d** - code for solving static density functional equations in 3D Cartesian mesh, generator of initial states for _pca_ code. (status: **up to date**)
* **st-wslda-2d** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ direction, generator of initial states for _pca_ and _cpca_ codes. (status: **up to date**)
* **st-wslda-1d** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ and _y_ directions, generator of initial states for _pca_ and _cpca_ and _ccpca_ codes. (status: **under construction**)


### Setting up and building
See instructios provided [here](http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/setting-up-calculations-and-compilation).  

### Wiki
For more details go to [Wiki Pages](http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/home)

### Authors
* Gabriel Wlazłowski, Warsaw University of Technology,<br/> Main developer
* Maciej Marchwiany, Interdisciplinary Centre for Mathematical and Computational Modelling (ICM), <br/> contribution to _pca_ code (2016-2018)
* Wojciech Pudełko, Warsaw University of Technology,<br/> Implmentation of Broyden algorithm (engineer thesis, 2020)