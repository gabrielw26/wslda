# W-SLDA Toolkit
Self-consistent solver of mathematical problems which have structure formally equivalent to Bogoliubov-de Gennes equations.

## [Official webpage](https://wslda.fizyka.pw.edu.pl/)

## Implemented functionalities
* Solvers: static equations and time dependent equations
* Dimensionality: 3D and 2D and 1D
* Functionals: BdG, (A)SLDA
* Zero temperature and finite temperature
* Spin balanced and spin imbalanced systems
* User defined external potential, external pairing potential, external velocity field
* Broyden mixing
* Quantum friction
* Integration with VisIt tool

## Codes
* **td-wslda-3d** - code for solving time-dependent density functional equations in 3D Cartesian mesh.
* **td-wslda-2d** - constrained-pca code for solving time-dependent density functional equations in 3D Cartesian mesh with imposed constraint that system is uniform in _z_ direction
* **td-wslda-1d** - constrained-constrained-pca code for solving time-dependent density functional equations in 3D Cartesian mesh with imposed constraint that system is uniform in _z_  and _y_ directions.
* **st-wslda-3d** - code for solving static density functional equations in 3D Cartesian mesh, generator of initial states for _td-wslda-3d_ code. 
* **st-wslda-2d** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ direction, generator of initial states for _td-wslda-3d_ and _td-wslda-2d_ codes. 
* **st-wslda-1d** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ and _y_ directions, generator of initial states for _td-wslda-3d_ and _td-wslda-2d_ and _td-wslda-1d_ codes. 

## Setting up and building
See instructions provided [here](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Installing%20the%20toolkit).

## Documentation
For documenation see [Wiki Pages](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/home)

## Developers
* Gabriel Wlazłowski, Warsaw University of Technology,<br/> Main developer
* Maciej Marchwiany, Interdisciplinary Centre for Mathematical and Computational Modelling (ICM), <br/> contribution to _td-wslda-3d_ code (2016-2018)
* Wojciech Pudełko, Warsaw University of Technology,<br/> Implementation of Broyden algorithm (engineer thesis, 2020)
