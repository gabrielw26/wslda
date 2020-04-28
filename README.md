# Polarized cold atoms - simulation toolkit

### Implemented functionalities
* Solvers: static equations and time dependent equations
* Dimensionality: 3D and 2D and 1D
* Functionals: BdG, (A)SLDA
* Zero temperature and finite temperature
* Spin balanced and spin imbalanced systems
* Computation in lab frame or rotating frame
* Quantum friction
* Integration with VisIt application

### Codes
* **pca** - code for solving time-dependent density functional equations in 3D Cartesian mesh. (status: **up to date**)
* **cpca** - constrained-pca code for solving time-dependent density functional equations in 3D Cartesian mesh with imposed constraint that system is uniform in _z_ direction (status: **up to date**)
* **ccpca** - constrained-constrained-pca code for solving time-dependent density functional equations in 3D Cartesian mesh with imposed constraint that system is uniform in _z_  and _y_ directions. (status: **up to date**)
* **s3dpca** - code for solving static density functional equations in 3D Cartesian mesh, generator of initial states for _pca_ code. (status: **up to date**)
* **s2dpca** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ direction, generator of initial states for _pca_ and _cpca_ codes. (status: **up to date**)
* **s1dpca** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ and _y_ directions, generator of initial states for _pca_ and _cpca_ and _ccpca_ codes. (status: **under construction**)

### Project structure 

* Directory **pca/** - folder with source files
* Directory **scripts/** - various useful scripts, mainly in python
* Directory **tex/** - files documenting formulas and numerical concepts 

### Building

In pca folder you will find make files of form Makefile._code_._machine_. In header of each make file are listed modules that you need to load before you compile the code. For example:

```bash
module load cuda
make -f Makefile.pca.summit
```

### Wiki
For more details go to [Wiki Pages](http://git2.if.pw.edu.pl/gabrielw/cold-atoms/wikis/home)

### Authors
* Gabriel Wlazłowski, Warsaw University of Technology<br/> Main developer
* Maciej Marchwiany, Interdisciplinary Centre for Mathematical and Computational Modelling (ICM), <br/> contribution to _pca_ code (2016-2018)
