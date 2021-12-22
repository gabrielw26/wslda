# W-SLDA Toolkit
Self-consistent solver of mathematical problems which have structure formally equivalent to Bogoliubov-de Gennes equations.  

The toolkit allows for simulating fermionic superfluids like ultracold atomic gases. Both static and time-depend phenomena can be investigated by means of W-SLDA. The software is optimized towards simulations of large systems, consisting of thousands of particles.   

## [Official webpage](https://wslda.fizyka.pw.edu.pl/)

## Implemented functionalities
* Solvers: static equations and time-dependent equations
* Dimensionality: 3D and 2D and 1D
* Functionals: BdG, (A)SLDA, SLDAE, custom functionals
* Zero temperature and finite temperature
* Spin balanced and spin imbalanced systems
* Mass imbalanced systems
* User-defined external potential, external pairing potential, external velocity field, time and position dependent scattering length
* Integration with visualization tool (VisIt)
* Extensions: templates for codes supporting data analysis, python lib for working with and manipulating data
* Speeding up of convergence: Broyden mixing, automatic interpolations, quantum friction
* Templates-based usage model.
* Results reporoducibility
* [W-data format](https://gitlab.fizyka.pw.edu.pl/wtools/wdata) for storing the results

## Codes
* **td-wslda-3d** - code for solving time-dependent density functional equations on 3D Cartesian mesh.
* **td-wslda-2d** - code for solving time-dependent density functional equations on 3D Cartesian mesh, with the imposed constraint that the system is uniform in _z_ direction
* **td-wslda-1d** - code for solving time-dependent density functional equations on 3D Cartesian mesh, with the imposed constraint that the system is uniform in _z_  and _y_ directions.
* **st-wslda-3d** - code for solving static density functional equations in 3D Cartesian mesh, generator of initial states for _td-wslda-3d_ code. 
* **st-wslda-2d** - code for solving static density functional equations in 3D Cartesian mesh, with the imposed constraint that system is uniform in _z_ direction, generator of initial states for _td-wslda-3d_ and _td-wslda-2d_ codes. 
* **st-wslda-1d** - code for solving static density functional equations in 3D Cartesian mesh, with imposed constraint that system is uniform in _z_ and _y_ directions, generator of initial states for _td-wslda-3d_ and _td-wslda-2d_ and _td-wslda-1d_ codes. 

For more info about the code types see [here](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Types%20of%20codes).

## Setting up and building
See instructions provided [here](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Installing%20the%20toolkit).

## Documentation
For documenation see [Wiki Pages](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/home)

## Developers
For complete list of developers and contributors see [here](https://wslda.fizyka.pw.edu.pl/index.php/Info2/Contributors).

## Bug reporting
To report bug:  
* use [Issues reporting system](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/issues)
* write email to wslda@fizyka.pw.edu.pl

