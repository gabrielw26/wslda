# W-SLDA Toolkit
The **W-SLDA Toolkit** is a self-consistent solver for mathematical problems whose structure is formally equivalent to the Bogoliubov–de Gennes equations.

The toolkit enables simulations of fermionic superfluids, with a primary focus on ultracold atomic gases. In addition, it can be used to model superconducting systems by solving static and time-dependent Bogoliubov–de Gennes equations within a unified numerical framework. Both static (ground-state) and time-dependent phenomena can be studied. The software is designed and optimized for large-scale simulations involving systems with thousands of particles and is suitable for execution on modern high-performance computing platforms.

## [Official webpage](https://wslda.fizyka.pw.edu.pl/)

## Implemented functionalities

* **Solvers**: static and time-dependent equations
* **Dimensionality**: 1D, 2D, and 3D geometries
* **Energy density functionals**: BdG, (A)SLDA, SLDAE, and user-defined custom functionals
* **Temperature regimes**: zero-temperature and finite-temperature calculations
* **Spin configurations**: spin-balanced and spin-imbalanced systems
* **Mass imbalance**: support for unequal particle masses
* **External fields**: user-defined external potentials, external pairing fields, external velocity fields, and time- and space-dependent scattering lengths
* **Visualization**: integration with the VisIt visualization tool
* **Extensions**: templates for data analysis workflows and a Python library for post-processing and data manipulation
* **Accelerated convergence**: Broyden mixing, automatic interpolations, and quantum friction techniques
* **Templates-based usage model**
* **Results reproducibility**: built-in mechanisms for reproducible simulations
* **Data format**: support for the [W-data format](https://gitlab.fizyka.pw.edu.pl/wtools/wdata) for storing simulation results

## Codes

* **td-wslda-3d** – solver for time-dependent density functional equations on a full 3D Cartesian mesh
* **td-wslda-2d** – solver for time-dependent density functional equations on a 3D Cartesian mesh with translational invariance imposed along the _z_ direction
* **td-wslda-1d** – solver for time-dependent density functional equations on a 3D Cartesian mesh with translational invariance imposed along the _y_ and _z_ directions
* **st-wslda-3d** – solver for static density functional equations on a full 3D Cartesian mesh; generates initial states for `td-wslda-3d`
* **st-wslda-2d** – solver for static density functional equations with translational invariance along the _z_ direction; generates initial states for `td-wslda-3d` and `td-wslda-2d`
* **st-wslda-1d** – solver for static density functional equations with translational invariance along the _y_ and _z_ directions; generates initial states for `td-wslda-3d`, `td-wslda-2d`, and `td-wslda-1d`

For a detailed description of the available code types, see
[Types of codes](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Types%20of%20codes).

## Setting up and building

Installation and build instructions are provided
[here](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/Installing%20the%20toolkit).

## Documentation

For full documentation, usage examples, and tutorials, see the
[Wiki Pages](https://gitlab.fizyka.pw.edu.pl/gabrielw/wslda/-/wikis/home).

## Developers

For the complete list of developers and contributors, see
[Contributors](https://wslda.fizyka.pw.edu.pl/index.php/Info2/Contributors).

## Bug reporting

To report bugs or issues:
* use the [Issue tracking system](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/issues)
* or contact the development team via email: **wslda@fizyka.pw.edu.pl**

## Contributing
If you would like to contribute to the W-SLDA Toolkit, please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines and further information.
