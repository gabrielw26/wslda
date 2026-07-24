# W-SLDA Toolkit
## High-Performance Simulation Platform for Fermionic Superfluids

The **W-SLDA Toolkit** is a large-scale, high-performance computing (HPC) framework for solving mathematical problems formally equivalent to the Bogoliubov–de Gennes (BdG) equations and density functional theories for superfluid fermionic systems.

It enables self-consistent static and time-dependent simulations of:

- ultracold atomic Fermi gases,
- superconducting systems,
- atomtronical devices,
- strongly interacting fermionic matter,
- vortex dynamics and quantum turbulence,
- Josephson effects,
- impurity dynamics,
- spin- and mass-imbalanced superfluids.

The toolkit is optimized for execution on modern HPC systems and is designed for simulations involving thousands of particles on large 3D grids.

## Scientific Scope

W-SLDA provides a unified numerical framework for:

- Bogoliubov–de Gennes (BdG) equations
- Superfluid Local Density Approximation (SLDA)
- Asymmetric SLDA (ASLDA)
- SLDA with extended terms (SLDAE)
- Custom user-defined energy density functionals

Both ground-state (static) and real-time (time-dependent) dynamics are supported in:

- 1D geometries,
- 2D geometries,
- full 3D Cartesian grids.

The code is particularly suited for studying emergent collective phenomena in strongly correlated fermionic systems.

## Architecture Philosophy

W-SLDA follows a **compile-time specialization model** to maximize performance.

Users define the physical problem in C (external potentials, pairing fields, constraints, etc.) and build a dedicated executable optimized for the specific geometry and lattice size.

This approach:

- maximizes numerical performance,
- enables strong compiler optimizations,
- ensures scalability on HPC systems,
- promotes reproducible workflows.

The toolkit is designed as an open, extensible research platform for large-scale density-functional simulations of fermionic superfluids.

## Reference
G. Wlazłowski, P. Magierski, M. M. Forbes, A. Bulgac,  
_W-SLDA Toolkit: A simulation platform for ultracold Fermi gases_,  
[[arXiv:2602.08982](https://arxiv.org/abs/2602.08982)].

If you use the W-SLDA Toolkit in your research, please cite the above work.

## [Official webpage](https://wslda.fizyka.pw.edu.pl/)

## Repositories:
* Main repository (authoritative source): [GitLab @ WUT](https://gitlab.fizyka.pw.edu.pl/wtools/wslda)
* Mirrors: [GitLab](https://gitlab.com/coldatoms/wslda), [GitHub](https://github.com/gabrielw26/wslda)  

## Implemented Functionalities

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
* **Template-based usage model**
* **Results reproducibility**: built-in mechanisms for reproducible simulations
* **Data format**: support for the [W-data format](https://gitlab.fizyka.pw.edu.pl/wtools/wdata) for storing simulation results

## Main Codes

Static solvers:
- st-wslda-1d
- st-wslda-2d
- st-wslda-3d

Time-dependent solvers:
- td-wslda-1d
- td-wslda-2d
- td-wslda-3d

For a detailed description of the available code types, see
[Types of codes](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Types-of-codes).

## System Requirements

W-SLDA is designed for HPC environments.

Minimum requirements:

- C compiler (GCC, Intel, or equivalent)
- MPI implementation (OpenMPI, MPICH, or vendor-specific)
- FFTW
- ScaLAPACK
- BLAS/LAPACK
- Optional: ELPA (recommended for large-scale diagonalizations)
- Optional: CUDA or HIP compiler (for time-dependent calculations)

Recommended:

- Multi-node cluster environment
- High-memory nodes for 3D simulations
- GPU-enabled systems (if using GPU-enabled builds)

Installation and build instructions are available
[here](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Installing-the-toolkit).

## Quick Start

Example: 1D static calculation.

1. Set environment variable:  
   `export WSLDA=/path/to/wslda`

2. Create a working directory from a template:   
   `cp -r $WSLDA/st-project-template my-project-name`  

   Then enter the project directory:  
   `cd my-project-name`

3. Edit the problem-definition files to configure your system:
   - `predefines.h`  
   - `problem-definition.h`  
   - `logger.h`
   - `input.txt` (edit before running the simulation)

4. Load required modules (example):  
   `source env.sh`

5. Compile:  
   `make 1d`

6. Run:  
   `mpirun -np 4 ./st-wslda-1d input.txt`

An example demonstrating the full workflow (Josephson junction dynamics) from configuration to visualization is available [here](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/Example-Josephson-junction). 

## Documentation

For full documentation, usage examples, and tutorials, see the
[Wiki Pages](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/wikis/home).  
Optionally, you can use wiki pages from mirror repositories: [GitLab wikis](https://gitlab.com/coldatoms/wslda/-/wikis/home), [GitHub wikis](https://github.com/gabrielw26/wslda/wiki).

## Reproducibility and Data

W-SLDA supports structured output via the W-data format.

Reproducibility packs and benchmark examples are provided through the repository and official webpage.

## License

W-SLDA Toolkit is distributed under the GNU General Public License v3 (GPLv3).

See the [COPYING](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/blob/public/COPYING) file for details.

## Developers

For the complete list of developers and contributors, see
[Contributors](https://wslda.fizyka.pw.edu.pl/index.php/Info2/Contributors).

## Bug Reporting

To report bugs or issues:
* use the [Issue tracking system](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/issues)
* or contact the development team via email: **wslda@fizyka.pw.edu.pl**

## Contributing
If you would like to contribute to the W-SLDA Toolkit, please see the [CONTRIBUTING.md](https://gitlab.fizyka.pw.edu.pl/wtools/wslda/-/blob/public/CONTRIBUTING.md) file for guidelines and further information.
