# W-data format

W-data format is designed to store and manipulate data defined on a spatial lattice.  
This format was originally derived from the [W-SLDA](https://gitlab.fizyka.pw.edu.pl/wtools/wslda) project.
The format is designed to be conceptually easy to understand, which turns out to be an important issue for academic applications where students are typically involved.  
Presently the format support lattices in 1D, 2D, and 3D. Real, complex, and vector variables can be stored. To learn more, see [Documenation](./doc/REDADME.md) in doc folder.

# C library compilation
The package contains C (standard C99) library and tools supporting data manipulation.  
Edit header of [Makefile](./Makefile) and execute
```bash
make
```
*Note*: to compile the tools you will need to install [W-interp library](https://gitlab.fizyka.pw.edu.pl/wtools/winterp) first. 
To compile only the library execute  
```bash
make lib
```

The C libary is stored in `libwdata.a` (static lib) and `libwdata.so` (dynamically linked lib) files.  

Once the library is compiled it is recommended to set system variable (preferably in `.bashrc`)
```bash
export WDATA=set-path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$WDATA
export PATH=$PATH:$WDATA/bin
```

# C library usage
To use the library you need to add in your C code header file
```c
#include "wdata.h"
```

To compile the code add to the compilation command
```bash
gcc ... -I$WDATA/c/ -L$WDATA -lwdata
```

See [here](https://gitlab.fizyka.pw.edu.pl/wtools/wdata/-/wikis/Examples/C-examples) for examples of usage.

# Python library
It is clone of repository:
```bash
git subtree pull --prefix python https://github.com/forbes-group/wdata branch/default --squash
```  

Follow instructions from [this repo](https://hg.iscimath.org/forbes-group/wdata) to learn how to use this module. 
The package can be found [here](https://pypi.org/project/wdata/).

# Developers
* Gabriel Wlazłowski, Warsaw University of Technology
* Andrzej Makowski, Warsaw University of Technology
* Michael McNeil Forbes, Washington State University

