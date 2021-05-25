Developer Notes
===============

This is a basic implementation of two classes: `wdata.io.Var` and
`wdata.io.WData`.  These implement the interfaces `wdata.io.IVar` and
`wdata.io.IWdata` respectively.  These provide a python interface to
the W-DATA format.

Thu Dec 17 2020
---------------
* Updated project to use [`poetry`] and [`nox`] for testing.  [Install
  poetry as](https://python-poetry.org/docs/#installation):
  
  ```bash
  curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
  ```

  Useful commands include:
  
  ```bash
  poetry run nox         # Run all tests
  poetry add [-D] <dep>  # Add [dev] dependence to pyproject.toml
  poetry remove <dep>    # Remove dependence to pyproject.toml
  poetry run nox         # Run all tests: requires python3.x installed
  ```
  
  To test against different versions of python, you will need to have
  those installed externally with either `pyenv`, `conda`, etc.  I did
  this:
  
  ```bash
  BINDIR=~/.local/bin/
  mkdir -p "${BINDIR}"
  for v in 3.6 3.7 3.8 3.9; do
    conda create -y -n py${v} python=${v}
    ln -fs $(conda run -n py${v} type -p python${v}) "${BINDIR}/"
  done
  conda clean -y --all    # Removes all downloaded files for conda.
  ```
  

[`poetry`]: https://python-poetry.org/
[`nox`]: https://nox.thea.codes/

Wed 23 Sept 2020
----------------
* Basic implementation complete with some tests.
* To test you should be able to run `pytest` or `nox`.  The latter
  will create an isolated environment.  Note: you will have to install
  the files first to use `pytest`:
  
      pip install --user -e .
      pytest

* Need some better documentation and more complete test coverage
  (especially all corner cases and error messages).
* I uploaded this to PyPi to reserver the `wdata` library name.

  https://pypi.org/project/wdata/0.1.0/

To Do:
* Implement WData.load to load a valid WData file when no infofile is
  specified.  The idea is that a folder with the following should be
  able to be loaded without any metadata:
  
      prefix_x.npy
      prefix_y.npy
      prefix_z.npy
      prefix_t.npy
      prefix_density.npy
      prefix_current.npy
      ...
      
  This should work with `io.WData.load('data_dir/prefix')` which would
  create a WData object with abscissa defined from the `x`, `y`, `z`
  and `t` files.  This will only work with `.npy` files which have the
  appropriate shape data inside.  Saving this data would then write
  the `prefix.wtxt` metadata file.
* Add methods to simplify adding data, extending data, and merging
  different datasets.
* Implement a C99 library for reading these.  Right now, we can use
  this library to convert everything to the `wdat` format (need a
  script for this in `bin/` and testing), but it would be better if we
  can embed Python in the C99 program so we can then use `numpy` to
  load data. 
* Implement a VisIt plugin and submit for inclusion with VisIt so we
  can have access to python data.  If we embed Python, then this
  should make life much easier, but I need to work with someone with
  this because I don't know VisIt well enough.
* Add proper version numbers, etc.
