affect
======

Tools for processing results of a finite element analysis.

License
-------

MIT. See [license file](https://github.com/kdcopps/affect/blob/master/LICENSE).

Dependencies
------------

`affect` requires Python 3.6+ with additional package dependencies. The easiest way to get going
 is by installing and using the [Anaconda](https://www.continuum.io/downloads) python environment.

* [numpy](htttps://www.numpy.org)
* [numexpr](https://github.com/pydata/numexpr)
* [netcdf4](https://pypi.python.org/pypi/netCDF4)

The netcdf4 package has additional dependencies including the hdf4 and hdf5 libraries.
Other dependencies which may be required for the development of affect are listed below.



Development
-----------

### Using the makefile
 
For building, testing, and documentation, the easiest way to work is using the `makefile` 
in the root directory of the project.

```
make clean; make build
make test
make docs
[...etc]
```

### Note on build warnings

Some warnings are hidden during the `make build` process, notably the warning about the 
deprecated numpy API that always occur with `cython` builds. These warnings are filtered using 
the `remove_warning.py` script.


The following type of warning occurs when building on Linux with gcc: 

```
cc1plus: warning: command line option "-Wstrict-prototypes" is valid for Ada/C/ObjC but not for C++
```

This is normal, and a side effect of using setuptools and distutils with the 
`Extension(language='c++')` parameter. These warnings may also be hidden during `make build`. 

### Development Dependencies

`affect` is developed with Cython.

* cython >= 0.24.1
* pytest >= 3.0.5
* pytest-cov >= 2.3.1
* pytest-runner >= 2.11.1
* sphinx >= 1.5.1
* sphinx_rtd_theme >= 0.1.9
* make

For building a conda recipe using `conda build` you may also need:

* jinja2

### Conda Recipe

To build the conda recipe, make sure you install the prerequisites:

`conda install conda-build`

The file `meta.yaml` in the root directory controls the building of the conda package.
Build the recipe using

`conda build <project-directory>`