affect
======

The affect Python package provides a API/library for processing meshes and field results of a finite element 
analysis.

License
-------

BSD. Specifically the "new/modified BSD" license, also called the "BSD three clause license".
See [license file](https://github.com/kdcopps/affect/blob/master/LICENSE).

Documentation
-------------

See the [affect documentation here](https://kevincopps.github.io/affect).

Basic Goals and Capabilities
----------------------------

* Designed for use on many-core shared memory workstations or servers using 
threaded parallelism.
* Easy to use within a Python (_Jupyter_) Notebook and targets modern Python
 (>3.5).
* Support for faster operations on very large analysis results on 
unstructured grids (billions of elements/nodes) using a compact storage scheme and CPU threaded parallelism. 
* Fast skinning of mesh/block boundaries.
* A clean Python Database class for reading mesh and fields from ExodusII files directly into Numpy arrays.
* Some optimized built in postprocessors, like the FRF (frequency response function).

Using the Library
----------------- 

### Dependencies

`affect` requires Python 3.6+ with additional package dependencies. The easiest way to get going
 is by installing and using the [Anaconda](https://www.continuum.io/downloads) python environment.

* [numpy](htttps://www.numpy.org)
* [numexpr](https://github.com/pydata/numexpr)
* [netcdf4](https://pypi.python.org/pypi/netCDF4)

The netcdf4 package has additional dependencies including the hdf4 and hdf5 libraries.
For a complete list of dependencies see the `requirements.txt` file.

### OpenMP and thread affinity

Setting OpenMP environment values for binding threads may improve performance. The matter of thread affinity becomes 
important on multi-socket nodes. Optimal settings of these environment variables will vary by platform.

Putting threads far apart, i.e. on different sockets

 * May improve the aggregated memory bandwidth available to your application
 * May improve the combined cache size available to your application
 * May decrease performance of synchronization constructs
 
Putting threads close together, i.e. on two adjacent cores which possibly shared some caches
 * May improve performance of synchronization constructs
 * May decrease the available memory bandwidth and cache size

See the OpenMP documentation for documentation on setting the two variables `OMP_PROC_BIND` and `OMP_PLACES`.
These can take on the following values:
 
One of `OMP_PROC_BIND={false|true|master|close|spread}` 

And one of `OMP_PLACES={sockets|cores|threads}` or an explicit list of places that are described by non-negative 
numbers.

For example,
```
export OMP_PLACES=threads
export OMP_PROC_BIND=close
```

Developing and Improving affect
-------------------------------

### Using the makefile
 
For building, testing, and documentation, the easiest way to work is using the `makefile` 
in the root directory of the project.

```
make clean; make build
make test
make docs
[...etc]
```

### Build warnings

Some warnings are hidden during the `make build` process, notably the warning about the 
deprecated numpy API that always occur with `cython` builds. These warnings are filtered using 
the `remove_warning.py` script.


The following type of warning occurs when building on Linux with gcc: 

```
cc1plus: warning: command line option "-Wstrict-prototypes" is valid for Ada/C/ObjC but not for C++
```

This is normal, and a side effect of using setuptools and distutils with the 
`Extension(language='c++')` parameter. These warnings may also be hidden during `make build`. 

### Setting up your environment

There are a few other dependencies required only for development purposes for 
building the code and documentation. `affect` is developed with Cython and the
documentation is built using sphinx. For building a conda recipe using 
`conda build` you may also need jinja2. 

* cython >= 0.24.1
* jinja2
* pytest >= 3.0.5
* pytest-cov >= 2.3.1
* pytest-runner >= 2.11.1
* sphinx >= 1.5.1
* sphinx_rtd_theme >= 0.1.9

Other command line tools required or that you may want to use on Mac or Linux 
are:

* ccache
* gcc6
* make

Using the Conda distribution is the easiest way to develop and manage the 
Python package dependencies. Each of these may be installed similar to:

```
conda install pytest-cov
```

If you find any of the dependencies are not in the default conda channel,
first try adding the conda-forge channel (you may also wish to disable SSL
verification behind firewalls or web proxies):  

```
conda config --set ssl_verify false
conda config --add channels conda-forge
```

If the package is still not found, fall back to using pip:

```
pip install sphinx-autodoc-typehints
```

### Conda Recipe

To build the conda recipe, make sure you install the prerequisites:

`conda install conda-build`

The file `meta.yaml` in the root directory controls the building of the conda package.
Build the recipe using

`conda build <project-directory>`

Immediate Plans
---------------

* Write meshes in the _ExodusII_ format (unstructured meshes)
* Read ExodusII meshes split across files, thus avoiding the need for 
executing the _Seacas_ `epu` tool. 
* Read and write _CGNS_ mesh format (structured meshes)
* Provide a local least squares recovery for mesh fields. This is a smoothing
algorithm that can be used to, for example, either promote a element field 
(piecewise constant) to a piecewise linear field, or provide a local error
indicator for spatial interpolation error on each element.

Future Plans
------------

* Fast visualization of mesh/block boundary.
* Divide mesh boundary (block skin) into separate non-intersecting 
surface regions at sharp edges. (This can be thought of as automatic 
creation of sidesets ExodusII terms).
* Fast search algorithm for nearest neighbor searches for element sides
and nodes on the boundary surface. 