affect
======

Read/Write and postprocess results of finite element analysis.


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

### Other Dependencies

* sphinx >= 1.5.1
* sphinx_rtd_theme >= 0.1.9

### Conda Builds

To build the conda package with Conda build, make sure you install the prerequisites:

`conda install conda-build`

The file `meta.yaml` in the root directory controls the building of the conda package.