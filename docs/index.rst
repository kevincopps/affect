.. affect documentation master file, created by
   sphinx-quickstart on Wed Dec 17 12:15:57 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Affect software package
===========================

.. toctree::
   :hidden:
   :glob:
   :maxdepth: 4

   ref/index

Exodus databases
================

Analysis data from unstructured finite element or finite volume models are accessed through the *Exodus II* library. The
:mod:`affect.exodus` module provides a Python interface for input and output of Exodus II databases.

    Exodus II is a model developed to store and retrieve data for finite element analyses. It is used for preprocessing
    (problem definition), postprocessing (results visualization), as well as code to code data transfer.
    An Exodus II data file is a random access, machine independent, binary file. The ExodusII file format and API is
    based on the NetCDF and HDF5 formats and API's, respectively.

    --- `EXODUS II: A Finite Element Data Model <http://gsjaardema.github.io/seacas/exodusII-new.pdf>`_, Gregory D.
    Sjaardema, et al. (Documentation for ExodusII database files, including the C/C++ and FORTRAN API.)

The :mod:`affect.exodus` module maintains compact representation of the array data accessed by
direct access through `Numpy array objects <https://docs.scipy.org/doc/numpy/reference/arrays.html>`_.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
   :maxdepth: 3

   ref/glossary
