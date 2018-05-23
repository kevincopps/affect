=================================================
``math`` --- Calculations on grid and mesh fields
=================================================

:mod:`affect.arithmetic`
=====================

.. automodule:: affect.arithmetic

.. currentmodule:: affect.arithmetic

Usage
-----

This module contains functions for making basic calculations on field values
associated with the nodes and cells of a mesh. For example, you may
want to compute the average of a node field in every cell.

It is not meant to contain every possible calculation that one would need in
practice. What distinguishes this module is that it includes highly optimized
and threaded functions for a few common operations.

.. autosummary::

   average_element_node_values

Functions
---------

.. autofunction:: average_element_node_values


Exceptions
----------

The exceptions thrown from this module are a part of the interface just as much as the functions and classes.
We define an Error root exception to allow you to insulate yourself from this API. All the exceptions
raised by this module inherit from it.

.. autoexception:: Error
   :show-inheritance:

.. autoexception:: MismatchedArrayShape
   :show-inheritance:

.. autoexception:: IllegalArgumentError
   :show-inheritance:

.. autoexception:: UnalignedArray
   :show-inheritance:

.. autoexception:: UnsupportedArrayType
   :show-inheritance:
