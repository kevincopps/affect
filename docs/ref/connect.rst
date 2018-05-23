======================================
``connect`` --- Connectivity Utilities
======================================

:mod:`affect.connect`
=====================

.. automodule:: affect.connect

.. currentmodule:: affect.connect

Usage
-----

Utilities for processing connectivity information. Often the only connectivity stored on disk is the element (cell) to
vertex, or element to node connectivity. The functions in this module determine element neighbor information and
the vertex to element connectivity.

.. autosummary::

   element_to_element
   vertex_to_element
   boundary_face_to_vertex
   convert_to_local_connectivity

Functions
---------

.. autofunction:: element_to_element
.. autofunction:: vertex_to_element
.. autofunction:: boundary_face_to_vertex
.. autofunction:: convert_to_local_connectivity

Cell and Element Topology
-------------------------

Standard cell to vertex connectivities are identified using the CellTopology enum.

.. autoclass:: CellTopology
   :members:
   :inherited-members:
   :show-inheritance:

   .. automethod:: max_vertex_per_face
   .. automethod:: max_node_per_face
   .. automethod:: num_edge
   .. automethod:: num_face
   .. automethod:: num_node
   .. automethod:: num_node_per_edge
   .. automethod:: num_vertex
   .. automethod:: spatial_dimension

   .. autoattribute:: HEX8
   .. autoattribute:: HEX
   .. autoattribute:: HEX20
   .. autoattribute:: HEX27
   .. autoattribute:: PYRAMID5
   .. autoattribute:: PYRAMID
   .. autoattribute:: PYRAMID13
   .. autoattribute:: QUAD4
   .. autoattribute:: QUAD
   .. autoattribute:: QUADRILATERAL4
   .. autoattribute:: QUADRILATERAL
   .. autoattribute:: QUAD8
   .. autoattribute:: QUADRILATERAL8
   .. autoattribute:: QUAD9
   .. autoattribute:: QUADRILATERAL9
   .. autoattribute:: QUADSHELL4
   .. autoattribute:: SHELL4
   .. autoattribute:: SHELL
   .. autoattribute:: QUADSHELL
   .. autoattribute:: QUADSHELL8
   .. autoattribute:: SHELL8
   .. autoattribute:: QUADSHELL9
   .. autoattribute:: SHELL9
   .. autoattribute:: TET4
   .. autoattribute:: TETRA4
   .. autoattribute:: TETRA
   .. autoattribute:: TET8
   .. autoattribute:: TETRA8
   .. autoattribute:: TET10
   .. autoattribute:: TETRA10
   .. autoattribute:: TRI3
   .. autoattribute:: TRIANGLE3
   .. autoattribute:: TRIANGLE
   .. autoattribute:: TRI
   .. autoattribute:: TRI4
   .. autoattribute:: TRIANGLE4
   .. autoattribute:: TRI6
   .. autoattribute:: TRIANGLE6
   .. autoattribute:: TRISHELL3
   .. autoattribute:: TRISHELL
   .. autoattribute:: TRIANGLESHELL3
   .. autoattribute:: TRIANGLESHELL
   .. autoattribute:: TRISHELL6
   .. autoattribute:: TRIANGLESHELL6
   .. autoattribute:: WEDGE6
   .. autoattribute:: WEDGE
   .. autoattribute:: WEDGE15
   .. autoattribute:: END_TOPOLOGY

.. autodata:: CELL_TOPOLOGY_ALL

.. autodata:: CELL_TOPOLOGY_2D

.. autodata:: CELL_TOPOLOGY_3D


Exceptions
----------

The exceptions thrown from this module are a part of the interface just as much as the functions and classes.
We define an Error root exception to allow you to insulate yourself from this API. All the exceptions
raised by this module inherit from it.

.. autoexception:: Error
   :show-inheritance:

.. autoexception:: InvalidFaceCount
   :show-inheritance:

.. autoexception:: UnknownCellTopology
   :show-inheritance:

