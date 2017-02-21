========================================
``exodus`` --- I/O of ExodusII Databases
========================================

:mod:`affect.exodus`
====================

.. automodule:: affect.exodus

.. currentmodule:: affect.exodus

Usage and Data Model
--------------------

The fundamental class for I/O is :class:`.Database`. Though it can be directly
used, there is a convenient :class:`.DatabaseFile` context manager.

Access to the stored values and the structure of the :class:`.Database` is through objects of the following classes:

.. autosummary::

    DatabaseFile
    Database
    Global
    Field
    FieldArray
    Fields
    Nodal
    Blocks
    Block
    Maps
    Map
    Sets
    Set

See also the Exodus data model :doc:`glossary` for more information.

Internal Memory Buffers
-----------------------

Some functions in this module that retrieve or modify data in an ExodusII database require allocating
a temporary memory buffer for working space. The buffer memory is released before the function returns. A note near the
documentation for these specific functions is provided calling out the maximum size in bytes of the required buffer.

This memory buffer may be required in order to rearrange the strings, integer or floating point
arrays in the correct order before supplying them or converting them to Python types or :class:`numpy.ndarray`. This
buffer is allocated on the C/C++ heap using malloc or equivalent function.


Database Objects
----------------

.. autoclass:: DatabaseFile
   :members:

.. autoclass:: Database
   :members:

.. autoclass:: Mode
   :members:
   :inherited-members:
   :show-inheritance:

   .. autoattribute:: Mode.READ_ONLY
   .. autoattribute:: Mode.READ_WRITE
   .. autoattribute:: Mode.CREATE
   .. autoattribute:: Mode.REPLACE

.. autoclass:: EntityType
   :members:
   :inherited-members:
   :show-inheritance:

   .. autoattribute:: EntityType.NODAL
      :annotation:

   .. autoattribute:: EntityType.NODE_SET
      :annotation:

   .. autoattribute:: EntityType.EDGE_BLOCK
      :annotation:

   .. autoattribute:: EntityType.EDGE_SET
      :annotation:

   .. autoattribute:: EntityType.FACE_BLOCK
      :annotation:

   .. autoattribute:: EntityType.FACE_SET
      :annotation:

   .. autoattribute:: EntityType.ELEM_BLOCK
      :annotation:

   .. autoattribute:: EntityType.ELEM_SET
      :annotation:

   .. autoattribute:: EntityType.SIDE_SET
      :annotation:

   .. autoattribute:: EntityType.ELEM_MAP
      :annotation:

   .. autoattribute:: EntityType.NODE_MAP
      :annotation:

   .. autoattribute:: EntityType.EDGE_MAP
      :annotation:

   .. autoattribute:: EntityType.FACE_MAP
      :annotation:

   .. autoattribute:: EntityType.GLOBAL
      :annotation:


Exceptions and Debug Messages
-----------------------------

The exceptions thrown from this module are a part of the interface just as much as the functions and classes.
We define an Error root exception to allow you to insulate yourself from this API. All the exceptions
raised by this module inherit from it.

.. autosummary::

    Error
    NoMemory
    FileError
    FileExists
    FileNotFound
    FileAccess
    ReadWriteError
    InternalError
    InvalidEntityType
    ArrayTypeError
    ArgumentTypeError
    InactiveComponent
    InactiveEntity
    InactiveField
    InvalidSpatialDimension
    NotYetImplemented
    RangeError

.. seealso:: example code in :class:`.Error`

.. |Error| replace:: :class:`.Error`
.. |NoMemory| replace:: :class:`.NoMemory`
.. |FileError| replace:: :class:`.FileError`
.. |FileExists| replace:: :class:`.FileExists`
.. |FileNotFound| replace:: :class:`.FileNotFound`
.. |FileAccess| replace:: :class:`.FileAccess`
.. |ReadWriteError| replace:: :class:`.ReadWriteError`
.. |InternalError| replace:: :class:`.InternalError`
.. |InvalidEntityType| replace:: :class:`.InvalidEntityType`
.. |ArrayTypeError| replace:: :class:`.ArrayTypeError`
.. |ArgumentTypeError| replace:: :class:`.ArgumentTypeError`
.. |InactiveComponent| replace:: :class:`.InactiveComponent`
.. |InactiveEntity| replace:: :class:`.InactiveEntity`
.. |InactiveField| replace:: :class:`.InactiveField`
.. |InvalidSpatialDimension| replace:: :class:`.InvalidSpatialDimension`
.. |NotYetImplemented| replace:: :class:`.NotYetImplemented`
.. |RangeError| replace:: :class:`.RangeError`

.. autofunction:: debug_messages

.. autoclass:: DebugMessages
   :members:
   :inherited-members:
   :show-inheritance:

Entry Types
-----------

.. autodata:: ELEMENT
   :annotation:

.. autodata:: NODE
   :annotation:

.. autodata:: FACE
   :annotation:

.. autodata:: EDGE
   :annotation:


Global Entitity
---------------

.. autoclass:: Global
   :members:
   :inherited-members:
   :show-inheritance:

Nodal Entity
------------

.. autoclass:: Nodal
   :members:
   :inherited-members:
   :show-inheritance:

Blocks Entities
---------------

.. autoclass:: Blocks
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Block
   :members:
   :inherited-members:
   :show-inheritance:

Sets Entities
-------------

.. autoclass:: Sets
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Set
   :members:
   :inherited-members:
   :show-inheritance:

Maps Entities
-------------

.. autoclass:: Maps
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Map
   :members:
   :inherited-members:
   :show-inheritance:

Fields
------

.. autoclass:: Field
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: Fields
   :members:
   :inherited-members:
   :show-inheritance:

.. autoclass:: FieldArray
   :members:
   :show-inheritance:

Exceptions
----------

.. autoclass:: Error
   :show-inheritance:

.. autoclass:: NoMemory
   :show-inheritance:

.. autoclass:: FileError
   :show-inheritance:

.. autoclass:: FileExists
   :show-inheritance:

.. autoclass:: FileNotFound
   :show-inheritance:

.. autoclass:: FileAccess
   :show-inheritance:

.. autoclass:: ReadWriteError
   :show-inheritance:

.. autoclass:: InternalError
   :show-inheritance:

.. autoclass:: InvalidEntityType
   :show-inheritance:

.. autoclass:: ArrayTypeError
   :show-inheritance:

.. autoclass:: ArgumentTypeError
   :show-inheritance:

.. autoclass:: InactiveComponent
   :show-inheritance:

.. autoclass:: InactiveEntity
   :show-inheritance:

.. autoclass:: InactiveField
   :show-inheritance:

.. autoclass:: InvalidSpatialDimension
   :show-inheritance:

.. autoclass:: NotYetImplemented
   :show-inheritance:

.. autoclass:: RangeError
   :show-inheritance:
