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

Array Data and Internal Memory Buffers
--------------------------------------

Wherever possible, native array data from the ExodusII C API is accessed directly without copying through a view on a
:class:`numpy.ndarray`. This maintains performance by eliminating copying, and it supplies
the arrays in a convenient form for computations with numpy or scipy functions.

Some methods in this module, however, require allocating a smaller temporary memory buffers for working space.
These buffers are small and the size is noted in the documentation for each method.
Typical examples of the temporary memory buffer include functions required to
translate Exodus C strings to the Python str type, or rearrange integer or floating point arrays in the correct order
before supplying them or converting them to :class:`numpy.ndarray`. Internal temporary buffers are allocated on the
C/C++ heap using malloc or equivalent function. The buffer memory is released before the function returns.
If an exception occurs, we are careful that the buffer is still released.


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

   .. autoattribute:: NODAL
   .. autoattribute:: NODE_SET
   .. autoattribute:: EDGE_BLOCK
   .. autoattribute:: EDGE_SET
   .. autoattribute:: FACE_BLOCK
   .. autoattribute:: FACE_SET
   .. autoattribute:: ELEM_BLOCK
   .. autoattribute:: ELEM_SET
   .. autoattribute:: SIDE_SET
   .. autoattribute:: ELEM_MAP
   .. autoattribute:: NODE_MAP
   .. autoattribute:: EDGE_MAP
   .. autoattribute:: FACE_MAP
   .. autoattribute:: GLOBAL

Entry Types
-----------

.. autoclass:: EntryType
   :members:
   :inherited-members:
   :show-inheritance:

   .. autoattribute:: COORDINATE
   .. autoattribute:: ELEMENT
   .. autoattribute:: NODE
   .. autoattribute:: FACE
   .. autoattribute:: EDGE

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

.. autoclass:: Messages
   :members:
   :inherited-members:
   :show-inheritance:

.. autofunction:: debug_messages

.. autoclass:: DebugMessages
   :members:
   :inherited-members:
   :show-inheritance:

Global Entity
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

.. autodata:: Component
   :annotation: = indicates either a str or int

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

Local Block Connectivity
------------------------

.. autoclass:: LocalConnectivity
   :members:
   :inherited-members:
   :show-inheritance:

Exceptions
----------

.. autoexception:: Error
   :show-inheritance:

.. autoexception:: NoMemory
   :show-inheritance:

.. autoexception:: FileError
   :show-inheritance:

.. autoexception:: FileExists
   :show-inheritance:

.. autoexception:: FileNotFound
   :show-inheritance:

.. autoexception:: FileAccess
   :show-inheritance:

.. autoexception:: ReadWriteError
   :show-inheritance:

.. autoexception:: InternalError
   :show-inheritance:

.. autoexception:: InvalidEntityType
   :show-inheritance:

.. autoexception:: ArrayTypeError
   :show-inheritance:

.. autoexception:: ArgumentTypeError
   :show-inheritance:

.. autoexception:: InactiveComponent
   :show-inheritance:

.. autoexception:: InactiveEntity
   :show-inheritance:

.. autoexception:: InactiveField
   :show-inheritance:

.. autoexception:: InvalidSpatialDimension
   :show-inheritance:

.. autoexception:: NotYetImplemented
   :show-inheritance:

.. autoexception:: RangeError
   :show-inheritance:


Utility Functions
-----------------

.. autofunction:: library_version

.. autofunction:: debug_messages


