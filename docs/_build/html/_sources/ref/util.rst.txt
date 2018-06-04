======================================
``util`` --- Utilities
======================================

:mod:`affect.util`
=====================

.. automodule:: affect.util

.. currentmodule:: affect.util

Usage
-----

This module provides for creating and performing operations with aligned and 
compressed arrays.

.. autosummary::

   CompressedArray
   empty_aligned
   zeros_aligned
   byte_align
   is_byte_aligned
   take
   compress
   decompress

It also provides some basic functions for debugging and testing.

.. autosummary::

   arrays_share_data
   get_array_base
   ConsoleCode
   print_blue
   print_bold
   print_green
   print_yellow
   print_function_starting
   print_array_info

Data alignment
--------------

Data alignment for arrays means putting the data at a memory address equal to
some multiple of the word size. This is done to increase efficiency of data 
loads and stores to and from the processor. Processors are designed to 
efficiently move data to and from memory addresses that are on specific byte 
boundaries.

In addition to creating the data on aligned boundaries (that 
aligns the base pointer), the compiler is able to make optimizations when 
the data access (including base-pointer plus index) is known to be aligned 
by 64 bytes. Special `SIMD <https://en.wikipedia.org/wiki/SIMD>`_
instructions can be utilized by the compiler for
certain platforms. For example, the compiler/platform may support the special
instructions on processors such as the Intel®
`AVX-512 <https://en.wikipedia.org/wiki/AVX-512>`_ instructions, which
enables processing of twice the number of data elements that
`AVX/AVX2 <https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#Advanced_Vector_Extensions_2>`_
can process with a single instruction and four times that of
`SSE <https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions>`_.

By default, the compiler cannot know nor assume data is aligned
inside loops without some help from the programmer. Thus, you must also inform 
the compiler of this alignment via a combination of pragmas (C/C++) or 
keywords, clauses, or attributes so that compilers can generate optimal code. 

For the
Intel® Many Integrated Core Architecture such (`Intel® Xeon Phi™ Coprocessor
<https://en.wikipedia.org/wiki/Xeon_Phi>`_), memory movement is optimal when
the data starting address lies on 64 byte boundaries. Thus, by default, at
least at the time of this writing it is optimal to create data objects with
starting addresses that are modulo 64 bytes. For slightly less ambitious
modern architectures, such as Intel® Skylake, 32 byte aligned addresses may be
recommended.


Aligned arrays
--------------

These functions create and perform other operations on
:class:`numpy.ndarray` objects. All arrays created by calls to `affect`,
and those used internally in `affect`, are aligned. The two main
functions used to created aligned arrays are :func:`.empty_aligned` and
:func:`.zeros_aligned` that behave similarly to :func:`numpy.empty` and
:func:`numpy.zeros`, respectively.

For now this module defaults to using a 64 byte boundary. To align the data
of :class:`numpy.ndarray` to the word boundaries, during allocation it may
be necessary to insert some unused bytes at the start of the block, this is
data padding.

.. autofunction:: empty_aligned
.. autofunction:: zeros_aligned
.. autofunction:: byte_align
.. autofunction:: is_byte_aligned
.. autofunction:: take


Array compression
-----------------

The method of array compression is multithreaded and fast and can usually
compress an array of integers to around a fourth of the original size. It does
not compress arrays of floating point values as efficiently.

.. autoclass:: CompressedArray
   :members:
   :inherited-members:
   :show-inheritance:

.. autofunction:: compress
.. autofunction:: decompress


Printing, Debugging, Testing
----------------------------

Most of these functions are used internally for testing, but you may find them
of value for regular use.

.. autoclass:: ConsoleCode
   :members:
   :inherited-members:
   :show-inheritance:

   .. autoattribute:: PURPLE
   .. autoattribute:: CYAN
   .. autoattribute:: DARK_CYAN
   .. autoattribute:: BLUE
   .. autoattribute:: GREEN
   .. autoattribute:: YELLOW
   .. autoattribute:: RED
   .. autoattribute:: UNDERLINE
   .. autoattribute:: END

.. autofunction:: arrays_share_data
.. autofunction:: get_array_base
.. autofunction:: print_blue
.. autofunction:: print_bold
.. autofunction:: print_green
.. autofunction:: print_yellow
.. autofunction:: print_function_starting
.. autofunction:: print_array_info


Exceptions
----------

The exceptions thrown from this module are a part of the interface just as
much as the functions and classes. We define an Error root exception to allow
you to insulate yourself from this API. All the exceptions raised by this
module inherit from it.

.. autosummary::

   Error
   IllegalArgumentError
   UnsupportedArrayType


