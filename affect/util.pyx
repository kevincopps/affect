"""
This module contains the functions and classes for basic operations.

.. currentmodule:: affect.utils

.. |uint32_1d| replace:: :obj:`ndarray[uint32_t, ndim=1, mode="c"]`
"""

cimport cython
cimport numpy
import numpy
from numpy cimport ndarray

import atexit
import threading
import sys
from typing import Iterable, Optional, Tuple, Union
import weakref

from libc.stdlib cimport free
from libc.stdint cimport intptr_t, uintptr_t, int32_t, int64_t, uint32_t, uint64_t, uint8_t
from libcpp.pair cimport pair

# noinspection PyPackageRequirements
import blosc  # from module python-blosc

from . cimport cutil

blosc.set_releasegil(True)  # release the Python global inter-lock during c-blosc compress and decompress for speed

TEST_FUNCTION_PREFIX = 'test_'

class Error(Exception):
    """Base class for all exceptions raised by this module."""

class IllegalArgumentError(Error, ValueError):
    """Exception raised if argument does not make sense for the current operation."""

class UnsupportedArrayType(Error, TypeError):
    """Exception raised if dtype of numpy array argument is not supported."""


# Align for AVX
#
# 32 bit alignment is recommended for AVX ans AVX2 intrinsics
#
# CPUs with AVX:
#
#       Sandy Bridge processor or later, Q1 2011
#       See https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX
#
# CPUs with AVX2:
#       Haswell processor or later, Q2 2013
#
# 64 bit alignment is recommended for AVX512
#
# CPUs with AVX-512:
#       Xeon Phi x200 (aka Knights Landing, either as a host processor or a coprocessor, 2016)
#       Skylake EP/EX Xeon "Purley" (Xeon E5-26xx V5) processors (expected in H2 2017)
#       Cannonlake (expected in 2017)
#
# Compilers supporting AVX-512
#   GCC 4.9 and newer
#   Clang 3.9 and newer
#   ICC 15.0.1 and newer
#

# initialize for use of C API functions
numpy.import_array()

# The optimum SIMD alignment for CPU using AVX and AVX2 is 32,
# for code that is not aware of this module, make it 64 anyway
cdef int _simd_alignment = 64

SUPPORTED_NUM_BYTES = (16, 32, 64, 128)

if cutil.can_use_intel_knl_features():
    _simd_alignment = 64

cdef inline size_t get_size(object shape):
    # get the size from an array shape, whether the shape argument object is either a single integer or a sequence
    # this method also alleviates possible issue with numpy.prod wrapping around on 32-bits on Windows 64-bit
    cdef size_t array_size
    if isinstance(shape, (int, numpy.integer)):
        array_size = <size_t> shape
    else:
        array_size = 1
        for each_dimension in shape:
            array_size *= each_dimension
    return array_size


class ConsoleCode:
    """
    Colors, bold, and underline codes for wrapping around printed strings to the terminal.

    Also, see the convenience functions :func:`.print_blue` :func:`.print_bold` etc.

    Example:

        Print bold and yellow text to stdout for display on the terminal.

        >>> print(ConsoleCode.BOLD + 'Hello, I and displayed in bold text' + ConsoleCode.END)
        >>> print(ConsoleCode.YELLOW + 'Hello, I and displayed in yellow text' + ConsoleCode.END)
    """
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARK_CYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def print_blue(*args: Iterable[object]):
    """
    Print object strings to console with a blue color.

    Args:
        *args: iterable sequence of any objects with a string representation
    """
    print(ConsoleCode.BLUE + ' '.join(map(str, args)) + ConsoleCode.END)


def print_bold(*args):
    """
    Print object strings to console with the console's representation of bold.

    Args:
        *args: iterable sequence of any objects with a string representation
    """
    print(ConsoleCode.BOLD + ' '.join(map(str, args)) + ConsoleCode.END)


def print_green(*args: Iterable[object]):
    """
    Print object strings to console with a blue color.

    Args:
        *args: iterable sequence of any objects with a string representation
    """
    print(ConsoleCode.GREEN + ' '.join(map(str, args)) + ConsoleCode.END)


def print_yellow(*args: Iterable[object]):
    """
    Print object strings to console with a blue color.

    Args:
        *args: iterable sequence of any objects with a string representation
    """
    print(ConsoleCode.YELLOW + ' '.join(map(str, args)) + ConsoleCode.END)


def print_function_starting():
    """
    Print a newline and then the current function name.

    Any underscore in the name is replaced with a space, and the prefix 'test' is removed if it exists.
    """
    # noinspection PyProtectedMember
    name = sys._getframe(0).f_code.co_name
    no_test_name = name[name.startswith(TEST_FUNCTION_PREFIX) and len(TEST_FUNCTION_PREFIX):]
    print_bold(f'\n{no_test_name.replace("_", " ")}:')


def print_array_info(name: str, x: numpy.ndarray):
    """
    Print an indented formatted summary of a numpy array to stdout.

    Example:

        Print a summary of a double precision array.

        >>> print(ConsoleCode.BOLD + 'Hello, I and displayed in bold text' + ConsoleCode.END)
        >>> print(ConsoleCode.YELLOW + 'Hello, I and displayed in yellow text' + ConsoleCode.END)

    Args:
        name (str): short name for the array parameter
        x (numpy.ndarray): the array
    """
    shape = (<object> x).shape
    print(f'    {name} array {type(x)} of {x.dtype} with shape {shape}:')
    print(f'{x}')


def get_array_base(array: numpy.ndarray) -> numpy.ndarray:
    """
    For a given Numpy array, finds the base array that 'owns' the actual data.

    Args:
        array (numpy.ndarray): any array

    Returns:
        numpy.ndarray: base - the underlying array holding the data
    """
    base = array
    while isinstance(base.base, numpy.ndarray):
        base = base.base
    return base


def arrays_share_data(x: numpy.ndarray, y: numpy.ndarray) -> bool:
    """
    Determine if arrays share the same underlying data.

    Args:
        x (numpy.ndarray): first array
        y (numpy.ndarray): second array

    Returns:
        bool: True if the two arrays share the same data, otherwise False.
    """
    return get_array_base(x) is get_array_base(y)


cpdef byte_align(a: numpy.ndarray, n=None, dtype=None):
    """
    If given numpy array is not aligned, return a copy that is.
    
    Args:
        a (numpy.ndarray): 
        n (int): alignment one of (16, 32, 64, 128)
        dtype (numpy.dtype): 

    Returns:
        numpy.ndarray: Original array if already aligned, or a new copy of the data.
    """
    if not isinstance(a, ndarray):
        raise UnsupportedArrayType('Invalid array: byte_align requires a subclass of numpy.ndarray')
    if n is None:
        n = _simd_alignment
    else:
        if n not in SUPPORTED_NUM_BYTES:
            raise IllegalArgumentError('n bytes was not a power of 2')
    update_dtype = False
    if dtype is not None:
        if not a.dtype == dtype:
            update_dtype = True
    else:
        dtype = a.dtype
        update_dtype = False

    # See if we're already n byte aligned. If so, do nothing.
    offset = <intptr_t>numpy.PyArray_DATA(a) %n

    if offset is not 0 or update_dtype:
        shape = (<object> a).shape
        _array_aligned = empty_aligned(shape, dtype, n=n)
        _array_aligned[:] = a
        a = _array_aligned.view(type=a.__class__)

    return a


cpdef is_byte_aligned(a: numpy.ndarray, n=None):
    """
    Check if a numpy array is aligned on an n-byte boundary.
    
    Args:
        a (numpy.ndarray): any array  
        n (int): byte alignment, (default=None which uses 64) 

    Returns:
        bool: True if array is aligned, False if not.
    """
    if not isinstance(a, ndarray):
        raise UnsupportedArrayType('Requires a subclass of ndarray')

    if n is None:
        n = _simd_alignment
    else:
        if n not in SUPPORTED_NUM_BYTES:
            raise IllegalArgumentError('n-byte was not a power of 2')

    offset = <intptr_t>numpy.PyArray_DATA(a) %n

    return not bool(offset)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef empty_aligned(shape: Union[int,Iterable[int]], dtype, order='C', n=None):
    """
    Create an empty numpy array that is n-byte aligned.
    
    Uses same arguments as per :func:`numpy.empty` with an additional argument for alignment boundary.
    
    Args:
        shape (Union[int,Iterable[int]]): shape of the array
        dtype (numpy.dtype): dtype of the array
        order (str): specify C or Fortran order
        n (int): byte alignment, (default=None which will use 32 or 64 dependent if CPU supports AVX512) 

    Returns:
        numpy.ndarray: a new ndarray with data aligned on the n-byte boundary.
    """
    if n is None:
        n = _simd_alignment
    else:
        if n not in SUPPORTED_NUM_BYTES:
            raise IllegalArgumentError('n-byte was not a power of 2')

    itemsize = numpy.dtype(dtype).itemsize

    array_size = get_size(shape)

    _array_aligned = numpy.empty(array_size * itemsize + n, dtype='int8')

    _array_aligned_offset = (n - <intptr_t>numpy.PyArray_DATA(_array_aligned))%n

    a = numpy.frombuffer(_array_aligned[_array_aligned_offset:_array_aligned_offset-n].data,
                          dtype=dtype).reshape(shape, order=order)
    return a


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef zeros_aligned(shape, dtype='float64', order='C', n=None):
    """
    Create a numpy array filled with zeros that is n-byte aligned.
    
    Uses same arguments as per :func:`numpy.zeros` with an additional argument for alignment boundary.
    
    Args:
        shape (Iterable[int]): shape of the array
        dtype (numpy.dtype): dtype of the array
        order (str): specify C or Fortran order
        n (int): byte alignment, (default=None which will use 32 or 64 dependent if CPU supports AVX512) 

    Returns:
        numpy.ndarray: a new ndarray with data aligned on the n-byte boundary filled with zero values
    """
    a = empty_aligned(shape, dtype=dtype, order=order, n=n)
    # call specialized template for char* array of bytes
    cdef size_t num_bytes = a.size * a.itemsize
    cdef uintptr_t ptr = <uintptr_t> a.__array_interface__['data'][0]
    with nogil:
        cutil.zero(<char *>ptr, num_bytes)
    return a


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef take(field: numpy.ndarray, indices: numpy.ndarray, shift = 0, boundscheck=False):
    """
    Copy selected values of a field at the given indices to a new byte aligned array.
    
    The index (indices[:] - shift) select values of the field along axis 0. Entries in indices array may be repeated.
    Avoid setting boundscheck=True except when you are debugging.
    
    This function is similar to :func:`numpy.take` except:
    
    * the return value is an array with byte aligned data,
    * performance is optimized, boundscheck is not performed by default and it uses threads and/or SIMD instructions 
    * there is no axis argument, the axis=0 always,
    * the field array is restricted to type of float32 or float64,
    * the indices array is restricted to type int32, int64, uint32, uint64
    * optionally allows a shift of the values in the indices array
    
    Example:
         
        Similar results could be obtained with either the numpy fancy indexing or the take function, as shown.
        For large array sizes this function may be up to twice as fast, for relatively smaller arrays it will be 
        slower because of thread setup and thread switching time.

        >>> local_coordinates = global_coordinates[local_connectivity.global_nodes, ...]         
        >>> local_coordinates = numpy.take(global_coordinates, local_connectivity.global_nodes, axis=0)
        >>> local_coordinates = util.take(global_coordinates, local_connectivity.global_nodes)

    Args:
        field (numpy.ndarray[order="C"]): array of any shape, with dtype in (numpy.float32, numpy.float64)
        indices (numpy.ndarray[shape=N]): 1D array of indices, with dtype in (int32, int64, uint32, uint64)
        shift (int): number to be subtracted from each entry of indices before being used as index into field
        boundscheck (bool): default is False, if any (indices[i] - shift) is outside the range (0, field.shape[0])
            and boundscheck is True an exception is raised, otherwise behavior is undefined and a segfault may
            occur

    Returns:
        numpy.ndarray[order="C"]: out - array of byte aligned data having the same dtype and order as `field`, 
        
        and having the same shape as `field` except the dimension of the first axis is replaced by the size of 
        `indices`, i.e., ``out.shape = (indices.size,) + field.shape[1:]``
                                     
    Raises:
        IndexError: if an (indices - shift) value is less than zero or greater than global_field.shape[0]
        UnsupportedArrayType: if field or indices do not meet requirements on dtype
    """
    cdef pair[int64_t, int64_t] min_max_int64
    cdef pair[int32_t, int32_t] min_max_int32
    cdef pair[uint64_t, uint64_t] min_max_uint64
    cdef pair[uint32_t, uint32_t] min_max_uint32

    cdef uintptr_t f_ptr
    cdef uintptr_t i_ptr
    cdef uintptr_t o_ptr

    cdef int64_t shift_int64
    cdef int32_t shift_int32
    cdef uint64_t shift_uint64
    cdef uint32_t shift_uint32

    cdef size_t num_components
    cdef size_t num_indices
    cdef size_t minimum, maximum

    cdef int out_of_bounds

    # shape of local values is same as the global one, but first dimension is the number of local nodes
    array_shape = list((<object> field).shape)
    array_shape[0] = indices.size

    # get space for the result, making sure it is aligned
    out = empty_aligned(array_shape, dtype=field.dtype)

    # for relatively small sizes of indices fall back to the numpy.take function which is faster than threaded
    if indices.size < 100000:
        if boundscheck:
            mode = 'raise'
        else:
            mode = 'wrap'
        numpy.take(field, indices, axis=0, out=out, mode=mode)
        return out

    num_components = numpy.prod(array_shape[1:])
    num_indices = array_shape[0]
    f_ptr   = field.__array_interface__['data'][0]
    i_ptr = indices.__array_interface__['data'][0]
    o_ptr     = out.__array_interface__['data'][0]
    
    if boundscheck:
        out_of_bounds = False
        minimum = shift
        maximum = shift + num_indices - 1
        if indices.dtype == numpy.dtype(numpy.int64):
            with nogil:
                out_of_bounds = cutil.is_index_out_of_range(<int64_t *> i_ptr, num_indices, minimum, maximum)
        elif indices.dtype == numpy.dtype(numpy.int32):
            with nogil:
                out_of_bounds = cutil.is_index_out_of_range(<int32_t *> i_ptr, num_indices, minimum, maximum)
        elif indices.dtype == numpy.dtype(numpy.uint64):
            with nogil:
                out_of_bounds = cutil.is_index_out_of_range(<uint64_t *> i_ptr, num_indices, minimum, maximum)
        elif indices.dtype == numpy.dtype(numpy.uint32):
            with nogil:
                out_of_bounds = cutil.is_index_out_of_range(<uint32_t *> i_ptr, num_indices, minimum, maximum)
        else:
            raise UnsupportedArrayType('aligned array take() for {} indices not implemented'.format(indices.dtype))
        if out_of_bounds:
            raise IndexError('(indices[i] - shift) out of range of (0, field.shape[0])')
        
    # take values along the first axis according to the indices specified in indices
    if field.dtype == numpy.dtype(numpy.float64):
        if indices.dtype == numpy.dtype(numpy.int64):
            shift_int64 = shift
            with nogil:
                cutil.take(<double *> f_ptr, num_components, <int64_t *> i_ptr, num_indices,
                           shift_int64, <double *> o_ptr)
        elif indices.dtype == numpy.dtype(numpy.int32):
            shift_int32 = shift
            with nogil:
                cutil.take(<double *> f_ptr, num_components, <int32_t *> i_ptr, num_indices,
                           shift_int32, <double *> o_ptr)
        elif indices.dtype == numpy.dtype(numpy.uint32):
            shift_uint32 = shift
            with nogil:
                cutil.take(<double *> f_ptr, num_components, <uint32_t *> i_ptr, num_indices,
                           shift_uint32, <double *> o_ptr)
        elif indices.dtype == numpy.dtype(numpy.uint64):
            shift_uint64 = shift
            with nogil:
                cutil.take(<double *> f_ptr, num_components, <uint64_t *> i_ptr, num_indices,
                           shift_uint64, <double *> o_ptr)
        else:
            raise UnsupportedArrayType('aligned array take() for {} indices not implemented'.format(indices.dtype))
    elif field.dtype == numpy.dtype(numpy.float32):
        if indices.dtype == numpy.dtype(numpy.int64):
            shift_int64 = shift
            with nogil:
                cutil.take(<float *> f_ptr, num_components, <int64_t *> i_ptr, num_indices,
                           shift_int64, <float *> o_ptr)
        elif indices.dtype == numpy.dtype(numpy.int32):
            shift_int32 = shift
            with nogil:
                cutil.take(<float *> f_ptr, num_components, <int32_t *> i_ptr, num_indices,
                           shift_int32, <float *> o_ptr)
        elif indices.dtype == numpy.dtype(numpy.uint32):
            shift_uint32 = shift
            with nogil:
                cutil.take(<float *> f_ptr, num_components, <uint32_t *> i_ptr, num_indices,
                           shift_uint32, <float *> o_ptr)
        elif indices.dtype == numpy.dtype(numpy.uint64):
            shift_uint64 = shift
            with nogil:
                cutil.take(<float *> f_ptr, num_components, <uint64_t *> i_ptr, num_indices,
                           shift_uint64, <float *> o_ptr)
        else:
            raise UnsupportedArrayType('aligned array take() for {} indices not implemented'.format(indices.dtype))
    else:
        raise UnsupportedArrayType('aligned array take() for {} field not implemented'.format(field.dtype))
    return out


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef compress(a: ndarray):
    """
    Compress a numpy array returning a bytes object that can later be decompressed using :func:`.decompress`.
    
    The CompressedArray class is more convenient and does not require saving the shape and dtype.
    
    Args:
        a (numpy.ndarray): an array of any size, dtype, and order 

    Returns:
        bytes: buffer - a Python bytes object
    """
    return blosc.compress_ptr(a.__array_interface__['data'][0], a.size, a.dtype.itemsize, clevel=5, cname='lz4')


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef decompress(buffer: bytes, shape: Union[int, Iterable[int]], dtype: numpy.dtype, order='C'):
    """
    Create a numpy array of the given shape and dtype by decompressing a buffer.
    
    The CompressedArray class is more convenient and does not require saving the shape and dtype.
    
    Args:
        buffer: the bytes object (created from a previous call to :func:`.compress`)
        shape (Iterable[int]): shape of the array
        dtype (numpy.dtype): dtype of the array
        order (str): specify C or Fortran order

    Returns:
        numpy.ndarray: an uncompressed aligned numpy array.
    """
    a = empty_aligned(shape, dtype=dtype, order=order)
    blosc.decompress_ptr(buffer, a.__array_interface__['data'][0])
    return a

# noinspection SpellCheckingInspection,SpellCheckingInspection
cdef class CompressedArray(object):

    # these are not class attributes because this is a cdef extension type

    cdef object _shape      # shape of the original array
    cdef object _dtype      # dtype of the original array
    cdef object _order      # order 'C' or 'F' column or row major order
    cdef object _buffer     # compressed array data, a python str / bytes object


    def __cinit__(self):
        self._shape = None
        self._dtype = None
        self._order = None
        self._buffer = None

    def __init__(self, array = None, shape = None, dtype = None, order = None, pointer = None):
        """
        A buffer storing compressed array data that can be later unpacked into a numpy array with aligned data.

        There are two ways to create the object: either by passing a numpy array as the array data; or
        by passing a shape, dtype, and pointer to the array data.

        Example:

            >>> # create a buffer from an existing array
            >>> b = CompressedArray(numpy.array([[1, 2.0], [3.0, 4]]))
            >>> c = b.unpack()
            >>> numpy.testing.assert_array_equal(b, c)
            >>>
            >>> # create a buffer from a shape, dtype, and a pointer
            >>> a = numpy.array([[1, 2.0], [3.0, 4]])
            >>> d = CompressedArray(shape=a.shape, dtype=a.dtype, pointer=a.__array_interface__['data'][0])
            >>> e = d.unpack()
            >>> numpy.testing.assert_array_equal(a, e)

        Args:
            array (ndarray): (optional) a numpy array (shape and dtype determined from the array)
            shape (tuple): (optional) shape for the uncompressed numpy array
            dtype (numpy.dtype): (optional) type of the elements for the uncompressed numpy array
            order (str): (optional, default='C') order of the array elements 'C' or 'F'
            pointer (uintptr_t): (optional) value of a void*, pointer to the data to be compressed and stored

        Raises:
            IllegalArgumentError: if one or more optional arguments are missing or the wrong type
        """
        cdef size_t length
        cdef uintptr_t p

        if array is not None:
            # we have an array only
            if shape is not None or dtype is not None or order is not None or pointer is not None:
                raise IllegalArgumentError('array must be the only argument')
            self._shape = array.shape
            self._dtype = array.dtype
            length = array.size
            p = array.__array_interface__['data'][0]
        else:
            # we need a shape, a dtype and pointer (a python int acting as a void pointer)
            # and possibly an order
            # raise any exceptions if these parameters are not present as expected with our own exception class
            if shape is None:
                raise IllegalArgumentError('parameter shape (Sequence) is required when a numpy array is not passed')
            self._shape = shape
            length = get_size(shape)
            if dtype is None:
                raise IllegalArgumentError('parameter dtype is required when a numpy array is not passed')
            try:
                self._dtype = numpy.dtype(dtype)
            except:
                raise IllegalArgumentError('parameter dtype cannot be converted to a numpy.dtype')
            if order is None:
                self._order = 'C'
            else:
                if order != 'F' or order != 'C':
                    raise IllegalArgumentError('parameter order must be "C" or "F"')
                self._order = order
            if pointer is None or pointer == 0:
                raise IllegalArgumentError('parameter pointer (non-zero) is required when a numpy array is not passed')
            p = pointer
        if self._dtype == numpy.dtype(numpy.float64):
            cname = 'zlib'
        else:
            cname = 'lz4'
        self._buffer = blosc.compress_ptr(p, length, self._dtype.itemsize, clevel=5, cname=cname)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef unpack(self):
        """
        Unpack (decompress) the compressed ArrayBuffer into an aligned NumPy array.

        Returns:	
            numpy.ndarray: array with aligned data and original shape, dtype, order 
        """
        out = empty_aligned(self._shape, dtype=self._dtype, order=self._order)
        blosc.decompress_ptr(self._buffer, out.__array_interface__['data'][0])
        return out


    @property
    def shape(self) -> Union[int, Iterable[int]]:
        """
        The shape of the stored array.

        Returns:
            shape - an `int` or Iterable[int]
        """
        return self._shape


    @property
    def dtype(self) -> numpy.dtype:
        """
        The scalar type of the stored array.

        Returns:
            numpy.dtype: type of the scalar entries
        """
        return self._dtype


    @property
    def size(self) -> int:
        """
        Number of elements in the stored array.

        Equivalent to np.prod(a.shape), i.e., the product of the array’s dimensions.

        Returns:
            int: total number of entries
        """
        return get_size(self._shape)


    @property
    def itemsize(self) -> int:
        """
        Length of one array element in bytes of the stored array.

        Returns:
            int: number of bytes
        """
        return self._dtype.itemsize


    @property
    def nbytes(self) -> int:
        """
        Total bytes consumed by the elements of the uncompressed array.

        Returns:
            int: number of total bytes of the entries
        """
        return self._dtype.itemsize * get_size(self._shape)


    def __len__(self) -> int:
        """
        Number of bytes used by compressed internal bytes buffer.

        Returns:
            int: length - size of internal python str / bytes buffer
        """
        return self._buffer.__len__()


    def __getstate__(self):
        # support for pickle
        return {'_shape': self._shape, '_dtype': self._dtype, '_buffer': self._buffer}


    def __setstate__(self, d):
        # support for pickle
        self._shape = d['_shape']
        self._dtype = d['_dtype']
        self._buffer = d['_buffer']
