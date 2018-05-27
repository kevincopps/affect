"""
This module contains the classes for reading/writing an ExodusII format mesh database.

.. currentmodule:: affect.exodus

.. |int8_1d| replace:: :obj:`ndarray[int8_t, ndim=1, mode="c"]`
.. |uint32_1d| replace:: :obj:`ndarray[uint32_t, ndim=1, mode="c"]`
.. |uint32_2d| replace:: :obj:`ndarray[uint32_t, ndim=2, mode="c"]`
"""

import collections
import os
import weakref
import logging
import sys
from abc import abstractmethod
from enum import IntEnum, IntFlag
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple, TypeVar, Union

cimport cython
cimport numpy
import numpy
from cpython cimport bool
from cython.operator cimport dereference
from libc.float cimport FLT_MAX
from libc.limits cimport LLONG_MIN
from libc.stdint cimport uintptr_t, int64_t, int32_t, int8_t, uint32_t, INT64_MAX, INT64_MIN
from libc.stdio cimport fflush, stderr
from libc.stdlib cimport malloc, free

from .cimport cexodus
from . import util

cdef extern from "unistd.h":
    enum: STDERR_FILENO


cexodus.omp_set_dynamic(0)      # Explicitly disable dynamic teams
cexodus.omp_set_num_threads(8)  # Use N threads for all consecutive parallel regions

# logger = logging.getLogger('affect')  # logger for debugging
# logger.setLevel(logging.DEBUG)
# ch = logging.StreamHandler(stream=sys.stdout)  # create console handler with a higher log level
# ch.setLevel(logging.DEBUG)
# formatter = logging.Formatter('%(asctime)s %(levelname)s  %(message)s')  # create formatter and add it to the handlers
# ch.setFormatter(formatter)
# logger.addHandler(ch)  # add the handlers to the logger

class _LoggerMessage(object):
    # internal class enabling {} style format to be used in log messages
    # example:
    #   logger.info(__('date = {}', date))
    def __init__(self, fmt, *args, **kwargs):
        self.fmt = fmt
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return self.fmt.format(*self.args, **self.kwargs)

__ = _LoggerMessage


class Error(Exception):
    """
    Base class for all exceptions raised by this module.

    Example:

        >>> import logging
        >>> try:
        >>>     with DatabaseFile('/tmp/myExodusFile.exo', Mode.READ_ONLY) as d:
        >>>         num_times = d.global.num_times()
        >>>         print(f'Number of time steps: {num_times}')
        >>> except Error as e:  # catch all the expected possible errors from the API
        >>>     logging.error(f'Bug in the calling code: {e}')
        >>> except Exception as e:  # catch unexpected errors inside the API
        >>>     logging.error(f'Bug in the exodus module: {e}')
        >>>     raise

        In the example above, the try/except statement on exodus.Error will not insulate you from unforeseen bugs in 
        this exodus module’s code itself, so you may also want to catch the root Python Exception. 
    
    Example:
        
        More specific exceptions can be caught if desired, as shown in the following example.

        >>> try:
        >>>     with DatabaseFile('/tmp/myExodusFile.exo', Mode.READ_ONLY) as d:
        >>>         coord = d.get_coordinates()
        >>> except FileError as e:
        >>>     logging.error(f'Problem with file access of database: {e}')
        >>> except ReadWriteError as e:
        >>>     logging.error(f'Problem reading data from file: {e}')
        >>> except InvalidSpatialDimension as e:
        >>>     logging.error(f'Bug in the calling code: {e}')
    """

class NoMemory(Error, MemoryError):
    """Exception raised during a memory error."""

class FileError(Error):
    """Exception related to opening or accessing a file."""

class FileExists(FileError):
    """Exception raised if existing file prevents writing."""

class FileNotFound(FileError):
    """Exception raised when file cannot be opened."""

class FileAccess(FileError):
    """Exception raised for an unexpected file mode."""

class ReadWriteError(Error, IOError):
    """Exception raised if an issue occurred reading or writing using the internal ExodusII library API."""

class InternalError(Error, RuntimeError):
    """Exception for unexpected internal errors in this module."""

class InvalidEntityType(Error, TypeError):
    """Exception raised if the type of a specified entity was inappropriate."""

class ArrayTypeError(Error, TypeError):
    """Exception raised if array or a field is unknown or a mismatched type or dimension."""

class ArgumentTypeError(Error, TypeError):
    """Exception raised if an argument to a function is not the expected type."""

class InactiveComponent(Error, IOError):
    """Exception raised if attempt is made to access inactive or non-existent component of a field."""

class InactiveEntity(Error, IOError):
    """Exception raised if attempt is made to access inactive or non-existent entity."""

class InactiveField(Error, IOError):
    """Exception raised if attempt is made to access inactive or non-existent field."""

class IndexRangeError(IndexError, IOError):
    """Exception raised if attempt is made to access array out of bounds."""

class InvalidSpatialDimension(Error, ValueError):
    """Exception raised if spatial dimension does not match or is invalid."""

class MaxLengthExceeded(Error, ValueError):
    """Exception raised if a 32 bit unsigned integer limit was exceeded."""

class NotYetImplemented(Error, NotImplementedError):
    """Exception raised if API function is not yet implemented."""

class RangeError(Error, ValueError):
    """Exception raised if given entry indices are out of range of existing data."""

class EntityKeyError(Error, KeyError):
    """Exception raised if an entity ID does not exist in the dictionary of that type."""

class InvalidType(Error, TypeError):
    """Exception raised if an argument to a function or method is the wrong type."""

class CompressionError(Error, TypeError):
    """Exception raised when data is the wrong type in the context of using compressed arrays"""

class Mode(IntEnum):
    """
    Database open modes. These cannot be combined.
    """

    READ_ONLY = 1
    """Mode to initialize a :class:`.Database` for reading an existing file, modifications not allowed."""

    READ_WRITE = 2
    """Mode to initialize a :class:`.Database` for reading and appending or modifying."""

    CREATE = 3
    """Mode to initialize a :class:`.Database` for writing safely, opening only if the file does not already exist."""

    REPLACE = 4
    """Mode to initialize a :class:`.Database` for writing, destroying if the file already exists."""


# Field components
#
# tuples (dimension, components, type)

_field_subscripts_union = ('x', 'y', 'z', 'xx', 'yy', 'zz', 'xy', 'yz', 'zx', 'yx', 'zy', 'xz')

cdef inline _is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

# internal utility functions
  
cdef inline unicode _to_unicode(char* s):
    """Convert C-string bytes to Python string object"""
    return s.decode('UTF-8', 'strict')

cdef inline unicode _to_unicode_with_length(char* s, size_t length):
    """Convert C-string with known length to Python string object."""
    return s[:length].decode('UTF-8', 'strict')

cdef inline bool _is_bulk_int64(int exodus_id):
    """Return True if bulk API is for int64_t, False if int32_t"""
    return (cexodus.ex_int64_status(exodus_id) & cexodus.EX_BULK_INT64_API) != 0

cdef inline bool _is_ids_int64(int exodus_id):
    """Return True if ids API is for int64_t, False if int32_t"""
    return (cexodus.ex_int64_status(exodus_id) & cexodus.EX_IDS_INT64_API) != 0

cdef inline bool _is_maps_int64(int exodus_id):
    """Return True if maps API is for int64_t, False if int32_t"""
    return (cexodus.ex_int64_status(exodus_id) & cexodus.EX_MAPS_INT64_API) != 0

cdef inline void* _allocate(size_t num_bytes):
    """Allocate aligned memory, raises NoMemory: if malloc fails."""
    cdef void* ptr = cexodus.ex_aligned_allocate(num_bytes)
    if ptr == NULL:
        message = f'error: memory allocation (size={num_bytes}) failed'
        raise NoMemory(message)
    # logger.debug(__('    _allocate({}) = {}', num_bytes, hex(<uintptr_t>ptr)))
    return ptr

cdef inline void* _unaligned_allocate(size_t num_bytes):
    """Allocate unaligned memory for c-strings and other types, raises NoMemory: if malloc fails."""
    cdef void * ptr = malloc(num_bytes)
    if ptr == NULL:
        message = f'error: memory allocation (size={num_bytes}) failed'
        raise NoMemory(message)
    # logger.debug(__('    _unaligned_allocate({}) = {}', num_bytes, hex(<uintptr_t>ptr)))
    return ptr

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

cdef inline object _create_array_data(object shape, object dtype, void** data, no_array=False):
    """
    Create data of size shape and type, return data pointer and numpy array (optional).
    
    If no_array=True, then caller must be responsible for calling `free(data)`.
    
    Args:
        shape(Sequence[int]): 
        dtype(numpy.dtype): type of entries of the array 
        no_array: if False (default) return the data and a wrapping numpy array, if True only return the data
        
'   Returns:
        array, data (numpy.ndarray, uintptr_t): the backing array and pointer to the data memory address.
        
    Raises:
        NoMemory: if memory could not be allocated
    """
    cdef uintptr_t ptr
    if no_array:
        array = None
        length = get_size(shape)
        data[0] = _allocate(sizeof(dtype.itemsize) * length)
        # logger.debug(__('    _create_array_data          ({}, {})', shape, hex(<uintptr_t> data[0])))
    else:
        array = util.empty_aligned(shape, dtype=dtype)
        ptr = array.__array_interface__['data'][0]
        data[0] = <void *> ptr
        return array

cdef inline _make_compressed_array(object shape, object dtype, void* pointer):
    """Compress memory buffer to be turned into numpy array at a later time."""
    return util.CompressedArray(shape=shape, dtype=dtype, pointer=<uintptr_t> pointer)

def _assert_is_file(path, exists) -> None:
    """
    Checks file existence or non-existence.
    
    Raises:
        FileNotFound: if parameter ``exists == True`` and the file does not exist
        FileExists: if parameter ``exists == False`` and the file exists
    """
    is_file = os.path.isfile(path)
    if exists:
        if not is_file:
            raise FileNotFound(f'File does not exist {path}.')
    else:
        if is_file:
            abs_path = os.path.abspath(path)
            raise FileExists(f'File already exists {abs_path}.')

def _assert_file_access(path: str, os_mode) -> None:
    """
    Checks OS file access permissions.
    
    Args:
        path: file path
        os_mode: os.F_OK, os.R_OK, os.W_OK, or os.X_OK
    
    Raises:
        FileAccess: if file is not one of the os_mode F_OK, R_OK, W_OK, or X_OK.
    """
    if not os.access(path, os_mode):
        abs_path = os.path.abspath(path)
        mode_type = {os.F_OK:'exists', os.R_OK:'read', os.W_OK:'write', os.X_OK:'execution'}
        raise FileAccess(f'file access "{mode_type[os_mode]}" is not valid for file {abs_path}.')

def _disambiguate_topology_name(name: str, nodes_per_entry: int, spatial_dim: int) -> str:
    """
    Determine an unambiguous element topology name, which is a upper case string. This is necessary because ExodusII
    does not enforce a strict policy on types of element names, so here we standardize common aliases for the same type
    of elements.

    The number of nodes per entry is typically from information about a :class:`Block` entity.
    
    Args:
        name: base name of the topology type
        nodes_per_entry: number of nodes per topology
        spatial_dim: spatial dimension of the topology

    Returns:
        A topology name, all uppercase.

    Raises:
        InternalError: if `nodes_per_entry` doesn't match in some cases of legacy element topology names.
    """
    topology = name.replace(' ','_').upper()
    if topology[:5] == 'SUPER':
        topology = 'SUPER' + str(nodes_per_entry)
    elif spatial_dim == 3:
        if topology[:3] == 'TRI':
            topology = 'TRISHELL' + str(nodes_per_entry)
    elif spatial_dim == 2:
        # normalize shell and bar/rod/truss names in 2D
        last_digit = topology[-1]
        if last_digit == '2':
            if topology == 'SHELL2':
                topology = 'SHELLLINE2D2'
                if nodes_per_entry != 2:
                    raise InternalError(f'Element type {topology} does not match the given 2 nodes_per_entry.')
            elif topology in ('ROD2', 'BAR2', 'TRUSS2'):
                topology = 'ROD2D2'
                if nodes_per_entry != 2:
                    raise InternalError(f'Element type {topology} does not match the given 2 nodes_per_entry.')
        elif last_digit == '3':
            if topology == 'SHELL3':
                topology = 'SHELLLINE2D3'
                if nodes_per_entry != 3:
                    raise InternalError(f'Element type {topology} does not match the given 3 nodes_per_entry.')
            elif topology in ('BAR3', 'ROD3', 'TRUSS3'):
                topology = 'ROD2D3'
                if nodes_per_entry != 3:
                    raise InternalError(f'Element type {topology} does not match the given 3 nodes_per_entry.')
    # append nodes per element if it does not already end with a number
    if not topology[-1].isdigit() and nodes_per_entry > 1:
        topology += str(nodes_per_entry)
    return topology


cdef _raise_io_error():
    """
    Internal function to raise a ReadWriteError signalling a problem accessing the Exodus II database file.
    
    Raises:
        ReadWriteError: always
    """
    cdef char* msg
    cdef char* func
    cdef int err_num
    message = 'Error reading or writing Exodus database file.'
    cexodus.ex_get_err(<const char**>&msg, <const char**>&func, &err_num)
    if err_num != 0 and msg[0] != 0:  # null character '\0 in ASCII'
        message += f'\n[{err_num}]: ({_to_unicode(func)}) {_to_unicode(msg)}'
    raise ReadWriteError(message)

cdef inline _inquire_int(int ex_id, int inquiry):
    """
    Internal function to return one of the metadata values on the Exodus database for integer values only.
    
    Args:
        ex_id: Exodus database ID.
        inquiry(int): variable requested, one of the cexodus.EX_INQ_* values.

    Returns:
        value of the inquired variable, a Python int
    """
    # use the internal exodus function for int only, always returns a signed int64_t
    # int64_t ex_inquire_int(int exoid, int req_info)
    cdef int64_t int64_val
    int64_val = cexodus.ex_inquire_int(ex_id, inquiry)
    if int64_val < 0:
        _raise_io_error()
    return int64_val

# noinspection SpellCheckingInspection
cdef inline _inquire_float(int ex_id, int inquiry):
    """
    Internal function to return one of the metadata values on the Exodus database for floating point
    values only.
    
    This function is probably not necessary. As of 2017-10-01 the only floating point cexodus.EX_INQ_* values are 
    the API, library and DB versions used. These are obtained through other functions elsewhere in the Database class.
    
    * EX_INQ_API_VERS
    * EX_INQ_DB_VERS
    * EX_INQ_LIB_VERS
    
    Args:
        ex_id: Exodus database ID.
        inquiry(cexodus.ex_inquiry): variable requested, one of the cexodus.EX_INQ_* values.

    Returns:
        value of the inquired variable, a Python float
    """
    # return either an integer or float value requested by the inquiry
    # init with some bad values so we will know if one changed
    cdef int64_t int64_val = LLONG_MIN
    cdef float float_val = -FLT_MAX
    cdef char* str_ptr = NULL # not used
    if 0 != cexodus.ex_inquire(ex_id, inquiry, &int64_val, &float_val, str_ptr):
        _raise_io_error()
    if int64_val != LLONG_MIN:
        raise InvalidType('An integer type was requested where a floating point type was expected.')
    elif float_val != -FLT_MAX:
        return float_val
    else:
        return 0

cdef inline _inquire_string(int ex_id, int inquiry):
    """
    Internal function to read one of the strings on an Exodus database.

    Uses an internal temporary allocated buffer of size up to ``cexodus.MAX_LINE_LENGTH + 1`` bytes.

    Args:
        ex_id: Exodus database ID.
        inquiry: variable requested, one of the `cexodus.EX_INQ_*` values.
        
    Returns:
        A UTF-8 encoded str

    Raises:
        InternalError: if length of this type of inquire string is unknown.
    """
    cdef int64_t int_val = LLONG_MIN
    cdef float float_val = -FLT_MAX
    cdef char* str_ptr = NULL
    cdef int len_str = -1
    cdef object unicode_str
    if inquiry == cexodus.EX_INQ_TITLE:
        len_str = cexodus.MAX_LINE_LENGTH
    elif inquiry == cexodus.EX_INQ_GROUP_NAME:
        len_str = cexodus.EX_INQ_GROUP_NAME_LEN
    elif inquiry == cexodus.EX_INQ_FULL_GROUP_NAME:
        len_str = cexodus.EX_INQ_GROUP_NAME_LEN
    else:
        raise InternalError(f'unknown string length for inquire = {inquiry}.')
    try:
        str_ptr = <char *> _unaligned_allocate(sizeof(char) * (len_str+1))
        if 0 != cexodus.ex_inquire(ex_id, inquiry, &int_val, &float_val, &str_ptr[0]):
            _raise_io_error()
        unicode_str = _to_unicode(str_ptr) # python string object
    finally:
        if str_ptr != NULL:
            # logger.debug(__('    _inquire_string ({}, {})', len_str+1, hex(<uintptr_t> str_ptr)))
            free(str_ptr)
    return unicode_str


cdef inline _get_db_max_read_name_length(int ex_id):
    cdef int len_str
    cdef int db_name_size
    len_str = _inquire_int(ex_id, cexodus.EX_INQ_MAX_READ_NAME_LENGTH)
    db_name_size = _inquire_int(ex_id, cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
    if db_name_size < len_str:
        len_str = db_name_size
    return len_str


@cython.boundscheck(False)
@cython.wraparound(False)
def _get_strings(max_string_length, num_strings, func) -> List[str]:
    """
    Internal member function to get strings from database using a given Exodus API function.

    Caller must supply a callback function taking a single argument, declared like: ``def func(uint_ptr str_ptr):``.

    Args:
        max_string_length: C-strings will be allocated with max_string_length+1
        num_strings: number of strings to allocate
        func: the wrapper around an Exodus API function that fills a char**

    Returns:
        strings - a list of python strings
    """
    cdef char** strings_ptr = NULL
    strings_list = []
    if num_strings > 0:
        try:
            strings_ptr = <char**> _unaligned_allocate(sizeof(char*) * num_strings)
            try:
                for i in range(num_strings):
                    strings_ptr[i] = <char *> _unaligned_allocate(sizeof(char) * (max_string_length+1))
                # call the wrapper of the ExodusII API function
                if 0 != func(<uintptr_t>strings_ptr):
                    _raise_io_error()
                for i in range(num_strings):
                    name = _to_unicode(strings_ptr[i])
                    strings_list.append(name)
            finally:
                for i in range(num_strings):
                    free(strings_ptr[i])
        finally:
            free(strings_ptr)
    return strings_list


cdef inline _set_param(int ex_id, cexodus.ex_entity_type set_type, cexodus.ex_entity_id set_id):
    """
    Internal function to return the pair of the (number of entries, number of distribution factors) which
    describe a single set.
    
    Args:
        ex_id(int): Exodus Database ID
        set_type: cexodus.EX_NODE_SET, .EX_SIDE_SET, .EX_EDGE_SET, .EX_FACE_SET, or .EX_ELEM_SET
        set_id(int): ID of the set
        
    Returns:
        num_entries, num_dist_face
    """
    cdef int64_t num_entries64
    cdef int64_t num_dist_fact64
    cdef int num_entries
    cdef int num_dist_fact
    if _is_bulk_int64(ex_id):
        if 0 != cexodus.ex_get_set_param(ex_id, set_type, set_id,
                                         &num_entries64, &num_dist_fact64):
            _raise_io_error()
        return num_entries64, num_dist_fact64
    else:
        if 0 != cexodus.ex_get_set_param(ex_id, set_type, set_id,
                                         &num_entries, &num_dist_fact):
            _raise_io_error()
        return num_entries, num_dist_fact

# global functions for the Exodus API

def library_version() -> str:
    """
    The ExodusII library API version number.
     
    This number reflects the version number of the C language implementation of the internal Exodus 
    library linked with this python module.
    
    Returns:
        Version string of the Exodus API library
    """
    return str(f'{cexodus.API_VERS:.6g}')


class Messages(IntFlag):
    """
    Flags for specifying what kind of messages are output by the internal ExodusII library.
    
    These can be combined using the bitwise operators ``(&, |, ^, ~)`` and the result is an :class:`enum.IntFlag` 
    member.
    
    =======  ============================================================================
    value    meaning
    =======  ============================================================================
    ABORT    Causes internal library fatal errors to force program exit. (Default is off)
    DEBUG    Certain informational messages will print for debug use. (Default is off)
    VERBOSE  All error messages will print. (Default is off)
    DEFAULT  No error messages will print. (Default is on)
    =======  ============================================================================
    """

    DEFAULT     = cexodus.EX_DEFAULT
    VERBOSE     = cexodus.EX_VERBOSE
    DEBUG       = cexodus.EX_DEBUG
    ABORT       = cexodus.EX_ABORT
    NULLVERBOSE = cexodus.EX_NULLVERBOSE


def debug_messages(options: int) -> int:
    """
    Change how internal messages are reported inside the internal Exodus Library API globally for any and all subsequent
    Database files, and Database operations.

    Any messages are written to the stderr device by the internal ExodusII library. Different options change how errors
    and warnings are written from all :class:`Database` objects to stderr. The options may be OR'ed together to provide 
    a combination of options.

    Consider the alternative context manager :class:`.DebugMessages`.

    Example:
        
        >>> debug_messages(Messages.DEBUG | Messages.VERBOSE)
        >>> with DatabaseFile('/tmp/myExodusFile.exo', Mode.READ_ONLY) as e:
        >>>     coord = e.get_coordinates()
        >>> debug_messages(Messages.DEFAULT)
        
    Args:
        options: one or more of :class:`Messages` values or'ed together.
        
    Return:
        The previous value of options internal to the Exodus II API library.
    """
    old_options = cexodus.ex_opts(options)
    return old_options


cdef class DebugMessages(object):
    """
    A context manager for capturing debug messages from the internal library in :class:`.Database` functions.

    Wrap your calls in this context if you are having trouble querying your database and wish to get more
    debug messages from the internal Exodus library. Use the `with` statement while instantiating this class to wrap
    calls to :class:`.Database` member functions. This context manager can be nested within the
    :class:`.DatabaseFile` context manager.

    Example:

        >>> with DatabaseFile('/tmp/myExodusFile.exo', Mode.READ_ONLY) as e:
        >>>     # append messages to any exceptions raised
        >>>     with DebugMessages as messages:
        >>>         coord = e.get_coordinates()
        >>>     # debugging turned off now
        >>> # file automatically now closed

    Upon entering this context, the internal Exodus library is set to DEBUG|VERBOSE options. And any Python
    exceptions raised as a result of Exodus API issues will have their messages appended with the contents of the Exodus
    debug and verbose output. Upon exiting the context, the previous values of the verbose or debug options are
    restored.

    An alternative to using this context manager---if messages to stderr are acceptable for your application---is to
    call :func:`.debug_messages` with the DEBUG and/or VERBOSE option.

    Note:

        This capability will only work on various Unix platforms. A side effect within this context manager is that any
        output to the C language stderr device---not only from the Exodus library but any other calls within the scope
        of the active context manager---is rerouted to a buffer. The stderr device is where the Exodus library writes
        its debug messages.
    """

    cdef int old_option
    cdef int saved_stderr
    cdef int err_pipe[2]
    cdef char* c_messages
    cdef bytearray py_messages

    def __init__(self, option = Messages.VERBOSE|Messages.DEBUG):
        """
        Initialize and activate buffering of Exodus API verbose or debug messages.

        Args:
            option: options to activate messages, either ``Messages.VERBOSE`` or 
                ``Messages.VERBOSE` | Messages.DEBUG``, etc.
        """
        self.py_messages = bytearray(cexodus.MAX_ERR_LENGTH)
        self.c_messages = self.py_messages # the c string shares memory with the python bytearray
        self.old_option = cexodus.ex_opts(option)

    def __enter__(self):
        self._redirect_stderr()
        return self.py_messages

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self._release_stderr()
        cexodus.ex_opts(self.old_option)
        if exc_type is not None:
            # an exception occurred, raise a new one with our message in it
            raise InternalError(f'{str(exc_value)} {self.py_messages}').with_traceback(exc_traceback)
        return False

    def _redirect_stderr(self):
        cdef extern from "unistd.h":
            int close(int fildes)
            int dup(int fildes)
            int dup2(int fildes, int fildes2)
            int pipe(int fildes[2])
        self.saved_stderr = dup(STDERR_FILENO) # save stderr for later
        if pipe(self.err_pipe) != 0:
            raise InternalError('unable to pipe error messages from ExodusII Library on stderr')
        dup2(self.err_pipe[1], STDERR_FILENO) # redirect stderr to the pipe
        close(self.err_pipe[1])

    def _release_stderr(self):
        cdef extern from "unistd.h":
            int dup2(int fildes, int fildes2)
            size_t read(int fildes, void *buf, size_t nbyte)
        # noinspection PyGlobalUndefined
        #global stderr
        fflush(stderr)  # make sure everything is now in the pipe
        read(self.err_pipe[0], self.c_messages, cexodus.MAX_ERR_LENGTH)  # get results from pipe into our buffer
        dup2(self.saved_stderr, STDERR_FILENO)  # reconnect stderr


class EntryType(IntEnum):
    """
    Types of entries that may make up a collection in a Database.
    """

    COORDINATE = cexodus.EX_COORDINATE
    """Denotes collection of coordinates entries."""

    ELEMENT = 65
    """Denotes entries of type element."""

    NODE    = 66
    """Denotes entries of type node."""

    FACE    = 67
    """Denotes entries of type face."""

    EDGE    = 68
    """Denotes entries of type edge."""


class EntityType(IntEnum):
    """
    Types of collections of entries for queries into Database.
    """

    NODAL      = cexodus.EX_NODAL
    NODE_BLOCK = cexodus.EX_NODE_BLOCK
    NODE_SET   = cexodus.EX_NODE_SET
    EDGE_BLOCK = cexodus.EX_EDGE_BLOCK
    EDGE_SET   = cexodus.EX_EDGE_SET
    FACE_BLOCK = cexodus.EX_FACE_BLOCK
    FACE_SET   = cexodus.EX_FACE_SET
    ELEM_BLOCK = cexodus.EX_ELEM_BLOCK
    ELEM_SET   = cexodus.EX_ELEM_SET
    SIDE_SET   = cexodus.EX_SIDE_SET
    ELEM_MAP   = cexodus.EX_ELEM_MAP
    NODE_MAP   = cexodus.EX_NODE_MAP
    EDGE_MAP   = cexodus.EX_EDGE_MAP
    FACE_MAP   = cexodus.EX_FACE_MAP
    GLOBAL     = cexodus.EX_GLOBAL

ALL_ENTITY_TYPES = frozenset(EntityType.__members__.values())
"""The set of all EntityType."""

BLOCK_ENTITY_TYPES = frozenset((EntityType.ELEM_BLOCK, EntityType.FACE_BLOCK, EntityType.EDGE_BLOCK))
"""The set of EntityType that are blocks."""

MAP_ENTITY_TYPES = frozenset((EntityType.ELEM_MAP, EntityType.NODE_MAP, EntityType.EDGE_MAP, EntityType.FACE_MAP))
"""The set of EntityType that are maps."""

SET_ENTITY_TYPES = frozenset((EntityType.SIDE_SET, EntityType.NODE_SET, EntityType.EDGE_SET, EntityType.FACE_SET,
                             EntityType.ELEM_SET))
"""The set of EntityType that are sets."""

ENTITY_TYPES_WITH_ATTRIBUTES = frozenset((EntityType.NODAL, EntityType.NODE_SET, EntityType.EDGE_BLOCK,
                                         EntityType.EDGE_SET, EntityType.FACE_BLOCK, EntityType.FACE_SET,
                                         EntityType.ELEM_BLOCK, EntityType.ELEM_SET, EntityType.SIDE_SET))
"""The set of EntityType with attributes."""


def _get_entity_type(value: int) -> str:
    """
    Return the name of the EntityType that is mapped to the given integer.
    
    Returns:
        name of the EntityType member.
    
    Raises:
        ValueError: if the value does not correspond to one of the EntityType.
    """
    return list(EntityType.__members__.keys())[list(EntityType.__members__.values()).index(value)]

def _raise_invalid_entity_type(string_or_set):
    """
    Raise an exception with the message saying the entity type is not a member of 'string_or_set'.
    
    Args:
        string_or_set:
         
    Raises:
        InvalidEntityType: always
    """
    if isinstance(string_or_set, frozenset):
        message = ', '.join(i.name for i in string_or_set)
    else:
        message = str(string_or_set)
    raise InvalidEntityType(f'The EntityType does not match {message}.')


Component = TypeVar('Component', int, str)  # Must be int or str


class Field(object):

    def __init__(self, name: str, components: Sequence[Component], variables = None, parent_field_dictionary = None):
        """
        A field name, number and name of components, and the ExodusII variable ID for each component.

        The method in which ExodusII stores a multi-component field is by storing multiple scalar variables — the
        scalar variable names associated with a field have the same starting prefix. This class makes it easy
        to access the data in multiple variables as one :class:`.FieldArray`.

        The array values of the field components are accessed through the member functions of the associated
        entity on which the field is active.

        Fields may be associated with an entity of type GLOBAL, NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
        FACE_SET, ELEM_BLOCK, ELEM_SET, or SIDE_SET.
        
        Args:
            name: base name of the field, without trailing subscripts
            components: a sequence containing names of the subscripts of the field components, for example,
                either containing all str ``('xx', 'yy', 'zz', 'xy', 'xz', 'yz')``,
                or containing all int ``[1, 2, 3, 4, 5, 6, 7]``.
            variables(Sequence[int]): a sequence of variable indices in the Exodus database corresponding to the 
                scalars storing the components with length equal to ``len(components)``
            parent_field_dictionary(:class:`.Fields`): the collection of fields that this :class:`Field` object belongs to
            
        Raises:
            ArrayTypeError: if components less than one, or components not equal to associated variable dimensions
        """
        self.name = name
        self.components = components
        self.variables = variables

        # if created from inside an Exodus database function, store the parent entity
        if parent_field_dictionary is not None:
            self.parent_field_dictionary = weakref.proxy(parent_field_dictionary)  # avoid circular reference

        if components is None or len(components) < 1:
            raise ArrayTypeError(f'Number of components of field {name} must be an integer greater than zero.')
        if variables and len(variables) != len(components):
            raise ArrayTypeError(f'Number of variable indices of field {name} != {len(components)} components.')

    def __str__(self):
        s = self.name
        if len(self.components) > 1:
            s += f'_{", ".join(str(n) for n in self.components)}'
        if self.variables is not None:
            s += f' [{", ".join(str(n) for n in self.variables)}]'
        return s


class FieldArray(numpy.ndarray):
    """
    A multidimensional array of double precision, which is a specialization of NumPy array.
    
    It holds a :class:`.Field` object for the field name and the names of its components.

    The constructor requires an already formed :class:`numpy.ndarray` instance (from any of the usual numpy 
    constructors) and optionally a :class:`Field` instance. The len(info.components) must equal 
    ``input_array.shape[-1]``. The ``input_array.dtype`` must be :obj:`numpy.double` or :obj:`numpy.float64`.

    Example: 
    
        >>> pos_info = Field('position', ('x', 'y', 'z'))
        >>> arr = numpy.empty(3, dtype=numpy.double)
        >>> obj = FieldArray(arr, pos_info)
        >>> type(obj)
        <class 'FieldArray'>
        >>> obj.info.name
        'position'
        >>> len(obj.info.components)
        3
        >>> v = obj[1:]
        >>> type(v)
        <class 'FieldArray'>
        >>> len(v.info.components)
        3

    You can use either the :meth:`FieldArray.empty` or :meth:`FieldArray.zeros` class factory to produce a 
    :class:`FieldArray` with given :class:`Field` information and a number of entries.

    Example:
        
        >>> my_coordinate_array = FieldArray.empty(Field('coordinates', ('x', 'y', 'z')), 256)
        >>> my_coordinate_array.shape
        (3, 256)
    """
    
    def __new__(cls, input_array, info=None):
        if input_array is not None:
            if input_array.dtype != numpy.float64:  # Input array is an already formed ndarray instance
                raise ArrayTypeError('input_array.dtype must be numpy.double or numpy.float64.')
        if info is not None:
            if not isinstance(info, Field):
                raise ArrayTypeError('info object must be an instance of Field or subclass of Field')
            if len(info.components) != input_array.shape[-1]:
                raise ArrayTypeError('expected len(info.components) == input_array.shape[-1]')
        obj = numpy.asarray(input_array).view(cls)  # We first cast to be our class type
        obj.info = info                          # add the new attribute to the created instance
        return obj                               # Finally, we must return the newly created object

    def __array_finalize__(self, obj):
        # ``self`` is a new object resulting from ndarray.__new__(FieldArray, ...), therefore it only has
        # attributes that the ndarray.__new__ constructor gave it, i.e. those of a standard ndarray.
        #
        # We could have got to the ndarray.__new__ call in 3 ways:
        # From an explicit constructor - e.g. FieldArray():
        #    obj is None
        #    (we're in the middle of the FieldArray.__new__
        #    constructor, and self.info will be set when we return to
        #    FieldArray.__new__)
        if obj is None: return
        # From view casting - e.g arr.view(FieldArray):
        #    obj is arr
        #    (type(obj) can be FieldArray)
        # From new-from-template - e.g infoarr[:3]
        #    type(obj) is FieldArray
        #
        # Note that it is here, rather than in the __new__ method, that we set the default value for 'info',
        # because this method sees all creation of default objects --- with the FieldArray.__new__ constructor,
        # but also with arr.view(FieldArray).
        self.info = getattr(obj, 'info', None)
        # We do not need to return anything

    @classmethod
    def empty(cls, field: Field, num_entries: int) -> 'FieldArray':
        """
        Construct a new field array with name and number of components and uninitialized values.
        
        The name and number of components are given by the :class:`.Field` object, and the number of entries given by 
        ``num_entries``.

        Args:
            field: object specifying name and names of components
            num_entries: number of entries
        
        Returns:
            Array of uninitialized values with shape ``(num_entries, len(field.components))``, or 
            
            :obj:`None` if ``num_entries < 1``.
        """
        cdef int64_t num_components = len(field.components)
        cdef int64_t n = num_entries
        if n < 1:
            return None
        shape = (num_entries, num_components)
        cdef numpy.ndarray[numpy.double_t, ndim=2] array = util.empty_aligned(shape, dtype=numpy.double)
        return FieldArray(array, field)

    @classmethod
    def zeros(cls, field: Field, num_entries: int) -> 'FieldArray':
        """
        Construct a new field array with name and number of components and values initialized to zero.

        The name and number of components are given by the :class:`.Field` object, and the number of entries given by 
        ``num_entries``.

        Args:
            field: object specifying name and names of components
            num_entries: number of entries
        
        Returns:
            Array of zero values with shape ``(num_entries, len(field.components))``, or 
            
            :obj:`None` if ``num_entries < 1``.
        """
        cdef int64_t num_components = len(field.components)
        cdef int64_t n = num_entries
        if n < 1:
            return None
        shape = (num_entries, num_components)
        cdef numpy.ndarray[numpy.double_t, ndim=2] array = util.zeros_aligned(shape, dtype=numpy.double)
        return FieldArray(array, field)


cdef class LocalConnectivity(object):

    # these are not class attributes because this is a cdef extension type
    cdef object _local_nodes
    cdef object _global_nodes
    cdef public size_t min_global
    cdef public size_t max_global
    cdef public size_t num_unique
    cdef bool _compress
    cdef bool _is_int64

    def __cinit__(self):
        self._local_nodes = None
        self._global_nodes = None
        self.min_global = 0
        self.max_global = 0
        self.num_unique = 0
        self._compress = False
        self._is_int64 = False

    def __init__(self,
                 uintptr_t entry_to_node_global_ptr,
                 object entry_to_node_global_shape,
                 bool is_int64,
                 uintptr_t global_to_local_ptr,
                 size_t num_global_nodes,
                 object nonzero_range=None,
                 bool compress=False):
        """
        Represents a local connectivity on a :class:`.Block`, for example, element-to_node connectivity.

        LocalConnectivity is a Python extension type; instances are only created internally and returned by other
        functions in the Database.

        It makes available the two data arrays:

        * local_nodes - an array storing the local node IDs for each entry
        * global_nodes - an array storing the global node ID corresponding to each local node

        The entries of the local_nodes array are zero-based. Thus, the following two conditions are True:
        ``local_nodes.min == 0`` and ``local_nodes.max == (global_nodes.size - 1)``.

        And it makes available three scalars with information about the content of the global_nodes array:

        * num_unique - the number of unique global nodes referenced in the connectivity, equal to global_nodes.size
        * min_global - the smallest global node ID used in the connectivity, equal to global_nodes.min
        * max_global - the largest global node ID used in the connectivity, equal to global_nodes.max

        (These are made available for convenient access without the need for decompressing the global_nodes array.)

        Args:
            entry_to_node_global_ptr: pointer to flattened array of entry-to-global-node connectivity
            entry_to_node_global_shape: shape of the entry_to_node_global_ptr data (num_entries, num_conn)
            is_int64: True if entry_to_node_global_ptr represents is int64 type
            global_to_local_ptr: pointer to working buffer of at least length (sizeof(uint32_t)*num_global_nodes)
            num_global_nodes: number of global nodes
            nonzero_range: if not None, the (start, end) range of working buffer that may be non-zero from previous use
            compress: if True, store the local_nodes and global_nodes data in compressed form

        Raises:
            MaxLengthExceeded: if number of maximum local index in connectivity exceeds unsigned 32 bit limit 2^32 - 1.
        """
        cdef size_t length_total, min_global, max_global, num_unique
        cdef uintptr_t ptr
        cdef uint32_t * entry_to_node_local_ptr = NULL
        cdef void * local_to_global_ptr = NULL
        cdef numpy.ndarray[numpy.uint32_t, ndim=2] entry_to_node_local = None
        cdef object local_to_global = None

        self._compress = compress
        self._is_int64 = is_int64

        # logger.debug(__('LocalConnectivity.__cinit__(is_int64={}, nonzero_range={}, compress={})',
        #                 is_int64, nonzero_range, compress))

        if nonzero_range is not None:
            if not nonzero_range[1] < num_global_nodes or not nonzero_range[0] < nonzero_range[1]:
                raise RangeError('range violates condition nonzero_range[0] < nonzero_range[1] < num_global_nodes')
            min_global = nonzero_range[0]
            max_global = nonzero_range[1]
        else:
            min_global = 0
            max_global = num_global_nodes

        # zero out space previous touched in our working buffer
        # logger.debug(__('    zero range      = {}-{}', min_global, max_global))
        cdef uint32_t * global_to_local_uint32_ptr = <uint32_t*> global_to_local_ptr
        with nogil:
            cexodus.ex_aligned_zero_uint32(global_to_local_uint32_ptr + min_global, max_global - min_global + 1)

        length_total = entry_to_node_global_shape[0] * entry_to_node_global_shape[1]

        try:
            entry_to_node_local = _create_array_data(entry_to_node_global_shape,
                                                     numpy.uint32,
                                                     <void **>&entry_to_node_local_ptr,
                                                     no_array=compress)
            if is_int64:
                with nogil:
                    num_unique = cexodus.compute_global_to_local64(length_total,
                                                                   <int64_t *> entry_to_node_global_ptr,
                                                                   &max_global,
                                                                   &min_global,
                                                                   entry_to_node_local_ptr,
                                                                   global_to_local_uint32_ptr)
                if num_unique == 0:
                        raise MaxLengthExceeded('Number of local indices exceed unsigned 32 bit limit 2^32 - 1')
                try:
                    local_to_global = _create_array_data(num_unique, numpy.int64,
                                                         &local_to_global_ptr, no_array=compress)
                    with nogil:
                        cexodus.fill_local_to_global64(max_global,
                                                       min_global,
                                                       global_to_local_uint32_ptr,
                                                       <int64_t *> local_to_global_ptr)
                    if compress:
                        # logger.debug('    _global_nodes _make_compressed_array')
                        self._global_nodes = _make_compressed_array(num_unique, numpy.int64,
                                                                    <int64_t*> local_to_global_ptr)
                    else:
                        self._global_nodes = local_to_global
                except:
                    raise
                finally:
                    if compress:
                        if local_to_global_ptr != NULL:
                            # logger.debug(__('    free(local_to_global_ptr {})', hex(<uintptr_t> local_to_global_ptr)))
                            free(local_to_global_ptr)
            else:
                with nogil:
                    num_unique = cexodus.compute_global_to_local32(length_total,
                                                                   <int32_t *> entry_to_node_global_ptr,
                                                                   &max_global,
                                                                   &min_global,
                                                                   entry_to_node_local_ptr,
                                                                   global_to_local_uint32_ptr)
                if num_unique == 0:
                        raise MaxLengthExceeded('Number of local indices exceed unsigned 32 bit limit 2^32 - 1')
                try:
                    local_to_global = _create_array_data(num_unique, numpy.int32,
                                                         &local_to_global_ptr, no_array=compress)
                    with nogil:
                        cexodus.fill_local_to_global32(max_global,
                                                       min_global,
                                                       global_to_local_uint32_ptr,
                                                       <int32_t *> local_to_global_ptr)
                    if compress:
                        # logger.debug('    _global_nodes _make_compressed_array')
                        self._global_nodes = _make_compressed_array(num_unique, numpy.int32,
                                                                    <int32_t*> local_to_global_ptr)
                    else:
                        self._global_nodes = local_to_global
                except:
                    raise
                finally:
                    if compress:
                        if local_to_global_ptr != NULL:
                            # logger.debug(__('    free(local_to_global_ptr {})', hex(<uintptr_t> local_to_global_ptr)))
                            free(local_to_global_ptr)
            if compress:
                # logger.debug('    _local_nodes _make_compressed_array')
                self._local_nodes = _make_compressed_array(entry_to_node_global_shape, numpy.uint32,
                                                           entry_to_node_local_ptr)
            else:
                self._local_nodes = entry_to_node_local
        except:
            raise
        finally:
            if compress:
                if entry_to_node_local_ptr != NULL:
                    # logger.debug(__('    free(entry_to_node_local_ptr {})', hex(<uintptr_t> entry_to_node_local_ptr)))
                    free(entry_to_node_local_ptr)
        self.min_global = min_global
        self.max_global = max_global
        self.num_unique = num_unique


    @property
    def local_nodes(self):
        """
        The local node connectivity, a 2D array where the [i,j] entry is the jth local node ID of the ith entry.

        Returns:
            local_nodes - numpy.ndarray of entry-to-local-node connectivity
        """
        if not self._compress:
            return self._local_nodes  # return the reference to the stored numpy.ndarray
        else:
            return self._local_nodes.unpack()  # unpack the util.CompressedArray object into appropriate numpy.ndarray

    @property
    def global_nodes(self):
        """
        The mapping of from local node ID to global node ID.

        Returns:
            global_nodes - numpy.ndarray of mapping from each local node in range (0, num_local_nodes) to a global node
                in the range (min_global, max_global).
        """
        if not self._compress:
            return self._global_nodes  # return the reference to the stored numpy.ndarray
        else:
            return self._global_nodes.unpack()  # unpack the util.CompressedArray object into appropriate numpy.ndarray
        
    def __getstate__(self):
        # support for pickle
        return {'_local_nodes':  self._local_nodes,
                '_global_nodes': self._global_nodes,
                'min_global':    self.min_global,
                'max_global':    self.max_global,
                'num_unique':    self.num_unique,
                '_compress':     self._compress,
                '_is_int64':     self._is_int64}

    def __setstate__(self, d):
        # support for pickle
        self._local_nodes = d['_local_nodes']
        self._global_nodes = d['_global_nodes']
        self.min_global = d['min_global']
        self.max_global = d['max_global']
        self.num_unique = d['num_unique']
        self._compress = d['_compress']
        self._is_int64 = d['_is_int64']


@cython.boundscheck(False)
@cython.wraparound(False)
cdef object _create_connectivity_local(object block,
                                       object entry_type,
                                       size_t num_entries,
                                       size_t num_conn_entries,
                                       size_t min_global_node,
                                       size_t max_global_node,
                                       size_t num_global_nodes,
                                       uint32_t * global_to_local_ptr,
                                       bool compress):
    """
    Internal factory method to return LocalConnectivity object using supplied temporary buffer.
    """
    cdef void * entry_to_node_global_ptr
    cdef bool is_int64 = _is_bulk_int64(block.ex_id)
    cdef size_t length = num_entries * num_conn_entries
    entry_to_node_global_shape = (num_entries, num_conn_entries)
    nonzero_range = (min_global_node, max_global_node)
    # logger.debug(__('_create_connectivity_local: entry_to_node_global_shape = {}', entry_to_node_global_shape))
    try:
        _create_connectivity(block,
                             entry_type,
                             entry_to_node_global_shape,
                             is_int64,
                             &entry_to_node_global_ptr,
                             zero_based=True,
                             no_array=True)
        # logger.debug('_create_connectivity_local: LocalConnectivity')
        local_connectivity = LocalConnectivity(<uintptr_t> entry_to_node_global_ptr,
                                               entry_to_node_global_shape,
                                               is_int64,
                                               <uintptr_t> global_to_local_ptr,
                                               num_global_nodes,
                                               nonzero_range,
                                               compress)
    except:
        raise
    finally:
        # logger.debug(__('    _create_connectivity_local: free(entry_to_node_global_ptr({}))',
        #                 hex(<uintptr_t> entry_to_node_global_ptr)))
        free(entry_to_node_global_ptr)
    return local_connectivity


class Entity(object):
    """
    An Entity is a generic association of variables, maps, or a collection of elements, faces, edges, or nodes.

    An Entity is abstract and is the base class of :class:`Block`, :class:`Set`, :class:`Map' classes.
    """

    def __init__(self, database_id=-1, entity_type=None, entity_id=-1, name=None, parent_entity_collection=None,
                 *args, **kwargs):
        """
        Initialize an entity from a specific Database.
        
        Args:
            database_id: the Exodus library ID of the database file
            entity_type: one of NODAL, GLOBAL;
                NODE_SET, SIDE_SET, ELEM_SET, EDGE_SET, FACE_SET;
                EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK;
                ELEM_MAP, NODE_MAP, EDGE_MAP, FACE_MAP.
            entity_id: the ID of the entity in the database
        
        Raises:
            InvalidEntityType: if the entity_type and parent_entity_collection.entity_type do not match
        """
        # logger.debug(__('Entity {}', self.__class__.__mro__))
        super(Entity, self).__init__(*args, **kwargs)
        if parent_entity_collection is not None and entity_type is not None:
            try:
                if entity_type != parent_entity_collection.entity_type:
                    _raise_invalid_entity_type('the type of collection of its parent')
            except AttributeError:
                # Global and Nodal instances may not have the entity_type attribute initialized yet
                pass
        self.ex_id = database_id
        self.entity_type = entity_type
        self.entity_id = entity_id
        self.parent_entity_collection = weakref.proxy(parent_entity_collection)  # avoid cyclic references
        self._name = name
        self._num_entries = -1


    def num_entries(self) -> int:
        """
        The number of entries contained within a single collection entity.

        For example return the number of elements in an element block, or return number of faces in a face set, etc.
        
        ======================================= ===================================================================
        entity_type                             range of num_entries
        ======================================= ===================================================================
        GLOBAL                                  1
        NODAL                                   number of global nodes
        EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK      number of edges, faces, or elements in the block, respectively
        NODE_SET, EDGE_SET, FACE_SET, ELEM_SET  number of nodes, edges, faces, or elements in the set, respectively
        SIDE_SET                                number of element sides in the side set
        NODE_MAP, EDGE_MAP, FACE_MAP, ELEM_MAP  number of global nodes, edges, faces, or elements, respectively
        ======================================= ===================================================================
        
        Returns:
            number of entries associated with the entity
            
        Raises:
            InvalidEntityType: if entity_type of this collection is invalid or does not support entries
        """
        cdef cexodus.ex_block block
        cdef int64_t num_entries = 0
        cdef int64_t num_entry_in_set
        cdef int64_t num_dist_fact_in_set
        if self._num_entries == -1:
            if self.entity_type == EntityType.NODAL:
                num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)
            elif self.entity_type == EntityType.GLOBAL:
                num_entries = 1
            elif self.entity_type in BLOCK_ENTITY_TYPES:
                block.id = self.entity_id
                block.type = self.entity_type
                if 0 != cexodus.ex_get_block_param(self.ex_id, &block):
                    _raise_io_error()
                num_entries = block.num_entry
            elif self.entity_type in SET_ENTITY_TYPES:
                num_entries, num_dist_fact = _set_param(self.ex_id, self.entity_type.value, self.entity_id)
            elif self.entity_type == EntityType.NODE_MAP:
                num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES) # returns regardless of map existing
            elif self.entity_type == EntityType.EDGE_MAP:
                num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE)  # returns regardless of map existing
            elif self.entity_type == EntityType.FACE_MAP:
                num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_FACE)  # returns regardless of map existing
            elif self.entity_type == EntityType.ELEM_MAP:
                num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM)  # returns regardless of map existing
            else:
                _raise_invalid_entity_type(ALL_ENTITY_TYPES)
            self._num_entries = num_entries
        return self._num_entries


class EntityWithProperty(Entity):

    def __init__(self, database_id=-1, entity_type=None, entity_id=-1, name=None, parent_entity_collection=None,
                 *args, **kwargs):
        """
        Initialize a type of entity which supports properties.

        :term:`properties` are named integer variables associated with every entity of a certain type in the database.
        
        Args:
            database_id(int): (optional) the ID of the associated Exodus Database
            entity_type(EntityType): one of NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK, FACE_SET, ELEM_BLOCK, ELEM_SET, 
                SIDE_SET, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP.
            entity_id(int): ID of the entity in the Database
            name(str): name of the entity
        """
        # logger.debug(__('EntityWithProperty {}', self.__class__.__mro__))
        super(EntityWithProperty, self).__init__(database_id=database_id, entity_type=entity_type, entity_id=entity_id,
                                                 name=name, parent_entity_collection=parent_entity_collection,
                                                 *args, **kwargs)

    def property(self, property_name: str) -> int:
        """
        The value of a property on a given entity collection.

        Properties are integer values associated with an entity.
        
        Args:
            property_name: name of the property
            
        Returns:
            An integer property value.
        """
        py_byte_string = property_name.encode('UTF-8')
        cdef char* c_string = py_byte_string
        cdef int value = 0
        if 0 != cexodus.ex_get_prop(self.ex_id, self.entity_type.value, self.entity_id, c_string, &value):
            _raise_io_error()
        return value


class EntityWithVariable(Entity):

    def __init__(self, database_id=-1, entity_type=None, entity_id=-1, name=None, parent_entity_collection=None,
                 *args, **kwargs):
        """
        Initialize an entity type with association of variables.
        
        The entity_type of ELEM_MAP, FACE_MAP, EDGE_MAP, and NODE_MAP do not have variables.

        A variable in Exodus is a scalar only. Multiple variables are required to make up either a vector or
        tensor field.

        Args:
            database_id(int): (optional) the ID of the associated Exodus Database
            entity_type(EntityType): one of NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK, FACE_SET, ELEM_BLOCK, ELEM_SET, 
                SIDE_SET, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP.
            entity_id(int): ID of the entity in the Database
            name(str): name of the entity
        
        Raises:
            InvalidEntityType: if the entity_type and parent_entity_collection.entity_type do not match
        """
        # logger.debug(__('EntityWithVariable {}', self.__class__.__mro__))
        super(EntityWithVariable, self).__init__(database_id=database_id, entity_type=entity_type, entity_id=entity_id,
                                                 name=name, parent_entity_collection=parent_entity_collection,
                                                 *args, **kwargs)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def is_field_active(self) -> numpy.ndarray:
        """
        Read the array of one or zero values indicating whether each entity field exists on this entity.

        All variables on :class:`Nodal` and :class:`Global` entities are always active. Uses an internal temporary
        allocated buffer up to the size ``num_vars * 8`` bytes.
        
        Returns:
            An array with shape ``len(num_vars,)`` of type :obj:`numpy.int32`, or 
            
            :obj:`None` if no variables are present.
            
        Raises:
            InactiveComponent: if components of one of the fields are not all active or all inactive.
        """
        cdef int * var_ptr = NULL
        cdef int num_var = 0
        fields_dict = self.parent_entity_collection.fields
        cdef int num_fields = len(fields_dict)
        if num_fields < 1:
            return None
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type.value, &num_var):
            _raise_io_error()
        if num_var != sum(len(f.components) for f in fields_dict.values()):  # consistency check
            raise InternalError('The sum of all the field components is not equal to the number of variables.')
        array_shape = num_fields
        cdef numpy.ndarray[numpy.int32_t, ndim=1] fields_table = util.empty_aligned(array_shape, dtype=numpy.int32)
        cdef int *fields_table_ptr = <int *> fields_table.data
        try:
            var_ptr = <int*> _allocate(sizeof(int) * num_var)
            if 0 != cexodus.ex_get_object_truth_vector(self.ex_id, self.entity_type.value, self.entity_id,
                                                       num_var, var_ptr):
                _raise_io_error()
            for index, k, in enumerate(fields_dict):
                v = fields_dict[k]
                variables = v.variables
                if var_ptr[variables[0]] == 1:
                    for i in range(1,len(variables)):  # ensure all the any other component variables are also active
                        if var_ptr[variables[i]] != 1:
                            raise InactiveComponent('Exodus database has a field with an inactive component:\n' 
                                f'  database {self.ex_id} field {k} component {i} of {len(variables)} ' 
                                f'on entity {self.entity_type.name}).')
                    fields_table_ptr[index] = 1
        finally:
            if var_ptr != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> var_ptr)))
                free(var_ptr)
        return fields_table

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def field(self, field: Field, time_step: int) -> FieldArray:
        """
        Read the array of field values of a single field at a time step for all entries.

        Uses a temporarily allocated buffer of size ``num_entries * 8`` bytes.
        
        Args:
            field: the desired field (name and component names)
            time_step: the time step index, at which the field component values are desired,
                the first time step is 0; a value of -1 will return the values at the last time step.

        Returns:
            a byte aligned array of :obj:`numpy.float64` values with shape ``(num_entries, len(field.components))``,
            or 
            
            :obj:`None` if ``num_entries == 0`` or field is not associated with ExodusII variable indices.
                
        Raises:
            InactiveField: if the given field does not exist on this entity.
        """
        cdef double * buffer = NULL
        cdef double * array_ptr = NULL
        cdef int64_t i
        cdef int64_t j
        cdef int64_t num_entries
        cdef int k, var_index, num_components
        cdef numpy.ndarray[numpy.double_t, ndim=2] array

        if not field.name in self.parent_entity_collection.fields:
            raise InactiveField(f'The given field {field.name} does not exist on this entity.')
        num_components = len(field.components)
        if num_components == 1:
            return self.variable(field.variables[0], time_step+1)

        num_entries = self.num_entries()
        if num_entries < 1:
            return None
        if field.variables is None:
            return None
        shape = (num_entries, num_components)
        array = util.empty_aligned(shape, dtype=numpy.double)
        array_ptr = <double *> array.data
        try:
            buffer = <double*> _allocate(sizeof(double) * num_entries)
            for k in range(0, num_components):
                var_index = field.variables[k] + 1
                if 0 != cexodus.ex_get_var(self.ex_id, time_step+1, self.entity_type.value, var_index, self.entity_id,
                                           num_entries, buffer):
                    _raise_io_error()
                with nogil:
                    cexodus.ex_aligned_copy_stride(buffer, num_entries, array_ptr + k, num_components)
        finally:
            if buffer != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> buffer)))
                free(buffer)
        return FieldArray(array, field)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def field_at_times(self, field: Field, entry_index: int, start_time_step: int,
                       stop_time_step: int) -> FieldArray:
        """
        Read an array of values of a field for a single entry (node, element, edge, face) on a range of time steps.

        The entry_index has the following meaning inside this method:
        
        =========== ===========================================
        entity_type entry_index
        =========== ===========================================
        NODAL       a global node index
        BLOCK       local entry index of the associated entries
        SET         local entry index of the associated entries
        GLOBAL      ignored
        =========== ===========================================
        
        Uses an internal temporary allocated buffer of size ``num_time * 8`` bytes.
        
        Args:
            field: the desired field (name and component names)
            entry_index: index of the desired entry in the :term:`entity`, in range ``[0, num_entries-1]``
            start_time_step: starting time step, in range ``[0, num_times]``
            stop_time_step: one past the end of the last time step, use -1 to indicate the final time step
            
        Returns:
            array of :obj:`numpy.float64` values with
            
            shape ``(stop_time_step - start_time_step, len(field.components))`` if ``len(field.components) > 1``, or
            
            shape ``(stop_time_step - start_time_step)`` if ``len(field.components) == 1``.
                
        Raises:
            InactiveField: if the given field does not exist on this entity.
        """
        if not field.name in self.parent_entity_collection.fields:
            raise InactiveField(f'The given field {field.name} does not exist on the entity.')
        cdef int num_components = len(field.components)
        if num_components == 1:
            return self.variable_at_times(field.variables[0], entry_index, start_time_step, stop_time_step)
        if field.variables is None:
            return None
        cdef int64_t end_time_step
        if stop_time_step == -1:
            end_time_step = _inquire_int(self.ex_id, cexodus.EX_INQ_TIME)
        else:
            end_time_step = stop_time_step
        cdef int num_time = end_time_step - start_time_step
        if num_time < 1:
            return None
        cdef double * array_ptr = NULL
        cdef double * var_ptr = NULL
        shape = (num_time, num_components)
        cdef numpy.ndarray[numpy.double_t, ndim=2] array = util.empty_aligned(shape, dtype=numpy.double)
        array_ptr = <double *> array.data
        cdef int i, j, k
        try:
            var_ptr = <double*> _allocate(sizeof(double) * num_time)
            k = 0
            while k < num_components:
                var_index = field.variables[k] + 1
                if 0 != cexodus.ex_get_var_time(self.ex_id, self.entity_type.value, var_index, entry_index+1,
                                                start_time_step+1, end_time_step, var_ptr):
                    _raise_io_error()
                i = 0
                j = k
                while i < num_time:
                    array_ptr[j] = var_ptr[i]
                    i += 1
                    j += num_components
                k += 1
        finally:
            if var_ptr != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> var_ptr)))
                free(var_ptr)
        return FieldArray(array, field)

    def field_partial(self, field: Field, time_step: int, start_entry: int, stop_entry: int) -> FieldArray:
        """
        Read an array of values of a single field on a subsection of the entries on an entity.

        Uses an internal temporary allocated buffer of size ``(stop_entry - start_entry) * 8`` bytes.
        
        Args:
            field: the field
            time_step: the time time_step index, which the object variable values are desired,
                the first time time_step is 0.
            start_entry: first entry at which field is read,
                in the range of ``[0, num_entry-1]``.
            stop_entry: one past the end of the last desired entry at which the field is read,
                in the range of ``[start_entry+1, num_entry]``
                
        Returns:
            array with shape ``(stop_entry - start_entry, len(field.components))`` of :obj:`numpy.double`
                
        Raises:
            InactiveField: if the given field does not exist on this entity.
        """
        if not field.name in self.parent_entity_collection.fields:
            raise InactiveField(f'The given field {field.name} does not exist on the entity.')
        if field.variables is None:
            return None
        num_components = len(field.components)
        if num_components == 1:
            return self.variable_partial(field.variables[0], time_step, start_entry, stop_entry)
        cdef int64_t i
        cdef int64_t j
        cdef int64_t num_entries = stop_entry - start_entry
        if num_entries < 1:
            return None
        shape = (num_entries, num_components)
        cdef numpy.ndarray[numpy.double_t, ndim=2] array = util.empty_aligned(shape, dtype=numpy.double)
        cdef double * array_ptr = <double *> array.data
        cdef double* buffer = NULL
        try:
            buffer = <double*> _allocate(sizeof(double) * num_entries)
            for k in range(0,num_components):
                var_index = field.variables[k] + 1
                if 0 != cexodus.ex_get_partial_var(self.ex_id, time_step, self.entity_type.value, var_index,
                                                   self.entity_id, start_entry+1, num_entries, buffer):
                    _raise_io_error()
                j = k
                for i in range(num_entries):
                    array_ptr[j] = buffer[i]
                    j += num_components
        finally:
            if buffer != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> buffer)))
                free(buffer)
        return FieldArray(array, field)

    def is_variable_active(self) -> numpy.ndarray:
        """
        Read an array of one or zero values indicating whether each variable exists or not, respectively,
        for a single associated entity.

        All variables on :term:`nodal` and :term:`global` entities are always active.
        
        Returns:
            array of :obj:`numpy.int32` with shape ``(num_vars,)``, or 
            
            :obj:`None` if no variables are present.
        """
        cdef int num_var = 0
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type.value, &num_var):
            _raise_io_error()
        if num_var < 1:
            return None
        array_shape = num_var
        cdef numpy.ndarray[numpy.int32_t, ndim=2] var_active = util.empty_aligned(array_shape, dtype=numpy.int32)
        cdef int *var_ptr = <int *> var_active.data
        if 0 != cexodus.ex_get_object_truth_vector (self.ex_id, self.entity_type.value, self.entity_id,
                                                    num_var, var_ptr):
            _raise_io_error()
        return var_active

    def variable(self, variable_index: int, time_step: int) -> numpy.ndarray:
        """
        Read the array of values of a single variable at a given time step.
        
        Args:
            variable_index: index of the desired variable, indices start at 0, in range ``[0, num_variables-1]``
            time_step: the time step index, which the object variable values are desired, the first time step is 0.

        Returns:
            array with shape ``(num_entries,)`` of type :obj:`numpy.float64`, or 
            
            :obj:`None` if there are no entries.
        """
        cdef int64_t var_length = self.num_entries()
        if var_length < 1:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] var = util.empty_aligned(var_length, dtype=numpy.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_var(self.ex_id, time_step+1, self.entity_type.value, variable_index+1, self.entity_id,
                                   var_length, var_ptr):
            _raise_io_error()
        return var

    def variable_at_times(self, variable_index: int, entry_index: int, start_time_step: int,
                          stop_time_step: int) -> numpy.ndarray:
        """
        Read the array of values of a single variable on a range of time steps.

        The entry_index has the following meaning, for a :class:`Nodal` entity it is a global node index,
        for :class:`Block` and :class:`Set` this is a the local entry index.

        The entry_index is ignored for the :class:`Global` entity.
        
        Args:
            variable_index: index of the desired variable, indices start at 0, in range ``[0, num_variables-1]``
            entry_index: index of the desired entry in the Entity, in range ``[0, num_entries-1]``
            start_time_step: starting time step, in range ``[0, num_times]``
            stop_time_step: one past the end of the last time step, use -1 to indicate the final time step

        Returns:
            array with shape ``(stop_time_step - start_time_step)`` of :obj:`numpy.float64` or 
            
            :obj:`None` if there are no entries
        """
        cdef int64_t num_entries = self.num_entries()
        cdef int64_t end_time_step
        if num_entries < 1:
            return None
        if stop_time_step == -1:
            end_time_step = _inquire_int(self.ex_id, cexodus.EX_INQ_TIME)
        else:
            end_time_step = stop_time_step
        cdef int64_t num_time = end_time_step - start_time_step
        if num_time < 1:
            return None
        array_shape = num_time
        cdef numpy.ndarray[numpy.double_t, ndim=1] var = util.empty_aligned(array_shape, dtype=numpy.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_var_time(self.ex_id, self.entity_type.value, variable_index+1, entry_index+1,
                                        start_time_step+1, end_time_step, var_ptr):
            _raise_io_error()
        return var

    def variable_partial(self, var_index: int, time_step: int, start_entry: int,
                         stop_entry: int) -> numpy.ndarray:
        """
        Read the array of values of a single variable on a subsection of the entries on an entity.
        
        Args:
            var_index: index of the desired variable, indices start at 0,
                in range ``[0, num_variables-1]``
            time_step: the time time_step index, which the object variable values are desired,
                the first time time_step is 0.
            start_entry: first entry at which variables are read,
                in the range of ``[0, num_entry-1]``.
            stop_entry:  one past the end of the last desired entry at which the variable is read,
                in the range of ``[start_entry+1, num_entry]``
                            
        Returns:
            array of type :obj:`numpy.double` with shape ``(stop_entry - start_entry)`` or
            
            :obj:`None` if there are no entries.
        """
        cdef int64_t num_entries = stop_entry - start_entry
        if num_entries < 1:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] var = util.empty_aligned(num_entries, dtype=numpy.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_partial_var(self.ex_id, time_step, self.entity_type.value, var_index,
                                           self.entity_id, start_entry+1, num_entries, var_ptr):
            _raise_io_error()
        return var


class EntityWithAttribute(EntityWithVariable):

    def __init__(self, database_id=-1, entity_type=None, entity_id=-1, name=None, parent_entity_collection=None,
                 *args, **kwargs):
        """
        Initialize an entity type with association of attributes.

        The entity_type of GLOBAL, ELEM_MAP, FACE_MAP, EDGE_MAP, and NODE_MAP do not have attributes.

        An attribute is a per entity scalar variable in Exodus.
        
        Args:
            database_id(int): (optional) the ID of the associated Exodus Database
            entity_type(EntityType): one of NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
                FACE_SET, ELEM_BLOCK, ELEM_SET, SIDE_SET.
            entity_id(int): ID of the entity in the Database
            name(str): name of the entity

        """
        # logger.debug(__('EntityWithAttribute {}', self.__class__.__mro__))
        super(EntityWithAttribute, self).__init__(database_id=database_id, entity_type=entity_type,
                                                  entity_id=entity_id, name=name,
                                                  parent_entity_collection=parent_entity_collection, *args, **kwargs)

    def num_attributes(self) -> int:
        """
        The number of attributes on a single collection entity.

        Only NODAL, EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK, NODE_SET, EDGE_SET, FACE_SET, SIDE_SET, ELEM_SET have
        attributes.

        Returns:
            Number of attributes
        """
        cdef int num_attribute = 0
        if self.entity_type in ENTITY_TYPES_WITH_ATTRIBUTES:
            if 0 != cexodus.ex_get_attr_param(self.ex_id, self.entity_type.value, self.entity_id, &num_attribute):
                _raise_io_error()
        return num_attribute

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def attribute_names(self) -> List[str]:
        """
        Read the attribute names on a single entity.

        Attribute names were added to the Exodus API library in version 4.26; databases written earlier than that
        will return empty strings.

        Uses an internal temporary allocated buffer of size
        ``num_attribute * (8 + cexodus.EX_INQ_MAX_READ_NAME_LENGTH + 1)`` bytes.
        
        Returns:
            list of string names
        """
        cdef int num_attribute = 0
        len_str = _get_db_max_read_name_length(self.ex_id)
        if self.entity_type in ENTITY_TYPES_WITH_ATTRIBUTES:
            if 0 != cexodus.ex_get_attr_param(self.ex_id, self.entity_type.value, self.entity_id, &num_attribute):
                _raise_io_error()

        def attribute_names_func(uintptr_t str_ptr):
            return cexodus.ex_get_attr_names(self.ex_id, self.entity_type.value, self.entity_id, <char **>str_ptr)

        return _get_strings(len_str, num_attribute, attribute_names_func)

    def attributes(self) -> numpy.ndarray:
        """
        Read all the attribute values for all the entries in the given collection.
        
        Returns:
            array of :obj:`numpy.double` with shape ``(num_entries, num_attribute)``
        """
        cdef int64_t num_entries = self.length()
        cdef int64_t num_attribute = self.num_attribute()
        if num_attribute < 1:
            return None
        array_shape = (num_entries, num_attribute)
        cdef numpy.ndarray[numpy.double_t, ndim=2] attributes = util.empty_aligned(array_shape, dtype=numpy.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_attr(self.ex_id, self.entity_type.value, self.entity_id, attributes_ptr):
                _raise_io_error()
        return attributes

    def attributes_partial(self, int start, int stop) -> numpy.ndarray:
        """
        Read all the attribute values for a subsection of the entries in the given collection.
        
        Args:
            start: first local entry index at which attributes are desired, in the range ``[0, num_entry - 1]``.
            stop:  one past the end of the last desired local entry index, in the range ``[start+1, num_entry]``.
            
        Returns:
            Array of attribute values of type :obj:`numpy.double` with shape ``(stop - start, num_attribute)`` or
            
            :obj:`None` if there are no attributes 
        """
        cdef int64_t num_entries = stop - start
        cdef int64_t num_attribute = self.num_attribute()
        if num_attribute < 1 or num_entries < 1:
            return None
        array_shape = (num_entries, num_attribute)
        cdef numpy.ndarray[numpy.double_t, ndim=2] attributes = util.empty_aligned(array_shape, dtype=numpy.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_partial_attr(self.ex_id, self.entity_type.value, self.entity_id, start+1,
                                            num_entries, attributes_ptr):
                _raise_io_error()
        return attributes

    def attribute(self, attrib_index: int) -> numpy.ndarray:
        """
        Read an array of values of a single attribute for all entries in the given collection.

        Args:
            attrib_index: desired attribute in range ``[0, num_attribute-1]``.
            
        Returns:
            Array of type :obj:`numpy.double` with shape ``(num_entries)``, or 
            
            :obj:`None` if this type of entity does not support attributes
        """
        cdef int64_t num_entries = self.length()
        if self.entity_type not in ENTITY_TYPES_WITH_ATTRIBUTES:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] attributes = util.empty_aligned(num_entries, dtype=numpy.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_one_attr(self.ex_id, self.entity_type.value, self.entity_id,
                                        attrib_index, attributes_ptr):
                _raise_io_error()
        return attributes

    def attribute_partial(self, start: int, stop: int, attrib_index: int) -> numpy.ndarray:
        """
        Read one entity attribute for a subsection of entries in the given collection.
        
        Args:
            start: first local entry index at which attributes are desired, in the range ``[0, num_entry - 1]``.
            stop:  one past the end of the last desired local entry index, in the range ``[start+1, num_entry]``.
            attrib_index: desired attribute in range ``[0, num_attribute-1]``
            
        Returns:
            Array of attribute values of type :obj:`numpy.double` with shape ``(stop - start)`` or
            
            :obj:`None` if there are no entries or this entity type does not support attributes
        """
        cdef int64_t num_entries = stop - start
        if num_entries < 1 or self.entity_type not in ENTITY_TYPES_WITH_ATTRIBUTES:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] attributes = util.empty_aligned(num_entries, dtype=numpy.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_partial_one_attr(self.ex_id, self.entity_type.value, self.entity_id, start+1,
                                                num_entries, attrib_index, attributes_ptr):
                _raise_io_error()
        return attributes


class BaseEntityCollection(object):
    """
    Base class of the Global, Nodal, Sets, Blocks, Maps classes acting as collections of entities.
    """

    def __init__(self, database_id=-1, entity_type=None, *args, **kwargs):
        super(BaseEntityCollection, self).__init__(*args, **kwargs)

        #: int: the ExodusII database ID
        self.ex_id = database_id

        #: EntityType: the type of entity held by this collection
        self.entity_type = entity_type

    def _num_entity(self) -> int:
        """
        The number of entities in this collection.
        
        Returns:
            number of entities
        """
        cdef int64_t num = 0
        if self.entity_type == EntityType.ELEM_BLOCK:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM_BLK)
        elif self.entity_type == EntityType.NODE_SET:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_NODE_SETS)
        elif self.entity_type == EntityType.SIDE_SET:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_SIDE_SETS)
        elif self.entity_type == EntityType.ELEM_MAP:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM_MAP)
        elif self.entity_type == EntityType.NODE_MAP:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_NODE_MAP)
        elif self.entity_type == EntityType.EDGE_BLOCK:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE_BLK)
        elif self.entity_type == EntityType.FACE_BLOCK:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_FACE_BLK)
        elif self.entity_type == EntityType.EDGE_SET:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE_SETS)
        elif self.entity_type == EntityType.FACE_SET:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_FACE_SETS)
        elif self.entity_type == EntityType.ELEM_SET:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM_SETS)
        elif self.entity_type == EntityType.EDGE_MAP:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE_MAP)
        elif self.entity_type == EntityType.FACE_MAP:
            num = _inquire_int(self.ex_id, cexodus.EX_INQ_FACE_MAP)
        elif self.entity_type == EntityType.GLOBAL:
            num = 1
        elif self.entity_type == EntityType.NODAL:
            num = 1
        else:
            _raise_invalid_entity_type(ALL_ENTITY_TYPES)
        return num

    def name(self, entity_id: int) -> str:
        """
        Read the name of one entity in an association or collection.

        Uses an internal temporary allocated buffer of up to size ``cexodus.EX_INQ_MAX_READ_NAME_LENGTH + 1`` bytes.
        
        Args:
            entity_id: the integer ID of an entity in the association or collection.
            
        Returns:
            Name of the :class:`Set`, :class:`Block`, or :class:`Map` object, given by the ID, or an empty string if 
            this is a :class:`Global` or :class:`Nodal`.
        """
        cdef char* name_ptr = NULL
        name = ''
        if self.entity_type != EntityType.GLOBAL and self.entity_type != EntityType.NODAL:
            len_str = _inquire_int(self.ex_id, cexodus.EX_INQ_MAX_READ_NAME_LENGTH)
            db_name_size = _inquire_int(self.ex_id, cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
            if db_name_size < len_str:
                len_str = db_name_size
            try:
                name_ptr = <char *> _unaligned_allocate(sizeof(char) * (len_str+1))
                if 0 != cexodus.ex_get_name(self.ex_id, self.entity_type.value, entity_id, &name_ptr[0]):
                    _raise_io_error()
                name = _to_unicode(name_ptr)
            finally:
                if name_ptr != NULL:
                    # logger.debug(__('    free({})', hex(<uintptr_t> name_ptr)))
                    free(name_ptr)
        return name

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def names(self) -> List[str]:
        """
        Read all the names of the entities in this association or collection.

        Uses an internal temporary allocated buffer of up to size 
        ``num_objects * (8 + cexodus.EX_INQ_MAX_READ_NAME_LENGTH + 1)`` bytes.

        Returns:
            List of names of all the :class:`Set`, :class:`Block`, or :class:`Map` objects in this collection, or an 
            empty list this is a :class:`Global` or :class:`Nodal` or 
            
            :obj:`None` if this is a :class:`Nodal` or :class:`Global`
        """
        if self.entity_type == EntityType.NODAL or self.entity_type == EntityType.GLOBAL:
            return None
        len_str = _get_db_max_read_name_length(self.ex_id)
        cdef int64_t num_objects = self.num_entity()

        def names_func(uintptr_t str_ptr):
            return cexodus.ex_get_names(self.ex_id, self.entity_type.value, <char **>str_ptr)

        return _get_strings(len_str, num_objects, names_func)

    def sum_entries(self) -> int:
        """
        Get the total count of all associated entries summed over all collections of entities of this type
        in the database.

        Returns:
            The count of *any* and *all* entries related to entities of this type.
            
        ============================ ==================================================================
        ``self.entity_type``         Return value
        ============================ ==================================================================
        :attr:`.EntityType.GLOBAL`   1 (always)
        :attr:`.EntityType.NODAL`    1 (always)
        :attr:`.EntityType.NODE_MAP` sum of count of nodes in all node maps 
        :attr:`.EntityType.EDGE_MAP` sum of count of edges in all edge maps
        :attr:`.EntityType.FACE_MAP` sum of count of faces in all face maps
        :attr:`.EntityType.ELEM_MAP` sum of count of elements in all element maps
        :attr:`.EntityType.NODE_SET` sum of count of nodes in all node sets
        :attr:`.EntityType.EDGE_SET` sum of count of edges in all edge sets
        :attr:`.EntityType.FACE_SET` sum of count of faces in all face sets
        :attr:`.EntityType.ELEM_SET` sum of count of elements in all element sets
        :attr:`.EntityType.SIDE_SET` (sum of count of elements, sum of count of nodes) in all side sets
        ============================ ==================================================================
        """
        cdef int64_t sum_entries = 0
        cdef int64_t num_nodes = 0
        if self.entity_type == EntityType.NODAL or self.entity_type == EntityType.GLOBAL:
            sum_entries = 1
        elif self.entity_type == EntityType.NODE_SET:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_NS_NODE_LEN)
        elif self.entity_type == EntityType.SIDE_SET:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_SS_ELEM_LEN)
            num_node = _inquire_int(self.ex_id, cexodus.EX_INQ_SS_NODE_LEN)
            return sum_entries, num_node
        elif self.entity_type == EntityType.NODE_MAP:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES) # returns regardless of map existing
        elif self.entity_type == EntityType.EDGE_MAP:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE)
        elif self.entity_type == EntityType.FACE_MAP:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_FACE)
        elif self.entity_type == EntityType.ELEM_MAP:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM)
        elif self.entity_type == EntityType.EDGE_SET:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_ES_LEN)
        elif self.entity_type == EntityType.FACE_SET:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_FS_LEN)
        elif self.entity_type == EntityType.ELEM_SET:
            sum_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_ELS_LEN)
        else:
            _raise_invalid_entity_type(ALL_ENTITY_TYPES)
        return sum_entries


class EntityDictionaryWithProperty(BaseEntityCollection, collections.MutableMapping):

    def __init__(self, database_id=-1, entity_type=None, *args, **kwargs):
        """
        Initialize a type of entity list which supports properties.
        
        If a property is associated with an entity, then the property is associated with *every* entity in that type 
        of collection in the :term:`database`.

        :term:`properties` can be associated with all EntityType except :attr:`EntityType.GLOBAL` and 
        :attr:`EntityType.NODAL.

        On creation, the dictionary is initialized with keys that are the IDs of Blocks, Maps, Sets,
        from the database and the values are initialized to None.

        Subclasses are responsible for implementing ``__getitem__``, ``__setitem__``, ``__delitem__``,
        ``__iter__``, and ``__len__``.
        
        Args:
            database_id: associated open database id
            entity_type: valid type of entity
        """
        # logger.debug(__('EntityDictionaryWithProperty {}', self.__class__.__mro__))
        super(EntityDictionaryWithProperty, self).__init__(database_id=database_id, entity_type=entity_type,
                                                           *args, **kwargs)
        if entity_type == EntityType.NODAL or entity_type == EntityType.GLOBAL:
            raise InternalError('entity_type NODAL or GLOBAL is not an EntityDictionaryWithProperty')
        else:
            # initialize dictionary with (entity ID, None) key, value pairs
            # logger.debug('EntityDictionaryWithProperty.__init__ calling _entity_ids()')
            self._store = collections.OrderedDict.fromkeys(self._entity_ids(), None)

    def __getitem__(self, key):
        value = self._store.get(key)  # _store has all the existing keys obtained at time of constructor
        if value is None:
            # this could be either because key does not exist in database, or just that value not created yet
            if key not in self._store:
                # logger.debug(__('EntityDictionaryWithProperty.__getitem__ key = {} not in self._store', key))
                if not key in self._entity_ids():
                    raise EntityKeyError(f'No {self.entity_type.name} with the id {key} exists in the database.')
                # we are OK the key was added to the database file since our constructor __init__
            # call derived class to create the value since it does not yet exist
            value = self._getitem_derived(key)
            # and put the value in the dictionary
            self._store[key] = value
        return value

    @abstractmethod
    def _getitem_derived(self, key):
        """
        Return the appropriate Block, Map, or Set object corresponding to the given ID.
        
        Derived class must implement.
        
        Args:
            key (int): ID of an entity_type 

        Returns:
            value - appropriate Block, Map or Set object
        """
        raise NotYetImplemented

    @abstractmethod
    def __setitem__(self, key, value):
        raise NotYetImplemented('Writing to database not yet supported.')
        #self._store[key] = value

    @abstractmethod
    def __delitem__(self, key):
        raise NotYetImplemented('Writing to database not yet supported.')
        #del self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return self._store.__len__()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def _entity_ids(self):
        """
        Read the list of positive entity IDs in the collection.

        Uses an internal temporary allocated buffer of up to size ``num_entity * 8`` bytes.
        
        Returns:
            entity_ids: list of positive integer IDs.
        """
        cdef int64_t num_entity
        cdef int64_t * entity_ids_ptr64 = NULL
        cdef int32_t * entity_ids_ptr32 = NULL
        ids = []
        if self.entity_type == EntityType.NODAL or self.entity_type == EntityType.GLOBAL:
            raise InternalError('NODAL or GLOBAL entity_type does not have IDs.')
        num_entity = self._num_entity()
        # logger.debug(__('_entity_ids {}', num_entity))
        if num_entity > 0:
            ids = [-1] * num_entity # fill with -1
            if _is_ids_int64(self.ex_id):
                try:
                    entity_ids_ptr64 = <int64_t*> _allocate(sizeof(int64_t) * num_entity)
                    if 0 != cexodus.ex_get_ids(self.ex_id, self.entity_type.value, entity_ids_ptr64):
                        _raise_io_error()
                    for i in range(num_entity):
                        ids[i] = entity_ids_ptr64[i]
                finally:
                    if entity_ids_ptr64 != NULL:
                        # logger.debug(__('    free({})', hex(<uintptr_t> entity_ids_ptr64)))
                        free(entity_ids_ptr64)
            else:
                try:
                    entity_ids_ptr32 = <int32_t*> _allocate(sizeof(int32_t) * num_entity)
                    if 0 != cexodus.ex_get_ids(self.ex_id, self.entity_type.value, entity_ids_ptr32):
                        _raise_io_error()
                    for i in range(num_entity):
                        ids[i] = entity_ids_ptr32[i]
                finally:
                    if entity_ids_ptr32 != NULL:
                        # logger.debug(__('    free({})', hex(<uintptr_t> entity_ids_ptr32)))
                        free(entity_ids_ptr32)
        return ids

    def name_ids(self) -> Dict[str, int]:
        """
        A dictionary of the mapping of entity name to entity id for every entity in this group.
        
        This method is a convenience for direct lookup and indexing entity ID's only by their name, which is often
        not needed. Using an iterator through he keys and items of this class is probably more useful.
        
        Example:

            >>> with DatabaseFile('myExodusFile.exo') as e:
            >>>     element_blocks = e.element_blocks
            >>>     name_ids = element_blocks.name_ids()
            >>>     block_id = name_ids['block_1']
            >>>     block = element_blocks[block_id]
            >>>     print(block)

            One of the element blocks is accessed by string name "block_1", which may be different depending on the
            individual ExodusII database file.
        
        Returns:
            a dictionary of `name` keys to `entity_id` values
        """
        # logger.debug('EntityDictionaryWithProperty.name_ids() calling _entity_ids()')
        return collections.OrderedDict(zip(self.names(), self._entity_ids()))

    def num_properties(self) -> int:
        """
        The number of properties for this type of collection entity.
        
        If an optional :term:`properties` is associated with an entity, it is associated with *every* entity of that 
        type of collection in the :term:`database`.
        
        Returns:
            number of integer properties
        """
        return cexodus.ex_get_num_props(self.ex_id, self.entity_type)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def property_names(self) -> List[str]:
        """
        Read all the property names present on this type of entity in the database.

        Uses an internal temporary allocated buffer of up to size 
        ``num_entity * (8 + cexodus.EX_INQ_MAX_READ_NAME_LENGTH + 1)`` bytes.
           
        Returns:
            a list of names
        """
        len_str = _get_db_max_read_name_length(self.ex_id)
        cdef int64_t num_entity = self.num_entity()

        def property_names_func(uintptr_t str_ptr):
            return cexodus.ex_get_prop_names(self.ex_id, self.entity_type.value, <char **>str_ptr)

        return _get_strings(len_str, num_entity, property_names_func)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def property(self, property_name: str) -> List[int]:
        """
        Read the values of an integer property on all entities of this type in the database.

        Uses an internal temporary allocated buffer of up to size ``num_properties * 8`` bytes.

        Args:
            property_name: name of the property
            
        Returns:    
            list of integer property values
        """
        cdef int num_properties = cexodus.ex_get_num_props(self.ex_id, self.entity_type.value)
        if num_properties < 1:
            return None
        py_byte_string = property_name.encode('UTF-8')
        cdef char* c_string = py_byte_string
        cdef int * val_ptr = NULL
        properties = []
        try:
            val_ptr = <int*> _unaligned_allocate(sizeof(int) * num_properties)
            if 0 != cexodus.ex_get_prop_array(self.ex_id, self.entity_type.value, c_string, val_ptr):
                _raise_io_error()
            for i in range(num_properties):
                properties.append(val_ptr[i])
        finally:
            if val_ptr != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> val_ptr)))
                free(val_ptr)
        return properties


class EntityCollectionWithVariable(BaseEntityCollection):

    def __init__(self, database_id=-1, entity_type=None, *args, **kwargs):
        """
        Initialize an entity list that has associated variables.

        The entity_type of ELEM_MAP, FACE_MAP, EDGE_MAP, and NODE_MAP do not have variables.

        A :term:`variable` in Exodus is a scalar only. Multiple variables are required to make up either a vector or
        tensor field.
        
        Args:
            database_id: associated open database id
            entity_type: valid type of entity, one of GLOBAL, NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
                FACE_SET, ELEM_BLOCK, ELEM_SET, SIDE_SET 
        """
        # logger.debug(__('EntityCollectionWithVariable {}', self.__class__.__mro__))
        super(EntityCollectionWithVariable, self).__init__(database_id=database_id, entity_type=entity_type,
                                                           *args, **kwargs)
        self._fields_dictionary = None

    @property
    def fields(self) -> Dict[str, Field]:
        """
        A dictionary of fields associated with this collection of entities.
         
        This object is a mapping of field names to :class:`Field` objects that are active on the collection of entities. 

        Each :term:`field` is a grouping of ExodusII :term:`variable` by common name prefix; the suffix of the
        :term:`variable` name (whatever follows the last underscore '_' in the name) becomes a component name.

        Example:

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     print(f'number of nodal fields: {len(e.nodal.fields)}')
            >>>     print(f'element field names: {e.globals.fields.keys()}')
            >>>     for f in e.element_blocks.fields.values():
            >>>         print(f)
            
        Returns:
            A :class:`collections.OrderedDict`, where `keys` of the dictionary are field names, and the 
            `values` of the dictionary are instances of :class:`.Field`.
        """
        if self._fields_dictionary is None:
            self._fields_dictionary = Fields(self)
        return self._fields_dictionary

    def num_variables(self) -> int:
        """
        The number of variables associated with all the entities of this type in the database.

        Returns:
            Number of variables active on this ``entity_type``.
        """
        cdef int num_vars = 0
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type.value, &num_vars):
            _raise_io_error()
        return num_vars

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def variable_names(self) -> List[str]:
        """
        Read the list of variable names associated with this type of entity.

        Uses an internal temporary allocated buffer of up to size
        ``num_variables * (8 + cexodus.MAX_STR_LENGTH + 1)`` bytes.

        Returns:
             List of variable names.
        """
        cdef int len_str = cexodus.MAX_STR_LENGTH + 1
        cdef int num_vars = 0
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type.value, &num_vars):
            _raise_io_error()

        def variable_names_func(uintptr_t str_ptr):
            return cexodus.ex_get_variable_names(self.ex_id, self.entity_type.value, num_vars, <char **>str_ptr)

        return _get_strings(len_str, num_vars, variable_names_func)

    def variable_name(self, variable_index: int) -> str:
        """
        The name of a variable associated with this entity_type.

        A variable is scalar only. Multiple variables make up a vector or tensor field.
        
        Args:
            variable_index: variable index in the range ``[0,num_variables - 1]``
            
        Returns:
            Name of the variable with the given index.
        """
        cdef char var_name[cexodus.MAX_STR_LENGTH+1]
        if 0 != cexodus.ex_get_variable_name(self.ex_id, self.entity_type.value, variable_index+1, &var_name[0]):
            _raise_io_error()
        variable_name = str(_to_unicode(&var_name[0]))
        return variable_name

    @property
    def variables(self) -> Dict[str, int]:
        """
        A dictionary of scalar variable names and indices associated with this collection of entities.
         
        This object is a mapping of scalar variable name -> index, that are active on the collection of entities. 

        Prefer using :meth:`fields` which is usually more convenient unless you only want one component, because 
        each :term:`variable` is only a scalar, and maybe a component value of a :term:`field` vector or tensor.

        Example:
            
            Get the Z component of a force variable at a single node at all times.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     node_vars = e.nodal.variables 
            >>>     fz = e.nodal.variable_at_times(node_vars['AForceZ'], 1, 0, e.globals.times().size)
            >>>     print(fz)
            
        Returns:
            A :class:`collections.OrderedDict`, where `keys` of the dictionary are variable names, and the 
            `values` of the dictionary are the integer index, or ID of the variable.
        """
        variable_dictionary = collections.OrderedDict()
        for index, name in enumerate(self.variable_names()):
            variable_dictionary[name] = index
        return variable_dictionary

    def is_variable_active_all(self) -> numpy.ndarray:
        """
        Read a table of indicating whether each of the variables on this type of entities is active on each entity. 
        
        A value of one in a table entry indicate a :term:`variable` exists, while a value of zero indicates the 
        variable does not exist on an entity.
        
        Returns:
            Array of type :obj:`numpy.int32`, with shape ``(num_entity, num_vars)``, with one or zero entries 
            indicating the variable is active on an entity.
            
            :obj:`None` if no variables are present on any entities of this type.
        """
        cdef int num_entity = self.num_entity()
        if num_entity < 1:
            return None
        cdef int num_var = self.num_variables()
        if num_var < 1:
            return None
        array_shape = (num_entity, num_var)
        cdef numpy.ndarray[numpy.int32_t, ndim=2] var_table = util.empty_aligned(array_shape, dtype=numpy.int32)
        cdef int *var_ptr = <int *> var_table.data
        if 0 != cexodus.ex_get_truth_table(self.ex_id, self.entity_type.value, num_entity, num_var, var_ptr):
            _raise_io_error()
        return var_table


class Fields(collections.MutableMapping):

    def __init__(self, parent_entity_collection: EntityCollectionWithVariable):
        """
        An accessor for fields, acting as an ordered dictionary with some extra member functions.

        For a database, this is an accessor for the fields on all the members of one type of entity in a
        database. The set of entities are one of the globals, the nodal, or one of
        element blocks, face blocks, or edge blocks; or side sets, node sets, element sets, face sets, or
        edge sets.

        The keys of the dictionary are the base names of fields.
        
        Args:
            parent_entity_collection: if reading/modifying an existing database, the owning collection of
                entities in an Exodus database, one of a :class:`.Global`, :class:`.Nodal`, 
                :class:`.Blocks`, or :class:`.Sets` object
                
        Raises:
            InvalidEntityType: if the parent collection of entities does not support fields
        """
        if not isinstance(parent_entity_collection, EntityCollectionWithVariable):
            _raise_invalid_entity_type('parent entity collection with support for fields')
        self.parent_entity_collection = weakref.proxy(parent_entity_collection)  # avoid cyclic references
        self._store = None

    def __getitem__(self, key):
        if self._store is None:
            self._store = self.__get_field_info_dictionary()
        return self._store[key]

    def __setitem__(self, key, value):
        #self._store[key] = value
        raise NotYetImplemented('Writing to database is not yet supported in this version of the module.')

    def __delitem__(self, key):
        #del self._store[key]
        raise NotYetImplemented('Writing to database is not yet supported in this version of the module.')

    def __iter__(self):
        if self._store is None:
            self._store = self.__get_field_info_dictionary()
        return iter(self._store)

    def __len__(self):
        if self._store is None:
            self._store = self.__get_field_info_dictionary()
        return self._store.__len__()

    def __get_field_info_dictionary(self):
        # Private member function constructing the store from an existing database
        #
        global _field_subscripts_union
        var_names = self.parent_entity_collection.variable_names()
        # build an intermediate dictionary of field names and list of pairs of
        # components and scalar variable ids. The result is an ordered dictionary of (k, v), where k is the field name,
        # v is a list of tuple pairs (i,j), where i is the component and j is the corresponding index for that component
        # in the list of scalar variables. For example,
        #   k, v = 'stress', [('xx', 7), ('yy', 8), ('zz', 9), ('xy', 10), ('yz', 11), ('zx', 12)]
        #   k, v = 'external_energy', [('', 6)]
        d = collections.OrderedDict()
        for index, name in enumerate(var_names):
            k,s,v = name.rpartition('_')                                # split string into key, _, subscript
            if k == '':
                d.setdefault(v, []).append(('', index))                 # no subscript here, just add it
            else:
                if v in _field_subscripts_union or _is_integer(v):      # look for an int, or x, y, z, xx, xy, etc.
                    d.setdefault(k, []).append((v, index))              # add the subscript to that key
                else:
                    if k in d:                          # we found a previous key, but subscript not recognized
                        for x in d[k]:                  # for each of the previous existing subscripts
                            d.setdefault(k + '_' + x[0], []).append(('',x[1]))  # add previous ones back
                        del(d[k])                       # delete the previous key and all the values
                    d.setdefault(name, []).append(('', index))  # add the name including subscript
        # we make a new dictionary of Field objects, from our intermediate data
        field_info_dictionary = collections.OrderedDict().fromkeys(d.keys())
        for k, v in d.items():
            components, variables = zip(*v)
            field_info_dictionary[k] = Field(k, components, variables, self)
        return field_info_dictionary

    def is_field_active_all(self) -> numpy.ndarray:
        """
        Read a table indicating whether each of the fields on this type of entities is active on each entity. 
        
        A value of one in a table entry indicate :term:`field` values exist, while a value of zero indicates the 
        field does not exist on an entity.
        
        Returns:
            Array of type numpy.int32, with shape ``(num_entity, num_fields)``, with one or zero entries indicating the
            field is active on an entity.
          
            :obj:`None` if no fields are present on any entities of this type.
            
        Raises:
            InactiveComponent: if there are inconsistencies if the related component :term:`variable` entries in the 
                ExodusII database corresponding to a :term:`field`.
        """
        fields_dict = self._store
        if fields_dict is None:
            return None
        parent_entity_collection = self.parent_entity_collection
        cdef int num_entity = parent_entity_collection.num_entity()
        if num_entity < 1:
            return None
        cdef int num_var = parent_entity_collection.num_variables()
        if num_var < 1:
            return None
        array_shape = (num_entity, num_var)
        cdef numpy.ndarray[numpy.int32_t, ndim=2] var_table = util.empty_aligned(array_shape, dtype=numpy.int32)
        cdef int *var_ptr = <int *> var_table.data
        if 0 != cexodus.ex_get_truth_table(parent_entity_collection.ex_id, parent_entity_collection.entity_type.value,
                                           num_entity, num_var, var_ptr):
            _raise_io_error()
        num_fields = len(fields_dict)
        array_shape = (num_entity, num_fields)
        cdef numpy.ndarray[numpy.int32_t, ndim=2] fields_table = util.empty_aligned(array_shape, dtype=numpy.int32)
        for entity in range(num_entity):
            for index, k, in enumerate(fields_dict):
                v = fields_dict[k]
                variables = v.variables
                if var_table[entity, variables[0]] == 1:
                    for i in range(1,len(variables)):  # ensure all the any other component variables are also active
                        if var_table[entity, variables[i]] != 1:
                            raise InactiveComponent('Exodus database has a field with an inactive component:\n'
                                f'  database {self.ex_id} field {k} component {i} of {len(v)} '
                                f'(entity id {entity} of type {self.entity_type.name}).')
                    fields_table[entity, index] = 1
        return fields_table


class Global(EntityCollectionWithVariable, EntityWithVariable):

    def __init__(self, database_id=-1, *args, **kwargs):
        """
        Initialize the global entity on a Database which acts as an BaseEntityCollection with only a single Entity
        member, which is itself.
        
        The number of entries in the :class:`Global` entity is also equal to 1.

        Args:
            database_id(int): ID of the Exodus database

        Attributes:
            ex_id (int): the ExodusII database ID
            entity_type(EntityType): always EntityType.GLOBAL, the type of entity held by this collection
        """
        # logger.debug(__('Global {}', self.__class__.__mro__))
        super(Global, self).__init__(database_id=database_id, entity_type=EntityType.GLOBAL, entity_id=-1, name=None,
                                     parent_entity_collection=self, *args, **kwargs)

    def dimension(self) -> int:
        """
        The spatial dimension of the coordinates in the database.
        
        Returns:
            integer dimension
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_DIM)

    def has_id_map(self, map_type: EntityType) -> bool:
        """
        Whether an ID map of the given type exists on the database.
        
        Args:
            map_type: type of map, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP
            
        Returns:
            True or False
            
        Raises:
            InvalidEntityType: if the map_type is invalid
        """
        cdef extern from "netcdf.h":
            int nc_inq_varid(int ncid, char *name, int *varidp)
        cdef int exoid
        cdef int mapid
        cdef char *vmap
        cdef bytes py_bytes
        if map_type == EntityType.NODE_MAP:
            py_bytes = "node_num_map".encode()
        elif map_type == EntityType.EDGE_MAP:
            py_bytes = "edge_num_map".encode()
        elif map_type == EntityType.FACE_MAP:
            py_bytes = "face_num_map".encode()
        elif map_type == EntityType.ELEM_MAP:
            py_bytes = "elem_num_map".encode()
        else:
            _raise_invalid_entity_type(MAP_ENTITY_TYPES)
        exoid = self.ex_id
        vmap = py_bytes
        # look for the netcdf map variable and check for error condition
        if nc_inq_varid(exoid, vmap, &mapid) != 0:
            return False
        else:
            return True

    def id_map(self, map_type: EntityType) -> numpy.ndarray:
        """
        Read the array of global integer IDs for all entries of the given type.

        This is the single map of the IDs of elements, nodes, edges, or faces (completely separate from the other list
        of maps. If no ID map exists, the default enumeration is returned, for example, [0, num_elements] when
        map_type == ELEM_MAP.
        
        Args:
            map_type: type of map, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP
        
        Returns:
            The map entries array with shape ``(num_elements,)``, ``(num_nodes,)``, ``(num_edges,)``, or 
            ``(num_faces,)``; or 
            
            :obj:`None` if no entries exist in the database.
            
        Raises:
            InvalidEntityType: if the `map_type` is invalid
        """
        cdef int64_t num_entries = 0
        cdef numpy.ndarray[numpy.int64_t, ndim=1] map_entries64
        cdef int64_t *map_ptr64
        cdef numpy.ndarray[numpy.int32_t, ndim=1] map_entries32
        cdef int *map_ptr32
        if map_type == EntityType.NODE_MAP:
            num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)
        elif map_type == EntityType.ELEM_MAP:
            num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM)
        elif map_type == EntityType.EDGE_MAP:
            num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE)
        elif map_type == EntityType.FACE_MAP:
            num_entries = _inquire_int(self.ex_id, cexodus.EX_INQ_FACE)
        else:
            _raise_invalid_entity_type(MAP_ENTITY_TYPES)
        if num_entries < 1:
            return None
        if _is_maps_int64(self.ex_id):
            map_entries64 = util.empty_aligned(num_entries, dtype=numpy.int64)
            map_ptr64 = <int64_t *> map_entries64.data
            if 0 != cexodus.ex_get_id_map(self.ex_id, map_type, map_ptr64):
                _raise_io_error()
            return map_entries64
        else:
            map_entries32 = util.empty_aligned(num_entries, dtype=numpy.int32)
            map_ptr32 = <int *> map_entries32.data
            if 0 != cexodus.ex_get_id_map(self.ex_id, map_type, map_ptr32):
                _raise_io_error()
            return map_entries32

    def num_elements(self) -> int:
        """
        The global number of element entries in the database.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_ELEM)

    def num_edges(self) -> int:
        """
        The global number of edge entries in the database.

        Most often the number of edges is zero even when num_nodes or num_elem is positive.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_EDGE)

    def num_faces(self) -> int:
        """
        The global number of edge entries in the database.

        Most often the number of faces is zero even when num_nodes or num_elem is positive.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_FACE)

    def num_nodes(self) -> int:
        """
        The global number of node entries in the database.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)

    def num_times(self) -> int:
        """
        The number of time steps on the database.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_TIME)

    def times(self) -> numpy.ndarray:
        """
        Read all time values.

        Returns:
            array of type :obj:`numpy.double` with shape ``(num_times)``, or 
            
            :obj:`None` if no time steps are present.
        """
        cdef int64_t array_shape = _inquire_int(self.ex_id, cexodus.EX_INQ_TIME)
        if array_shape < 1:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] times = util.empty_aligned(array_shape, dtype=numpy.double)
        cdef double *times_ptr = <double *> times.data
        if 0 != cexodus.ex_get_all_times(self.ex_id, times_ptr):
            _raise_io_error()
        return times

    def variables_at_time(self, step: int) -> numpy.ndarray:
        """
        Read an array of values of all the global variables at a single time step.
        
        Args:
            step: the time step index, which the variable values are desired, the first time step is 0.
            
        Returns:
            array of :obj:`numpy.double` with shape ``(num_variables)``, or 
            
            :obj:`None` if there are no variables.
        """
        cdef int64_t var_length = self.num_variables()
        if var_length < 1:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] var = util.empty_aligned(var_length, dtype=numpy.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_glob_vars(self.ex_id, step, var_length, var_ptr):
            _raise_io_error()
        return var


# noinspection SpellCheckingInspection
class Nodal(EntityCollectionWithVariable, EntityWithAttribute):

    def __init__(self, database_id=-1, *args, **kwargs):
        """
        Initialize the Nodal entity on a Database which acts a both an entity, and a list of entities, whose sole
        member is itself.
        
        Args:
            database_id: ID of the Exodus database.
        """
        # logger.debug(__('Nodal {}', self.__class__.__mro__))
        super(Nodal, self).__init__(database_id=database_id, entity_type=EntityType.NODAL, entity_id=-1, name=None,
                                     parent_entity_collection=self, *args, **kwargs)
        self._displacement_name = None

    def num_coordinate_frames(self) -> int:
        """
        Get the number of coordinate frames
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_COORD_FRAMES)

    def coordinate_names(self) -> List[str]:
        """
        Read names of the coordinate axes.

        Uses an internal temporary allocated buffer of up to size ``num_dimensions * cexodus.MAX_STR_LENGTH``.
        
        Returns:
            list of strings of the names of the coordinate axes
            
        Raises:
            InvalidSpatialDimension: if the stored dimension is not in the range ``[0,3]``.
        """
        cdef int64_t ndim = _inquire_int(self.ex_id, cexodus.EX_INQ_DIM)
        cdef char *coord_names_ptr[3];
        if ndim < 0 or ndim > 3:
            raise InvalidSpatialDimension(f'unexpected spatial dimension = {ndim}')
        coord_names = []
        try:
            for i in range(ndim):
                coord_names_ptr[i] = <char*> _unaligned_allocate(sizeof(char) * cexodus.MAX_STR_LENGTH)
            if 0 != cexodus.ex_get_coord_names(self.ex_id, coord_names_ptr):
                _raise_io_error()
            for i in range(ndim):
                coord_names.append(_to_unicode(coord_names_ptr[i]))
        finally:
            for i in range(ndim):
                free(coord_names_ptr[i])
        return coord_names

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def coordinates(self) -> numpy.ndarray:
        """
        Read all the global node coordinates.

        Uses an internal temporary allocated buffer of size ``num_nodes * 8`` bytes.

        Returns:
            array of :obj:`numpy.double` with shape ``(num_nodes, num_dim)``.
            
        Raises:
            InvalidSpatialDimension: if the dimension is not in the range ``[0,3]``
        """
        cdef int64_t i = 0
        cdef int64_t j = 0
        array_shape = None
        cdef int64_t nnodes = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)
        cdef int64_t ndim = _inquire_int(self.ex_id, cexodus.EX_INQ_DIM)
        if ndim == 1:
            array_shape = nnodes
        elif ndim == 2:
            array_shape = (nnodes,2)
        elif ndim == 3:
            array_shape = (nnodes,3)
        else:
            raise InvalidSpatialDimension(f'Unexpected number of dimensions = {ndim}')
        cdef numpy.ndarray[numpy.double_t, ndim=2] coords = util.empty_aligned(array_shape, dtype=numpy.double)
        cdef double * coords_ptr = <double *> coords.data
        cdef double * coords_buffer = NULL
        try:
            coords_buffer = <double*> _allocate(sizeof(double) * nnodes)

            # x-coord
            if 0 != cexodus.ex_get_coord(self.ex_id, coords_buffer, NULL, NULL):
                _raise_io_error()
            with nogil:
                cexodus.ex_aligned_copy_stride(coords_buffer, nnodes, coords_ptr, ndim)

            if ndim > 1:  # y-coord
                if 0 != cexodus.ex_get_coord(self.ex_id, NULL, coords_buffer, NULL):
                    _raise_io_error()
                with nogil:
                    cexodus.ex_aligned_copy_stride(coords_buffer, nnodes, coords_ptr + 1, ndim)

                if ndim > 2:  # z-coord
                    if 0 != cexodus.ex_get_coord(self.ex_id, NULL, NULL, coords_buffer):
                        _raise_io_error()
                    with nogil:
                        cexodus.ex_aligned_copy_stride(coords_buffer, nnodes, coords_ptr + 2, ndim)
        finally:
            if coords_buffer != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> coords_buffer)))
                free(coords_buffer)
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def coordinates_partial(self, start: int, stop: int) -> numpy.ndarray:
        """
        Read a range of the global node coordinates in an array.

        Uses an internal temporary allocated buffer of size ``(stop - start) * 8`` bytes.
        
        Args:
            start: first desired node, in the range ``[0, num_nodes-1]``.
            stop:  one past the end of the last desired node, in the range ``[start+1, num_nodes]``
            
        Returns:
            array of :obj:`numpy.double` with shape ``(stop - start, num_dim)``
            
        Raises:
            InvalidSpatialDimension: if the dimension is not in the range ``[0,3]``
        """
        cdef int64_t i = 0
        cdef int64_t j = 0
        array_shape = None
        cdef int64_t nnodes = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)
        cdef int64_t ndim = _inquire_int(self.ex_id, cexodus.EX_INQ_DIM)
        cdef int64_t start_node = start + 1
        cdef int64_t len_nodes = stop - start
        if len_nodes < 1:
            return None
        if ndim == 1:
            array_shape = len_nodes
        elif ndim == 2:
            array_shape = (len_nodes,2)
        elif ndim == 3:
            array_shape = (len_nodes,3)
        else:
            raise InvalidSpatialDimension(f'Unexpected number of dimensions = {ndim}')
        cdef numpy.ndarray[numpy.double_t, ndim=2] coords = util.empty_aligned(array_shape, dtype=numpy.double)
        cdef double * coords_ptr = <double *> coords.data
        cdef double * coords_buffer = NULL
        try:
            coords_buffer = <double*> _allocate(sizeof(double) * len_nodes)

            # x-coord
            if 0 != cexodus.ex_get_partial_coord(self.ex_id, start_node, len_nodes, coords_buffer, NULL, NULL):
                _raise_io_error()
            with nogil:
                cexodus.ex_aligned_copy_stride(coords_buffer, len_nodes, coords_ptr, ndim)

            if ndim > 1:  # y-coord
                if 0 != cexodus.ex_get_partial_coord(self.ex_id, start_node, len_nodes, NULL, coords_buffer, NULL):
                    _raise_io_error()
                with nogil:
                    cexodus.ex_aligned_copy_stride(coords_buffer, len_nodes, coords_ptr + 1, ndim)

                if ndim > 2:  # z-coord
                    if 0 != cexodus.ex_get_partial_coord(self.ex_id, start_node, len_nodes, NULL, NULL, coords_buffer):
                        _raise_io_error()
                    with nogil:
                        cexodus.ex_aligned_copy_stride(coords_buffer, len_nodes, coords_ptr + 2, ndim)
        finally:
            if coords_buffer != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> coords_buffer)))
                free(coords_buffer)
        return coords

    def coordinates_local(self, local: LocalConnectivity) -> numpy.ndarray:
        """
        Read the coordinates specified by the local connectivity of a :class:`.Block`.

        Examples:

            If you require local coordinates for *more than one* block, it is almost always more efficient to read
            *all* the global coordinates in a single read operation and copy them block by block.

            Let's compute the mean coordinates of every element block in the database.

            We use the :func:`numpy.take` function to index local coordinates from the global coordinates. Use the array
            of global node IDs from the block's local connectivity as indices.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     global_coordinates = e.nodal.coordinates()
            >>>     block_iterator = e.element_blocks.connectivity_local_all()
            >>>     for key, block, local in block_iterator:
            >>>         local_coordinates = global_coordinates.take(local.global_nodes, axis=0)  # use numpy.take()
            >>>         print(numpy.mean(local_coordinates, axis=0))

            Otherwise, if you do not expect to read any other coordinates outside than those of a single block, use this
            function. The read operation uses less memory and is probably faster since it reads only the range of
            data required for the block.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     block = e.element_blocks[7]
            >>>     local_connectivity = block.connectivity_local()
            >>>     local_coordinates = e.nodal.coordinates_local(local_connectivity)
            >>>     print(numpy.mean(local_coordinates, axis=0))

        Args:
            local: structure specifying the range of local connectivity of a block

        Returns:
            a byte aligned coordinates array (numpy.ndarray[numpy.float64, order="C"]) local to a block, where
            shape is (local.global_nodes.size, ndim)
        """
        partial_coordinates = self.coordinates_partial(start=local.min_global, stop=local.max_global+1)
        return util.take(partial_coordinates, local.global_nodes, shift=local.min_global)


    @property
    def displacement_name(self) -> str:
        """
        Name of the nodal field acting as the displacements.

        Returns:
            string name, or 
            
            :obj:`None` if no suitable field for displacements was found
        """
        if self._displacement_name is None:
            ndim = -1
            for c in ('displacement', 'displace', 'displ', 'disp'):
                if c in self.fields:
                    if ndim == -1:
                        ndim = _inquire_int(self.ex_id, cexodus.EX_INQ_DIM)
                    if len(self.fields[c].components) == ndim:
                        self._displacement_name = c
        return self._displacement_name

    def coordinates_displaced(self, time_step: int) -> numpy.ndarray:
        """
        The sum of the coordinates and displacements at a time step.
        
        This is the coordinate vector (returned by :meth:`coordinates`) added to the displacement vector 
        (the :class:`.FieldArray` corresponding to :meth:`displacement_name`) at all nodes at a single time step.

        Uses an internal temporary allocated buffer for coordinate and displacement arrays of size
        ``2 * num_nodes * num_dim * 8`` bytes.
        
        Args:
            time_step: the time step index, at which the nodal displaced coordinate values are desired;             
            the first time step is 0; a value of -1 will return the values at the last time step.
        
        Returns:
            array with shape ``(num_nodes, num_dim)``, or 
            
            :obj:`None` if no displacements exist
        """
        if self.displacement_name is None:
            return None
        else:
            coordinates = self.coordinates()
            displacement_info = self.fields[self.displacement_name]
            displacements = self.field(displacement_info, time_step)
            return numpy.add(coordinates, displacements)

    def coordinates_displaced_partial(self, time_step, start_entry, stop_entry) -> numpy.ndarray:
        """
        The sum of the coordinates and displacements on a range of nodes.
        
        This is the coordinate vector (returned by :meth:`coordinates`) added to the displacement vector 
        (the :class:`.FieldArray` corresponding to :meth:`displacement_name`) on the given range nodes at a single 
        time step.

        Uses an internal temporary allocated buffer for coordinate and displacement arrays of size
        ``2 * (stop_entry - start_entry) * num_dim * 8`` bytes.
        
        Args:
            time_step: the time step index, at which the nodal displaced coordinate values are desired,
                the first time step is 0; a value of -1 will return the values at the last time step.
            start_entry: first desired node, in the range ``[0, num_nodes-1].``
            stop_entry:  one past the end of the last desired node, in the range ``[start+1, num_nodes]``
            
        Returns:
             array of :obj:`numpy.double` with shape ``(stop_entry - start_entry, num_dim)``, or 
             
             :obj:`None` if no displacements exist.
        """
        if self.displacement_name is None:
            return None
        else:
            coordinates = self.coordinates_partial(time_step, start_entry, stop_entry)
            displacement_info = self.fields[self.displacement_name]
            displacements = self.field_partial(displacement_info, time_step, start_entry, stop_entry)
            return numpy.add(coordinates, displacements)

    def coordinates_displaced_local(self, time_step, local: LocalConnectivity) -> numpy.ndarray:
        """
        The sum of the coordinates and displacements specified by the local connectivity of a :class:`.Block` at a
        time step.

        Examples:

            See :meth:`coordinates_local`.

        Args:
            time_step: the time step index, at which the nodal displaced coordinate values are desired,
                the first time step is 0; a value of -1 will return the values at the last time step.
            local: structure specifying the range of local connectivity of a block

        Returns:
            a byte aligned displaced coordinates array (numpy.ndarray[numpy.float64, order="C"]) local to a block, where
            shape is (local.global_nodes.size, ndim)
        """
        if self.displacement_name is None:
            return None
        else:
            coordinates = self.coordinates_local(local)
            displacement_info = self.fields[self.displacement_name]
            partial_displacements = self.field_partial(displacement_info, time_step,
                                                       local.min_global, local.max_global+1)
            displacements = util.take(partial_displacements, local.global_nodes, shift=local.min_global)
            return numpy.add(coordinates, displacements)

    def field_local(self, field: Field, time_step, local: LocalConnectivity) -> FieldArray:
        """
        The array of field values of a single field at a time step for all nodes specified by the local
        connectivity of a :class:`.Block`.

        Note:

            If you require local fields for *more than one* block, it is almost always more efficient to read
            the global field for all nodes in a single read operation and copy them block by block.

            Use the :func:`numpy.take` function to index local indices from the global nodes. Use the array
            of global node IDs from the block's :class:`.LocalConnectivity`.

        Args:
            field: the desired field (name and component names)
            time_step: the time step index, at which the field component values are desired,
                the first time step is 0; a value of -1 will return the values at the last time step.
            local: structure specifying the range of local connectivity of a block

        Returns:
            a byte aligned array of :obj:`numpy.float64` values with shape
            ``(local.global_nodes.size, len(field.components))``, or

            :obj:`None` if ``num_entries == 0`` or field is not associated with ExodusII variable indices.

        Raises:
            InactiveField: if the given field does not exist on this entity.
        """
        partial_field = self.field_partial(field, time_step, local.min_global, local.max_global+1)
        return util.take(partial_field, local.global_nodes, shift=local.min_global)


# noinspection PyProtectedMember
cdef inline uint32_t _num_connectivity_entries(block, entry_type):
    """
    Internal function, second dimension of the connectivity array (returned by :meth:`connectivity` method) 
    for a given type of entry type.
    
    Returns:
        number of nodes, edges or faces per entry given the entry type
    
    Raises:
        EntryTypeError: if entry_type is not NODE, EDGE, or FACE.
    """
    if entry_type == EntryType.NODE:
        num_conn_entries = block._num_nodes_per_entry
    elif entry_type == EntryType.EDGE:
        num_conn_entries = block._num_edges_per_entry
    elif entry_type == EntryType.FACE:
        num_conn_entries = block._num_faces_per_entry
    else:
        block._raise_entry_type_error()
    return num_conn_entries

cdef inline object _create_connectivity(block, entry_type, shape, bool is_int64, void** data,
                                        bool zero_based, bool no_array):
    """
    Internal function to allocate and read block connectivity data.
    If `no_array=True` the returned pointer must be freed by the caller.
    
    Args:
        block(Block): the block instance
        entry_type(EntryType): type of entry
        length(int): num_entries * num_connected_entries
        is_int64(bool): True if block is using int64 array data on disk
        zero_based(bool): True if connected entries ID's start at zero, False they start at one
        no_array(bool): True causes only pointer to data to be returned with no array 

    Returns:
        array, data (numpy.ndarray, uintptr_t):
        
        * array - an ndarray[shape=(num_entries, num_connected_entries)] or None if argument no_array=True
        * data - a raw C pointer to array data to int32 or int64 data
    
    Raises:
        InvalidEntityError: when entry type is not NODE, EDGE, or FACE.
    """
    cdef int error = -1
    cdef size_t length
    # logger.debug(__('    _create_connectivity array shape = {} is_int64 = {}', shape, is_int64))
    if is_int64:
        array = _create_array_data(shape, numpy.int64, data, no_array)
    else:
        array = _create_array_data(shape, numpy.int32, data, no_array)
    if entry_type == EntryType.NODE:
        error = cexodus.ex_get_conn(block.ex_id, block.entity_type.value, block.entity_id, data[0], NULL, NULL)
    elif entry_type == EntryType.EDGE:
        error = cexodus.ex_get_conn(block.ex_id, block.entity_type.value, block.entity_id, NULL, data[0], NULL)
    elif entry_type == EntryType.FACE:
        error = cexodus.ex_get_conn(block.ex_id, block.entity_type.value, block.entity_id, NULL, NULL, data[0])
    else:
        # noinspection PyProtectedMember
        block._raise_entry_type_error()
    if error != 0:
        _raise_io_error()
    if zero_based:
        length = get_size(shape)
        if is_int64:
            cexodus.ex_aligned_to_zero_based_int64(<int64_t *> data[0], length)
        else:
            cexodus.ex_aligned_to_zero_based_int32(<int32_t *> data[0], length)
    return array


class Blocks(EntityDictionaryWithProperty, EntityCollectionWithVariable):
    
    def __init__(self, database_id=-1, block_type=None, *args, **kwargs):
        """
        The collection of all Block entities of a specific type on a Database.

        Args:
            database_id: ID of the Exodus database.
            block_type: the type of the entity, one of ELEM_BLOCK, EDGE_BLOCK, or FACE_BLOCK
            
        Raises:
            InvalidEntityType: if the `block_type` is invalid.
        """
        # logger.debug(__('Blocks {}', self.__class__.__mro__))
        if block_type not in BLOCK_ENTITY_TYPES:
            _raise_invalid_entity_type(BLOCK_ENTITY_TYPES)
        super(Blocks, self).__init__(database_id=database_id, entity_type=block_type, *args, **kwargs)

    def _getitem_derived(self, key):
        """
        Construct a block from the database given the block ID.
        
        Args:
            key: the ID of the block, one of those returned by :meth:`entity_ids`

        Returns:
            A :class:'Block' initialized with name, topology, sizes of entries and attributes, etc.
            
        Raises:
            EntityKeyError: if the key is not present in the database.
        """
        cdef cexodus.ex_block block
        # logger.debug('Blocks._getitem_derived()')
        block.id = key
        block.type = self.entity_type
        if 0 != cexodus.ex_get_block_param(self.ex_id, &block):
            _raise_io_error()
        ndim = _inquire_int(self.ex_id, cexodus.EX_INQ_DIM)
        topology_name = _disambiguate_topology_name(_to_unicode(block.topology), block.num_nodes_per_entry, ndim)
        name = self.name(key)
        return Block(self.ex_id, self.entity_type, key, name, topology_name, block.num_entry,
                     block.num_nodes_per_entry, block.num_edges_per_entry, block.num_faces_per_entry,
                     parent_entity_collection=self)

    def __setitem__(self, key, value):
        EntityDictionaryWithProperty.__setitem__(self, key, value)

    def __delitem__(self, key):
        EntityDictionaryWithProperty.__delitem__(self, key)

    # noinspection PyProtectedMember
    def connectivity_local_all(self, keys=None, entry_type=EntryType.NODE,
                               compress=False) -> Iterator[Tuple[int, Block, LocalConnectivity]]:
        """
        Construct a generator that can be used to iterate over the local connectivity of each block.

        The next() of the iterator returns a tuple (key, block, local_connectivity).

        If a block has no connected entries, the local_connectivity returned will be :obj:`None`.

        Using this iterator avoids repeatedly allocating, zeroing, and freeing some temporary memory used to
        compute local connectivity arrays. This method is more efficient than calling
        :meth:`.Block.connectivity_local` for each block separately.

        See :meth:`.Block.connectivity_local` for more.

        Args:
            keys(Iterable[int]): an ordered sequence of Block keys to use, which is a subset of all the Blocks
                (as opposed to all the Blocks which is the default)
            entry_type: (default EntryType.NODE) type of entries returned in the connectivity array, one of
                NODE, EDGE, FACE
            compress: if False (default) return a numpy.ndarray, if True return a compressed array in
                :class:`affect.util.CompressedArray`

        Returns:
            connectivity_local_iterator - an iterator over tuples of (key, block, local_connectivity), where key
                is the block ID, block is a :class:`.Block`, and local_connectivity is a :class:`.LocalConnectivity`

        Raises:
            EntityKeyError: if keys is not None and one of the keys is not present in the database.
            InvalidEntityType: if the connect entry_type does not make sense
            MaxLengthExceeded: if the number of unique local nodes in this block exceeds UINT32_MAX, 2^32 - 1,
        """
        cdef uint32_t * global_to_local_ptr = NULL
        cdef void * entry_to_node_global_ptr = NULL
        cdef size_t length, min_global_node, max_global_node
        cdef size_t num_entries, num_conn_entries, num_global_nodes

        # get iterator on the blocks we are accessing
        if keys is None:
            blocks = ((k, block) for k, block in self.items())
        else:
            blocks = ((k, self[k]) for k in keys)  # will eventually raise an EntityKeyError below if k is not in Blocks
        try:
            # allocate our working buffer to be used for all blocks
            num_global_nodes = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)
            global_to_local_ptr = <uint32_t*> _allocate(sizeof(uint32_t) * num_global_nodes)
            # and initial range of possible non-zero entries in buffer
            min_global_node = 0  # initial range of possible non-zero entries in buffer
            max_global_node = num_global_nodes - 1
            for key, block in blocks:
                if block._num_entries < 1:
                    yield key, block, None
                num_entries = block._num_entries
                num_conn_entries = _num_connectivity_entries(block, entry_type)
                if num_conn_entries == 0:
                    yield key, block, None
                local_connectivity = _create_connectivity_local(block,
                                                                entry_type,
                                                                num_entries,
                                                                num_conn_entries,
                                                                min_global_node,
                                                                max_global_node,
                                                                num_global_nodes,
                                                                global_to_local_ptr,
                                                                compress)
                # save range of dirty global entries for next
                min_global_node = local_connectivity.min_global
                max_global_node = local_connectivity.max_global
                # logger.debug(__('    min-max range was = {}-{}', min_global_node, max_global_node))
                yield key, block, local_connectivity
        except:
            raise
        finally:
            if global_to_local_ptr != NULL:
                # logger.debug(__('    generator free(global_to_local_ptr {})', hex(<uintptr_t> global_to_local_ptr)))
                free(global_to_local_ptr)  # clean up after the generator, temporary buffer


class Block(EntityWithAttribute, EntityWithProperty):

    def __init__(self, database_id=-1, block_type=None, block_id=-1, name=None, topology_name=None,
                 num_entries=0, num_nodes_per_entry=0, num_edges_per_entry=0, num_faces_per_entry=0,
                 attribute_names=None, parent_entity_collection=None, *args, **kwargs):
        """
        A block of entries — elements, faces, or edges of the same topology — with connectivity information.
        
        Args:
            block_type(int): entity type of the block, one of EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK
            block_id(int): unique integer identifier
            name(str): name of the block
            topology_name(str): topology_name of the entries in the block
            num_entries(int): number of the entries in this block: nodes, edges, faces, or elements
            num_nodes_per_entry(int): number of nodes per topology_name
            num_edges_per_entry(int): number of edges per topology_name
            num_faces_per_entry(int): number of edges per topology_name
            attribute_names(Iterable(str)): iterable collection of attribute names
            parent_entity_collection(:class:`.Blocks`): the containing blocks object
        
        Raises:
            InvalidEntityType: if entity_type is not one of EDGE_BLOCK, FACE_BOCK or ELEM_BLOCK
            ArgumentTypeError: if the attribute_names is not :obj:`None` and a list could not be created from them
        """
        # logger.debug(__('Block {}', self.__class__.__mro__))
        if block_type not in BLOCK_ENTITY_TYPES:
            _raise_invalid_entity_type(BLOCK_ENTITY_TYPES)
        super(Block, self).__init__(database_id=database_id, entity_type=block_type, entity_id=block_id, name=name,
                                    parent_entity_collection=parent_entity_collection, *args, **kwargs)
        self._num_entries = num_entries
        self._num_nodes_per_entry = num_nodes_per_entry
        self._num_edges_per_entry = num_edges_per_entry
        self._num_faces_per_entry = num_faces_per_entry
        self._topology_name = topology_name
        self._attribute_names = None
        if attribute_names:
            try:
                self._attribute_names = list(attribute_names)
            except:
                raise ArgumentTypeError('A list could not be created from the attribute names.')

    def __str__(self):
        """
        A human readable representation of the block.
        """
        lookup = { EntityType.EDGE_BLOCK:'edges', EntityType.FACE_BLOCK:'faces', EntityType.ELEM_BLOCK:'elements' }
        s = f'block {self.entity_id} '
        if self._name != '':
            s += f'{self._name} '
        s += f'of {self._num_entries} {self._topology_name} {lookup[self.entity_type]}'
        if self._num_nodes_per_entry > 0:
            s += f' {self._num_nodes_per_entry} nodes'
        if self._num_edges_per_entry > 0:
            s += f' {self._num_edges_per_entry} edges'
        if self._num_faces_per_entry > 0:
            s += f' {self._num_faces_per_entry} faces'
        if self._attribute_names is not None:
            s += f' attributes {str(self._attribute_names)}'
        return s

    def _raise_entry_type_error(self):
        _raise_invalid_entity_type('one with connected entries of type NODE, EDGE, or FACE')

    @property
    def topology_name(self) -> str:
        """
        The name of the uniform topology of all elements in the block.

        Returns:
            Unambiguous name of the element topology as an uppercase string.
        """
        return self._topology_name

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def connectivity(self, entry_type=EntryType.NODE, zero_based=True, compress=False) -> Union[numpy.ndarray, util.CompressedArray]:
        """
        Read an array of the entries-to-node connectivity for a given block.

        The enumeration of the connected entries are zero based, unless zero_based is set to False, in which case
        the entries begin at one.
        
        Args:
            entry_type: (default EntryType.NODE) type of entries returned in the connectivity array, one of
                NODE, EDGE, FACE
            zero_based: if True (default) the enumeration of the connected entries begins at zero
            compress: if False (default) return a numpy.ndarray, if True return a compressed array in
                :class:`affect.util.CompressedArray`
            
        Returns:
            Array of integer type with shape ``(num_entry, num_conn_entries)``, for example,
            when entry_type is NODE, ``conn[i,j]`` is the jth node of the ith entry in the block; or
            
            :obj:`None` if there are no entries
                            
        Raises:
            InvalidEntityType: if the connect entry_type does not make sense
        """
        cdef size_t num_entries, num_conn_entries
        cdef void * data = NULL
        if self._num_entries < 1:
            return None
        num_entries = self._num_entries
        num_conn_entries = _num_connectivity_entries(self, entry_type)
        if num_conn_entries == 0:
            return None
        shape = (num_entries, num_conn_entries)
        is_int64 = _is_bulk_int64(self.ex_id)
        array = None
        try:
            array = _create_connectivity(self, entry_type, shape, is_int64, &data, zero_based, compress)
            if compress:
                # array was None from above, just overwrite
                if is_int64:
                    array = _make_compressed_array(shape, numpy.int64, <int64_t *> data)
                else:
                    array = _make_compressed_array(shape, numpy.int32, <int32_t *> data)
        except:
            raise
        finally:
            if compress:
                if data != NULL:
                    # logger.debug(__('    free({})', hex(<uintptr_t> data)))
                    free(data)
        return array

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def connectivity_local(self, entry_type=EntryType.NODE, compress=False) -> LocalConnectivity:
        """
        An array of the entries-to-local-node connectivity for a given block.

        This returns an enumeration of the node IDs connected to entries in the range (0, num_local_nodes),
        where num_local_nodes in the number of unique global nodes connected to entries
        in this block only. The local nodes are always enumerated starting with zero.

        The local node connectivity array uses unsigned 32 bit integers, limiting the unique local nodes in the block
        to 4,294,967,295. For a rough estimate, for a cube shaped block of HEX8 elements with this many unique nodes,
        reaching this limit would mean the element-to-local-node array requires 137 GB of memory, and the
        local-to-global mapping from (local to global node ID) requires 34 GB of memory.

        Uses an internal temporary allocated buffer up to the size ``num_global_nodes * 8`` bytes. For more than
        four billion global nodes this requires more than 32 GB.

        Args:
            entry_type: (default EntryType.NODE) type of entries returned in the connectivity array, one of
                NODE, EDGE, FACE
            compress: if False (default) do not store the data of the LocalConnectivity instance in compressed form,
                      if True, internal storage of the instance of LocalConnectivity uses
                      :class:`affect.util.CompressedArray`

        Returns:
            local_connectivity - instance of :class:`.LocalConnectivity`, with an array of local_nodes of
            shape ``(num_entry, num_conn_entries)`` and global_nodes of shape ``(num_local_nodes,)``

        Raises:
            InvalidEntityType: if the connect entry_type does not make sense
            MaxLengthExceeded: if the number of unique local nodes in this block exceeds UINT32_MAX, 2^32 - 1,
        """
        cdef uint32_t * global_to_local_ptr = NULL
        if self._num_entries < 1:
            return None
        num_entries = self._num_entries
        num_conn_entries = _num_connectivity_entries(self, entry_type)
        if num_conn_entries == 0:
            return None
        num_global_nodes = _inquire_int(self.ex_id, cexodus.EX_INQ_NODES)
        min_global_node = 0  # initial range of possible non-zero entries in buffer
        max_global_node = num_global_nodes - 1
        try:
            # allocate our working buffer
            global_to_local_ptr = <uint32_t*> _allocate(sizeof(uint32_t) * num_global_nodes)
            local_connectivity = _create_connectivity_local(self,
                                                            entry_type,
                                                            num_entries,
                                                            num_conn_entries,
                                                            min_global_node,
                                                            max_global_node,
                                                            num_global_nodes,
                                                            global_to_local_ptr,
                                                            compress)
        except:
            raise
        finally:
            if global_to_local_ptr != NULL:
                # logger.debug(__('    free({})', hex(<uintptr_t> global_to_local_ptr)))
                free(global_to_local_ptr)  # clean up the temporary buffer
        return local_connectivity

    def connectivity_partial(self, start: int, stop: int, entry_type=EntryType.NODE, zero_based=True) -> numpy.ndarray:
        """
        Read an array of a subsection of the entries-to-node connectivity for a given block.

        The enumeration of the connected entries are zero based, unless zero_based is set to False, in which case
        the entries begin at one.
        
        Args:
            start: first local element index at which nodes are desired, in the range [0, num_entry - 1].
            stop:  one past the end of the last desired local entry index, in the range [start+1, num_entry]
            entry_type(EntryType): (default EntryType.NODE) type of entries returned in the connectivity array, one of
                            EntryType.NODE, EntryType.EDGE, EntryType.FACE
            zero_based(bool): if True (default) the enumeration of the connected entries in the array begins at 0,
                if False, entries start at 1.
            
        Returns:
             An array of :obj:`numpy.int` with shape ``(num_entry, num_conn_entries)``, for example,
             when entry_type is NODE, array entry [i,j] is the jth node of the ith entry in the block; or
             
             :obj:`None` if there are no entries
             
        Raises:
            RangeError: if the start, stop range is invalid
        """
        cdef int64_t num_entries
        cdef numpy.ndarray[numpy.int64_t, ndim=2] conn64
        cdef int64_t *conn_ptr64
        cdef numpy.ndarray[numpy.int32_t, ndim=2] conn32
        cdef int32_t *conn_ptr32
        cdef error = -1
        cdef int64_t index64
        cdef int32_t index32
        cdef size_t length_total
        num_entries = self._num_entries
        if num_entries < 1:
            return None
        if start < 0 or stop <= start or stop > num_entries:
            raise RangeError(f'(start, stop) arguments are out of range [0, {self._num_entries}]')
        num_entries = stop - start
        if entry_type == EntryType.NODE:
            num_conn_entries = self._num_nodes_per_entry
        elif entry_type == EntryType.EDGE:
            num_conn_entries = self._num_edges_per_entry
        elif entry_type == EntryType.FACE:
            num_conn_entries = self._num_faces_per_entry
        else:
            self._raise_entry_type_error()
        if num_conn_entries == 0:
            return None
        array_shape = (num_entries, num_conn_entries)
        if _is_bulk_int64(self.ex_id):
            conn64 = util.empty_aligned(array_shape, dtype=numpy.int64)
            conn_ptr64 = <int64_t *> conn64.data
            if entry_type == EntryType.NODE:
                error = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type.value, self.entity_id, start+1, stop,
                                                     conn_ptr64, NULL, NULL)
            elif entry_type == EntryType.EDGE:
                error = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type.value, self.entity_id, start+1, stop,
                                                     NULL, conn_ptr64, NULL)
            elif entry_type == EntryType.FACE:
                error = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type.value, self.entity_id, start+1, stop,
                                                     NULL, NULL, conn_ptr64)
            else:
                self._raise_entry_type_error()
            if error != 0:
                _raise_io_error()
            if zero_based:
                length_total = num_entries * num_conn_entries
                with nogil:
                    cexodus.ex_aligned_to_zero_based_int64(conn_ptr64, length_total)
            return conn64
        else:
            conn32 = util.empty_aligned(array_shape, dtype=numpy.int32)
            conn_ptr32 = <int32_t *> conn32.data
            if entry_type == EntryType.NODE:
                error = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type.value, self.entity_id, start+1, stop,
                                                     conn_ptr32, NULL, NULL)
            elif entry_type == EntryType.EDGE:
                error = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type.value, self.entity_id, start+1, stop,
                                                     NULL, conn_ptr32, NULL)
            elif entry_type == EntryType.FACE:
                error = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type.value, self.entity_id, start+1, stop,
                                                     NULL, NULL, conn_ptr32)
            else:
                self._raise_entry_type_error()
            if error != 0:
                _raise_io_error()
            if zero_based:
                length_total = num_entries * num_conn_entries
                with nogil:
                    cexodus.ex_aligned_to_zero_based_int32(conn_ptr32, length_total)
            return conn32


class Maps(EntityDictionaryWithProperty):

    def __init__(self, database_id=-1, map_type=None, *args, **kwargs):
        """
        The collection of all Maps of a specific type in a Database.

        An ordered dictionary of all ( map_id, :class:`.Map` ) pairs, based on the :term:`map` concept.
        
        Args:   
            database_id(int): ID of the Exodus database.
            map_type(EntityType): (default :obj:`None` the type of the entities collected, one of EntityType.ELEM_MAP, 
                EntityType.NODE_MAP, EntityType.EDGE_MAP, or EntityType.FACE_MAP
                
        Raises:
            InvalidEntityType: if the map_type is invalid
        """
        # logger.debug(__('Maps {}', self.__class__.__mro__))
        if map_type not in MAP_ENTITY_TYPES:
            _raise_invalid_entity_type(MAP_ENTITY_TYPES)
        super(Maps, self).__init__(database_id=database_id, entity_type=map_type, *args, **kwargs)

    def _getitem_derived(self, key):
        return Map(self.ex_id, self.entity_type, key, parent_entity_collection=self)

    def __setitem__(self, key, value):
        EntityDictionaryWithProperty.__setitem__(self, key, value)

    def __delitem__(self, key):
        EntityDictionaryWithProperty.__delitem__(self, key)


class Map(EntityWithProperty):

    def __init__(self, database_id=-1, map_type=None, map_id=-1, name=None, parent_entity_collection=None,
                *args, **kwargs):
        """
        A map entity.
        
        A :term:`map` has :term:`properties` and a renumbering the entries of a certain type.

        Args:
            database_id(int): ID of the Exodus database.
            map_type(EntityType): (default :obj:`None`) the type of the entities collected, one of ELEM_MAP, NODE_MAP, 
                EDGE_MAP, or FACE_MAP
            map_id(int): ID of the Map
            name(str): name of the Map
            parent_entity_collection(Maps): (default :obj:`None`) the parent :class:`Maps` collection
                
        Raises:
            InvalidEntityType: if the map_type is invalid
        """
        if map_type not in MAP_ENTITY_TYPES:
            _raise_invalid_entity_type(MAP_ENTITY_TYPES)
        super(Map, self).__init__(database_id=database_id, entity_type=map_type, entity_id=map_id, name=name,
                                  parent_entity_collection=parent_entity_collection, *args, **kwargs)

    def entries(self, zero_based=True) -> numpy.ndarray:
        """
        Read the integer entries of the map.
        
        Args:
            zero_based(bool): (default=True) or False if entry ID's beginning at 1 are desired.
            
        Returns:
            An array of :obj:`numpy.int` with shape ``(num_entries,)``; or
            
            :obj:`None` if there are no entries.
        """
        cdef numpy.ndarray[numpy.int64_t, ndim=1] map_entries64
        cdef int64_t *map_ptr64
        cdef numpy.ndarray[numpy.int32_t, ndim=1] map_entries32
        cdef int *map_ptr32
        cdef int64_t num_entries = self.num_entries()
        cdef int64_t i
        if num_entries < 1:
            return None
        if _is_maps_int64(self.ex_id):
            map_entries64 = util.empty_aligned(num_entries, dtype=numpy.int64)
            map_ptr64 = <int64_t *> map_entries64.data
            if 0 != cexodus.ex_get_num_map(self.ex_id, self.entity_type.value, self.entity_id, map_ptr64):
                _raise_io_error()
            if zero_based:
                with nogil:
                    cexodus.ex_aligned_to_zero_based_int64(map_ptr64, num_entries)
            return map_entries64
        else:
            map_entries32 = util.empty_aligned(num_entries, dtype=numpy.int32)
            map_ptr32 = <int *> map_entries32.data
            if 0 != cexodus.ex_get_num_map(self.ex_id, self.entity_type.value, self.entity_id, map_ptr32):
                _raise_io_error()
            if zero_based:
                for i in range(num_entries):
                    map_ptr32[i] -= 1
            return map_entries32

    def entries_partial(self, start: int, stop: int, zero_based=True) -> numpy.ndarray:
        """
        Read a section of the the integer entries of the map.
        
        Args:
            start: first entry, in range ``[0, len_entity(map_type)-1]``
            stop: one past the last desired entry, in range ``[start+1, len_entity(map_type)]``
            zero_based(bool): (default=True) or False if entry ID's beginning at 1 are desired.
        
        Returns:
            An array of :obj:`numpy.int` with shape ``(stop - start,)``; or 
            
            :obj:`None` if there are no entries
        """
        cdef int64_t num_entries
        cdef numpy.ndarray[numpy.int64_t, ndim=1] map_entries64
        cdef int64_t * map_ptr64
        cdef numpy.ndarray[numpy.int32_t, ndim=1] map_entries32
        cdef int32_t * map_ptr32
        num_entries = stop - start
        if num_entries < 1:
            return None
        if _is_maps_int64(self.ex_id):
            map_entries64 = util.empty_aligned(num_entries, dtype=numpy.int64)
            map_ptr64 = <int64_t *> map_entries64.data
            if 0 != cexodus.ex_get_partial_num_map(self.ex_id, self.entity_type.value, self.entity_id,
                                                   start+1, stop, map_ptr64):
                _raise_io_error()
            if zero_based:
                with nogil:
                    cexodus.ex_aligned_to_zero_based_int64(map_ptr64, num_entries)
            return map_entries64
        else:
            map_entries32 = util.empty_aligned(num_entries, dtype=numpy.int32)
            map_ptr32 = <int32_t *> map_entries32.data
            if 0 != cexodus.ex_get_partial_num_map(self.ex_id, self.entity_type.value, self.entity_id,
                                                   start+1, stop, map_ptr32):
                _raise_io_error()
            if zero_based:
                with nogil:
                    cexodus.ex_aligned_to_zero_based_int32(map_ptr32, num_entries)
            return map_entries32


class Sets(EntityDictionaryWithProperty, EntityCollectionWithVariable):

    def __init__(self, database_id=-1, set_type=None, *args, **kwargs):
        """
        A collection of all the Set entities of a specific type on a Database.
        
        An ordered dictionary of all ( set_id, :class:`.Set` ) pairs, based on the :term:`set` concept.
         
        Args:   
            database_id(int): ID of the Exodus database.
            set_type(EntityType): (default :obj:`None`) the type of the contained sets, one of ELEM_SET, NODE_SET, 
                SIDE_SET, EDGE_SET, or FACE_SET
                
        Raises:
            InvalidEntityType: if the set_type is invalid
        """
        # logger.debug(__('Sets {}', self.__class__.__mro__))
        if set_type not in SET_ENTITY_TYPES:
            _raise_invalid_entity_type(SET_ENTITY_TYPES)
        super(Sets, self).__init__(database_id=database_id, entity_type=set_type, *args, **kwargs)

    def __setitem__(self, key, value):
        EntityDictionaryWithProperty.__setitem__(self, key, value)

    def __delitem__(self, key):
        EntityDictionaryWithProperty.__delitem__(self, key)

    def _getitem_derived(self, key):
        """
        Construct a set from the database given the set ID.
        
        Args:
            key: the ID of the Set, one of those returned by :meth:`entity_ids`.
            
        Returns:
            a :class:`Set` initialized with name, sizes of entries and distribution factors, etc.
        """
        num_entries, num_dist_fact = _set_param(self.ex_id, self.entity_type.value, key)
        name = self.name(key)
        return Set(self.ex_id, self.entity_type, key, name, num_entries, num_dist_fact,
                   parent_entity_collection=self)

    def num_distribution_factors_all(self) -> int:
        """
        The total length of distribution factors over all existing sets of this type in this collection.
        
        Returns:
            the length of the concatenated set distribution factors
        """
        cdef int64_t count = 0
        if self.entity_type == EntityType.SIDE_SET:
            count  = _inquire_int(self.ex_id, cexodus.EX_INQ_SS_DF_LEN)
        elif self.entity_type == EntityType.NODE_SET:
            count  = _inquire_int(self.ex_id, cexodus.EX_INQ_NS_DF_LEN)
        elif self.entity_type == EntityType.ELEM_SET:
            count  = _inquire_int(self.ex_id, cexodus.EX_INQ_ELS_DF_LEN)
        elif self.entity_type == EntityType.EDGE_SET:
            count  = _inquire_int(self.ex_id, cexodus.EX_INQ_ES_DF_LEN)
        elif self.entity_type == EntityType.FACE_SET:
            count = _inquire_int(self.ex_id, cexodus.EX_INQ_FS_DF_LEN)
        else:
            _raise_invalid_entity_type(SET_ENTITY_TYPES)
        return count


class Set(EntityWithAttribute, EntityWithProperty):

    def __init__(self, database_id=1, set_type=None, set_id=-1, name=None, num_entries=0, num_dist_fact=0,
                 attribute_names=None, parent_entity_collection=None, *args, **kwargs):
        """
        A set of entries in a Database.
        
        The :term:`set` may contain :term:`attributes`, :term:`properties` and :term:`distribution factors`.
        
        Args:
            database_id(int): ID of the Exodus database.
            set_type(EntityType): (default :obj:`None`) the type of the entities collected, one of NODE_SET, SIDE_SET, 
                EDGE_SET, FACE_SET, ELEM_SET
            set_id(int): ID of the set
            name(str): name of the set
            num_entries(int): (default 0) number of the entries in this set: nodes, sides, edges, faces, or elements
            num_dist_fact(int): (default 0) number of distribution factors on this set
            attribute_names(Iterable[str]): names of the attributes on this set
            parent_entity_collection(Sets): (default :obj:`None`) the parent :class:`Sets` collection
                
        Raises:
            InvalidEntityType: if the set_type is invalid
            ArgumentTypeError: if the attribute_names is not :obj:`None` and a list could not be created from them
        """
        # logger.debug(__('Set {}', self.__class__.__mro__))
        if set_type not in SET_ENTITY_TYPES:
            _raise_invalid_entity_type(SET_ENTITY_TYPES)
        super(Set, self).__init__(database_id=database_id, entity_type=set_type, entity_id=set_id, name=name,
                                  parent_entity_collection=parent_entity_collection, *args, **kwargs)
        self._num_entries = num_entries
        self._num_dist_fact = num_dist_fact
        self._attribute_names = None
        if attribute_names:
            try:
                self._attribute_names = list(attribute_names)
            except:
                raise ArgumentTypeError('A list could not be created from the attribute names.')

    def entries(self, zero_based=True):
        """
        Read an array of numbers that identifying the members of the set.

        The return value is different when it is a :attr:`EntityType.SIDE_SET`.

        ======== ============================================================================================
        set_type return value
        ======== ============================================================================================
        SIDE_SET two arrays of size num_entries, [i] is the the ith element number, [j] is the jth local side
        NODE_SET array of size num_entries, [i] is the ith node number in the set
        ELEM_SET array of size num_entries, [i] is the ith element number in the set
        FACE_SET array of size num_entries, [i] is the ith face number in the set
        EDGE_SET array of size num_entries, [i] is the ith edge number in the set
        ======== ============================================================================================
        
        May use temporary allocated buffers up to size ``num_entries * 8`` bytes when ``set_type == SIDE_SET``.
        
        Args:
            zero_based: (default=True) False if entry ID's beginning at 1 are desired.
            
        Returns:
            set_entries(|uint32_1d|) - array of entry IDs of size ``(num_entries)``, or

            set_entries(|uint32_1d|), set_sides(|int8_1d|) - array of element entries IDs and array of local element
                sides, both size ``(num_entries)``
        """
        cdef int64_t * set_entries_ptr64 = NULL
        cdef int64_t * set_extra_ptr64 = NULL
        cdef int32_t * set_entries_ptr32 = NULL
        cdef int32_t * set_extra_ptr32 = NULL
        cdef int8_t * set_sides_ptr8 = NULL
        cdef numpy.ndarray[numpy.int64_t, ndim=1] set_entries64
        cdef numpy.ndarray[numpy.int32_t, ndim=1] set_entries32
        cdef numpy.ndarray[numpy.int8_t, ndim=1] set_sides8
        cdef int64_t num_entries = self.num_entries()
        if num_entries < 1:
            return None
        if _is_bulk_int64(self.ex_id):
            set_entries64 = util.empty_aligned(num_entries, dtype=numpy.int64)
            set_entries_ptr64 = <int64_t*> set_entries64.data
            if self.entity_type == EntityType.SIDE_SET:
                try:
                    set_extra_ptr64 = <int64_t*> _allocate(sizeof(int64_t) * num_entries)
                    if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type.value, self.entity_id,
                                               set_entries_ptr64, set_extra_ptr64):
                        raise _raise_io_error()
                    set_sides8 = util.empty_aligned(num_entries, dtype=numpy.int8)
                    set_sides_ptr8 = <int8_t *> set_sides8.data
                    if zero_based:
                        with nogil:
                            cexodus.ex_aligned_to_zero_based_int64(set_entries_ptr64, num_entries)
                        for i in range(num_entries):
                            set_sides_ptr8[i] = <int8_t> (set_extra_ptr64[i] - 1)
                    else:
                        for i in range(num_entries):
                            set_sides_ptr8[i] = <int8_t> set_extra_ptr64[i]
                except:
                    raise
                finally:
                    if set_extra_ptr64 != NULL:
                        # logger.debug(__('    free({})', hex(<uintptr_t> set_extra_ptr64)))
                        free(set_extra_ptr64)
                return set_entries64, set_sides8
            else:
                if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type.value, self.entity_id, set_entries_ptr64, NULL):
                    raise _raise_io_error()
                if zero_based:
                    with nogil:
                        cexodus.ex_aligned_to_zero_based_int64(set_entries_ptr64, num_entries)
                return set_entries64
        else: # 32 bit integer version
            set_entries32 = util.empty_aligned(num_entries, dtype=numpy.int32)
            set_entries_ptr32 = <int32_t *> set_entries32.data
            if self.entity_type == EntityType.SIDE_SET:
                try:
                    set_extra_ptr32 = <int32_t *> _allocate(sizeof(int32_t) * num_entries)
                    if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type.value, self.entity_id,
                                               set_entries_ptr32, set_extra_ptr32):
                        raise _raise_io_error()
                    set_sides8 = util.empty_aligned(num_entries, dtype=numpy.int8)
                    set_sides_ptr8 = <int8_t *> set_sides8.data
                    if zero_based:
                        with nogil:
                            cexodus.ex_aligned_to_zero_based_int32(set_entries_ptr32, num_entries)
                        for i in range(num_entries):
                            set_sides_ptr8[i] = <int8_t> (set_extra_ptr32[i] - 1)
                    else:
                        for i in range(num_entries):
                            set_sides_ptr8[i] = <int8_t> set_extra_ptr32[i]
                except:
                    raise
                finally:
                    if set_extra_ptr32 != NULL:
                        # logger.debug(__('    free({})', hex(<uintptr_t> set_extra_ptr32)))
                        free(set_extra_ptr32)
                return set_entries32, set_sides8
            else:
                if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type.value, self.entity_id, set_entries_ptr32, NULL):
                    raise _raise_io_error()
                if zero_based:
                    with nogil:
                        cexodus.ex_aligned_to_zero_based_int32(set_entries_ptr32, num_entries)
                return set_entries32

    def num_distribution_factors(self) -> int:
        """
        The count of distribution factors on this set.
        
        Returns:
            number of distribution factors
        """
        return self._num_dist_fact

    def distribution_factors(self) -> numpy.ndarray:
        """
        Read the distribution factors for the entries in this set.

        A special scalar variable, :term:`distribution factors` can be associated with a set collection type.
        There is either zero or one distribution factors on a given set.
        
        Returns:
             Array of type :obj:`numpy.double` with shape ``(num_entries,)``, or 
             
             :obj:`None` if no distribution factors exist.
        """
        cdef int64_t num_entry_in_set = self._num_entries
        cdef int num_dist_fact_in_set = self._num_dist_fact
        if num_dist_fact_in_set < 1 or num_entry_in_set < 1:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] factors = util.empty_aligned(num_entry_in_set, dtype=numpy.double)
        cdef double *factors_ptr = <double *> factors.data
        if 0 != cexodus.ex_get_set_dist_fact(self.ex_id, self.entity_type.value, self.entity_id, factors_ptr):
            _raise_io_error()
        return factors

    def distribution_factors_partial(self, start: int, stop: int) -> numpy.ndarray:
        """
        Read the distribution factors for a subsection of the entries in a given set.
        
        Args:
            start: first entry at which variables are read, must be in the range of ``[0, num_entry-1]``.
            stop:  one past the end of the last desired entry at which the factor is read, must be in the
                            range of ``[start+1, num_entry]``.
                            
        Returns:
            Array of type :obj:`numpy.double` with shape ``(stop - start,)``, or 
            
            :obj:`None` if no distribution factors exists
        """
        cdef int64_t num_entry_in_set = self._num_entries
        cdef int num_dist_fact_in_set = self._num_dist_fact
        cdef int64_t len_dist_fact = stop - start
        if num_dist_fact_in_set < 1 or num_entry_in_set < 1 or len_dist_fact < 1:
            return None
        cdef numpy.ndarray[numpy.double_t, ndim=1] factors = util.empty_aligned(len_dist_fact, dtype=numpy.double)
        cdef double *factors_ptr = <double *> factors.data
        if 0 != cexodus.ex_get_partial_set_dist_fact(self.ex_id, self.entity_type.value, self.entity_id,
                                                     start+1, stop, factors_ptr):
            _raise_io_error()
        return factors


class DatabaseFile(object):

    def __init__(self, path: str, mode = Mode.READ_ONLY):
        """A context manager for opening and using a Database.

        Use the *with* statement while instantiating this class to easily connect to a :class:`.Database`, operate on 
        it, and ensure the :class:`Database` is closed properly after all read/writes are complete. There is no need
        to surround your code with a try/finally block in order to ensure the :class:`Database` is closed after an
        error occurs.

        Example:

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     coord = e.nodal.coordinates()
            >>> # file automatically now closed

        Args:
            path: The path to the file to open or create.
            mode: (default :attr:`Mode.READ_ONLY`) the mode on which to open the database

        The mode is one of four exclusive values (which may not be OR'ed together).

        ======================= ===============================================
        mode                    allowed actions on the database
        ======================= ===============================================
        :attr:`Mode.READ_ONLY`  read if it exists but not modify
        :attr:`Mode.READ_WRITE` read and/or append and modify if it exists
        :attr:`Mode.CREATE`     write only if it does not already exist
        :attr:`Mode.REPLACE`    write and destroy if it already exists
        ======================= ===============================================
            
        Raises:
            IOError: if the file given by path cannot be accessed
        """
        self.database = Database(path, mode)

    def __enter__(self):
        return self.database
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.database.close()
        return False


cdef class Database:
    """
    An Exodus database for reading/wring Exodus format files.
    """
    
    # these are not class attributes because this is a cdef extension type
    cdef int ex_id
    cdef int _mode
    cdef float _version

    cdef public object path  # str
    cdef object _globals
    cdef object _nodal
    cdef object _blocks
    cdef object _maps
    cdef object _sets

    def __cinit__(self, path: str, mode = Mode.READ_ONLY):
        # __cinit__ is called prior to any __init__ method.
        # Be careful what you do in the __cinit__() method, because the object
        # may not yet be fully valid Python object when it is called. Therefore,
        # you should be careful invoking any Python operations which might touch
        # the object; in particular, its methods.
        self.ex_id = 0
        # check file existence and read-ability/write-ability in python before we
        # enter the Exodus library so we can raise more informative exceptions
        if mode == Mode.READ_ONLY:
            _assert_is_file(path, True)
            _assert_file_access(path, os.R_OK)
        elif mode == Mode.READ_WRITE:
            _assert_is_file(path, True)
            _assert_file_access(path, os.R_OK)
            _assert_file_access(path, os.W_OK)
        elif mode == Mode.CREATE:
            _assert_is_file(path, False)
            _assert_file_access(path, os.W_OK)
        elif mode == Mode.REPLACE:
            _assert_file_access(path, os.W_OK)
        else:
            abs_path = path
            try:
                abs_path = os.path.abspath(path)
            except OSError:
                pass
            finally:
                raise FileAccess(f'Unknown mode for file {abs_path}.')

        # Size/precision of floating point numbers in the database:
        # Almost all platforms map Python floats to IEEE-754 “double precision”.
        # Try to force reading all floating point data as doubles.
        cdef int io_ws = 0     # for reading we just let the word size be whatever it is inside the file
        cdef int comp_ws = 8   # but we convert it to double precision during queries

        cdef float db_version
        py_byte_string = path.encode('UTF-8')  # where path is a py unicode string (Python 3 and later)
        cdef char* c_path = py_byte_string     # takes the pointer to the byte buffer of the Python byte string
        self.ex_id = cexodus.ex_open_int(c_path, mode.value, &comp_ws, &io_ws, &db_version, cexodus.EX_API_VERS_NODOT)
        if self.ex_id < 0:
            _raise_io_error()

        self.path = path
        self._mode = mode.value
        self._version = db_version
        self._blocks = {}
        self._maps = {}
        self._sets = {}

        # do not do this, instead use the 32 bit wherever we can to save memory unless the Exodus library says 64
        # self._int_size = 32
        # cdef int status = cexodus.ex_int64_status(self.ex_id)
        # if  (status & cexodus.EX_MAPS_INT64_DB) or \
        #     (status & cexodus.EX_IDS_INT64_DB) or  \
        #     (status & cexodus.EX_BULK_INT64_DB):
        #     cexodus.ex_set_int64_status(self.ex_id, cexodus.EX_ALL_INT64_API)
        #     self._int_size = 64


    def __init__(self, path: str, mode = Mode.READ_ONLY):
        """
        Opens a connection to a Database object with the given path and mode.
        
        Args:
            path(str): The path to the file to open or create.
            mode(:class:`Mode`): (default :attr:`Mode.READ_ONLY`) mode with which to act on the database
            
        Raises:
            FileNotFound: reading a file that does not exist
            FileExists: if trying to write to existing file without replace mode
            FileAccess: if the mode is invalid
        """
        # (we now have a valid ex_id, since __cinit__ would have been already called, also
        # arguments were only used by the __cinit__ method)
        self._globals = None
        self._nodal = None


    def __del__(self):
        """
        Cleanup any resources and close the ExodusII database.
        """
        # the __del__ function does not often actually get called automatically by Python, but just in case.
        self.close()

    def __str__(self):
        """Return a human readable representation of the object.
        """
        return f'{{ id: {self.ex_id}, filename: {self.path}, version: {self._version:.6g}}}'

    @property
    def mode(self) -> Mode:
        """
        Get the mode with which the ExodusII database was opened.

        Returns:
            Mode with which the :class:`Database` object was opened.
        """
        return Mode(self._mode)
    
    def summary(self) -> str:
        """
        Get a human readable string representation of the global parameters of the Database.
        
        This is a shortcut instead of using multiple calls to the :class:`.Global` object. 
        
        Returns:
            String representation of the global parameters with line breaks.
        """
        cdef cexodus.ex_init_params param
        if 0 != cexodus.ex_get_init_ext(self.ex_id, &param):
            _raise_io_error()
        return str(f'{{title:         {_to_unicode(param.title)},\n' 
                f' num_dim:       {param.num_dim},\n'
                f' num_node:      {param.num_nodes},\n'
                f' num_elem:      {param.num_elem},\n'
                f' num_face:      {param.num_face},\n'
                f' num_edge:      {param.num_edge},\n'
                f' num_elem_blk:  {param.num_elem_blk},\n'
                f' num_face_blk:  {param.num_face_blk},\n'
                f' num_edge_blk:  {param.num_edge_blk},\n'
                f' num_node_sets: {param.num_node_sets},\n'
                f' num_side_sets: {param.num_side_sets},\n'
                f' num_elem_sets: {param.num_elem_sets},\n'
                f' num_face_sets: {param.num_face_sets},\n'
                f' num_edge_sets: {param.num_edge_sets},\n'
                f' num_node_maps: {param.num_node_maps},\n'
                f' num_elem_maps: {param.num_elem_maps},\n'
                f' num_face_maps: {param.num_face_maps},\n'
                f' num_edge_maps: {param.num_edge_maps}}}')

    def _len_max_names(self):
        cdef int max_all_name_length = <int> cexodus.ex_inquire_int(self.ex_id,
                                                                    cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
        if max_all_name_length < 0:
            _raise_io_error()
        cdef int max_use_name_length = <int> cexodus.ex_inquire_int(self.ex_id, cexodus.EX_INQ_DB_MAX_USED_NAME_LENGTH)
        if max_use_name_length < 0:
            _raise_io_error()
        # logger.debug(__('File {} can use at most {}-character names', self.path, max_all_name_length))
        # logger.debug(__('File {} used no more than {}-character names', self.path, max_use_name_length))
        max_name_length = max_use_name_length
        if 0 != cexodus.ex_set_max_name_length(self.ex_id, max_use_name_length):
            _raise_io_error()
        return max_use_name_length

    def close(self) -> None:
        """
        Close an existing Exodus file database.

        You should *never* have to call close because it is called automatically when the :class:`Database` object goes out of 
        scope, that is, when the :meth:`.Database.__del__` method is automatically called.
        """
        if 0 != cexodus.ex_close(self.ex_id):
            _raise_io_error()
        # noinspection PyAttributeOutsideInit
        self.ex_id = -1
        # noinspection PyAttributeOutsideInit
        self.path = None

    def version(self) -> float:
        """
        The version used to create this ExodusII database.

        The database version number reflects the version of the ExodusII API library that was used to create the file.

        Returns:
            Floating point number representation of an ExodusII API library version.
        """
        return self._version

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @property
    def info_records(self) -> List[str]:
        """
        A property that returns the information records as a list of text strings.

        Uses temporary allocated buffer of size ``num_records * (8 + cexodus.MAX_LINE_LENGTH + 1)`` bytes.

        Returns:
            list of :term:`Information Data` strings
        """
        cdef int len_str = cexodus.MAX_LINE_LENGTH
        cdef int64_t num_inf_rec = _inquire_int(self.ex_id, cexodus.EX_INQ_INFO)

        def info_records_func(uintptr_t str_ptr):
            return cexodus.ex_get_info(self.ex_id, <char **>str_ptr)

        return _get_strings(len_str, num_inf_rec, info_records_func)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @property
    def qa_records(self) -> List[str]:
        """
        A property that returns the list of of the quality assurance records.
        
        The :term:`quality assurance records` are returned as one concatenated list of strings, 
        in multiples of four, sequentially:
        
        1. Code name, an application code that has created or modified the database.
        2. Code QA descriptor, usually the version of the application code.
        3. Date the application code touched the database, in the format MM/DD/YY.
        4. Time the application code touched the database, in the 24 hour format HH:MM:SS.

        Example:

            Print all the quality assurance records in a database.
    
            >>> with DatabaseFile('/tmp/myExodusFile.exo', Mode.READ_ONLY) as e:
            >>>     records = e.qa_records
            >>>     if len(records):
            >>>         print(f'len(qa_records): {len(records)}')
            >>>         print('qa records:')
            >>>         for record in records:
            >>>             print('  {record}')
    
        Uses temporary allocated buffer of size ``4 * num_qa_records * (8 + cexodus.MAX_STR_LENGTH + 1)`` bytes.

        Returns:
            list with length `4 * num_qa_records` of quality assurance strings
        """
        cdef char** qa_record_ptr = NULL
        cdef int len_str = cexodus.MAX_STR_LENGTH
        cdef int64_t num_qa_rec = _inquire_int(self.ex_id, cexodus.EX_INQ_QA)
        cdef int64_t length_qa = 4 * num_qa_rec
        records = []
        if num_qa_rec > 0:
            try:
                qa_record_ptr = <char**> _unaligned_allocate(sizeof(char*) * length_qa)
                try:
                    for i in range(length_qa):
                        qa_record_ptr[i] = <char *> _unaligned_allocate(sizeof(char) * (len_str+1))
                    if 0 != cexodus.ex_get_qa(self.ex_id, <char *(*)[4]> &qa_record_ptr[0]):
                        _raise_io_error()
                    for i in range(length_qa):
                        record = _to_unicode(qa_record_ptr[i])
                        records.append(record)
                finally:
                    for i in range(length_qa):
                        free(qa_record_ptr[i])
            finally:
                free(qa_record_ptr)
        return records

    @property
    def globals(self) -> Global:
        """
        A property that returns the :class:`.Global` object for this database.
        
        The :class:`.Global` includes the number and value of time planes, any global variables, and the total number of 
        entries across all blocks such as the global number of elements and nodes, and global id maps.

        Example:

            Print the number of time steps and global number of nodes of a database.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     num_time_steps = e.globals.num_times()
            >>>     print(f'num_times = {num_time_steps}')
            >>>     num_nodes = e.globals.num_nodes()
            >>>     print(f'num_nodes = {num_nodes}')
            
        Returns:
            The singleton :class:`.Global` object for this Database.
        """
        if self._globals is None:
            self._globals = Global(self.ex_id)
        return self._globals

    @property
    def nodal(self) -> Nodal:
        """
        A property that returns the :class:`.Nodal` object on this database for accessing global nodal coordinates, 
        nodal variables, and nodal attributes.

        Example:

            Get a bounding box (minimum, maximum) over all the node coordinates in the database.

            >>> import numpy
            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     coord = e.nodal.coordinates()
            >>>     min_coord  = numpy.min(coord, axis=0)
            >>>     max_coord = numpy.max(coord, axis=0)
        """
        if self._nodal is None:
            self._nodal = Nodal(self.ex_id)
        return self._nodal

    def _get_blocks(self, block_type):
        value = self._blocks.get(block_type, None)
        if value is None:
            value = Blocks(self.ex_id, block_type)
            self._blocks[block_type] = value
        return value
    
    @property
    def element_blocks(self) -> Blocks:
        """
        A property that returns the :class:`.Blocks` dictionary for accessing the collection of entities of
        ``Block(entity_type = ELEM_BLOCK)``.

        Example:

            Iterate over the elements blocks in a database and print the block ID's and number of elements in each 
            block.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     print(f'database has {len(e.element_blocks)} element blocks')
            >>>     for block_id, block in e.element_blocks.items():
            >>>         print(f'block_{block_id} num_elements = {block.num_entries}')
        """
        return self._get_blocks(EntityType.ELEM_BLOCK)

    @property
    def face_blocks(self) -> Blocks:
        """
        A property that returns the :class:`.Blocks` dictionary for accessing the collection of entities of
        ``Block(entity_type = EntityType.FACE_BLOCK)``.
        """
        return self._get_blocks(EntityType.FACE_BLOCK)

    @property
    def edge_blocks(self) -> Blocks:
        """
        A property that returns the :class:`.Blocks` dictionary for accessing the collection of entities of
        ``Block(entity_type = EntityType.EDGE_BLOCK)``.
        """
        return self._get_blocks(EntityType.EDGE_BLOCK)

    def _get_maps(self, map_type):
        value = self._maps.get(map_type, None)
        if value is None:
            value = Maps(self.ex_id, map_type)
            self._maps[map_type] = value
        return value
    
    @property
    def element_maps(self) -> Maps:
        """
        A property that returns the :class:`.Maps` dictionary for accessing the collection of entities of
        ``Map(entity_type = EntityType.ELEM_MAP)``.
        """
        return self._get_maps(EntityType.ELEM_MAP)

    @property
    def face_maps(self) -> Maps:
        """
        A property that returns the :class:`.Maps` dictionary for accessing the collection of entities of
        ``Map(entity_type = EntityType.FACE_MAP)``.
        """
        return self._get_maps(EntityType.FACE_MAP)

    @property
    def edge_maps(self) -> Maps:
        """
        A property that returns the :class:`.Maps` dictionary for accessing the collection of entities of
        ``Map(entity_type = EntityType.EDGE_MAP)``.
        """
        return self._get_maps(EntityType.EDGE_MAP)

    @property
    def node_maps(self) -> Maps:
        """
        A property that returns the :class:`.Maps` dictionary for accessing the collection of entities of
        ``Map(entity_type = EntityType.NODE_MAP)``.
        """
        return self._get_maps(EntityType.NODE_MAP)

    def _get_sets(self, set_type):
        value = self._sets.get(set_type, None)
        if value is None:
            value = Sets(self.ex_id, set_type)
            self._sets[set_type] = value
        return value

    @property
    def element_sets(self) -> Sets:
        """
        A property that returns the :class:`.Sets` dictionary for accessing the collection of entities of
        ``Set(entity_type = EntityType.ELEM_SET)``.
        """
        return self._get_sets(EntityType.ELEM_SET)

    @property
    def face_sets(self) -> Sets:
        """
        A property that returns the :class:`.Sets` dictionary for accessing the collection of entities of
        ``Set(entity_type = EntityType.FACE_SET)``.
        """
        return self._get_sets(EntityType.FACE_SET)

    @property
    def edge_sets(self) -> Sets:
        """
        A property that returns the :class:`.Sets` dictionary for accessing the collection of entities of
        ``Set(entity_type = EntityType.EDGE_SET)``.
        """
        return self._get_sets(EntityType.EDGE_SET)

    @property
    def node_sets(self) -> Sets:
        """
        A property that returns the :class:`.Sets` dictionary for accessing the collection of entities of
        ``Set(entity_type = EntityType.NODE_SET)``.
        """
        return self._get_sets(EntityType.NODE_SET)

    @property
    def side_sets(self) -> Sets:
        """
        A property that returns the :class:`.Sets` dictionary for accessing the collection of entities of
        ``Set(entity_type = EntityType.SIDE_SET)``.
        """
        return self._get_sets(EntityType.SIDE_SET)

    #
    # Groups
    #

    def num_children(self):
        """
        Get the number of groups contained in this group (self.ex_id).
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_NUM_CHILD_GROUPS)

    def get_parent(self):
        """
        Return the id of parent of this (self.ex_id) group; returns (self.ex_id) if at root.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_GROUP_PARENT)

    def get_root(self):
        """
        Return the _id of root group of this (self.ex_id) group; returns self.ex_id if at root.
        """
        return _inquire_int(self.ex_id, cexodus.EX_INQ_GROUP_ROOT)

    def get_group_name(self):
        """
        Return the name of this group self.ex_id.
        """
        return _inquire_string(self.ex_id, cexodus.EX_INQ_GROUP_NAME)

    def get_full_group_name(self):
        """
        Return the full -separated path name of this (self.ex_id) group.
        """
        return _inquire_string(self.ex_id, cexodus.EX_INQ_FULL_GROUP_NAME)
