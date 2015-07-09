#!/usr/bin/env python

# -*- coding: utf-8 -*-

# cython: embedsignature=True

"""
This module contains the classes for reading/wring an ExodusII format mesh database.

The fundamental class for I/O is :class:`~exodus.Database`. Though it can be directly
used, there is a convenient :class:`~exodus.DatabaseFile` context manager.

The main classes for database access are

.. autosummary::

    ~exodus.DatabaseFile
    ~exodus.Database
    ~exodus.Global
    ~exodus.Field
    ~exodus.FieldArray
    ~exodus.Fields
    ~exodus.Nodal
    ~exodus.Blocks
    ~exodus.Block
    ~exodus.Maps
    ~exodus.Map
    ~exodus.Sets
    ~exodus.Set
"""

cimport cexodus
cimport cconnect
cimport cython

import os
import sys
from contextlib import contextmanager
import collections
from abc import abstractmethod

import numpy as np
cimport numpy as np

from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t, int32_t
from libc.float cimport FLT_MAX
from libc.limits cimport LLONG_MIN, LONG_MIN
from libc.stdio cimport fflush, stderr

cdef extern from "unistd.h":
    enum: STDERR_FILENO

#-----------------------------------------------------------------------------------------------------------------------
# Exodus File open modes

OPEN_READ_ONLY = 1
"""Mode to initialize a :class:`~exodus.Database` for reading an existing file, modifications not allowed."""

OPEN_READ_WRITE = 2
"""Mode to initialize a :class:`~exodus.Database` for reading and appending or modifying."""

OPEN_CREATE = 3
"""Mode to initialize a :class:`~exodus.Database` for writing safely, opening only if the file does not already exist."""

OPEN_REPLACE = 4
"""Mode to initialize a :class:`~exodus.Database` for writing, destroying if the file already exists."""

#-----------------------------------------------------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------------------------------------------------
# internal utility functions
  
cdef inline unicode _to_unicode(char* s):
    """Convert C-string bytes to Python string object"""
    return s.decode('UTF-8', 'strict')

cdef inline unicode _to_unicode_with_length(char* s, size_t length):
    """Convert C-string with known length to Python string object"""
    return s[:length].decode('UTF-8', 'strict')

cdef inline void* _allocate(size_t num_bytes):
    """Allocate memory with C malloc, raise a MemoryError"""
    cdef void* ptr = malloc(num_bytes)
    if ptr == NULL:
        message = 'malloc: (size={}) failed'.format(num_bytes)
        raise MemoryError(message)
    return ptr

def _assert_is_file(path, exists):
    """
    Raise an IOError exception if (exists == True and file does not exist) or if (exists == False and file exists).
    """
    is_file = os.path.isfile(path)
    if exists:
        if not is_file:
            message = 'File does not exist {}.'.format(path)
            raise IOError(message)
    else:
        if is_file:
            abs_path = os.path.abspath(path)
            message = 'File already exists {}.'.format(abs_path)
            raise IOError(message)

def _assert_file_access(path, os_mode):
    """
    Raise an IOError exception if file is not one of the os_mode F_OK, R_OK, W_OK, or X_OK.
    """
    if not os.access(path, os_mode):
        abs_path = os.path.abspath(path)
        mode_type = {os.F_OK:'Exists', os.R_OK:'Read', os.W_OK:'Write', os.X_OK:'Execution'}
        message = '{} is not allowed for file {}.'.format(mode_type[os_mode], abs_path)
        raise IOError(message)

def _topology_name(name, nodes_per_entry, spatial_dim):
    """
    Construct an unambiguous topology type, which is a string name.

    This is typically from information about a Block entity.

    :param name: base name of the topology type
    :type name: str
    :param nodes_per_entry: number of nodes per topology
    :type nodes_per_entry: int
    :param spatial_dim: spatial dimension of the topology
    :type spatial_dim: int
    :return: type -- str
    """
    topo_name = name.replace(" ","_").lower()

    # append nodes per element if it does not end with a number
    if not topo_name[-1].isdigit() and nodes_per_entry > 1:
        topo_name += str(nodes_per_entry)

    if topo_name[:5] == "super":
        topo_name = "super" + str(nodes_per_entry)
    elif spatial_dim == 3:
        if topo_name[:3] == "tri":
            topo_name = "trishell" + str(nodes_per_entry)
    elif spatial_dim == 2:
        if topo_name == "shell2":
            topo_name = "shellline2d2"
        elif topo_name == "rod2" or topo_name == "bar2" or topo_name == "truss2":
            topo_name = "rod2d2"
        elif topo_name == "shell3":
            topo_name = "shellline2d3"
        elif topo_name == "bar3" or topo_name == "rod3" or topo_name == "truss3":
            topo_name = "rod2d3"
    return topo_name

cdef _raise_io_error():
    cdef char* msg
    cdef char* func
    cdef int err_num
    message = 'Error reading or writing Exodus database file.'
    cexodus.ex_get_err(<const char**>&msg, <const char**>&func, &err_num)
    if err_num != 0 and msg[0] != '\0':
        message += '\n[{}]: ({}) {}'.format(err_num,_to_unicode(func),_to_unicode(msg))
    raise IOError(message)

cdef inline _inquire_value(ex_id, inquiry):
    """
    Return one of the metadata values on the Exodus database, for integer and floating point values.

    :param ex_id: Exodus database ID.
    :param inquiry: variable requested, one of the cexodus.EX_INQ_* values.
    :type cexodus.ex_inquiry:
    :return: val
    :rtype: Python int or float
    """
    # return either an integer or float value requested by the inquiry
    # init with some bad values so we will know if one changed
    cdef int64_t int64_val = LLONG_MIN
    cdef float float_val = -FLT_MAX
    cdef char* str_ptr = NULL # not used
    if 0 != cexodus.ex_inquire(ex_id, inquiry, &int64_val, &float_val, str_ptr):
        _raise_io_error()
    if int64_val != LLONG_MIN:
        return int64_val
    elif float_val != -FLT_MAX:
        return float_val
    else:
        return 0

cdef inline _inquire_string(ex_id, inquiry):
    cdef int64_t int_val = LLONG_MIN
    cdef float float_val = -FLT_MAX
    cdef char* str_ptr = NULL
    cdef int len_str = -1
    if inquiry == cexodus.EX_INQ_TITLE:
        len_str = cexodus.MAX_LINE_LENGTH
    elif inquiry == cexodus.EX_INQ_GROUP_NAME:
        len_str = cexodus.EX_INQ_GROUP_NAME_LEN
    elif inquiry == cexodus.EX_INQ_FULL_GROUP_NAME:
        len_str = cexodus.EX_INQ_GROUP_NAME_LEN
    else:
        raise RuntimeError(
            "unknown string length for inquire = {}".format(inquiry))
    try:
        str_ptr = <char *> _allocate(sizeof(char) * (len_str+1))
        if 0 != cexodus.ex_inquire(ex_id, inquiry, &int_val, &float_val, &str_ptr[0]):
            _raise_io_error()
        return _to_unicode(str_ptr) # python string object
    finally:
        free(str_ptr)

cdef inline _set_param(ex_id, set_type, set_id):
    """
    Internal function to return the pair of the number of entries and the number of distribution factors which
    describe a single set.

    :param ex_id: Exodus Database ID
    :param set_type: one of NODE_SET, SIDE_SET, EDGE_SET, FACE_SET, or ELEM_SET
    :param set_id: ID of the set
    :return: num_entries, num_dist_face
    :rtype num_entries: Python int
    :rtype num_dist_fact: Python int
    """
    cdef int64_t num_entries64
    cdef int64_t num_dist_fact64
    cdef int num_entries
    cdef int num_dist_fact
    if cexodus.ex_int64_status(ex_id) & cexodus.EX_BULK_INT64_API:
        if 0 != cexodus.ex_get_set_param(ex_id, set_type, set_id,
                                         &num_entries64, &num_dist_fact64):
            _raise_io_error()
        return num_entries64, num_dist_fact64
    else:
        if 0 != cexodus.ex_get_set_param(ex_id, set_type, set_id,
                                         &num_entries, &num_dist_fact):
            _raise_io_error()
        return num_entries, num_dist_fact

#-------------------------------------------------------------------------------
# global functions for the Exodus API

def library_version():
    """
    Return the ExodusII library API version number which reflects the version number of the C language implementation
    of the Exodus library linked with this module.

    :rtype : floating point
    """
    return '{:.6g}'.format(cexodus.API_VERS)

#-----------------------------------------------------------------------------------------------------------------------
# Error capturing of Exodus library messages

DEFAULT     = cexodus.EX_DEFAULT
VERBOSE     = cexodus.EX_VERBOSE
DEBUG       = cexodus.EX_DEBUG
ABORT       = cexodus.EX_ABORT
NULLVERBOSE = cexodus.EX_NULLVERBOSE

def debug_messages(options):
    """
    Set message reporting options in the Exodus Library that control error
    and informational messages when calling :class:`~exodus.Database` functions.

    These control how errors and warnings are written to stderr. Values of options may be OR'ed together to provide any
    combination of these capabilities.

    Consider using the context manager :class:`~exodus.DebugMessages`.

    =======  ============================================================================
    option   meaning
    =======  ============================================================================
    ABORT    Causes internal library fatal errors to force program exit. (Default is off)
    DEBUG    Certain messages print for debug use. (Default is off)
    VERBOSE  All error messages print. (Default is off)
    DEFAULT  Does not print error messages. (Default is on)
    =======  ============================================================================

    :param options: one or more of ABORT, DEBUG, VERBOSE, DEFAULT
    :type options: integer
    :return: previous options value
    """
    old_options = cexodus.ex_opts(options)
    return old_options

cdef class DebugMessages(object):
    """
    A context manager for capturing debug messages from the internal Exodus library when calling
    :class:`~exodus.Database` functions.

    Within this context, the internal Exodus library is set to the DEBUG|VERBOSE options, active for all calls to the
    API. Also within this context, any Python exceptions raised will have their messages appended with the contents of
    this Exodus debug and verbose output. Upon exiting the context, the previous values of the verbose or debug
    options are restored.

    A side effect within this context manager is that any output to the C language stderr device from not only the
    Exodus library but any other C library will be rerouted to a buffer. The stderr device is where the
    Exodus library writes its debug messages.

    An alternative to using this context manager is just calling :func:`~exodus.debug_messages` with the
    DEBUG or VERBOSE options and letting it write to stderr.

    Use the `with` statement while instantiating this class to wrap calls to :class:`~exodus.Database`
    member functions. This context manager can be nested within the :class:`~exodus.DatabaseFile` context
    manager.

    :Example:

        >>> with DatabaseFile('/tmp/myExodusFile.exo', OPEN_READ_ONLY) as e:
        >>>     # messages will be appended to any exceptions raised
        >>>     with DebugMessages as messages:
        >>>         coord = e.get_coordinates()
        >>>     # debugging turned off now
        >>> # file automatically now closed
    """

    cdef int old_option
    cdef int saved_stderr
    cdef int err_pipe[2]
    cdef char* c_messages
    cdef bytearray py_messages

    def __init__(self, option=VERBOSE|DEBUG):
        """
        __init__(option=VERBOSE|DEBUG)
        Initialize and activate buffering of Exodus API verbose or debug messages.

        :param option: Either exodus.VERBOSE or exodus.VERBOSE|exodus.DEBUG etc.
        :type option: int
        :returns messages: a bytearray object used to buffer messages from stderr
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
            raise type(exc_type)(str(exc_value) + ' ' + self.py_messages).with_traceback(exc_traceback)
        return False

    def _redirect_stderr(self):
        cdef extern from "unistd.h":
            int close(int fildes)
            int dup(int fildes)
            int dup2(int fildes, int fildes2)
            int pipe(int fildes[2])
        self.saved_stderr = dup(STDERR_FILENO) # save stderr for later
        if pipe(self.err_pipe) != 0:
            raise IOError("unable to pipe error messages from ExodusII Library on stderr")
        dup2(self.err_pipe[1], STDERR_FILENO) # redirect stderr to the pipe
        close(self.err_pipe[1])

    def _release_stderr(self):
        cdef extern from "unistd.h":
            int dup2(int fildes, int fildes2)
            size_t read(int fildes, void *buf, size_t nbyte)
        global stderr
        fflush(stderr) # make sure everything is now in the pipe
        read(self.err_pipe[0], self.c_messages, cexodus.MAX_ERR_LENGTH) # get results from pipe into our buffer
        dup2(self.saved_stderr, STDERR_FILENO)  # reconnect stderr

#-----------------------------------------------------------------------------------------------------------------------
# Entity and other identifiers

NODAL      = cexodus.EX_NODAL
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
COORDINATE = cexodus.EX_COORDINATE

ELEMENT = 65
"""Denotes entries of type element."""

NODE    = 66
"""Denotes entries of type node."""

FACE    = 67
"""Denotes entries of type face."""

EDGE    = 68
"""Denotes entries of type edge."""

#-----------------------------------------------------------------------------------------------------------------------
class FieldArray(np.ndarray):
    """
    A specialization of NumPy array holding additional metadata in a Field.

    The constructor takes an already formed ndarray (from the usual numpy calls to np.array) and a given Field
    instance and returns an object.

    The len(info.components) must equal input_array.shape[-1]. The input_array.dtype must be np.double or np.float64.

    :Example:
        >>> pos_info = Field('position', ('x', 'y', 'z'))
        >>> arr = np.empty(3, dtype=np.double)
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
    """
    def __new__(cls, input_array, info=None):
        if input_array is not None:
            if input_array.dtype != np.float64:  # Input array is an already formed ndarray instance
                raise ValueError("The dtype of the input_array must be numpy.double or numpy.float64.")
        if info is not None:
            if not isinstance(info, Field):
                raise TypeError('info object must be an instance of Field or subclass of Field')
            if len(info.components) != input_array.shape[-1]:
                raise ValueError("The len(info.components) does not equal input_array.shape[-1].")
        obj = np.asarray(input_array).view(cls)  # We first cast to be our class type
        obj.info = info                          # add the new attribute to the created instance
        return obj                               # Finally, we must return the newly created object

    def __array_finalize__(self, obj):
        # ``self`` is a new object resulting from
        # ndarray.__new__(FieldArray, ...), therefore it only has
        # attributes that the ndarray.__new__ constructor gave it -
        # i.e. those of a standard ndarray.
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
        # Note that it is here, rather than in the __new__ method,
        # that we set the default value for 'info', because this
        # method sees all creation of default objects - with the
        # FieldArray.__new__ constructor, but also with
        # arr.view(FieldArray).
        self.info = getattr(obj, 'info', None)
        # We do not need to return anything

#-----------------------------------------------------------------------------------------------------------------------
class Field(object):

    def __init__(self, name, components, variables = None, parent_field_dictionary = None):
        """
        A Field represents a multi-component field stored in a database, including it's name, the number and name of it's
        components, the internal Exodus variables in which it's components are stored.

        Array data of a field must be accessed through the member functions of the associated entity on which it is active.

        Fields may be associated with an entity of type GLOBAL, NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
                                FACE_SET, ELEM_BLOCK, ELEM_SET, or SIDE_SET.

        :param name: base name of the field, without trailing subscripts
        :type name: str
        :param components: names of the subscripts of the field components, a sequence
                            either containing all str ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'],
                            or containing all int [1, 2, 3, 4, 5, 6, 7]
        :param variables: a sequence of variable indices in the Exodus database corresponding to the scalars storing the
                          components with length equal to len(components)
        :param parent_field_dictionary: the collection of fields that this Field object belongs to
        :type parent_field_dictionary: :class:`~exodus.Fields`
        """
        self.name = name
        self.components = components
        self.variables = variables
        self.parent_field_dictionary = parent_field_dictionary

        if components is None or len(components) < 1:
            raise ValueError("Length of components must be one or greater.")
        if variables and len(variables) != len(components):
            raise ValueError("Length of variable indices of field {} do not equal len(components).".format(name))
        
    def __str__(self):
        s = self.name
        if len(self.components) > 1:
            s += '_{}'.format(','.join(str(n) for n in self.components))
        if self.variables is not None:
            s += ' [{}]'.format(', '.join(str(n) for n in self.variables))
        return s

#-----------------------------------------------------------------------------------------------------------------------
class Fields(collections.MutableMapping):

    def __init__(self, parent_entity_dictionary = None):
        """
        An accessor for fields, acting as an ordered dictionary with some extra member functions.

        For a database, this is an accessor for the fields on all the members of one type of entity in a
        database. The set of entities are one of the globals, the nodal, or one of
        element blocks, face blocks, or edge blocks; or side sets, node sets, element sets, face sets, or
        edge sets.

        The keys of the dictionary are the base names of fields.

        :param parent_entity_dictionary: if reading/modifying an existing database, the owning collection of
                                        entities in an Exodus database
        :type parent_entity_dictionary: one of a Global, Nodal, Blocks, or Sets object
        :raises TypeError: if the parent collection of entities does not support fields
        """
        self.parent_entity_dictionary = parent_entity_dictionary
        self._store = None
        if not isinstance(parent_entity_dictionary, EntityDictionaryWithVariable):
            raise TypeError("Parent collection of entities does not support fields.")

    def __getitem__(self, key):
        if self._store is None:
            self._store = self.__get_field_info_dictionary()
        return self._store[key]

    def __setitem__(self, key, value):
        #self._store[key] = value
        raise NotImplementedError("Writes are not yet supported.")

    def __delitem__(self, key):
        del self._store[key]

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
        var_names = self.parent_entity_dictionary.variable_names()
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
        field_info_dictionary = collections.OrderedDict().fromkeys(d.iterkeys())
        for k, v in d.items():
            components, variables = zip(*v)
            field_info_dictionary[k] = Field(k, components, variables, self)
        return field_info_dictionary

    def is_field_active_all(self):
        """
        Returns a table of one or zero values indicating whether each field exists, or not,
        on all entities of this type on the database.

        :return: fields_table -- numpy.ndarray of type np.int32, with shape (num_entity, num_fields), or None if no
                                fields are present.
        """
        fields_dict = self._store
        if fields_dict is None:
            return None
        parent_entity_dictionary = self.parent_entity_dictionary
        cdef int num_entity = parent_entity_dictionary.num_entity()
        if num_entity < 1:
            return None
        cdef int num_var = parent_entity_dictionary.num_variables()
        if num_var < 1:
            return None
        array_shape = (num_entity, num_var)
        cdef np.ndarray[np.int32_t, ndim=2] var_table = np.empty(array_shape, dtype=np.int32)
        cdef int *var_ptr = <int *> var_table.data
        if 0 != cexodus.ex_get_truth_table(parent_entity_dictionary.ex_id, parent_entity_dictionary.entity_type,
                                           num_entity, num_var, var_ptr):
            _raise_io_error()
        num_fields = len(fields_dict)
        array_shape = (num_entity, num_fields)
        cdef np.ndarray[np.int32_t, ndim=2] fields_table = np.empty(array_shape, dtype=np.int32)
        for entity in range(num_entity):
            for index, k, in enumerate(fields_dict):
                v = fields_dict[k]
                variables = v.variables
                if var_table[entity, variables[0]] == 1:
                    for i in range(1,len(variables)):  # ensure all the any other component variables are also active
                        if var_table[entity, variables[i]] != 1:
                            msg = "ERROR: Exodus database has a field with an inactive component:\n"
                            msg += "  database {} field {} component {} of {} (entity id {} of type {})."
                            msg.format(self.ex_id, k, i, len(v), entity, self.entity_type)
                            raise AssertionError(msg)
                    fields_table[entity, index] = 1

#-----------------------------------------------------------------------------------------------------------------------
class Entity(object):
    """
    An Entity is a generic association of variables, maps, or a collection of elements, faces, edges, or nodes.

    An Entity is abstract and is the base of the Block, Set, Map class.
    """

    def __init__(self, database_id, entity_type, entity_id, name = None, parent_entity_dictionary = None):
        """
        Initialize an entity from a specific Database.

        :param database_id: the Exodus library ID of the database file
        :param entity_type: one of NODAL, GLOBAL;
                                    NODE_SET, SIDE_SET, ELEM_SET, EDGE_SET, FACE_SET;
                                    EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK;
                                    ELEM_MAP, NODE_MAP, EDGE_MAP, FACE_MAP.
        :param entity_id: the ID of the entity in the database
        :raises ValueError: if the entity_type and parent_entity_dictionary.entity_type do not match
        """
        super(Entity, self).__init__()
        self.ex_id = database_id
        self.entity_type = entity_type
        self.entity_id = entity_id
        self.parent_entity_dictionary = parent_entity_dictionary
        self._name = name
        self._num_entries = -1
        if parent_entity_dictionary is not None:
            if entity_type != parent_entity_dictionary.entity_type:
                raise ValueError("entity_type does not match parent_entity_dictionary.entity_type")

    def num_entries(self):
        """
        Return the number of entries contained within a single collection entity.

        For example return the number of elements in an element block, or return number of faces in a face set, etc.

        For entity_type GLOBAL, the num_entries is always 1.
        For entity_type (NODE_MAP, EDGE_MAP, FACE_MAP and ELEM_MAP) num_entries is
        the number of global (nodes, edges, faces, and elements), respectively.
        For entity_type SIDE_SET, the num_entries is the number of element sides.

        :return: num_entries -- number of entries in the entity.
        """
        cdef cexodus.ex_block block
        cdef int64_t num_entries
        cdef int64_t num_entry_in_set
        cdef int64_t num_dist_fact_in_set
        if self._num_entries == -1:
            if self.entity_type == NODAL:
                num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_NODES)
            elif self.entity_type == GLOBAL:
                num_entries = 1
            elif self.entity_type in (EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK):
                block.id = self.entity_id
                block.type = self.entity_type
                if 0 != cexodus.ex_get_block_param(self.ex_id, &block):
                    _raise_io_error()
                num_entries = block.num_entry
            elif self.entity_type in (NODE_SET, SIDE_SET, ELEM_SET, EDGE_SET, FACE_SET):
                num_entries, num_dist_fact = _set_param(self.ex_id, self.entity_type, self.entity_id)
            elif self.entity_type == NODE_MAP:
                num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_NODES) # returns regardless of map existing
            elif self.entity_type == EDGE_MAP:
                num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE)  # returns regardless of map existing
            elif self.entity_type == FACE_MAP:
                num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_FACE)  # returns regardless of map existing
            elif self.entity_type == ELEM_MAP:
                num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM)  # returns regardless of map existing
            else:
                raise RuntimeError("Invalid entity type.")
            self._num_entries = num_entries
        return self._num_entries

#-----------------------------------------------------------------------------------------------------------------------
class EntityWithProperty(Entity):

    def __init__(self, database_id, entity_type, entity_id, name=None, parent_entity_dictionary=None):
        """
        Initialize a type of entity which supports properties.

        Properties are named integer variables associated with every entity of a certain type in the database.

        Properties can be associated with entity_type NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK, FACE_SET, ELEM_BLOCK,
        ELEM_SET, SIDE_SET, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP.
        """
        super(EntityWithProperty, self).__init__(database_id, entity_type, entity_id, name, parent_entity_dictionary)

    def property(self, property_name):
        """
        Return the value of a property on a given entity collection.

        Properties are integer values associated with an entity.

        :type property_name: str
        :return: property_value -- int
        """
        py_byte_string = property_name.encode('UTF-8')
        cdef char* c_string = py_byte_string
        cdef int value = 0
        if 0 != cexodus.ex_get_prop(self.ex_id, self.entity_type, self.entity_id, c_string, &value):
            _raise_io_error()
        return value

#-----------------------------------------------------------------------------------------------------------------------
class EntityWithVariable(Entity):

    def __init__(self, database_id, entity_type, entity_id, name=None, parent_entity_dictionary=None):
        """
        Initialize an entity type with association of variables.

        The entity_type is one of GLOBAL, NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
        FACE_SET, ELEM_BLOCK, ELEM_SET, SIDE_SET.

        The entity_type of ELEM_MAP, FACE_MAP, EDGE_MAP, and NODE_MAP do not have variables.

        A variable in Exodus is a scalar only. Multiple variables are required to make up either a vector or
        tensor field.
        """
        super(EntityWithVariable, self).__init__(database_id, entity_type, entity_id, name, parent_entity_dictionary)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def is_field_active(self):
        """
        Returns an array of one or zero values indicating whether each entity field exists on this entity.

        All variables on Nodal and Global entities are always active.

        :return: field_table -- numpy.ndarray of type np.int32_t, with shape len(num_vars,), or None if no
                                variables are present.
        """
        cdef int *var_ptr
        cdef int num_var = 0
        fields_dict = self.parent_entity_dictionary.fields
        cdef int num_fields = len(fields_dict)
        if num_fields < 1:
            return None
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type, &num_var):
            _raise_io_error()
        if num_var != sum(len(f.components) for f in fields_dict.itervalues()):  # consistency check
            AssertionError("the sum of all the field components is not equal to the number of variables")
        array_shape = num_fields
        cdef np.ndarray[np.int32_t, ndim=1] fields_table = np.empty(array_shape, dtype=np.int32)
        cdef int *fields_table_ptr = <int *> fields_table.data
        try:
            var_ptr = <int*> _allocate(sizeof(int) * num_var)
            if 0 != cexodus.ex_get_object_truth_vector(self.ex_id, self.entity_type, self.entity_id, num_var, var_ptr):
                _raise_io_error()
            for index, k, in enumerate(fields_dict):
                v = fields_dict[k]
                variables = v.variables
                if var_ptr[variables[0]] == 1:
                    for i in range(1,len(variables)):  # ensure all the any other component variables are also active
                        if var_ptr[variables[i]] != 1:
                            msg = "ERROR: Exodus database has a field with an inactive component:\n"
                            msg += "  database {} field {} component {} of {} (entity of type {})."
                            msg.format(self.ex_id, k, i, len(variables), self.entity_type)
                            raise AssertionError(msg)
                    fields_table_ptr[index] = 1
        finally:
            free(var_ptr)
        return fields_table

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def field(self, field, time_step):
        """
        Return array of field values of a single field at a time step for all entries.

        :param field: index of the desired field, indices start at 0, in range [0, num_fields-1]
        :type field: Field
        :param time_step: the time step index, at which the object field values are desired, the first time step is 0.
        :return: field_array -- array of values with shape (num_entries, len(field.components)), or None if
                                num_entries == 0 or field is not associated with ExodusII variable indices.
        :rtype: :class:`~FieldArray`
        """
        cdef int64_t i
        cdef int64_t j
        cdef int k
        if not field.name in self.parent_entity_dictionary.fields:
            raise KeyError("The given field {} does not exist on this entity.".format(field.name))
        cdef int num_components = len(field.components)
        if num_components == 1:
            return self.variable(field.variables[0], time_step)
        cdef int64_t num_entries = self.num_entries()
        if num_entries < 1:
            return None
        if field.variables is None:
            return None
        shape = (num_entries, num_components)
        cdef np.ndarray[np.double_t, ndim=2] array = np.empty(shape, dtype=np.double)
        cdef double *array_ptr = <double *> array.data
        cdef double* buffer
        try:
            buffer = <double*> _allocate(sizeof(double) * num_entries)
            k = 0
            while k < num_components:
                var_index = field.variables[k] + 1
                if 0 != cexodus.ex_get_var(self.ex_id, time_step+1, self.entity_type, var_index, self.entity_id,
                                           num_entries, buffer):
                    _raise_io_error()
                i = 0
                j = k
                while i < num_entries:
                    array_ptr[j] = buffer[i]
                    i += 1
                    j += num_components
                k += 1
        finally:
            free(buffer)
        return FieldArray(array, field)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def field_at_times(self, field, entry_index, start_time_step, stop_time_step):
        """
        Return array of values of a field for a single entry (node, element, edge, face)
        on a range of time steps.

        The entry_index has the following meaning, for a Nodal entity it is a global node index,
        for Block and Set this is a the local entry index.

        The entry_index is ignored for the Global entity.

        :param field: the field
        :type field: :class:`~exodus.Field`
        :param entry_index: index of the desired entry in the Entity, in range [0, num_entries-1]
        :param start_time_step: starting time step, in range [0, num_times]
        :type start_time_step: int
        :param stop_time_step: one past the end of the last time step, use -1 to indicate the final time step
        :type stop_time_step: int
        :return: var -- numpy.ndarray with
                        shape (stop_time_step - start_time_step, len(field.components)) if len(field.components) > 1,
                    and shape (stop_time_step - start_time_step) if len(field.components) == 1.
        """
        if not field.name in self.parent_entity_dictionary.fields:
            raise KeyError("The given field {} does not exist on this entity.".format(field.name))
        cdef int num_components = len(field.components)
        if num_components == 1:
            return self.variable_at_times(field.variables[0], entry_index, start_time_step, stop_time_step)
        if field.variables is None:
            return None
        cdef int64_t end_time_step
        if stop_time_step == -1:
            end_time_step = _inquire_value(self.ex_id, cexodus.EX_INQ_TIME)
        else:
            end_time_step = stop_time_step
        cdef int num_time = end_time_step - start_time_step
        if num_time < 1:
            return None
        cdef double *array_ptr
        cdef double* var_ptr
        shape = (num_time, num_components)
        cdef np.ndarray[np.double_t, ndim=2] array = np.empty(shape, dtype=np.double)
        array_ptr = <double *> array.data
        cdef int i, j, k
        try:
            var_ptr = <double*> _allocate(sizeof(double) * num_time)
            k = 0
            while k < num_components:
                var_index = field.variables[k] + 1
                if 0 != cexodus.ex_get_var_time(self.ex_id, self.entity_type, var_index, entry_index+1,
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
            free(var_ptr)
        return FieldArray(array, field)

    def partial_field(self, field, time_step, start_entry, stop_entry):
        """
        Return array of values of a single field on a subsection of the entries on an entity.

        :param field: the field
        :type field: :class:`~exodus.Field`
        :param time_step: the time time_step index, which the object variable values are desired,
                            the first time time_step is 0.
        :param start_entry: first entry at which field is read,
                            in the range of [0, num_entry-1].
        :param stop_entry:  one past the end of the last desired entry at which the field is read,
                            in the range of [start_entry+1, num_entry]
        :return: array -- numpy.ndarray of type np.double with shape (stop_entry - start_entry, len(field.components))
        """
        if not field.name in self.parent_entity_dictionary.fields:
            raise KeyError("The given field {} does not exist on this entity.".format(field.name))
        if field.variables is None:
            return None
        num_components = len(field.components)
        if num_components == 1:
            return self.partial_variable(field.variables[0], time_step, start_entry, stop_entry)
        cdef int64_t i
        cdef int64_t j
        cdef int64_t num_entries = stop_entry - start_entry
        if num_entries < 1:
            return None
        shape = (num_entries, num_components)
        cdef np.ndarray[np.double_t, ndim=2] array = np.empty(shape, dtype=np.double)
        cdef double *array_ptr = <double *> array.data
        cdef double* buffer
        try:
            buffer = <double*> _allocate(sizeof(double) * num_entries)
            for k in range(0,num_components):
                var_index = field.variables[k] + 1
                if 0 != cexodus.ex_get_partial_var(self.ex_id, time_step, self.entity_type, var_index,
                                                   self.entity_id, start_entry+1, num_entries, buffer):
                    _raise_io_error()
                j = k
                for i in range(num_entries):
                    array_ptr[j] = buffer[i]
                    j += num_components
        finally:
            free(buffer)
        return FieldArray(array, field)

    def is_variable_active(self):
        """
        Returns an array of one or zero values indicating whether each variable exists or not, respectively,
        for a single associated entity.

        All variables on Nodal and Global entities are always active.

        :return: var_table -- numpy.ndarray of type np.int32_t, with shape (num_vars,), or None if no
                            variables are present.
        """
        cdef int num_var = 0
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type, &num_var):
            _raise_io_error()
        if num_var < 1:
            return None
        array_shape = num_var
        cdef np.ndarray[np.int32_t, ndim=2] var_active = np.empty(array_shape, dtype=np.int32)
        cdef int *var_ptr = <int *> var_active.data
        if 0 != cexodus.ex_get_object_truth_vector (self.ex_id, self.entity_type, self.entity_id, num_var, var_ptr):
            _raise_io_error()
        return var_active

    def variable(self, variable_index, time_step):
        """
        Return array of values of a single variable at a given time step.

        :param variable_index: index of the desired variable, indices start at 0,
                                in range [0, num_variables-1]
        :param time_step: the time step index, which the object variable values are desired, the first time step is 0.
        :return: var -- numpy.ndarray with shape (num_entries,) or None.
        """
        cdef int64_t var_length = self.num_entries()
        if var_length < 1:
            return None
        cdef np.ndarray[np.double_t, ndim=1] var = np.empty(var_length, dtype=np.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_var(self.ex_id, time_step+1, self.entity_type, variable_index+1, self.entity_id,
                                   var_length, var_ptr):
            _raise_io_error()
        return var

    def variable_at_times(self, variable_index, entry_index, start_time_step, stop_time_step):
        """
        Return array of values of a single variable on a range of time steps.

        The entry_index has the following meaning, for a Nodal entity it is a global node index,
        for Block and Set this is a the local entry index.

        The entry_index is ignored for the Global entity.

        :param variable_index: index of the desired variable, indices start at 0,
                                in range [0, num_variables-1]
        :param entry_index: index of the desired entry in the Entity, in range [0, num_entries-1]
        :param start_time_step: starting time step, in range [0, num_times]
        :type start_time_step: int
        :param stop_time_step: one past the end of the last time step, use -1 to indicate the final time step
        :type stop_time_step: int
        :return: var -- numpy.ndarray with shape (stop_time_step - start_time_step) or None.
        """
        cdef int64_t num_entries = self.num_entries()
        cdef int64_t end_time_step
        if num_entries < 1:
            return None
        if stop_time_step == -1:
            end_time_step = _inquire_value(self.ex_id, cexodus.EX_INQ_TIME)
        else:
            end_time_step = stop_time_step
        cdef int64_t num_time = end_time_step - start_time_step
        if num_time < 1:
            return None
        array_shape = num_time
        cdef np.ndarray[np.double_t, ndim=1] var = np.empty(array_shape, dtype=np.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_var_time(self.ex_id, self.entity_type, variable_index+1, entry_index+1,
                                        start_time_step+1, end_time_step, var_ptr):
            _raise_io_error()
        return var

    def partial_variable(self, var_index, time_step, start_entry, stop_entry):
        """
        Return array of values of a single variable on a subsection of the entries on an entity.

        :param var_index: index of the desired variable, indices start at 0,
                            in range [0, num_variables-1]
        :param time_step: the time time_step index, which the object variable values are desired,
                            the first time time_step is 0.
        :param start_entry: first entry at which variables are read,
                            in the range of [0, num_entry-1].
        :param stop_entry:  one past the end of the last desired entry at which the variable is read,
                            in the range of [start_entry+1, num_entry]
        :return: var -- numpy.ndarray of type np.double with shape (stop_entry - start_entry)
        """
        cdef int64_t num_entries = stop_entry - start_entry
        if num_entries < 1:
            return None
        cdef np.ndarray[np.double_t, ndim=1] var = np.empty(num_entries, dtype=np.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_partial_var(self.ex_id, time_step, self.entity_type, var_index, self.entity_id, start_entry+1,
                                           num_entries, var_ptr):
            _raise_io_error()
        return var

#-----------------------------------------------------------------------------------------------------------------------
class EntityWithAttribute(EntityWithVariable):

    def __init__(self, database_id, entity_type, entity_id, name=None, parent_entity_dictionary=None):
        """
        Initialize an entity type with association of attributes.

        The entity_type is one of NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
            FACE_SET, ELEM_BLOCK, ELEM_SET, SIDE_SET.

        The entity_type of GLOBAL, ELEM_MAP, FACE_MAP, EDGE_MAP, and NODE_MAP do not have attributes.

        An attribute is a per entity scalar variable in Exodus.
        """
        super(EntityWithAttribute, self).__init__(database_id, entity_type, entity_id, name, parent_entity_dictionary)

    def num_attributes(self):
        """
        Return the number of attributes on a single collection entity.

        Only NODAL, EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK, NODE_SET, EDGE_SET, FACE_SET, SIDE_SET, ELEM_SET have
        attributes.

        :return: num_attribute -- integer
        """
        cdef int num_attribute = 0
        if self.entity_type not in (GLOBAL, ELEM_MAP, FACE_MAP, EDGE_MAP, NODE_MAP):
            if 0 != cexodus.ex_get_attr_param(self.ex_id, self.entity_type, self.entity_id, &num_attribute):
                _raise_io_error()
        return num_attribute

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def attribute_names(self):
        """
        Return the attribute names on a single collection entity.

        Attribute names were added to the Exodus library in version 4.26; earlier databases will return empty strings.

        :return: names -- list of strings
        """
        cdef char** name_ptr
        cdef int num_attribute = 0
        len_str = _inquire_value(self.ex_id, cexodus.EX_INQ_MAX_READ_NAME_LENGTH)
        db_name_size = _inquire_value(self.ex_id, cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
        if db_name_size < len_str:
            len_str = db_name_size
        if self.entity_type not in (GLOBAL, ELEM_MAP, FACE_MAP, EDGE_MAP, NODE_MAP):
            if 0 != cexodus.ex_get_attr_param(self.ex_id, self.entity_type, self.entity_id, &num_attribute):
                _raise_io_error()
        names = []
        if num_attribute > 0:
            try:
                name_ptr = <char**> _allocate(sizeof(char*) * num_attribute)
                try:
                    for i in range(num_attribute):
                        name_ptr[i] = <char *> _allocate(sizeof(char) * (len_str+1))
                    if 0 != cexodus.ex_get_attr_names(self.ex_id, self.entity_type, self.entity_id, name_ptr):
                        _raise_io_error()
                    for i in range(num_attribute):
                        name = _to_unicode(name_ptr[i])
                        names.append(name)
                finally:
                    for i in range(num_attribute):
                        free(name_ptr[i])
            finally:
                free(name_ptr)
        return names

    def attributes(self):
        """
        Return all the attribute values for all the entries in the given collection.

        :return: attributes -- numpy.ndarray of type double with shape (num_entries, num_attribute)
        """
        cdef int64_t num_entries = self.length()
        cdef int64_t num_attribute = self.num_attribute()
        if num_attribute < 1:
            return None
        array_shape = (num_entries, num_attribute)
        cdef np.ndarray[np.double_t, ndim=2] attributes = np.empty(array_shape, dtype=np.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_attr(self.ex_id, self.entity_type, self.entity_id, attributes_ptr):
                _raise_io_error()
        return attributes

    def partial_attributes(self, start, stop):
        """
        Return all the attribute values for a subsection of the entries in the given collection.

        :param start: first local entry index at which attributes are desired, in the range [0, num_entry - 1].
        :param stop:  one past the end of the last desired local entry index, in the range [start+1, num_entry]
        :return: attributes -- numpy.ndarray of type double with shape (stop - start, num_attribute)
        """
        cdef int64_t num_entries = stop - start
        cdef int64_t num_attribute = self.num_attribute()
        if num_attribute < 1 or num_entries < 1:
            return None
        array_shape = (num_entries, num_attribute)
        cdef np.ndarray[np.double_t, ndim=2] attributes = np.empty(array_shape, dtype=np.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_partial_attr(self.ex_id, self.entity_type, self.entity_id, start+1,
                                            num_entries, attributes_ptr):
                _raise_io_error()
        return attributes

    def attribute(self, attrib_index):
        """
        Return an array of values of a single attribute for all entries in the given collection.

        :param attrib_index: desired attribute in range [0, num_attribute-1]
        :type attrib_index: int
        :return: attributes -- numpy.ndarray of type double with shape num_entries
        """
        cdef int64_t num_entries = self.length()
        if self.entity_type in (GLOBAL, ELEM_MAP, FACE_MAP, EDGE_MAP, NODE_MAP):
            return None
        cdef np.ndarray[np.double_t, ndim=1] attributes = np.empty(num_entries, dtype=np.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_one_attr(self.ex_id, self.entity_type, self.entity_id, attrib_index, attributes_ptr):
                _raise_io_error()
        return attributes

    def partial_attribute(self, start, stop, attrib_index):
        """
        Read one entity attribute for a subsection of entries in the given collection.

        :param start: first local entry index at which attributes are desired, in the range [0, num_entry - 1].
        :param stop:  one past the end of the last desired local entry index, in the range [start+1, num_entry]
        :param attrib_index: desired attribute in range [0, num_attribute-1]
        :type attrib_index: int
        :return: attributes -- numpy.ndarray of type double with shape (stop - start)
        """
        cdef int64_t num_entries = stop - start
        if num_entries < 1 or self.entity_type in (GLOBAL, ELEM_MAP, FACE_MAP, EDGE_MAP, NODE_MAP):
            return None
        cdef np.ndarray[np.double_t, ndim=1] attributes = np.empty(num_entries, dtype=np.double)
        cdef double *attributes_ptr = <double *> attributes.data
        if 0 != cexodus.ex_get_partial_one_attr(self.ex_id, self.entity_type, self.entity_id, start+1, num_entries,
                                                attrib_index, attributes_ptr):
                _raise_io_error()
        return attributes

#-----------------------------------------------------------------------------------------------------------------------
class EntityDictionary(collections.MutableMapping):

    def __init__(self, database_id, entity_type):
        """
        Ordered dictionary of entity ID -> Entity.

        On creation, the dictionary is initialized with keys that are the IDs of Blocks, Maps, Sets,
        from the database and the values are initialized to None.

        No keys or values exist when entity_type is Nodal or Global.

        Subclasses are responsible for implementing __getitem__, __setitem__, __delitem__,
        __iter__.
        """
        self.ex_id = database_id
        self.entity_type = entity_type
        if entity_type == NODAL or entity_type == GLOBAL:
            self._store = None
        else:
            # initialize dictionary with (entity ID, None) key, value pairs
            self._store = collections.OrderedDict.fromkeys(self._entity_ids(), None)

    @abstractmethod
    def __getitem__(self, key):
        raise KeyError

    @abstractmethod
    def __setitem__(self, key, value):
        raise KeyError

    @abstractmethod
    def __delitem__(self, key):
        raise KeyError

    @abstractmethod
    def __iter__(self):
        raise KeyError

    @abstractmethod
    def __len__(self):
        raise KeyError

    def _num_entity(self):
        cdef int64_t num = 0
        if self.entity_type == ELEM_BLOCK:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM_BLK)
        elif self.entity_type == NODE_SET:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_NODE_SETS)
        elif self.entity_type == SIDE_SET:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_SIDE_SETS)
        elif self.entity_type == ELEM_MAP:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM_MAP)
        elif self.entity_type == NODE_MAP:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_NODE_MAP)
        elif self.entity_type == EDGE_BLOCK:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE_BLK)
        elif self.entity_type == FACE_BLOCK:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_FACE_BLK)
        elif self.entity_type == EDGE_SET:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE_SETS)
        elif self.entity_type == FACE_SET:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_FACE_SETS)
        elif self.entity_type == ELEM_SET:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM_SETS)
        elif self.entity_type == EDGE_MAP:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE_MAP)
        elif self.entity_type == FACE_MAP:
            num = _inquire_value(self.ex_id, cexodus.EX_INQ_FACE_MAP)
        elif self.entity_type == GLOBAL:
            num = 1
        elif self.entity_type == NODAL:
            num = 1
        else:
            raise ValueError("Invalid entity type.")
        return num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def _entity_ids(self):
        """
        Return list of entity IDs.

        :return: entity_ids -- list of positive integer IDs, or empty list if entity_type is GLOBAL or NODAL.
        """
        cdef int64_t num_entity
        cdef int64_t *entity_ids_ptr64
        cdef int *entity_ids_ptr32
        ids = []
        if self.entity_type == NODAL or self.entity_type == GLOBAL:
            return ids
        num_entity = self._num_entity()
        if num_entity > 0:
            ids = [-1] * num_entity # fill with -1
            if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_IDS_INT64_API:
                try:
                    entity_ids_ptr64 = <int64_t*> _allocate(sizeof(int64_t) * num_entity)
                    if 0 != cexodus.ex_get_ids(self.ex_id, self.entity_type, entity_ids_ptr64):
                        _raise_io_error()
                    for i in range(num_entity):
                        ids[i] = entity_ids_ptr64[i]
                finally:
                    free(entity_ids_ptr64)
            else:
                try:
                    entity_ids_ptr32 = <int*> _allocate(sizeof(int) * num_entity)
                    if 0 != cexodus.ex_get_ids(self.ex_id, self.entity_type, entity_ids_ptr32):
                        _raise_io_error()
                    for i in range(num_entity):
                        ids[i] = entity_ids_ptr32[i]
                finally:
                    free(entity_ids_ptr32)
        return ids

    def name(self, entity_id):
        """
        Return the name of the association or collection.

        :return: name -- str
        """
        cdef char* name_ptr
        name = ""
        if self.entity_type != GLOBAL and self.entity_type != NODAL:
            len_str = _inquire_value(self.ex_id, cexodus.EX_INQ_MAX_READ_NAME_LENGTH)
            db_name_size = _inquire_value(self.ex_id, cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
            if db_name_size < len_str:
                len_str = db_name_size
            try:
                name_ptr = <char *> _allocate(sizeof(char) * (len_str+1))
                if 0 != cexodus.ex_get_name(self.ex_id, self.entity_type, entity_id, &name_ptr[0]):
                    _raise_io_error()
                name = _to_unicode(name_ptr)
            finally:
                free(name_ptr)
        return name

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def names(self):
        """
        Return all the names of a specific type of collection entity, such as block, set, or map.

        :return: names -- a list of strings, or None if entity_type is GLOBAL or NODAL
        """
        cdef char** name_ptr
        if self.entity_type == NODAL or self.entity_type == GLOBAL:
            return None
        len_str = _inquire_value(self.ex_id, cexodus.EX_INQ_MAX_READ_NAME_LENGTH)
        db_name_size = _inquire_value(self.ex_id, cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
        if db_name_size < len_str:
            len_str = db_name_size
        cdef int64_t num_objects = self.num_entity()
        names = []
        if num_objects > 0:
            try:
                name_ptr = <char**> _allocate(sizeof(char*) * num_objects)
                try:
                    for i in range(num_objects):
                        name_ptr[i] = <char *> _allocate(sizeof(char) * (len_str+1))
                    if 0 != cexodus.ex_get_names(self.ex_id, self.entity_type, name_ptr):
                        _raise_io_error()
                    for i in range(num_objects):
                        name = _to_unicode(name_ptr[i])
                        names.append(name)
                finally:
                    for i in range(num_objects):
                        free(name_ptr[i])
            finally:
                free(name_ptr)
        return names

    def sum_entries(self):
        """
        Return the total count of entries (elements, faces, edges, nodes) summed over all
        collection entities of this type.

        For example, return the sum of the number of elements included in all element blocks, or
        return the sum of the number of faces in all face sets, etc.

        For entity_type GLOBAL or NODAL, the return value is always 1.

        For entity_type (NODE_MAP, EDGE_MAP, FACE_MAP and ELEM_MAP) the return value is
        always the number of global (nodes, edges, faces, and elements), respectively.

        For entity_type SIDE_SET there a tuple of size two is returned (num_elem, num_node), the number of element
        sides in the side set and the number of nodes, respectively.

        :return: sum_entries -- integer, the number of the concatenated entries in that type, or for SIDE_SET
                                (sum_elem, sum_node) the concatenated number of elements and nodes
        """
        cdef int64_t sum_entries = 0
        cdef int64_t num_nodes = 0
        if self.entity_type == NODAL or self.entity_type == GLOBAL:
            sum_entries = 1
        elif self.entity_type == NODE_SET:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_NS_NODE_LEN)
        elif self.entity_type == SIDE_SET:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_SS_ELEM_LEN)
            num_node = _inquire_value(self.ex_id, cexodus.EX_INQ_SS_NODE_LEN)
            return sum_entries, num_node
        elif self.entity_type == NODE_MAP:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_NODES) # returns regardless of map existing
        elif self.entity_type == EDGE_MAP:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE)
        elif self.entity_type == FACE_MAP:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_FACE)
        elif self.entity_type == ELEM_MAP:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM)
        elif self.entity_type == EDGE_SET:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_ES_LEN)
        elif self.entity_type == FACE_SET:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_FS_LEN)
        elif self.entity_type == ELEM_SET:
            sum_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_ELS_LEN)
        else:
            raise ValueError("Invalid entity type.")
        return sum_entries

#-----------------------------------------------------------------------------------------------------------------------
class EntityDictionaryWithProperty(EntityDictionary):

    def __init__(self, database_id, entity_type):
        """
        Initialize a type of entity list which supports properties. Properties are integer variables associated
        with every entity in the list.

        Properties are named integer variables associated with every entity of a certain type in the database.

        Properties can be associated with entity_type NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK, FACE_SET, ELEM_BLOCK,
        ELEM_SET, SIDE_SET, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP.
        Properties are not associated with entity_type GLOBAL, NODAL.

        :param database_id:
        :param entity_type:
        :return:
        """
        super(EntityDictionaryWithProperty, self).__init__(database_id, entity_type)

    # Note: Nodal and Global do not subclass from this---they don't have properties---therefore we
    # can implement some of the MutableMapping members for this case. Subclasses must still implement
    # __getitem__

    def __setitem__(self, key, value):
        self._store[key] = value

    def __delitem__(self, key):
        del self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return self._store.__len__()

    def num_properties(self):
        """
        Return the number of properties for a type of collection entity.

        :return: num -- int, the number of integer properties
        """
        return cexodus.ex_get_num_props(self.ex_id, self.entity_type)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def property_names(self):
        """
        Return all the property names present on this entity_type in the database.

        :return: names -- a list of strings
        """
        cdef char** name_ptr
        len_str = _inquire_value(self.ex_id, cexodus.EX_INQ_MAX_READ_NAME_LENGTH)
        db_name_size = _inquire_value(self.ex_id, cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
        if db_name_size < len_str:
            len_str = db_name_size
        cdef int64_t num_entity = self.num_entity()
        names = []
        if num_entity > 0:
            try:
                name_ptr = <char**> _allocate(sizeof(char*) * num_entity)
                try:
                    for i in range(num_entity):
                        name_ptr[i] = <char *> _allocate(sizeof(char) * (len_str+1))
                    if 0 != cexodus.ex_get_prop_names(self.ex_id, self.entity_type, name_ptr):
                        _raise_io_error()
                    for i in range(num_entity):
                        name = _to_unicode(name_ptr[i])
                        names.append(name)
                finally:
                    for i in range(num_entity):
                        free(name_ptr[i])
            finally:
                free(name_ptr)
        return names

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def property(self, property_name):
        """
        Return the values of an integer property on all entities of this type in the database.

        :param property_name: name of the property
        :type property_name: str
        :return: properties -- list of integer values
        """
        cdef int num_properties = cexodus.ex_get_num_props(self.ex_id, self.entity_type)
        if num_properties < 1:
            return None
        py_byte_string = property_name.encode('UTF-8')
        cdef char* c_string = py_byte_string
        cdef int* val_ptr
        properties = []
        try:
            val_ptr = <int*> _allocate(sizeof(int) * num_properties)
            if 0 != cexodus.ex_get_prop_array(self.ex_id, self.entity_type, c_string, val_ptr):
                _raise_io_error()
            for i in range(num_properties):
                properties.append(val_ptr[i])
        finally:
            free(val_ptr)
        return properties

#-----------------------------------------------------------------------------------------------------------------------
class EntityDictionaryWithVariable(EntityDictionary):

    def __init__(self, database_id, entity_type):
        """
        Initialize an entity list that has associated variables.

        The entity_type is one of GLOBAL, NODAL, NODE_SET, EDGE_BLOCK, EDGE_SET, FACE_BLOCK,
        FACE_SET, ELEM_BLOCK, ELEM_SET, SIDE_SET.

        The entity_type of ELEM_MAP, FACE_MAP, EDGE_MAP, and NODE_MAP do not have variables.

        A variable in Exodus is a scalar only. Multiple variables are required to make up either a vector or
        tensor field.
        """
        super(EntityDictionaryWithVariable, self).__init__(database_id, entity_type)
        self._fields_dictionary = None

    @property
    def fields(self):
        """
        A dictionary of :class:`~exodus.Field` objects describing the current multi-component fields on this
        collection of Entity's.

        Fields are a grouping of the variables by common name prefix; whatever follows the last '_' underscore in the
        variable names implicitly becomes the component names.

        The keys of this dictionary are field names, and the values of the dictionary are instances of
        :class:`~exodus.Field`.

        :Example:

            >>> print 'number of nodal fields: {}'.format(len(exo.nodal.fields))
            >>> print 'element field names: {}'.format(exo.globals.fields.keys())
            >>> for f in exo.element_blocks.fields.itervalues():
            >>>     print f

        :return: fields_dict -- collections.OrderedDict of Field items, with the name of the field as the key.
        """
        if self._fields_dictionary is None:
            self._fields_dictionary = Fields(self)
        return self._fields_dictionary

    def num_variables(self):
        """
        Return the number of variables associated with all the entities of this type in the database.

        :return: num_vars -- number of variables on this entity_type
        """
        cdef int num_vars = 0
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type, &num_vars):
            _raise_io_error()
        return num_vars

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def variable_names(self):
        """
        Return a list of variable names associated with a given type of entity.

        :return: variable_names -- list of strings
        """
        cdef char** var_name_ptr
        cdef int len_str = cexodus.MAX_STR_LENGTH + 1
        cdef int num_vars = 0
        if 0 != cexodus.ex_get_variable_param(self.ex_id, self.entity_type, &num_vars):
            _raise_io_error()
        variable_names = []
        if num_vars > 0:
            try:
                var_name_ptr = <char**> _allocate(sizeof(char*) * num_vars)
                try:
                    for i in range(num_vars):
                        var_name_ptr[i] = <char *> _allocate(sizeof(char) * len_str)
                    if 0 != cexodus.ex_get_variable_names(self.ex_id, self.entity_type, num_vars, var_name_ptr):
                        _raise_io_error()
                    for i in range(num_vars):
                        variable_name = _to_unicode(var_name_ptr[i])
                        variable_names.append(variable_name)
                finally:
                    for i in range(num_vars):
                        free(var_name_ptr[i])
            finally:
                free(var_name_ptr)
        return variable_names

    def variable_name(self, variable_index):
        """
        Return the name for a variable associated with this entity_type.

        A variable is scalar only. Multiple variables make up a vector or tensor field.

        :param variable_index: variable index in the range [0,num_variables - 1]
        :return: variable_name -- str
        """
        cdef char var_name[cexodus.MAX_STR_LENGTH+1]
        if 0 != cexodus.ex_get_variable_name(self.ex_id, self.entity_type, variable_index+1, &var_name[0]):
            _raise_io_error()
        variable_name = _to_unicode(&var_name[0])
        return variable_name

    def is_variable_active_all(self):
        """
        Returns a table of one or zero values indicating whether each variable exists, or not,
        on all entities of this type on the database.

        :return: var_table -- numpy.ndarray of type np.int32, with shape (num_entity, num_vars), or None if no
                            variables are present.
        """
        cdef int num_entity = self.num_entity()
        if num_entity < 1:
            return None
        cdef int num_var = self.num_variables()
        if num_var < 1:
            return None
        array_shape = (num_entity, num_var)
        cdef np.ndarray[np.int32_t, ndim=2] var_table = np.empty(array_shape, dtype=np.int32)
        cdef int *var_ptr = <int *> var_table.data
        if 0 != cexodus.ex_get_truth_table(self.ex_id, self.entity_type, num_entity, num_var, var_ptr):
            _raise_io_error()
        return var_table

#-----------------------------------------------------------------------------------------------------------------------
class Global(EntityDictionaryWithVariable, EntityWithVariable):

    def __init__(self, database_id):
        """
        Initialize the global entity on a Database which acts as an EntityDictionary with only a single Entity member, which
        is itself. The number of entries in the Global entity is also equal to 1.

        :param database_id: ID of the Exodus database
        """
        EntityDictionaryWithVariable.__init__(self, database_id, GLOBAL)
        EntityWithVariable.__init__(self, database_id, GLOBAL, -1, name=None, parent_entity_dictionary=self)

    @classmethod
    def _key_error(cls):
        raise KeyError("{} does not have an ID (key) and there can be only one instance per database."
                       .format(cls.__name__))

    def __getitem__(self, key):
        return self

    def __setitem__(self, key, value):
        self.__class__._key_error()

    def __delitem__(self, key):
        self.__class__._key_error()

    def __iter__(self):
        self.__class__._key_error()

    def __len__(self):
        return 1

    def dimension(self):
        """
        Return the spatial dimension of the coordinates in the database.
        :return: ndim -- integer dimension
        """
        return _inquire_value(self.ex_id, cexodus.EX_INQ_DIM)

    def has_id_map(self, map_type):
        """
        Return whether an ID map of the given type exists on the database.

        :param map_type: type of map, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP
        :return: True or False
        :raises ValueError: if the map_type is invalid
        """
        cdef extern from "netcdf.h":
            int nc_inq_varid(int ncid, char *name, int *varidp)
        cdef int exoid
        cdef int mapid
        cdef char *vmap
        cdef bytes py_bytes
        if map_type == NODE_MAP:
            py_bytes = "node_num_map".encode()
        elif map_type == EDGE_MAP:
            py_bytes = "edge_num_map".encode()
        elif map_type == FACE_MAP:
            py_bytes = "face_num_map".encode()
        elif map_type == ELEM_MAP:
            py_bytes = "elem_num_map".encode()
        else:
            raise ValueError("Bad map_type is not one of ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP.")
        exoid = self.ex_id
        vmap = py_bytes
        # look for the netcdf map variable and check for error condition
        if nc_inq_varid(exoid, vmap, &mapid) != 0:
            return False
        else:
            return True

    def id_map(self, map_type):
        """
        Return the array of global integer IDs for all entries of the given type.

        This is the single map of the IDs of elements, nodes, edges, or faces (completely separate from the other list
        of maps. If no ID map exists, the default enumeration is returned, for example, [0, num_elements] when
        map_type == ELEM_MAP.

        :param map_type: type of map, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP
        :return: map_entries -- numpy.ndarray with shape num_elements, num_nodes, num_edges, or num_faces; or None
                                if no entries exist in the database
        :raises ValueError: if the map_type is invalid
        """
        cdef int64_t num_entries
        cdef np.ndarray[np.int64_t, ndim=1] map_entries64
        cdef int64_t *map_ptr64
        cdef np.ndarray[np.int32_t, ndim=1] map_entries32
        cdef int *map_ptr32
        if map_type == NODE_MAP:
            num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_NODES)
        elif map_type == ELEM_MAP:
            num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM)
        elif map_type == EDGE_MAP:
            num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE)
        elif map_type == FACE_MAP:
            num_entries = _inquire_value(self.ex_id, cexodus.EX_INQ_FACE)
        else:
            raise ValueError("Invalid map type.")
        if num_entries < 1:
            return None
        if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_MAPS_INT64_API:
            map_entries64 = np.empty(num_entries, dtype=np.int64)
            map_ptr64 = <int64_t *> map_entries64.data
            if 0 != cexodus.ex_get_id_map(self.ex_id, map_type, map_ptr64):
                _raise_io_error()
            return map_entries64
        else:
            map_entries32 = np.empty(num_entries, dtype=np.int32)
            map_ptr32 = <int *> map_entries32.data
            if 0 != cexodus.ex_get_id_map(self.ex_id, map_type, map_ptr32):
                _raise_io_error()
            return map_entries32

    def num_elements(self):
        """Return the global number of element entries in the database."""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_ELEM)

    def num_edges(self):
        """
        Return the global number of edge entries in the database.

        Most often the number of edges is zero even when num_nodes or num_elem is positive.
        """
        return _inquire_value(self.ex_id, cexodus.EX_INQ_EDGE)

    def num_faces(self):
        """
        Return the global number of edge entries in the database.

        Most often the number of faces is zero even when num_nodes or num_elem is positive.
        """
        return _inquire_value(self.ex_id, cexodus.EX_INQ_FACE)

    def num_nodes(self):
        """Return the global number of node entries in the database."""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_NODES)

    def num_times(self):
        """Return the number of time steps on the database"""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_TIME)

    def times(self):
        """
        Return all time values.
        :return: times -- numpy.ndarray of type double_t with shape (num_times), or None if no time steps are present.
        """
        cdef int64_t array_shape = _inquire_value(self.ex_id, cexodus.EX_INQ_TIME)
        if array_shape < 1:
            return None
        cdef np.ndarray[np.double_t, ndim=1] times = np.empty(array_shape, dtype=np.double)
        cdef double *times_ptr = <double *> times.data
        if 0 != cexodus.ex_get_all_times(self.ex_id, times_ptr):
            _raise_io_error()
        return times

    def variables(self, step):
        """
        Return an array of values of all the global variables at a single time step.

        :param step: the time step index, which the variable values are desired, the first time step is 0.
        :return: var -- numpy.ndarray with shape (num_variables), or None if there are no variables.
        """
        cdef int64_t var_length = self.num_variables()
        if var_length < 1:
            return None
        cdef np.ndarray[np.double_t, ndim=1] var = np.empty(var_length, dtype=np.double)
        cdef double *var_ptr = <double *> var.data
        if 0 != cexodus.ex_get_glob_vars(self.ex_id, step, var_length, var_ptr):
            _raise_io_error()
        return var


#-----------------------------------------------------------------------------------------------------------------------
class Nodal(EntityDictionaryWithVariable, EntityWithAttribute):

    def __init__(self, database_id):
        """
        Initialize the Nodal entity on a Database which acts a both an entity, and a list of entities, whose sole
        member is itself.

        :param database_id: ID of the Exodus database.
        """
        EntityDictionaryWithVariable.__init__(self, database_id, NODAL)
        EntityWithAttribute.__init__(self, database_id, NODAL, entity_id=-1, name=None, parent_entity_dictionary=self)

    @classmethod
    def _key_error(cls):
        raise KeyError("{} does not have an ID (key) and there can be only one instance per database."
                       .format(cls.__name__))

    def __getitem__(self, key):
        return self

    def __setitem__(self, key, value):
        self.__class__._key_error()

    def __delitem__(self, key):
        self.__class__._key_error()

    def __iter__(self):
        self.__class__._key_error()

    def __len__(self):
        return 1

    def num_coordinate_frames(self):
        """Return the number of coordinate frames."""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_COORD_FRAMES)

    def coordinate_names(self):
        """
        Return names of the coordinate axes.

        :return: list of strings of the names of the coordinate axes
        :raises RuntimeError: if the dimension is not in the range [0,3]
        """
        cdef int64_t ndim = _inquire_value(self.ex_id, cexodus.EX_INQ_DIM)
        cdef char *coord_names_ptr[3];
        if ndim < 0 or ndim > 3:
            raise RuntimeError("unexpected spatial dimension = {}".format(ndim))
        coord_names = []
        try:
            for i in range(ndim):
                coord_names_ptr[i] = <char*> _allocate(sizeof(char) * cexodus.MAX_STR_LENGTH)
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
    def coordinates(self):
        """
        Return the global node coordinates in an array.

        :return: coords -- numpy.ndarray with shape (num_nodes, num_dim).
        :raises RuntimeError: if the dimension is not in the range [0,3]
        """
        cdef int64_t i = 0
        cdef int64_t j = 0
        array_shape = None
        cdef int64_t nnodes = _inquire_value(self.ex_id, cexodus.EX_INQ_NODES)
        cdef int64_t ndim = _inquire_value(self.ex_id, cexodus.EX_INQ_DIM)
        if ndim == 1:
            array_shape = nnodes
        elif ndim == 2:
            array_shape = (nnodes,2)
        elif ndim == 3:
            array_shape = (nnodes,3)
        else:
            raise RuntimeError('Unexpected num_dim = {}'.format(ndim))
        cdef np.ndarray[np.double_t, ndim=2] coords = np.empty(array_shape, dtype=np.double)
        cdef double *coords_ptr = <double *> coords.data
        cdef double* coords_buffer
        try:
            coords_buffer = <double*> _allocate(sizeof(double) * nnodes)
            # x-coord
            if 0 != cexodus.ex_get_coord(self.ex_id, coords_buffer, NULL, NULL):
                _raise_io_error()
            j = 0
            i = 0
            while i < nnodes:
                coords_ptr[j] = coords_buffer[i]
                i += 1
                j += ndim
            if ndim > 1:  # y-coord
                if 0 != cexodus.ex_get_coord(self.ex_id, NULL, coords_buffer, NULL):
                    _raise_io_error()
                j = 1
                i = 0
                while i < nnodes:
                    coords_ptr[j] = coords_buffer[i]
                    i += 1
                    j += ndim
                if ndim > 2:  # z-coord
                    if 0 != cexodus.ex_get_coord(self.ex_id, NULL, NULL, coords_buffer):
                        _raise_io_error()
                    j = 2
                    i = 0
                    while i < nnodes:
                        coords_ptr[j] = coords_buffer[i]
                        i += 1
                        j += ndim
        finally:
            free(coords_buffer)
        return coords

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def partial_coordinates(self, start, stop):
        """
        Return a range of the global node coordinates in an array.

        :param start: first node at which coordinates are read, in the range [0, num_nodes-1].
        :param stop:  one past the end of the last desired node, in the range [start+1, num_nodes]
        :return: coords -- numpy.ndarray of type np.double with shape (stop - start, num_dim)
        :raises RuntimeError: if the dimension is not in the range [0,3]
        """
        cdef int64_t i = 0
        cdef int64_t j = 0
        array_shape = None
        cdef int64_t nnodes = _inquire_value(self.ex_id, cexodus.EX_INQ_NODES)
        cdef int64_t ndim = _inquire_value(self.ex_id, cexodus.EX_INQ_DIM)
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
            raise RuntimeError('Unexpected num_dim = {}'.format(ndim))
        cdef np.ndarray[np.double_t, ndim=2] coords = np.empty(array_shape, dtype=np.double)
        cdef double *coords_ptr = <double *> coords.data
        cdef double* coords_buffer
        try:
            coords_buffer = <double*> _allocate(sizeof(double) * len_nodes)
            # x-coord
            if 0 != cexodus.ex_get_partial_coord(self.ex_id, start_node, len_nodes, coords_buffer, NULL, NULL):
                _raise_io_error()
            j = 0
            for i in xrange(len_nodes):
                coords_ptr[j] = coords_buffer[i]
                j += ndim
            if ndim > 1:  # y-coord
                if 0 != cexodus.ex_get_partial_coord(self.ex_id, start_node, len_nodes, NULL, coords_buffer, NULL):
                    _raise_io_error()
                j = 1
                for i in xrange(nnodes):
                    coords_ptr[j] = coords_buffer[i]
                    j += ndim
                if ndim > 2:  # z-coord
                    if 0 != cexodus.ex_get_partial_coord(self.ex_id, start_node, len_nodes, NULL, NULL, coords_buffer):
                        _raise_io_error()
                    j = 2
                    for i in xrange(nnodes):
                        coords_ptr[j] = coords_buffer[i]
                        j += ndim
        finally:
            free(coords_buffer)
        return coords

    def displacement_name(self):
        """
        Return the name of the nodal field acting as the displacements.

        :return: displacements_name -- str, or None if no suitable field for displacements was found
        """
        ndim = -1
        for c in ('displacement', 'displace', 'displ', 'disp'):
            if c in self.fields:
                if ndim == -1:
                    ndim = _inquire_value(self.ex_id, cexodus.EX_INQ_DIM)
                if len(self.fields[c].components) == ndim:
                    return c
        return None

    def displacement(self, time_step):
        """
        Return the array of displaced coordinates at a time step for all nodes.

        This is equivalent to adding the array returned by :meth:`coordinates` to the :class:`~exodus.FieldArray`
        corresponding to the name returned by :meth:`displacement_name`.

        :param time_step: the time step index, at which the object displacements values are desired,
                    the first time step is 0; a value of -1 will return the displacements at the last time step.
        :return: array -- numpy.ndarray with shape (num_nodes, num_components) or None if no displacements exist.
        """
        # TODO: implement Nodal.displacement
        raise NotImplementedError()

    def partial_displacement(self, name, time_step, start_entry, stop_entry):
        # TODO: implement Nodal.partial_displacement
        raise NotImplementedError()

#-----------------------------------------------------------------------------------------------------------------------
class Blocks(EntityDictionaryWithVariable, EntityDictionaryWithProperty):
    
    def __init__(self, database_id, block_type):
        """
        The collection of all :class:`~exodus.Block` entities of a specific type on a Database.

        :param database_id: ID of the Exodus database.
        :param block_type: the type of the entity, one of ELEM_BLOCK, EDGE_BLOCK, or FACE_BLOCK
        :raises ValueError: if the block_type is invalid
        """
        super(Blocks, self).__init__(database_id, block_type)
        if block_type not in (EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK):
            raise ValueError("The block_type is invalid.")

    def __getitem__(self, key):
        value = self._store.get(key)        # throws KeyError if the ID is incorrect
        if value is None:
            value = self._get_block(key)    # create the value if it does not yet exist
            self._store[key] = value        # put the value in the dictionary
        return value

    def _get_block(self, block_id):
        """
        Construct a block from the database given the block ID.

        :param block_id: the ID of the block, one of those returned by entity_ids()
        :return: block -- Block initialized with name, topology, sizes of entries and attributes, etc.
        """
        cdef cexodus.ex_block block
        block.id = block_id
        block.type = self.entity_type
        if 0 != cexodus.ex_get_block_param(self.ex_id, &block):
            _raise_io_error()
        ndim = _inquire_value(self.ex_id, cexodus.EX_INQ_DIM)
        topology_name = _topology_name(_to_unicode(block.topology), block.num_nodes_per_entry, ndim)
        name = self.name(block_id)
        return Block(self.ex_id, self.entity_type, block_id, name, topology_name, block.num_entry,
                     block.num_nodes_per_entry, block.num_edges_per_entry, block.num_faces_per_entry,
                     parent_entity_dictionary=self)

#-----------------------------------------------------------------------------------------------------------------------
class Block(EntityWithAttribute, EntityWithProperty):

    def __init__(self, database_id, block_type, block_id, name, topology_name,
                 num_entries, num_nodes_per_entry, num_edges_per_entry, num_faces_per_entry, attribute_names=None,
                 parent_entity_dictionary=None):
        """
        A block of entries of all the same topology with connectivity information.

        :param block_type: entity type of the block, one of EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK
        :type block_id: int
        :param block_id: unique integer identifier
        :type block_id: int
        :param name: name of the block
        :type name: str
        :param topology_name: topology_name of the entries in the block
        :type topology_name: str
        :param num_entries: number of the entries in this block: nodes, edges, faces, or elements
        :type num_entries: int
        :param num_nodes_per_entry: number of nodes per topology_name
        :type num_nodes_per_entry: int
        :param num_edges_per_entry: number of edges per topology_name
        :type num_edges_per_entry: int
        :param num_faces_per_entry: number of edges per topology_name
        :type num_faces_per_entry: int
        :param attribute_names: list or tuple of attribute names
        :type attribute_names: list or tuple of strings
        :raises ValueError: if the block_type is invalid
        """
        super(Block, self).__init__(database_id, block_type, block_id, name, parent_entity_dictionary)
        if block_type not in (EDGE_BLOCK, FACE_BLOCK, ELEM_BLOCK):
            raise ValueError("The block_type is invalid.")
        self._num_entries = num_entries
        self._num_nodes_per_entry = num_nodes_per_entry
        self._num_edges_per_entry = num_edges_per_entry
        self._num_faces_per_entry = num_faces_per_entry
        self._topology_name = topology_name
        self._attribute_names = None
        if attribute_names:
            assert isinstance(attribute_names, (list, tuple))
            self._attribute_names = list(attribute_names)

    def __str__(self):
        """Return a human readable representation of the block."""
        lookup = { EDGE_BLOCK:'edges', FACE_BLOCK:'faces', ELEM_BLOCK:'elements' }
        s = '{{block {} '.format(self.entity_id)
        if self._name != '':
            s += '{} '.format(self._name)
        s += 'with {} {} {}'.format(self._num_entries, self._topology_name, lookup[self.entity_type])
        if self._num_nodes_per_entry > 0:
            s += ' with {} nodes'.format(self._num_nodes_per_entry)
        if self._num_edges_per_entry > 0:
            s += ' {} edges'.format(self._num_edges_per_entry)
        if self._num_faces_per_entry > 0:
            s += ' {} edges'.format(self._num_faces_per_entry)
        if self._attribute_names is not None:
            s += ' attributes {}'.format(str(self._attribute_names))
        s += '}'
        return s

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def connectivity(self, entry_type=NODE, zero_based=True):
        """
        Return an array of the entries-to-node connectivity for a given block.

        The enumeration of the connected entries are zero based, unless zero_based is set to False, in which case
        the entries begin at one.

        :param entry_type: (default NODE) type of entries returned in the connectivity array, one of
                            NODE, EDGE, FACE
        :param zero_based: if True (default) the enumeration of the connected entries begins at zero
        :return: conn -- numpy.ndarray of integer type with shape (num_entry, num_conn_entries), for example,
                            when entry_type is NODE, conn[i,j] is the jth node of the ith entry in the block
        :raises ValueError: if the connect entry_type does not make sense
        """
        cdef int64_t num_entries
        cdef np.ndarray[np.int64_t, ndim=2] conn64
        cdef int64_t *conn_ptr64
        cdef np.ndarray[np.int32_t, ndim=2] conn32
        cdef int *conn_ptr32
        cdef retval = -1
        cdef int64_t index
        if self._num_entries < 1:
            return None
        num_entries = self._num_entries
        if entry_type == NODE:
            num_conn_entries = self._num_nodes_per_entry
        elif entry_type == EDGE:
            num_conn_entries = self._num_edges_per_entry
        elif entry_type == FACE:
            num_conn_entries = self._num_faces_per_entry
        else:
            raise ValueError("The connected entry_type is not one of NODE, EDGE, FACE.")
        if num_conn_entries == 0:
            return None
        array_shape = (num_entries, num_conn_entries)
        if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_MAPS_INT64_API:
            conn64 = np.empty(array_shape, dtype=np.int64)
            conn_ptr64 = <int64_t *> conn64.data
            if entry_type == NODE:
                retval = cexodus.ex_get_conn(self.ex_id, self.entity_type, self.entity_id, conn_ptr64, NULL, NULL)
            elif entry_type == EDGE:
                retval = cexodus.ex_get_conn(self.ex_id, self.entity_type, self.entity_id, NULL, conn_ptr64, NULL)
            elif entry_type == FACE:
                retval = cexodus.ex_get_conn(self.ex_id, self.entity_type, self.entity_id, NULL, NULL, conn_ptr64)
            else:
                raise ValueError("The connected entry_type is not one of NODE, EDGE, FACE.")
            if retval != 0:
                _raise_io_error()
            if zero_based:
                index = num_entries * num_conn_entries
                while index:
                    index -= 1
                    conn_ptr64[index] -= 1
            return conn64
        else:
            conn32 = np.empty(array_shape, dtype=np.int32)
            conn_ptr32 = <int *> conn32.data
            if entry_type == NODE:
                retval = cexodus.ex_get_conn(self.ex_id, self.entity_type, self.entity_id, conn_ptr32, NULL, NULL)
            elif entry_type == EDGE:
                retval = cexodus.ex_get_conn(self.ex_id, self.entity_type, self.entity_id, NULL, conn_ptr32, NULL)
            elif entry_type == FACE:
                retval = cexodus.ex_get_conn(self.ex_id, self.entity_type, self.entity_id, NULL, NULL, conn_ptr32)
            else:
                raise ValueError("The connected entry_type is not one of NODE, EDGE, FACE.")
            if retval != 0:
                _raise_io_error()
            if zero_based:
                for i in xrange(num_entries * num_conn_entries):
                    conn_ptr32[i] -= 1
            return conn32

    def partial_connectivity(self, start, stop, entry_type=NODE, zero_based=True):
        """
        Return an array of a subsection of the entries-to-node connectivity for a given block.

        The enumeration of the connected entries are zero based, unless zero_based is set to False, in which case
        the entries begin at one.

        :param start: first local element index at which nodes are desired, in the range [0, num_entry - 1].
        :param stop:  one past the end of the last desired local entry index, in the range [start+1, num_entry]
        :param entry_type: (default NODE) type of entries returned in the connectivity array, one of
                            NODE, EDGE, FACE
        :param zero_based: if True (default) the enumeration of the connected entries begins at zero
        :return: conn -- numpy.ndarray of integer type with shape (num_entry, num_conn_entries), for example,
                            when entry_type is NODE, conn[i,j] is the jth node of the ith entry in the block
        :raises ValueError: if the start, stop range is invalid
        :raises ValueError: if the connected entry_type does not make sense
        """
        cdef int64_t num_entries
        cdef np.ndarray[np.int64_t, ndim=1] conn64
        cdef int64_t *conn_ptr64
        cdef np.ndarray[np.int32_t, ndim=1] conn32
        cdef int *conn_ptr32
        cdef retval = -1
        cdef int64_t index
        num_entries = self._num_entries
        if num_entries < 1:
            return None
        if start < 0 or stop > num_entries:
            raise ValueError("The (start, stop) parameters are out of range [{}, {}]".format(0,self._num_entries))
        num_entries = stop - start
        if entry_type == NODE:
            num_conn_entries = self._num_nodes_per_entry
        elif entry_type == EDGE:
            num_conn_entries = self._num_edges_per_entry
        elif entry_type == FACE:
            num_conn_entries = self._num_faces_per_entry
        else:
            raise ValueError("The connected entry_type is not one of NODE, EDGE, FACE.")
        if num_conn_entries == 0:
            return None
        array_shape = (num_entries, num_conn_entries)
        if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_MAPS_INT64_API:
            conn64 = np.empty(array_shape, dtype=np.int64)
            conn_ptr64 = <int64_t *> conn64.data
            if entry_type == NODE:
                retval = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type, self.entity_id, start+1, stop,
                                                     conn_ptr64, NULL, NULL)
            elif entry_type == EDGE:
                retval = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type, self.entity_id, start+1, stop,
                                                     NULL, conn_ptr64, NULL)
            elif entry_type == FACE:
                retval = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type, self.entity_id, start+1, stop,
                                                     NULL, NULL, conn_ptr64)
            else:
                raise ValueError("The connected entry_type is not one of NODE, EDGE, FACE.")
            if retval != 0:
                _raise_io_error()
            if zero_based:
                index = num_entries * num_conn_entries
                while index:
                    index -= 1
                    conn_ptr64[index] -= 1
            return conn64
        else:
            conn32 = np.empty(array_shape, dtype=np.int32)
            conn_ptr32 = <int *> conn32.data
            if entry_type == NODE:
                retval = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type, self.entity_id, start+1, stop,
                                                     conn_ptr32, NULL, NULL)
            elif entry_type == EDGE:
                retval = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type, self.entity_id, start+1, stop,
                                                     NULL, conn_ptr32, NULL)
            elif entry_type == FACE:
                retval = cexodus.ex_get_partial_conn(self.ex_id, self.entity_type, self.entity_id, start+1, stop,
                                                     NULL, NULL, conn_ptr32)
            else:
                raise ValueError("The connected entry_type is not one of NODE, EDGE, FACE.")
            if retval != 0:
                _raise_io_error()
            if zero_based:
                for i in range(num_entries * num_conn_entries):
                    conn_ptr32[i] -= 1
            return conn32

#-----------------------------------------------------------------------------------------------------------------------
class Maps(EntityDictionaryWithProperty):

    def __init__(self, database_id, map_type):
        """
        The collection of all :class:`~exodus.Map` of entities of a specific type on a Database.

        An ordered dictionary of all ( map_id, :class:`~exodus.Map` ) pairs of a specific map_type.

        :param database_id: ID of the Exodus database.
        :param map_type: the type of the entity list, one of ELEM_BLOCK, EDGE_BLOCK, or FACE_BLOCK
        :raises ValueError: if the map_type is invalid
        """
        super(Maps, self).__init__(database_id, map_type)

        if map_type not in (ELEM_MAP, NODE_MAP, EDGE_MAP, FACE_MAP):
            raise ValueError("The map_type is invalid.")

    def __getitem__(self, key):
        value = self._store.get(key)                        # throws KeyError if the ID is incorrect
        if value is None:
            value = Map(self.ex_id, self.entity_type, key, parent_entity_dictionary=self)
            self._store[key] = value                        # put the value in the dictionary
        return value

#-----------------------------------------------------------------------------------------------------------------------
class Map(EntityWithProperty):

    def __int__(self, database_id, map_type, map_id, name=None, parent_entity_dictionary=None):
        """
        A map entity renumbering the entries of a certain type.

        :param database_id: ID of the Exodus database.
        :param map_type: the type of the map entity, one of type of map, ELEM_MAP, NODE_MAP, EDGE_MAP, or FACE_MAP
        :raises ValueError: if the map_type is invalid
        """
        super(Map, self).__init__(database_id, map_type, map_id, name, parent_entity_dictionary)

        if map_type not in (ELEM_MAP, FACE_MAP, EDGE_MAP, NODE_MAP):
            raise ValueError("The map_type is invalid.")

    def entries(self, zero_based=True):
        """
        Return the integer entries of a given map.

        :return: map -- numpy.ndarray with shape num_entries().
        """
        cdef np.ndarray[np.int64_t, ndim=1] map_entries64
        cdef int64_t *map_ptr64
        cdef np.ndarray[np.int32_t, ndim=1] map_entries32
        cdef int *map_ptr32
        cdef int64_t num_entries = self.num_entries()
        cdef int64_t i
        if num_entries < 1:
            return None
        if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_MAPS_INT64_API:
            map_entries64 = np.empty(num_entries, dtype=np.int64)
            map_ptr64 = <int64_t *> map_entries64.data
            if 0 != cexodus.ex_get_num_map(self.ex_id, self.entity_type, self.entity_id, map_ptr64):
                _raise_io_error()
            if zero_based:
                i = 0
                while i < num_entries:
                    map_ptr64[i] -= 1
                    i += 1
            return map_entries64
        else:
            map_entries32 = np.empty(num_entries, dtype=np.int32)
            map_ptr32 = <int *> map_entries32.data
            if 0 != cexodus.ex_get_num_map(self.ex_id, self.entity_type, self.entity_id, map_ptr32):
                _raise_io_error()
            if zero_based:
                for i in range(num_entries):
                    map_ptr32[i] -= 1
            return map_entries32

    def partial_entries(self, start, stop, zero_based=True):
        """
        Return a section of the the integer entries of a given map.

        :param start: first entry, in range [0, len_entity(map_type)-1]
        :param stop: one past the last desired entry, in range [start+1, len_entity(map_type)]
        :return: map_entries -- numpy.ndarray with shape (stop - start).
        """
        cdef int64_t num_entries
        cdef np.ndarray[np.int64_t, ndim=1] map_entries64
        cdef int64_t *map_ptr64
        cdef np.ndarray[np.int32_t, ndim=1] map_entries32
        cdef int *map_ptr32
        cdef int64_t i
        num_entries = stop - start
        if num_entries < 1:
            return None
        if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_MAPS_INT64_API:
            map_entries64 = np.empty(num_entries, dtype=np.int64)
            map_ptr64 = <int64_t *> map_entries64.data
            if 0 != cexodus.ex_get_partial_num_map(self.ex_id, self.entity_type, self.entity_id,
                                                   start+1, stop, map_ptr64):
                _raise_io_error()
            if zero_based:
                i = 0
                while i < num_entries:
                    map_ptr64[i] -= 1
                    i += 1
            return map_entries64
        else:
            map_entries32 = np.empty(num_entries, dtype=np.int32)
            map_ptr32 = <int *> map_entries32.data
            if 0 != cexodus.ex_get_partial_num_map(self.ex_id, self.entity_type, self.entity_id,
                                                   start+1, stop, map_ptr32):
                _raise_io_error()
            if zero_based:
                for i in range(num_entries):
                    map_ptr32[i] -= 1
            return map_entries32

#-----------------------------------------------------------------------------------------------------------------------
class Sets(EntityDictionaryWithVariable, EntityDictionaryWithProperty):

    def __init__(self, database_id, set_type):
        """
        The collection of all :class:`~exodus.Set` entities of a specific type on a Database.

        :param database_id: ID of the Exodus database.
        :param set_type: the type of the entity list, one of ELEM_SET, NODE_SET, SIDE_SET, EDGE_SET, or FACE_SET
        :raises ValueError: if the set_type is invalid
        """
        super(Sets, self).__init__(database_id, set_type)
        if set_type not in (SIDE_SET, NODE_SET, ELEM_SET, EDGE_SET, FACE_SET):
            raise ValueError("The set_type is invalid.")

    def __getitem__(self, key):
        value = self._store.get(key)    # throws KeyError if the ID is incorrect
        if value is None:
            value = self._get_set(key)  # create the value if it does not yet exist
            self._store[key] = value    # put the value in the dictionary
        return value

    def _get_set(self, set_id):
        """
        Construct a set from the database given the set ID.

        :param set_id: the ID of the block, one of those returned by entity_ids()
        :return: set -- Set initialized with name, sizes of entries and distribution factors, etc.
        """
        num_entries, num_dist_fact = _set_param(self.ex_id, self.entity_type, set_id)
        name = self.name(set_id)
        return Set(self.ex_id, self.entity_type, set_id, name, num_entries, num_dist_fact,
                   parent_entity_dictionary=self)

    def num_distribution_factors_all(self):
        """
        Return the the total length of distribution factors over all existing sets of a given type.

        :return: count -- integer, the length of the concatenated set distribution factors
        """
        cdef int64_t count = 0
        if self.entity_type == SIDE_SET:
            count  = _inquire_value(self.ex_id, cexodus.EX_INQ_SS_DF_LEN)
        elif self.entity_type == NODE_SET:
            count  = _inquire_value(self.ex_id, cexodus.EX_INQ_NS_DF_LEN)
        elif self.entity_type == ELEM_SET:
            count  = _inquire_value(self.ex_id, cexodus.EX_INQ_ELS_DF_LEN)
        elif self.entity_type == EDGE_SET:
            count  = _inquire_value(self.ex_id, cexodus.EX_INQ_ES_DF_LEN)
        elif self.entity_type == FACE_SET:
            count = _inquire_value(self.ex_id, cexodus.EX_INQ_FS_DF_LEN)
        else:
            raise RuntimeError("Invalid type of set, not one of SIDE_SET, NODE_SET, ELEMENT_SET, EDGE_SET, FACE_SET")
        return count

#-----------------------------------------------------------------------------------------------------------------------
class Set(EntityWithAttribute, EntityWithProperty):

    def __init__(self, database_id, set_type, set_id, name, num_entries, num_dist_fact, attribute_names=None,
                 parent_entity_dictionary=None):
        """
        Initialize a set of entries from the database.

        :param set_type: type of set, one of NODE_SET, SIDE_SET, EDGE_SET, FACE_SET, ELEM_SET
        :param set_id: the ID of the set
        :param name: name of this set
        :type name: str
        :param num_entries: number of the entries in this set: nodes, sides, edges, faces, or elements
        :type num_entries: int
        :param num_dist_fact: number of distribution factors in this set
        :type num_dist_fact: int
        :param attribute_names: list or tuple of attribute names
        :type attribute_names: list or tuple of strings
        """
        super(Set, self).__init__(database_id, set_type, set_id, name, parent_entity_dictionary)
        self._num_entries = num_entries
        self._num_dist_fact = num_dist_fact
        if attribute_names:
            assert isinstance(attribute_names, (list, tuple))
            self._attribute_names = list(attribute_names)

    def entries(self, zero_based=True):
        """
        Return the list of entries, an array of numbers that describe the members of the set.

        The shape and content of the array that is returned is different when it is a SIDE_SET.

        ======== =============== ================================================================
        set_type shape           content
        ======== =============== ================================================================
        SIDE_SET (num_entries,2) [i,1] is the the ith element number, [i,2] is the jth local side
        NODE_SET num_entries     [i] is the ith node number in the set
        ELEM_SET num_entries     [i] is the ith element number in the set
        FACE_SET num_entries     [i] is the ith face number in the set
        EDGE_SET num_entries     [i] is the ith edge number in the set
        ======== =============== ================================================================

        :return: set_entries -- numpy.ndarray of integers, with shape num_entries for NODE_SET and ELEM_SET, and
                                FACE_SET, and EDGE_SET, and shape (num_entries,2) for SIDE_SET.
        """
        cdef int64_t *set_entry_ptr64
        cdef int64_t *set_extra_ptr64
        cdef int64_t *set_entries_and_extra_ptr64
        cdef int *set_entry_ptr32
        cdef int *set_extra_ptr32
        cdef int *set_entries_and_extra_ptr32
        cdef np.ndarray[np.int64_t, ndim=1] set_entries64
        cdef np.ndarray[np.int32_t, ndim=1] set_entries32
        cdef np.ndarray[np.int64_t, ndim=2] set_entries_and_extra64
        cdef np.ndarray[np.int32_t, ndim=2] set_entries_and_extra32
        cdef int64_t num_entries = self.num_entries()
        if num_entries < 1:
            return None
        if cexodus.ex_int64_status(self.ex_id) & cexodus.EX_BULK_INT64_API:
            if self.entity_type == SIDE_SET:
                try:
                    set_entry_ptr64 = <int64_t*> _allocate(sizeof(int64_t) * num_entries)
                    set_extra_ptr64 = <int64_t*> _allocate(sizeof(int64_t) * num_entries)
                    if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type, self.entity_id,
                                               set_entry_ptr64, set_extra_ptr64):
                        raise _raise_io_error()
                    set_entries_and_extra64 = np.empty((num_entries,2), dtype=np.int64)
                    set_entries_and_extra_ptr64 = <int64_t *> set_entries_and_extra64.data
                    if zero_based:
                        for i in range(num_entries):
                            k = 2*i
                            set_entries_and_extra_ptr64[k]   = set_entry_ptr64[i] - 1
                            set_entries_and_extra_ptr64[k+1] = set_extra_ptr64[i] - 1
                    else:
                        for i in range(num_entries):
                            k = 2*i
                            set_entries_and_extra_ptr64[k]   = set_entry_ptr64[i]
                            set_entries_and_extra_ptr64[k+1] = set_extra_ptr64[i]
                finally:
                    free(set_entry_ptr64)
                    free(set_extra_ptr64)
                return set_entries_and_extra64
            else:
                set_entries64 = np.empty(num_entries, dtype=np.int64)
                set_entry_ptr64 = <int64_t *> set_entries64.data
                if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type, self.entity_id, set_entry_ptr64, NULL):
                    raise _raise_io_error()
                if zero_based:
                    for i in range(num_entries):
                        set_entry_ptr64[i] -= 1
                return set_entries64
        else: # 32 bit integer version
            if self.entity_type == SIDE_SET:
                try:
                    set_entry_ptr32 = <int *> _allocate(sizeof(int) * num_entries)
                    set_extra_ptr32 = <int *> _allocate(sizeof(int) * num_entries)
                    if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type, self.entity_id,
                                               set_entry_ptr32, set_extra_ptr32):
                        raise _raise_io_error()
                    set_entries_and_extra32 = np.empty((num_entries,2), dtype=np.int32)
                    set_entries_and_extra_ptr32 = <int *> set_entries_and_extra32.data
                    if zero_based:
                        for i in range(num_entries):
                            k = 2*i
                            set_entries_and_extra_ptr32[k]   = set_entry_ptr32[i] - 1
                            set_entries_and_extra_ptr32[k+1] = set_extra_ptr32[i] - 1
                    else:
                        for i in range(num_entries):
                            k = 2*i
                            set_entries_and_extra_ptr32[k]   = set_entry_ptr32[i]
                            set_entries_and_extra_ptr32[k+1] = set_extra_ptr32[i]
                finally:
                    free(set_entry_ptr32)
                    free(set_extra_ptr32)
                return set_entries_and_extra32
            else:
                set_entries32 = np.empty(num_entries, dtype=np.int32)
                set_entry_ptr32 = <int *> set_entries32.data
                if 0 != cexodus.ex_get_set(self.ex_id, self.entity_type, self.entity_id, set_entry_ptr32, NULL):
                    raise _raise_io_error()
                if zero_based:
                    for i in range(num_entries):
                        set_entry_ptr32[i] -= 1
                return set_entries32

    def num_distribution_factors(self):
        """
        Return the count of distribution factors on this set.

        :return: num_dist_fact -- int
        """
        return self._num_dist_fact

    def distribution_factors(self):
        """
        Return the distribution factors for the entries in this set.

        Distribution factors are a special scalar variable that can be associated with a set collection type.
        There is either zero or one distribution factors on a given set.

        :return: factors -- numpy.ndarray of type np.double with shape num_entries(), or None if no
                            distribution factors exist
        """
        cdef int64_t num_entry_in_set = self._num_entries
        cdef int num_dist_fact_in_set = self._num_dist_fact
        if num_dist_fact_in_set < 1 or num_entry_in_set < 1:
            return None
        cdef np.ndarray[np.double_t, ndim=1] factors = np.empty(num_entry_in_set, dtype=np.double)
        cdef double *factors_ptr = <double *> factors.data
        if 0 != cexodus.ex_get_set_dist_fact(self.ex_id, self.entity_type, self.entity_id, factors_ptr):
            _raise_io_error()
        return factors

    def partial_distribution_factors(self, start, stop):
        """
        Return the distribution factors for a subsection of the entries in a given set.

        :param start: first entry at which variables are read, must be in the range of [0, num_entry-1].
        :param stop:  one past the end of the last desired entry at which the factor is read, must be in the
                            range of [start+1, num_entry]
        :return: var -- numpy.ndarray of type np.double with shape (stop - start), or None if no distribution factor
                        exists
        """
        cdef int64_t num_entry_in_set = self._num_entries
        cdef int num_dist_fact_in_set = self._num_dist_fact
        cdef int64_t len_dist_fact = stop - start
        if num_dist_fact_in_set < 1 or num_entry_in_set < 1 or len_dist_fact < 1:
            return None
        cdef np.ndarray[np.double_t, ndim=1] factors = np.empty(len_dist_fact, dtype=np.double)
        cdef double *factors_ptr = <double *> factors.data
        if 0 != cexodus.ex_get_partial_set_dist_fact(self.ex_id, self.entity_type, self.entity_id,
                                                     start+1, stop, factors_ptr):
            _raise_io_error()
        return factors

#-----------------------------------------------------------------------------------------------------------------------
class DatabaseFile(object):

    def __init__(self, file_path, mode = OPEN_READ_ONLY):
        """A context manager for opening and using :class:`~exodus.Database`.

        Use the `with` statement while instantiating this class to easily open a Database, operate on it, and ensure the
        Database is closed properly after all read/writes are complete. There is no need to surround your code with a
        try/finally block in order to ensure the Database is closed after an error occurs.

        :Example:

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     coord = e.nodal.coordinates()
            >>> # file automatically now closed

        :param file_path: The path to the file to open or create.
        :type file_path: str
        :param mode: Either exodus.READ or exodus.WRITE.
        :type mode: int
        :raises: IOError
        """
        self.database = Database(file_path, mode)

    def __enter__(self):
        return self.database
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.database.close()
        return False

#-----------------------------------------------------------------------------------------------------------------------
cdef class Database:
    """
    An Exodus database for reading/wring Exodus format files.
    """
    
    # these are not class attributes because this is a cdef extension type
    cdef int ex_id
    cdef public object path
    cdef public int mode
    cdef float _version
    cdef public int _max_name_length
    cdef object _globals
    cdef object _nodal
    cdef object _blocks
    cdef object _maps
    cdef object _sets

    def __cinit__(self, path, mode = OPEN_READ_ONLY):
        # __cinit__ is called prior to any __init__ method.
        # Be careful what you do in the __cinit__() method, because the object
        # may not yet be fully valid Python object when it is called. Therefore,
        # you should be careful invoking any Python operations which might touch
        # the object; in particular, its methods.
        self.ex_id = 0
        # check existence and read-ability/write-ability in python in order catch things before we
        # enter the Exodus library so that we can raise more informative exceptions
        if mode == OPEN_READ_ONLY:
            _assert_is_file(path, True)
            _assert_file_access(path, os.R_OK)
        elif mode == OPEN_READ_WRITE:
            _assert_is_file(path, True)
            _assert_file_access(path, os.R_OK)
            _assert_file_access(path, os.W_OK)
        elif mode == OPEN_CREATE:
            _assert_is_file(path, False)
            _assert_file_access(path, os.W_OK)
        elif mode == OPEN_REPLACE:
            _assert_file_access(path, os.W_OK)
        else:
            abs_path = path
            try:
                abs_path = os.path.abspath(path)
            except OSError:
                pass
            finally:
                message = 'Unknown mode for file {}.'.format(abs_path)
                raise ValueError(message)
        cdef int io_ws = 8
        cdef int comp_ws = 8
        cdef float db_version
        self.ex_id = cexodus.ex_open_int(path, mode,
                                       &io_ws, &comp_ws,
                                       &db_version,
                                       cexodus.EX_API_VERS_NODOT)
        if self.ex_id < 0:
            _raise_io_error()
        self.path = path
        self.mode = mode
        self._version = db_version
        self._globals = None
        self._nodal = None
        self._blocks = {}
        self._maps = {}
        self._sets = {}
        # set what 64 bit integers in the inquire API by default, however, not the
        #       inquire values, EX_INQ_INT64_API
        #       maps (id, order, ...), EX_MAPS_INT64_API
        #       entity ids (sets, blocks, maps), EX_IDS_INT64_API
        #       or bulk data (local indices, counts, maps) EX_BULK_INT64_API
        #       or all EX_ALL_INT64_API
        cexodus.ex_set_int64_status(self.ex_id, cexodus.EX_ALL_INT64_API)

    def __init__(self, path, mode = OPEN_READ_ONLY):
        """
        __init__(path, mode = OPEN_READ_ONLY)
        Initialize a Database object with the given path and mode.

        The mode is one of four exclusive values (which may not be OR'ed together).

        =============== ===============================================
        mode            allowed actions on the database
        =============== ===============================================
        OPEN_READ_ONLY  read if it exists but not modify
        OPEN_READ_WRITE read and/or append and modify if it exists
        OPEN_CREATE     write only if it does not already exist
        OPEN_REPLACE    write and destroy if it already exists
        =============== ===============================================

        :param path: The path to the file to open or create.
        :type path: str
        :param mode: one of the four values of modes, default value is OPEN_READ_ONLY if none are given
        :type mode: int
        
        :raises: IOError depending on values of mode and whether the file exists or not
        :raises ValueError: if the mode is invalid
        """
        # we now have a valid id, since __cinit__ would have been already called
        # arguments were only used by the __cinit__ method
        pass

    def __del__(self):
        """Cleanup any resources and close the ExodusII database"""
        self.close()

    def __str__(self):
        """Return a human readable representation of the object."""
        return '{{ id: {}, filename: {}, version: {:.6g}}}'.format(self.ex_id, self.path, self._version)
    
    def summary(self):
        """
        Return a human readable string representation of the global parameters.
        """
        cdef cexodus.ex_init_params param
        if 0 != cexodus.ex_get_init_ext(self.ex_id, &param):
            _raise_io_error()
        return ('{{  title:         {},\n   num_dim:       {},\n'
                '   num_node:      {},\n   num_elem:      {},\n   num_face:      {},\n   num_edge:      {},\n'
                '   num_elem_blk:  {},\n   num_face_blk:  {},\n   num_edge_blk:  {},\n'
                '   num_node_sets: {},\n   num_side_sets: {},\n   num_elem_sets: {},\n   num_face_sets: {},\n'
                '   num_edge_sets: {},\n'
                '   num_node_maps: {},\n   num_elem_maps: {},\n   num_face_maps: {},\n   num_edge_maps: {}}}'
               ).format(_to_unicode(param.title), param.num_dim,
                        param.num_nodes, param.num_elem, param.num_face, param.num_edge,
                        param.num_elem_blk, param.num_face_blk, param.num_edge_blk,
                        param.num_node_sets, param.num_side_sets, param.num_elem_sets, param.num_face_sets,
                        param.num_edge_sets,
                        param.num_node_maps, param.num_elem_maps, param.num_face_maps, param.num_edge_maps
                       )

    def _len_max_names(self):
        cdef int max_all_name_length = <int> cexodus.ex_inquire_int(self.ex_id,
                                                                    cexodus.EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH)
        if max_all_name_length < 0:
            _raise_io_error()
        cdef int max_use_name_length = <int> cexodus.ex_inquire_int(self.ex_id, cexodus.EX_INQ_DB_MAX_USED_NAME_LENGTH)
        if max_use_name_length < 0:
            _raise_io_error()
        #print "File {} can use at most {}-character names".format(self.path, max_all_name_length)
        #print "File {} used no more than {}-character names".format(self.path, max_use_name_length)
        max_name_length = max_use_name_length
        if 0 != cexodus.ex_set_max_name_length(self.ex_id, max_use_name_length):
            _raise_io_error()
        return max_use_name_length

    def close(self):
        """
        Close an existing Exodus file database.

        You should never have to call this as it is called automatically when the object goes out of scope,
        i.e., the __del__() method is automatically called.
        """
        if 0 != cexodus.ex_close(self.ex_id):
            _raise_io_error()
        self.ex_id = -1
        self.path = None

    def get_file_type(self):
        """Return ExodusII file type"""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_FILE_TYPE)
    
    def version(self):
        """
        Return ExodusII database version number.

        The database version number reflects the version of the library that was used to write the file.

        :return: database_version -- str
        """
        return self._version

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @property
    def info_records(self):
        """
        Information records---a list of text strings stored in the database.

        :return: records -- list of string
        """
        cdef char** inf_record_ptr
        cdef int len_str = cexodus.MAX_LINE_LENGTH
        cdef int64_t num_inf_rec = _inquire_value(self.ex_id, cexodus.EX_INQ_INFO)
        records = []
        if num_inf_rec > 0:
            try:
                inf_record_ptr = <char**> _allocate(sizeof(char*) * num_inf_rec)
                try:
                    for i in range(num_inf_rec):
                        inf_record_ptr[i] = <char *> _allocate(sizeof(char) * (len_str+1))
                    if 0 != cexodus.ex_get_info(self.ex_id, inf_record_ptr):
                        _raise_io_error()
                    for i in range(num_inf_rec):
                        record = _to_unicode(inf_record_ptr[i])
                        records.append(record)
                finally:
                    for i in range(num_inf_rec):
                        free(inf_record_ptr[i])
            finally:
                free(inf_record_ptr)
        return records

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @property
    def qa_records(self):
        """
        List of lists of the quality assurance records.
        
        Records are returned as one concatenated list of strings, in multiples of four sequentially:
        
        1. Code name, an application code that has created or modified the database.
        2. Code QA descriptor, a location for a version identifier of the
           application code.
        3. Date in the format 01/25/93.
        4. Time in the 24 hour format hours:minutes:seconds, such as 16:30:15.

        :return: records -- list of string with len = 4 * num_qa_records
        """
        cdef char** qa_record_ptr
        cdef int len_str = cexodus.MAX_STR_LENGTH
        cdef int64_t num_qa_rec = _inquire_value(self.ex_id, cexodus.EX_INQ_QA)
        cdef int64_t length_qa = 4 * num_qa_rec
        records = []
        if num_qa_rec > 0:
            try:
                qa_record_ptr = <char**> _allocate(sizeof(char*) * length_qa)
                try:
                    for i in range(length_qa):
                        qa_record_ptr[i] = <char *> _allocate(sizeof(char) * (len_str+1))
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
    def globals(self):
        """
        A property that is the :class:`~exodus.Global` object for accessing global coordinates, time planes, global
        variables, total number of elements and nodes, and global id maps.
        """
        if self._globals is None:
            self._globals = Global(self.ex_id)
        return self._globals

    @property
    def nodal(self):
        """
        A property that is the :class:`~exodus.Nodal` object for accessing global nodal coordinates, nodal variables,
        and nodal attributes.

        :Example:

        Get a bounding box (minimum, maximum) over all the coordinates in the database.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     coord = e.nodal.coordinates()
            >>>     mins  = numpy.min(coord, axis=0)
            >>>     maxes = numpy.max(coord, axis=0)
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
    def element_blocks(self):
        """
        A property that is the :class:`~exodus.Blocks` object for accessing the collection of entities of
        Block(entity_type = ELEM_BLOCK).

        :Example:

        Iterate over the elements blocks in a database and print the block ID's and number of elements in the block.

            >>> with DatabaseFile('/tmp/myExodusFile.exo') as e:
            >>>     print 'database has {} element blocks'.format(len(e.element_blocks))
            >>>     for id, block in e.element_blocks.iteritems():
            >>>         print 'block_{} num_elements = {}'.format(id, block.num_entries)
        """
        return self._get_blocks(ELEM_BLOCK)

    @property
    def face_blocks(self):
        """
        A property that is the :class:`~exodus.Blocks` object for accessing the collection of entities of
        Block(entity_type = FACE_BLOCK).
        """
        return self._get_blocks(FACE_BLOCK)

    @property
    def edge_blocks(self):
        """
        A property that is the :class:`~exodus.Blocks` object for accessing the collection of entities of
        Block(entity_type = EDGE_BLOCK).
        """
        return self._get_blocks(EDGE_BLOCK)

    def _get_maps(self, map_type):
        value = self._maps.get(map_type, None)
        if value is None:
            value = Maps(self.ex_id, map_type)
            self._maps[map_type] = value
        return value
    
    @property
    def element_maps(self):
        """
        A property that is the :class:`~exodus.Maps` object for accessing the collection of entities of
        Map(entity_type = ELEM_MAP).
        """
        return self._get_maps(ELEM_MAP)

    @property
    def face_maps(self):
        """
        A property that is the :class:`~exodus.Maps` object for accessing the collection of entities of
        Map(entity_type = FACE_MAP).
        """
        return self._get_maps(FACE_MAP)

    @property
    def edge_maps(self):
        """
        A property that is the :class:`~exodus.Maps` object for accessing the collection of entities of
        Map(entity_type = EDGE_MAP).
        """
        return self._get_maps(EDGE_MAP)

    @property
    def node_maps(self):
        """
        A property that is the :class:`~exodus.Maps` object for accessing the collection of entities of
        Map(entity_type = NODE_MAP).
        """
        return self._get_maps(NODE_MAP)

    def _get_sets(self, set_type):
        value = self._sets.get(set_type, None)
        if value is None:
            value = Sets(self.ex_id, set_type)
            self._sets[set_type] = value
        return value

    @property
    def element_sets(self):
        """
        A property that is the :class:`~exodus.Sets` object for accessing the collection of entities of
        Set(entity_type = ELEM_SET).
        """
        return self._get_sets(ELEM_SET)

    @property
    def face_sets(self):
        """
        A property that is the :class:`~exodus.Sets` object for accessing the collection of entities of
        Set(entity_type = FACE_SET).
        """
        return self._get_sets(FACE_SET)

    @property
    def edge_sets(self):
        """
        A property that is the :class:`~exodus.Sets` object for accessing the collection of entities of
        Set(entity_type = EDGE_SET).
        """
        return self._get_sets(EDGE_SET)

    @property
    def node_sets(self):
        """
        A property that is the :class:`~exodus.Sets` object for accessing the collection of entities of
        Set(entity_type = NODE_SET).
        """
        return self._get_sets(NODE_SET)

    @property
    def side_sets(self):
        """
        A property that is the :class:`~exodus.Sets` object for accessing the collection of entities of
        Set(entity_type = SIDE_SET).
        """
        return self._get_sets(SIDE_SET)

    #-------------------------------------------------------------------------------------------------------------------
    # Groups

    def num_children(self):
        """Return the number of groups contained in this group (self.ex_id)."""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_NUM_CHILD_GROUPS)

    def get_parent(self):
        """Return the id of parent of this (self.ex_id) group;
           returns (self.ex_id) if at root."""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_GROUP_PARENT)

    def get_root(self):
        """Return the _id of root group of this (self.ex_id) group;
        returns self.ex_id if at root."""
        return _inquire_value(self.ex_id, cexodus.EX_INQ_GROUP_ROOT)

    def get_group_name(self):
        """Return the name of this group self.ex_id."""
        return _inquire_string(self.ex_id, cexodus.EX_INQ_GROUP_NAME)

    def get_full_group_name(self):
        """Return the full -separated path name of this (self.ex_id) group."""
        return _inquire_string(self.ex_id, cexodus.EX_INQ_FULL_GROUP_NAME)
