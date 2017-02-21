#!/usr/bin/env python

# -*- coding: utf-8 -*-

# cython: embedsignature=True

"""
This module contains methods for operations on the connectivity of unstructured meshes.

We assume local enumerations for element-to-node, element-to-face, face-to-node, etc. We use enumerations commonly
used in the ExodusII mesh database format.

"""

cimport cython
import numpy as np
cimport numpy as np
from . cimport cconnect
from libc.stdlib cimport malloc, free
from libc.stdint cimport int64_t, int32_t

#-----------------------------------------------------------------------------------------------------------------------

HEX8 =          cconnect.HEX8
HEX20 =         cconnect.HEX20
HEX27 =         cconnect.HEX27
PYRAMID5 =      cconnect.PYRAMID5
PYRAMID13 =     cconnect.PYRAMID13
QUAD4 =         cconnect.QUAD4
QUAD8 =         cconnect.QUAD8
QUAD9 =         cconnect.QUAD9
QUADSHELL4 =    cconnect.QUADSHELL4
QUADSHELL8 =    cconnect.QUADSHELL8
QUADSHELL9 =    cconnect.QUADSHELL9
TET4 =          cconnect.TET4
TET8 =          cconnect.TET8
TET10 =         cconnect.TET10
TRI3 =          cconnect.TRI3
TRI4 =          cconnect.TRI4
TRI6 =          cconnect.TRI6
TRISHELL3 =     cconnect.TRISHELL3
TRISHELL6 =     cconnect.TRISHELL6
WEDGE6 =        cconnect.WEDGE6
WEDGE15 =       cconnect.WEDGE15
END_TOPOLOGY =  cconnect.END_TOPOLOGY

TOPOLOGIES_2D = [QUAD4, QUAD8, QUAD9, TRI3, TRI4, TRI6]
"""
List of 2D topologies.
"""

TOPOLOGIES_3D = [HEX8, HEX20, HEX27, PYRAMID5, PYRAMID13, QUAD4, QUAD8, QUAD9, QUADSHELL4, QUADSHELL8,
                 QUADSHELL9, TET4, TET8, TET10, TRI3, TRI4, TRI6, TRISHELL3, TRISHELL6, WEDGE6, WEDGE15]
"""
List of 3D topologies.
"""

def raise_topology_not_found(topology):
    if topology < 0 or topology >= END_TOPOLOGY:
        raise ValueError('Element topology is unknown ({})'.format(topology))

def get_topology(topology_name):
    """
    Given a string name, return an integer representing the topology.

    :param topology_name: name of a topology, for example "TET4", or "hex20"
    :type topology_name: string
    :return: integer -- enumerated value for the topology, in the range [0, END_TOPOLOGY-1]
    :raises ValueError: if the topology_name does not match any known topology
    """
    t = cconnect.get_topology(topology_name)
    if t < 0:
        raise ValueError('No topology matches ' + topology_name + '.')
    return t

def is_topology(topology_name, topology):
    """
    Return if the element described by the string name is one of the given topology.

    :param topology_name: name of a topology, for example "TET4", or "hex20"
    :type topology_name: string
    :param topology: one of HEX8, HEX20, HEX27, TET4, TET8, TRI3, QUAD4, PYRAMID5, etc.
    :type topology: cconnect.topology_type
    :return: True or False
    """
    if cconnect.is_element_name(topology_name, cconnect.aliases[topology]):
        return True
    else:
        return False

#-----------------------------------------------------------------------------------------------------------------------

def vertex_to_element(
        num_global_vertices,
        topology,
        num_elements,
        np.ndarray[np.int64_t, ndim=2, mode="c"] element_to_vertex not None):
    """
    For a uniform block of elements, construct an array which
    contains the sums of counts of the number of elements connected to each vertex.
    The purpose of this array is to be used as a lookup index into the
    array holding the set of elements connected to each vertex.

    Note: the length of the input vertexToElementCount array must be one longer
    than the number of vertices in the block.

    :param num_global_vertices: number of global vertices
    :param num_elements: number of elements in this block
    :param topology: one of HEX8, TET4, WEDGE6, etc.
    :param element_to_vertex: array of element-to-vertex connectivity, length num_elements * num_vertex_per_element
    :type element_to_vertex: numpy.ndarray of type int64_t
    :return: vertex_to_element_begin[:], vertex_to_element[:], of type numpy.ndarray
    """
    raise_topology_not_found(topology)
    num_vertex_per_element = cconnect.num_vertex[topology]

    cdef np.ndarray[np.int64_t, ndim=1] vertex_to_element_begin = np.empty(num_global_vertices + 1, dtype=np.int64)
    cdef int64_t *vertex_to_element_begin_ptr = <int64_t *> vertex_to_element_begin.data

    cconnect.count_vertex_to_element(num_global_vertices,
                                     num_vertex_per_element,
                                     num_elements,
                                     <const int64_t *> &element_to_vertex[0,0],
                                     vertex_to_element_begin_ptr)

    return vertex_to_element_begin

    # connect_vertex_to_element(
    #   int64_t numVertex,
    #   int64_t numVertexPerElement,
    #   int64_t numElement,
    #   const int64_t * elementToVertex,
    #   int64_t * vertexToElementCount,
    #   int64_t * vertexToElement