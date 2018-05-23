"""
This module contains classes and methods for computing geometry calculations on meshes.

.. currentmodule:: affect.geomtry

.. |uint32_1d| replace:: :obj:`ndarray[uint32_t, ndim=1, mode="c"]`
.. |uint32_2d| replace:: :obj:`ndarray[uint32_t, ndim=2, mode="c"]`
.. |double_2d| replace:: :obj:`ndarray[double, ndim=2, mode="c"]`
"""


cimport cython
cimport numpy
import numpy
from libc.stdint cimport uintptr_t, int32_t, int64_t, uint32_t, uint64_t

from . cimport carithmetic
from . import util


class Error(Exception):
    """Base class for all exceptions raised by this module."""

class MismatchedArrayShape(Error, ValueError):
    """Exception raised if number of array is wrong shape or length to support the operation."""

class IllegalArgumentError(Error, ValueError):
    """Exception raised if argument does not make sense for the current operation."""

class UnalignedArray(Error, ValueError):
    """Exception raised if an array argument is not aligned on a (64) byte boundary"""

class UnsupportedArrayType(Error, TypeError):
    """Exception raised if dtype of numpy array argument is not supported."""

    
def average_element_node_values(
        element_to_node: numpy.ndarray,
        node_values: numpy.ndarray,
        vertices = None) -> numpy.ndarray:
    """
    Compute average value of nodal vector field for each element in a block of elements.
    Alternatively, for higher order elements, average the nodal field only from the *vertices* of the element.

    Some relevant sizes and dimensions are determined from the shapes of the array arguments, these include:

    * num_elements, the number of elements in the block, ``element_to_nodes.shape[0]``
    * nodes_per_element, number of nodes used by each element, ``element_to_nodes.shape[1]``
    * num_components, length of the vector values at nodes, ``node_values.shape[1]``

    If vertices is not None, it must be an integer `N` <= `nodes_per_element`, and the average on each element
    will be computed from the first `N` nodal values of the field components.

    Example:

        Compute the geometric centroid of the elements by passing in the node coordinates. The centroid in
        n-dimensional space is the arithmetic mean ("average") position of all the points in the shape.

        >>> import affect.exodus as exodus
        >>> import affect.arithmetic as arithmetic
        >>>
        >>> with exodus.DatabaseFile('/tmp/myExodusFile.exo') as e:
        >>>
        >>>     coordinates = e.nodal.coordinates()
        >>>
        >>>     for block in e.element_blocks.values():
        >>>
        >>>         num_vertices = connect.CellTopology[block.topology_name].
        >>>
        >>>         connectivity = block.connectivity()  # element to global node connectivity
        >>>
        >>>         # compute the centroid (average of node values) restricted to vertices
        >>>         centroids = arithmetic.average_element_node_values(connectivity, coordinates, vertices=True)

    Args:
        element_to_node(numpy.ndarray): array of dtype=integer, shape(num_elements, num_nodes_per_element)
        node_values(numpy.ndarray): array of dtype=float64, with shape(num_vertex_local, num_components)
        vertices(int): if value is not None (default), compute using only the values at this number of nodes,
            This ignores any higher order edge, face, interior nodes.
            
    Returns:
        element_average - array of |double_2d| with shape(num_elements, num_components).
        
    Raises:
        MismatchedArrayShape: if ``element_to_node[1]`` is not equal to the number of nodes in the cell topology.    
    """
    cdef uintptr_t element_to_node_ptr
    cdef uintptr_t node_values_temp_uintptr
    cdef double * node_values_ptr
    cdef double * element_average_ptr
    cdef numpy.ndarray[numpy.double_t, ndim=2] element_average

    cdef size_t num_elements = element_to_node.shape[0]
    cdef uint32_t num_node_per_element = element_to_node.shape[1]
    cdef uint32_t num_components = node_values.shape[1]
    cdef uint32_t num_vertex_per_element

    element_to_node_ptr = element_to_node.__array_interface__['data'][0]
    if node_values.dtype != numpy.dtype(numpy.float64):
        raise UnsupportedArrayType('node_values.dtype {} is not a double type'.format(node_values.dtype))
    node_values_temp_uintptr = node_values.__array_interface__['data'][0]
    node_values_ptr = <double *> node_values_temp_uintptr

    array_shape = (num_elements, num_components)
    element_average = util.empty_aligned(array_shape, dtype=numpy.double)
    element_average_ptr = <double *> element_average.data
    
    if vertices is not None:
        if not isinstance(vertices, int) or vertices < 1 or vertices > num_node_per_element:
            raise IllegalArgumentError('vertices not an integer in range [1, {}]'.format(element_to_node.shape[1]))
        num_vertex_per_element = vertices
    else:
        num_vertex_per_element = num_node_per_element

    if element_to_node.dtype == numpy.dtype(numpy.int64):
        with nogil:
            carithmetic.average_element_vertex(
                num_vertex_per_element,
                num_elements,
                num_node_per_element,
                <int64_t *> element_to_node_ptr,
                num_components,
                node_values_ptr,
                element_average_ptr)
    elif element_to_node.dtype == numpy.dtype(numpy.int32):
        with nogil:
            carithmetic.average_element_vertex(
                num_vertex_per_element,
                num_elements,
                num_node_per_element,
                <int32_t *> element_to_node_ptr,
                num_components,
                node_values_ptr,
                element_average_ptr)
    elif element_to_node.dtype == numpy.dtype(numpy.uint32):
        with nogil:
            carithmetic.average_element_vertex(
                num_vertex_per_element,
                num_elements,
                num_node_per_element,
                <uint32_t *> element_to_node_ptr,
                num_components,
                node_values_ptr,
                element_average_ptr)
    elif element_to_node.dtype == numpy.dtype(numpy.uint64):
        with nogil:
            carithmetic.average_element_vertex(
                num_vertex_per_element,
                num_elements,
                num_node_per_element,
                <uint64_t *> element_to_node_ptr,
                num_components,
                node_values_ptr,
                element_average_ptr)
    else:
        raise UnsupportedArrayType('aligned array take() for {} element_to_node not implemented'.format(element_to_node.dtype))

    return element_average
