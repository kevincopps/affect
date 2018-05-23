"""
This module contains methods for operations on the connectivity of unstructured meshes.

We assume local enumerations for element-to-node, element-to-face, face-to-node, etc.

.. currentmodule:: affect.connect

.. |int64_1d| replace:: :obj:`ndarray[int64_t, ndim=1]`
.. |int64_2d| replace:: :obj:`ndarray[int64_t, ndim=2, mode="c"]`
.. |int32_1d| replace:: :obj:`ndarray[int32_t, ndim=1]`
.. |int32_2d| replace:: :obj:`ndarray[int32_t, ndim=2, mode="c"]`
.. |uint32_1d| replace:: :obj:`ndarray[uint32_t, ndim=1, mode="c"]`
.. |uint32_2d| replace:: :obj:`ndarray[uint32_t, ndim=2, mode="c"]`
.. |(uint32_2d,uint32_1d)| replace:: :obj:`Tuple[ndarray[uint32_t, ndim=2, mode="c"], ndarray[numpy.uint32_t,ndim=1]]`
.. |(uint32_2d,int64_1d)| replace:: :obj:`Tuple[ndarray[uint32_t, ndim=2, mode="c"], ndarray[numpy.int64_t,ndim=1]]`
.. |(int,uint32_1d,uint32_1d)| replace:: :obj:`Tuple[int, ndarray[numpy.uint32_t, ndim=1], ndarray[numpy.uint32_t, ndim=1]]`
.. |(int64_1d,int,int)| replace:: :obj:`Tuple[ndarray[numpy.int64_t, ndim=1], int, int]`
.. |(int32_1d,int,int)| replace:: :obj:`Tuple[ndarray[numpy.int32_t, ndim=1], int, int]`
"""

import time
from enum import IntEnum
from typing import Tuple

cimport cython
import cython
cimport numpy
import numpy
from libc.stdint cimport int32_t, int64_t, uint32_t, int8_t, uintptr_t, UINT32_MAX
from numpy cimport ndarray

from .cimport cconnect
from .cimport cutil
from . import util

# initialize for use of Numpy array C API
numpy.import_array()

cutil.omp_set_dynamic(0)      # Explicitly disable dynamic teams
cutil.omp_set_num_threads(8)  # Use N threads for all consecutive parallel regions


class Error(Exception):
    """Base class for all exceptions raised by this module."""

class InvalidFaceCount(Error, ValueError):
    """Exception raised if number of triangle or quadrilateral faces does not make sense for a cell topology."""

class MaxLengthExceeded(Error, ValueError):
    """Exception raised if a 32 bit unsigned integer limit was exceeded."""

class UnknownCellTopology(Error, ValueError):
    """Exception raised if the string name of a cell topology was unknown."""

class UnalignedArray(Error, ValueError):
    """Exception raised if an array argument is not aligned on a (64) byte boundary"""


class CellTopology(IntEnum):
    """
    Types of element or cell connectivities and aliases.
    
    Examples:
        
        If you want to access members by name (passing the name of an alias will return the original),

        >>> CellTopology['HEX8']
        <CellTopology.HEX8: 0>
        >>> CellTopology['HEX']
        <CellTopology.HEX8: 0>
        >>> CellTopology['TriangleShell3'.upper()]
        <CellTopology.TRISHELL3: 18>
        
        If you want to access members by value,
        
        >>> CellTopology(0)
        <CellTopology.HEX8: 0>
        >>> CellTopology(5)
        <CellTopology.QUAD4: 5>
        
        If you have an instance and need its name or value:
        
        >>> topology = CellTopology.TET4
        >>> topology.name
        'TET4'
        >>> topology.value
        11

        If you have an instance and need a property, for example the number of connected nodes, vertices, faces, etc.
        of the WEDGE15 element:

        >>> topology = CellTopology.WEDGE15
        >>> topology.spatial_dimension
        3
        >>> topology.max_edge_per_face
        4
        >>> topology.max_node_per_face
        8
        >>> topology.max_vertex_per_face
        4
        >>> topology.num_edge
        9
        >>> topology.num_face
        5
        >>> topology.num_node
        15
        >>> topology.num_node_per_edge
        3
        >>> topology.num_vertex
        6

        You need the local number of vertices, nodes, or edges on each face, for example, the WEDGE15 element:

        >>> topology = CellTopology.WEDGE15
        >>> print(topology.edge_per_face)
        [4 4 4 3 3]
        >>> print(topology.node_per_face)
        [8 8 8 6 6]
        >>> print(topology.vertex_per_face)
        [4 4 4 3 3]

        The local order (enumeration) of vertices or nodes on the edges. The shape of the arrays is
        (num_edge,num_vertex) and (num_edge,num_node), respectively. For example, the WEDGE15 element:

        >>> topology = CellTopology.WEDGE15
        >>> print(topology.edge_vertex_order)
        [[0 1]
         [1 2]
         [2 0]
         [3 4]
         [4 5]
         [5 3]
         [0 3]
         [1 4]
         [2 5]]
        >>> print(topology.edge_node_order)
        [[ 0  1  6]
         [ 1  2  7]
         [ 2  0  8]
         [ 3  4 12]
         [ 4  5 13]
         [ 5  3 14]
         [ 0  3  9]
         [ 1  4 10]
         [ 2  5 11]]

        The local order (enumeration) of vertices, nodes, or edge on faces.

        Note that WEDGE and PYRAMID topologies
        require slightly different treatment because the number of entries is different on different faces. So, for
        indices that exceed the local number of entries on that face, the value is UINT32_MAX. The second dimension of
        the shape of the arrays is max_edge_per_face, max_node_per_face, or max_vertex_per_face, respectively.
        For example, the WEDGE15 element:

        >>> topology = CellTopology.WEDGE15
        >>> print(topology.face_edge_order)
        [[         0          7          3          6]
         [         1          8          4          7]
         [         6          5          8          2]
         [         2          1          0 4294967295]
         [         3          4          5 4294967295]]
        >>> print(topology.face_node_order)
        [[         0          1          4          3          6         10          12          9]
         [         1          2          5          4          7         11          13         10]
         [         0          3          5          2          9         14          11          8]
         [         0          2          1          8          7          6  4294967295 4294967295]
         [         3          4          5         12         13         14  4294967295 4294967295]]
        >>> print(topology.face_vertex_order)
        [[         0          1          4          3]
         [         1          2          5          4]
         [         0          3          5          2]
         [         0          2          1 4294967295]
         [         3          4          5 4294967295]]

    The array returned by any property of a CellTopology instance are immutable and will raise an Exception if you
    attempt to modify them.

        >>> topology = CellTopology.HEX8
        >>> a = topology.face_edge_order
        >>> a[0] = 15
        ValueError: assignment destination is read-only

    Pre-constructed sets of the 2D, 3D, and *all* the cell topologies are available:

        >>> CELL_TOPOLOGY_2D
        frozenset({<CellTopology.QUAD4: 5>,
                   <CellTopology.QUAD8: 6>,
                   <CellTopology.QUAD9: 7>,
                   <CellTopology.TRI3: 14>,
                   <CellTopology.TRI4: 15>,
                   <CellTopology.TRI6: 16>})

        >>> CELL_TOPOLOGY_3D
        frozenset({<CellTopology.HEX8: 0>,
                   <CellTopology.HEX20: 1>,
                   <CellTopology.HEX27: 2>,
                   <CellTopology.PYRAMID5: 3>,
                   <CellTopology.PYRAMID13: 4>,
                   <CellTopology.QUAD4: 5>,
                   [...]
                   <CellTopology.WEDGE15: 20>})

        >>> CELL_TOPOLOGY_ALL
        frozenset({<CellTopology.HEX8: 0>,
                   <CellTopology.HEX20: 1>,
                   <CellTopology.HEX27: 2>,
                   <CellTopology.PYRAMID5: 3>,
                   <CellTopology.PYRAMID13: 4>,
                   <CellTopology.QUAD4: 5>,
                   [...]
                   <CellTopology.END_TOPOLOGY: 21>})
    """

    HEX8           = cconnect.HEX8
    HEX            = cconnect.HEX8
    HEX20          = cconnect.HEX20
    HEX27          = cconnect.HEX27
    PYRAMID5       = cconnect.PYRAMID5
    PYRAMID        = cconnect.PYRAMID5
    PYRAMID13      = cconnect.PYRAMID13
    QUAD4          = cconnect.QUAD4
    QUAD           = cconnect.QUAD4
    QUADRILATERAL4 = cconnect.QUAD4
    QUADRILATERAL  = cconnect.QUAD4
    QUAD8          = cconnect.QUAD8
    QUADRILATERAL8 = cconnect.QUAD8
    QUAD9          = cconnect.QUAD9
    QUADRILATERAL9 = cconnect.QUAD9
    QUADSHELL4     = cconnect.QUADSHELL4
    SHELL4         = cconnect.QUADSHELL4
    SHELL          = cconnect.QUADSHELL4
    QUADSHELL      = cconnect.QUADSHELL4
    QUADSHELL8     = cconnect.QUADSHELL8
    SHELL8         = cconnect.QUADSHELL8
    QUADSHELL9     = cconnect.QUADSHELL9
    SHELL9         = cconnect.QUADSHELL9
    TET4           = cconnect.TET4
    TETRA4         = cconnect.TET4
    TETRA          = cconnect.TET4
    TET8           = cconnect.TET8
    TETRA8         = cconnect.TET8
    TET10          = cconnect.TET10
    TETRA10        = cconnect.TET10
    TRI3           = cconnect.TRI3
    TRIANGLE3      = cconnect.TRI3
    TRIANGLE       = cconnect.TRI3
    TRI            = cconnect.TRI3
    TRI4           = cconnect.TRI4
    TRIANGLE4      = cconnect.TRI4
    TRI6           = cconnect.TRI6
    TRIANGLE6      = cconnect.TRI6
    TRISHELL3      = cconnect.TRISHELL3
    TRISHELL       = cconnect.TRISHELL3
    TRIANGLESHELL3 = cconnect.TRISHELL3
    TRIANGLESHELL  = cconnect.TRISHELL3
    TRISHELL6      = cconnect.TRISHELL6
    TRIANGLESHELL6 = cconnect.TRISHELL6
    WEDGE6         = cconnect.WEDGE6
    WEDGE          = cconnect.WEDGE6
    WEDGE15        = cconnect.WEDGE15
    END_TOPOLOGY   = cconnect.END_TOPOLOGY

    @property
    def max_edge_per_face(self) -> int:
        """
        A property that is the maximum number of edges per face of the element.

        Returns:
            maximum edges per face
        """
        return cconnect.max_edge_per_face[self.value]

    @property
    def max_node_per_face(self) -> int:
        """
        A property that is the maximum number of nodes per face of the element.

        Returns:
            maximum nodes per face
        """
        return cconnect.max_node_per_face[self.value]

    @property
    def max_vertex_per_face(self) -> int:
        """
        A property that is the maximum number of vertices per face of the element.
        
        Returns:
            maximum vertices per face
        """
        return cconnect.max_vertex_per_face[self.value]

    @property
    def num_edge(self) -> int:
        """
        A property that is the number of edges per element.
        
        Returns:
            number of edges
        """
        return cconnect.num_edge[self.value]
    
    
    @property
    def num_face(self) -> int:
        """
        A property that is the number of faces per element.
        
        Returns:
            number of faces
        """
        return cconnect.num_face[self.value]

    @property
    def num_node(self) -> int:
        """
        A property that is the number of nodes per element.
        
        Returns:
            number of nodes
        """
        return cconnect.num_node[self.value]

    @property
    def num_node_per_edge(self) -> int:
        """
        A property that is the number of nodes per edge of the element.
        
        Returns:
            nodes per edge
        """
        return cconnect.num_node_per_edge[self.value]

    @property
    def num_vertex(self) -> int:
        """
        A property that is the number of vertices per element.

        Returns:
            number of vertices
        """
        return cconnect.num_vertex[self.value]

    @property
    def spatial_dimension(self) -> int:
        """
        A property that is the number of spatial dimensions spanned by the element, in range [0,3].
        
        Returns:
            number of spatial dimensions
        """
        return cconnect.num_spatial_dim[self.value]

    @property
    def edge_per_face(self) -> numpy.ndarray:
        """
        A property that is an array of the number of edges per face of the element.
        
        Returns:
            array of number of edges per face
        """
        cdef uint32_t num_face = cconnect.num_face[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_face]> cconnect.edge_per_face[self.value])  # get a memory view
        a.flags['WRITEABLE'] = False
        return a

    @property
    def vertex_per_face(self) -> numpy.ndarray:
        """
        A property that is an array of the number of vertices per face of the element.
        
        Returns:
            array of number of vertices per face
        """
        cdef uint32_t num_face = cconnect.num_face[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_face]> cconnect.vertex_per_face[self.value])
        a.flags['WRITEABLE'] = False
        return a

    @property
    def node_per_face(self) -> numpy.ndarray:
        """
        A property that is an array of the number of nodes per face of the element.
        
        Returns:
            array of number of nodes per face
        """
        cdef uint32_t num_face = cconnect.num_face[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_face]> cconnect.node_per_face[self.value])
        a.flags['WRITEABLE'] = False
        return a

    @property
    def edge_vertex_order(self) -> numpy.ndarray:
        """
        A property that is an array of local vertex enumeration for each edge.

        Returns:
            array of local edge, local vertex numbers, with shape(num_edge, 2)
        """
        cdef uint32_t num_edge = cconnect.num_edge[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_edge*2]> cconnect.edge_vertex_order[self.value])
        r = numpy.reshape(a, (num_edge, 2))
        r.flags['WRITEABLE'] = False
        return r

    @property
    def edge_node_order(self) -> numpy.ndarray:
        """
        A property that is an array of local node enumeration for each edge.

        Returns:
            array of local edge, local node numbers, with shape(num_edge, num_node_per_edge)
        """
        cdef uint32_t num_edge = cconnect.num_edge[self.value]
        cdef uint32_t num_node_per_edge = cconnect.num_node_per_edge[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_edge*num_node_per_edge]> cconnect.edge_node_order[self.value])
        r = numpy.reshape(a, (num_edge, num_node_per_edge))
        r.flags['WRITEABLE'] = False
        return r

    @property
    def face_vertex_order(self) -> numpy.ndarray:
        """
        A property that is an array of local vertex enumeration for each face.
        
        Entries corresponding to indices greater than the number of vertices on a single face
        have values UINT32_MAX. This case occurs for WEDGE and PYRAMID topologies, where some faces have 3 vertices and
        some have 4 vertices.

        Returns:
            array of local face, local vertex numbers, with shape(num_face, max_vertex_per_face)
        """
        cdef uint32_t num_face = cconnect.num_face[self.value]
        cdef uint32_t max_vertex_per_face = cconnect.max_vertex_per_face[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_face*max_vertex_per_face]> cconnect.face_vertex_order[self.value])
        r = numpy.reshape(a, (num_face, max_vertex_per_face))
        r.flags['WRITEABLE'] = False
        return r

    @property
    def face_node_order(self) -> numpy.ndarray:
        """
        A property that is an array of local node enumeration for each face.
        
        Entries corresponding to indices greater than the number of nodes on a single face
        have values UINT32_MAX. This case occurs for WEDGE and PYRAMID topologies, where faces have a different number
        of nodes.

        Returns:
            array of local face, local node numbers, with shape(num_face, max_node_per_face)
        """
        cdef uint32_t num_face = cconnect.num_face[self.value]
        cdef uint32_t max_node_per_face = cconnect.max_node_per_face[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_face*max_node_per_face]> cconnect.face_node_order[self.value])
        r = numpy.reshape(a, (num_face, max_node_per_face))
        r.flags['WRITEABLE'] = False
        return r

    @property
    def face_edge_order(self) -> numpy.ndarray:
        """
        A property that is an array of local edge enumeration for each face.
        
        Entries corresponding to indices greater than the number of edges on a single face
        have values UINT32_MAX. This case occurs for WEDGE and PYRAMID topologies, where faces have a different number
        of edges.

        Returns:
            array of local face, local edge numbers, with shape(num_face, max_edge_per_face)
        """
        cdef uint32_t num_face = cconnect.num_face[self.value]
        cdef uint32_t max_edge_per_face = cconnect.max_edge_per_face[self.value]
        a = numpy.asarray(<numpy.uint32_t[:num_face*max_edge_per_face]> cconnect.face_edge_order[self.value])
        r = numpy.reshape(a, (num_face, max_edge_per_face))
        r.flags['WRITEABLE'] = False
        return r



CELL_TOPOLOGY_ALL = frozenset(CellTopology.__members__.values())
"""The set of all CellTopology values."""

CELL_TOPOLOGY_2D = frozenset((CellTopology.QUAD4,
                           CellTopology.QUAD8,
                           CellTopology.QUAD9,
                           CellTopology.TRI3,
                           CellTopology.TRI4,
                           CellTopology.TRI6
                           ))
"""
Set of 2D topologies.
"""

CELL_TOPOLOGY_3D = frozenset((CellTopology.HEX8,
                           CellTopology.HEX20,
                           CellTopology.HEX27,
                           CellTopology.PYRAMID5,
                           CellTopology.PYRAMID13,
                           CellTopology.QUAD4,
                           CellTopology.QUAD8,
                           CellTopology.QUAD9,
                           CellTopology.QUADSHELL4,
                           CellTopology.QUADSHELL8,
                           CellTopology.QUADSHELL9,
                           CellTopology.TET4,
                           CellTopology.TET8,
                           CellTopology.TET10,
                           CellTopology.TRI3,
                           CellTopology.TRI4,
                           CellTopology.TRI6,
                           CellTopology.TRISHELL3,
                           CellTopology.TRISHELL6,
                           CellTopology.WEDGE6,
                           CellTopology.WEDGE15
                           ))
"""
Set of 3D topologies.
"""







@cython.boundscheck(False)
def vertex_to_element(upper_bound_vertex: int,
                      ndarray[numpy.uint32_t, ndim=2, mode="c"] element_to_vertex not None
                      ) -> Tuple[int, ndarray, ndarray]:
    """
    For a uniform block of elements construct the vertex to element connectivity.

    The elements connected to the ith vertex are located at ``vertices_to_elements[begin:end]`` where
    ``begin = vertex_to_element_begin[i]`` and ``end = vertex_to_element_begin[i+1]``

    Number of elements in the given block is element_to_vertex.shape[0].
    Number of vertices per element is element_to_vertex.shape[1].

    Uses unsigned 32 bit integers, so upper_bound_vertex must be equal or less than 4,294,967,295 which equals 2^32 − 1.

    Args:
        upper_bound_vertex (int): upper limit on largest vertex entry that may possibly appear in `element_to_vertex`,
            i.e., at least one more than the max value appearing in `element_to_vertex`.
            (supplying a larger `upper_bound_vertex` is fine but uses more memory)
        element_to_vertex(|uint32_2d|): array of element to vertex connectivity with
            shape(num_elements, num_vertex_per_element)
    Returns:
        |(int,uint32_1d,uint32_1d)|:
        (max_element_per_vertex, vertex_to_element_begin, vertices_to_elements) -

        * maximum number of elements occurring at a single vertex,
        * an array with shape(upper_bound_vertex+1),
        * an array shape(indices[-1])
    """
    cdef uint32_t _upper_bound_vertex
    cdef uint32_t size
    cdef uint32_t num_vertex_per_element
    cdef uint32_t num_elements
    cdef cutil.uint32_ptr element_to_vertex_ptr
    cdef cutil.uint32_ptr vertices_to_elements_ptr
    cdef cutil.uint32_ptr indices_ptr
    cdef ndarray[numpy.uint32_t, ndim=1] indices
    cdef ndarray[numpy.uint32_t, ndim=1] vertices_to_elements

    _upper_bound_vertex = upper_bound_vertex  # coercion to C type before we enter nogil

    if not util.is_byte_aligned(element_to_vertex):
        raise UnalignedArray("Argument array element_to_vertex is not properly aligned on byte boundary")

    element_to_vertex_ptr = <cutil.uint32_ptr> element_to_vertex.data
    indices = util.empty_aligned(_upper_bound_vertex + 1, dtype=numpy.uint32)
    indices_ptr = <cutil.uint32_ptr> indices.data
    if element_to_vertex.shape[0] > UINT32_MAX:
        raise MaxLengthExceeded('Number of local elements greater than 2^32 - 1 (UINT32_MAX).')
    num_elements = <uint32_t>element_to_vertex.shape[0]
    num_vertex_per_element = <uint32_t>element_to_vertex.shape[1]

    with nogil:
        max_element_per_vertex = cconnect.count_vertex_to_element(_upper_bound_vertex,
                                                                  num_vertex_per_element,
                                                                  num_elements,
                                                                  element_to_vertex_ptr,
                                                                  indices_ptr)
    #util.print_array_info('indices', indices)
    size = indices[-1]
    vertices_to_elements = util.empty_aligned(size, dtype=numpy.uint32)
    vertices_to_elements_ptr = <cutil.uint32_ptr> vertices_to_elements.data

    with nogil:
        cconnect.connect_vertex_to_element(_upper_bound_vertex,
                                           num_vertex_per_element,
                                           num_elements,
                                           element_to_vertex_ptr,
                                           indices_ptr,
                                           vertices_to_elements_ptr)

    return max_element_per_vertex, indices, vertices_to_elements


@cython.boundscheck(False)
@cython.wraparound(False)
def element_to_element(topology: CellTopology,
                       max_elements_per_vertex: int,
                       ndarray[numpy.uint32_t, ndim=2, mode="c"] element_to_vertex not None,
                       ndarray[numpy.uint32_t, ndim=1] vertex_to_element_begin not None,
                       ndarray[numpy.uint32_t, ndim=1] vertex_to_element not None) -> Tuple[ndarray, int, int]:
    """
    For a single block of elements with 3D topology, return an array of 
    the element-to-element connectivity.
    
    Uses temporary space of ``sizeof(uint32_t) * 2 * max_elements_per_vertex``
    
    Args:
        topology(CellTopology): type of cells in this block
        max_elements_per_vertex(int): maximum count of number of cells connected to a single vertex
        element_to_vertex(|uint32_2d|): array of shape(num_elements, num_vertex_per_element)
        vertex_to_element_begin(|uint32_1d|): array with length (num_vertices+1), where the ith 
            entry is the starting index into the vertex_to_element array for the ith element
        vertex_to_element(|uint32_1d|): array of vertex-to-element connectivity
  
    Returns:
        |(int64_1d,int,int)|: (neighbor, num_boundary_quad_face, num_boundary_tri_face) - 
        
        * the array with an entry for each element side, containing either the ID of the element neighbor or 
          -1 if no neighbor exists, shape(num_face_per_element * num_element),
        * the number of quadrilateral element faces without a neighbor
        * the number of triangular element faces without a neighbor
    """
    cdef uint32_t num_boundary_quad_face = 0
    cdef uint32_t num_boundary_tri_face = 0
    cdef uint32_t num_elements
    cdef uint32_t num_face_per_element
    cdef uint32_t *element_to_vertex_ptr
    cdef uint32_t *vertex_to_element_begin_ptr
    cdef uint32_t *vertex_to_element_ptr
    cdef uint32_t *element_set_ptr
    cdef int64_t *neighbor_ptr
    cdef ndarray[numpy.uint32_t, ndim=1] element_set
    cdef ndarray[numpy.int64_t, ndim=1] neighbor
    cdef size_t neighbor_length

    if element_to_vertex.shape[0] > UINT32_MAX:
        raise MaxLengthExceeded('Number of local elements greater than 2^32 - 1 (UINT32_MAX).')
    num_elements = <uint32_t>element_to_vertex.shape[0]
    if num_elements < 1 or  max_elements_per_vertex < 1:
        return None, 0, 0

    element_to_vertex_ptr = <uint32_t *> element_to_vertex.data
    vertex_to_element_begin_ptr = <uint32_t *> vertex_to_element_begin.data
    vertex_to_element_ptr = <uint32_t *> vertex_to_element.data

    num_face_per_element = cconnect.num_face[topology]

    # element_set = numpy.empty(2 * max_elements_per_vertex, dtype=numpy.uint32)
    element_set = util.empty_aligned(2 * max_elements_per_vertex, dtype=numpy.uint32)
    element_set_ptr = <uint32_t *> element_set.data

    neighbor_length = num_elements * num_face_per_element
    neighbor = util.empty_aligned(neighbor_length, dtype=numpy.int64)
    neighbor_ptr = <int64_t *> neighbor.data
    cutil.fill[int64_t](neighbor_ptr, neighbor_length, <int64_t> -2)

    #start_time = time.process_time()
    if topology == CellTopology.HEX8:
        num_boundary_quad_face = cconnect.neighbor_hex(num_elements,
                                                       max_elements_per_vertex,
                                                       element_set_ptr,
                                                       element_to_vertex_ptr,
                                                       vertex_to_element_begin_ptr,
                                                       vertex_to_element_ptr,
                                                       neighbor_ptr)
        stop_time = time.process_time()
    elif topology == CellTopology.TET4:
        num_boundary_tri_face = cconnect.neighbor_tet(num_elements,
                                                      max_elements_per_vertex,
                                                      element_set_ptr,
                                                      element_to_vertex_ptr,
                                                      vertex_to_element_begin_ptr,
                                                      vertex_to_element_ptr,
                                                      neighbor_ptr)
    elif topology == CellTopology.WEDGE6:
        cconnect.neighbor_wedge(num_elements,
                                max_elements_per_vertex,
                                element_set_ptr,
                                element_to_vertex_ptr,
                                vertex_to_element_begin_ptr,
                                vertex_to_element_ptr,
                                neighbor_ptr,
                                &num_boundary_quad_face,
                                &num_boundary_tri_face)
    else:
        raise UnknownCellTopology('Unsupported cell topology, only HEX8, TET4, and WEDGE5 are supported.')
    #stop_time = time.process_time()
    #print('    element_to_element process time: {:f}s'.format(stop_time - start_time))
    return neighbor, num_boundary_quad_face, num_boundary_tri_face


@cython.boundscheck(False)
@cython.wraparound(False)
def boundary_face_to_vertex(topology: CellTopology,
                            num_boundary_quad_face: int,
                            num_boundary_tri_face: int,
                            ndarray[numpy.uint32_t, ndim=2, mode="c"] element_to_vertex not None,
                            ndarray[numpy.int64_t, ndim=1] neighbor not None
                            ) -> ndarray:
    """
    For a uniform block of elements construct the boundary face to vertex connectivity.
    
    Args:
        topology(CellTopology): type of cells in this block
        num_boundary_quad_face(int): the number of quad faces on the boundary
        num_boundary_tri_face(int): the number of triangle faces on the boundary 
        element_to_vertex(|uint32_2d|): the element to vertex connectivity
        neighbor(|int64_1d|): the element to element connectivity, where -1 indicates an element 
            face on the boundary

    Returns:
        |uint32_1d|: boundary_face_to_vertex - array of shape ``(4*num_boundary_quad_face + 3*num_boundary_tri_face,)``
    """
    cdef uint32_t num_elements, num_boundary_vertices
    cdef uint32_t *element_to_vertex_ptr
    cdef uint32_t *boundary_face_to_vertex_ptr
    cdef int64_t  *neighbor_ptr
    cdef ndarray[numpy.uint32_t, ndim=1] boundary_face_to_vertex

    element_to_vertex_ptr = <uint32_t *> element_to_vertex.data
    neighbor_ptr = <int64_t *> neighbor.data
    num_boundary_vertices = 4 * num_boundary_quad_face + 3 * num_boundary_tri_face

    if element_to_vertex.shape[0] > UINT32_MAX:
        raise MaxLengthExceeded('Number of local elements greater than 2^32 - 1 (UINT32_MAX).')
    num_elements = <uint32_t>element_to_vertex.shape[0]

    #boundary_face_to_vertex = numpy.empty(num_boundary_vertices, dtype=numpy.uint32)
    boundary_face_to_vertex = util.empty_aligned(num_boundary_vertices, dtype=numpy.uint32)

    boundary_face_to_vertex_ptr = <uint32_t *> boundary_face_to_vertex.data

    if topology == CellTopology.HEX8:
        if num_boundary_tri_face != 0:
            raise InvalidFaceCount('num_boundary_tri_face was not zero for HEX8')
        if num_boundary_quad_face < 1:
            raise InvalidFaceCount('num_boundary_quad_face was zero for HEX8')
        cconnect.create_boundary_faces_hex(
            num_elements,
            element_to_vertex_ptr,
            neighbor_ptr,
            boundary_face_to_vertex_ptr)
    elif topology == CellTopology.TET4:
        if num_boundary_quad_face != 0:
            raise InvalidFaceCount('num_boundary_quad_face was not zero for TET4')
        if num_boundary_tri_face < 1:
            raise InvalidFaceCount('num_boundary_tri_face was zero for TET4')
        cconnect.create_boundary_faces_tet(
            num_elements,
            element_to_vertex_ptr,
            neighbor_ptr,
            boundary_face_to_vertex_ptr)

    elif topology == CellTopology.WEDGE6:
        if num_boundary_quad_face == 0 or num_boundary_tri_face == 0:
            raise InvalidFaceCount('Number of tri faces and quad faces was zero for WEDGE6')
        cconnect.create_boundary_faces_wedge(
            num_elements,
            element_to_vertex_ptr,
            neighbor_ptr,
            boundary_face_to_vertex_ptr,
            &boundary_face_to_vertex_ptr[num_boundary_quad_face * 4])
    else:
        raise UnknownCellTopology('Unsupported cell topology, only HEX8, TET4, and WEDGE5 are supported.')
    return boundary_face_to_vertex


# cpdef array_str_line(input_array):
#
#     cdef uint64_t i, count
#
#     if input_array is not None:
#         if type(input_array) is numpy.ndarray:
#             cdef uint64_t count = input_array.size
#             if input_array.dtype == numpy.float64:
#                 cdef double *a_ptr = <double *> a.data
#                 ## write all the element to a strstream
#                 for i in range(count):
#                     s << a_ptr[i]
#
#
#             elif input_array.dtype == numpy.uint32:
#
#             else:
#                 raise TypeError('input_array.dtype must be numpy.double or numpy.float64.')
#
#
#         #obj = numpy.asarray(input_array).view(cls)  # We first cast to be our class type
#         #obj.info = info                          # add the new attribute to the created instance
#         #return obj


@cython.boundscheck(False)
@cython.wraparound(False)
def element_neighbors(topology: CellTopology,
                      num_vertices: int,
                      ndarray[numpy.uint32_t, ndim=2, mode="c"]
                      element_to_vertex not None) -> Tuple[ndarray, ndarray]:
    """
    Construct array of neighbor element IDs and array of local neighbor face IDs for each element face using the sibling
    half-facet algorithm.

    Args:
        num_vertices:
        topology(CellTopology): type of cells in this block
        element_to_vertex(|uint32_2d|): array of shape(num_elements, num_vertex_per_element)

    Returns:
        neighbor_faces(|uint32_2d|): array of neighbor element ID for each element face
        neighbor_faces(|int8_2d|): array of neighbor local face ID for each element face
    """
    cdef uint32_t * element_to_vertex_ptr
    cdef uint32_t * neighbor_elements_ptr
    cdef int8_t * neighbor_faces_ptr
    #cdef cconnect.facet_struct * neighbor_facets_ptr
    cdef uint32_t num_elements = <uint32_t> element_to_vertex.shape[0]
    cdef uint32_t num_face_per_element = cconnect.num_face[topology.value]
    cdef numpy.ndarray[numpy.uint32_t, ndim=2] neighbor_elements
    cdef numpy.ndarray[numpy.int8_t, ndim=2] neighbor_faces
    cdef uint32_t _num_vertices = num_vertices
    cdef int topology_value = topology.value

    array_shape = (num_elements, num_face_per_element)
    element_to_vertex_ptr = <uint32_t *> element_to_vertex.data
    neighbor_elements = util.zeros_aligned(array_shape, dtype=numpy.uint32)
    neighbor_elements_ptr = <uint32_t *> neighbor_elements.data
    neighbor_faces = util.empty_aligned(array_shape, dtype=numpy.int8)
    neighbor_faces_ptr = <int8_t *> neighbor_faces.data

    with nogil:
        cconnect.connect_element_neighbors(topology_value,
                                           num_elements,
                                           _num_vertices,
                                           element_to_vertex_ptr,
                                           neighbor_elements_ptr,
                                           neighbor_faces_ptr)

    return neighbor_elements, neighbor_faces
    # return neighbor_facets  # when using alternate struct storage


@cython.boundscheck(False)
@cython.wraparound(False)
def connect_vertex_to_element_face(
        topology: CellTopology,
        num_vertices: int,
        ndarray[numpy.uint32_t, ndim=2, mode="c"] element_to_vertex not None,
        ndarray[numpy.int8_t, ndim=2, mode="c"] neighbor_faces not None) -> Tuple[ndarray, ndarray]:
    """
    Construct arrays for each vertex to one incident element and the local element face. Boundary faces are prioritized
    so that if the vertex is on the boundary, the incident element face is one that is on the boundary.

    Note: Storage of the result of this function (and the given half-facets (neighbor, and neighbor face) allows any
    adjacency query on the mesh in constant time.

    Args:
        topology(CellTopology): type of cells in this block
        num_vertices: number of vertices in this block
        element_to_vertex(|uint32_2d|): array of shape(num_elements, num_vertex_per_element)
        neighbor_faces(|int8_1d|): array of neighbor local face ID for each element face

    Returns:
        vertex_facet_element(|uint32_1d|): array of one incident element for each vertex
        vertex_facet_face(|int8_1d|): array of local incident element face for each vertex
    """
    cdef cutil.uint32_ptr element_to_vertex_ptr
    cdef int8_t * neighbor_faces_ptr
    cdef uint32_t * vertex_facet_element_ptr
    cdef int8_t * vertex_facet_face_ptr

    cdef numpy.ndarray[numpy.uint32_t, ndim=1] vertex_facet_element
    cdef numpy.ndarray[numpy.int8_t, ndim=1] vertex_facet_face
    cdef uint32_t num_elements = <uint32_t> element_to_vertex.shape[0]
    cdef uint32_t _num_vertices = num_vertices
    cdef int topology_value = topology.value

    element_to_vertex_ptr = <cutil.uint32_ptr> element_to_vertex.data
    neighbor_faces_ptr = <int8_t *> neighbor_faces.data

    vertex_facet_element = util.empty_aligned(_num_vertices, dtype=numpy.uint32)
    vertex_facet_element_ptr = <uint32_t *> vertex_facet_element.data
    vertex_facet_face = util.empty_aligned(_num_vertices, dtype=numpy.int8)
    vertex_facet_face_ptr = <int8_t *> vertex_facet_face.data

    with nogil:
        cconnect.connect_vertex_to_element_face(
            topology_value,
            num_elements,
            _num_vertices,
            element_to_vertex_ptr,
            neighbor_faces_ptr,
            vertex_facet_element_ptr,
            vertex_facet_face_ptr)

    return vertex_facet_element, vertex_facet_face


@cython.boundscheck(False)
@cython.wraparound(False)
def convert_to_local_connectivity_buffer(
        element_to_vertex_global: numpy.ndarray,
        max_global: int,
        ndarray[numpy.uint32_t, ndim=1, mode="c"] global_to_local_buffer) -> Tuple[ndarray, ndarray]:
    """
    Convert array of element to global vertex to array of element to local vertex.

    Creates a map of local to global vertex indices and a new array of element to local vertex connectivity.

    In the process, this also converts the global array of 64 integers to the local array which uses
    unsigned 32 bit integers, thus the maximum number of unique local vertex indices must be less than or
    equal to 4,294,967,295, which equals 2^32 − 1.

    Args:
        max_global(int): highest value of global ID that may be present int element_to_vertex_global
        element_to_vertex_global(|int64_2d| or |int32_2d|): global element_to_vertex array
        global_to_local_buffer((|uint_1d|): buffer of size (max_global) used internally for sorting

    Returns:
        |(uint32_2d,int64_1d)|: (element_to_vertex_local, local_to_global_vertex) -

        * the local element to vertex connectivity and,
        * the array of the global vertex for each local vertex.
    """
    cdef int32_t * element_to_vertex_global_ptr32
    cdef int32_t * local_to_global_ptr32
    cdef ndarray[numpy.int32_t, ndim=1] local_to_global32
    cdef int32_t _max_global32
    cdef int32_t max_global_index32
    cdef int32_t min_global_index32

    cdef int64_t * element_to_vertex_global_ptr64
    cdef int64_t * local_to_global_ptr64
    cdef ndarray[numpy.int64_t, ndim=1] local_to_global64
    cdef int64_t _max_global64
    cdef int64_t max_global_index64
    cdef int64_t min_global_index64

    cdef uint32_t * global_to_local_ptr
    cdef uint32_t * element_to_vertex_local_ptr
    cdef ndarray[numpy.uint32_t, ndim=2] element_to_vertex_local
    cdef size_t num_entry
    cdef uint32_t num_unique
    cdef uintptr_t ptr

    if max_global < 1:
        raise ValueError('Maximum number of global indices was less than one.')
    if not isinstance(element_to_vertex_global, numpy.ndarray):
        raise TypeError('element_to_vertex_global is not a numpy.ndarray')
    num_entry = element_to_vertex_global.size
    if num_entry < 1:
        raise ValueError('Connectivity is empty.')
    array_shape = (element_to_vertex_global.shape[0], element_to_vertex_global.shape[1])

    element_to_vertex_local = util.empty_aligned(array_shape, dtype=numpy.uint32)
    element_to_vertex_local_ptr = <uint32_t *> element_to_vertex_local.data

    if global_to_local_buffer.size < max_global:
        raise ValueError('global_to_local_buffer.size is less than max_global')
    global_to_local_ptr = <uint32_t *> global_to_local_buffer.data

    if element_to_vertex_global.dtype.type is numpy.int32:

        _max_global32 = max_global
        ptr = <uintptr_t> element_to_vertex_global.__array_interface__['data'][0]
        element_to_vertex_global_ptr32 = <int32_t *> ptr
        with nogil:
            num_unique = cconnect.compute_global_to_local(num_entry,
                                                          _max_global32,
                                                          element_to_vertex_global_ptr32,
                                                          &max_global_index32,
                                                          &min_global_index32,
                                                          element_to_vertex_local_ptr,
                                                          global_to_local_ptr)
        local_to_global32 = util.empty_aligned(num_unique, dtype=numpy.int32)
        ptr = <uintptr_t> local_to_global32.__array_interface__['data'][0]
        local_to_global_ptr32 = <int32_t *> ptr
        with nogil:
            cconnect.fill_local_to_global(max_global_index32,
                                          min_global_index32,
                                          global_to_local_ptr,
                                          local_to_global_ptr32)
        return element_to_vertex_local, local_to_global32

    elif element_to_vertex_global.dtype.type is numpy.int64:

        _max_global64 = max_global
        ptr = <uintptr_t> element_to_vertex_global.__array_interface__['data'][0]
        element_to_vertex_global_ptr64 = <int64_t *> ptr
        with nogil:
            num_unique = cconnect.compute_global_to_local(num_entry,
                                                          _max_global64,
                                                          element_to_vertex_global_ptr64,
                                                          &max_global_index64,
                                                          &min_global_index64,
                                                          element_to_vertex_local_ptr,
                                                          global_to_local_ptr)
        #print('    max_global_index64 = {}, min_global_index64 = {}'.format(max_global_index64, min_global_index64))
        if num_unique == 0:
            raise MaxLengthExceeded('Number of local indices exceed unsigned 32 bit limit: 2^32 − 1')
        local_to_global64 = util.empty_aligned(num_unique, dtype=numpy.int64)
        ptr = <uintptr_t> local_to_global64.__array_interface__['data'][0]
        local_to_global_ptr64 = <int64_t *> ptr
        with nogil:
            cconnect.fill_local_to_global(max_global_index64,
                                          min_global_index64,
                                          global_to_local_ptr,
                                          local_to_global_ptr64)
        return element_to_vertex_local, local_to_global64
    else:
        raise TypeError('element_to_vertex_global.dtype is not numpy.int32 and is not numpy.int64')

    # NOTE: A slower method that conserves memory using using numpy.unique:
    #   Collect the unique set of global vertices that are actually used in this connectivity and a map.
    #   We will get the sorted unique entries of the array, and
    #   the indices to reconstruct the (flattened) original array from the unique array.
    #   The indices array will have the same dtype as the input to unique
    # local_to_global_vertex64, indices = numpy.unique(element_to_vertex_global, return_inverse=True)  #, max_entry)
    # if len(indices) > UINT32_MAX:
    #     raise MaxLengthExceeded('Number of unique vertices in local connectivity exceeded UINT32_MAX')
    # # convert to uint32 and reshape 1D back to 2D array
    # # local_to_global_vertex32 = local_to_global_vertex64.astype(numpy.uint32)
    # array_shape = (element_to_vertex_global.shape[0], element_to_vertex_global.shape[1])
    # element_to_vertex_local32 = indices.reshape(array_shape).astype(numpy.uint32)
    # return element_to_vertex_local32, local_to_global_vertex64


@cython.boundscheck(False)
@cython.wraparound(False)
def convert_to_local_connectivity(
        element_to_vertex_global: numpy.ndarray,
        num_global_vertices: int) -> Tuple[ndarray, ndarray]:
    """
    Convert array of element to global vertex to array of element to local vertex.

    Creates a map of local to global vertex indices and a new array of element to local vertex connectivity.

    In the process, this also converts the global array of 64 integers to the local array which uses
    unsigned 32 bit integers, thus the maximum number of unique local vertex indices must be less than or
    equal to 4,294,967,295, which equals 2^32 − 1.

    Args:
        element_to_vertex_global(|int64_2d| or |int32_2d|): global element_to_vertex array
        num_global_vertices(int): upper bound of global vertex ID

    Returns:
        |(uint32_2d,int64_1d)|: (element_to_vertex_local, local_to_global_vertex) -

        * the local element to vertex connectivity and,
        * the array of the global vertex for each local vertex.
    """

    cdef ndarray[numpy.uint32_t, ndim=1, mode="c"] global_to_local_buffer
    global_to_local_buffer = util.empty_aligned(num_global_vertices, dtype=numpy.uint32)
    return convert_to_local_connectivity_buffer(element_to_vertex_global,
                                                num_global_vertices,
                                                global_to_local_buffer)
