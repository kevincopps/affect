#!/usr/bin/env python

from libc.stdint cimport uint32_t
from libc.stdint cimport int64_t
from libc.stdint cimport int8_t

from libcpp cimport bool

cimport cutil

cdef extern from "connect/element_topology.hpp":

    ctypedef enum topology_type:
        HEX8
        HEX20
        HEX27
        PYRAMID5
        PYRAMID13
        QUAD4
        QUAD8
        QUAD9
        QUADSHELL4
        QUADSHELL8
        QUADSHELL9
        TET4
        TET8
        TET10
        TRI3
        TRI4
        TRI6
        TRISHELL3
        TRISHELL6
        WEDGE6
        WEDGE15
        END_TOPOLOGY

    bool is_element_name(const char *topology_name, const char *topology_aliases)

    int get_topology(const char* topology_name)

    const char *names[]

    const char *aliases[]

    const uint32_t num_spatial_dim[]

    const uint32_t num_vertex[]

    const uint32_t num_node[]

    const uint32_t num_edge[]

    const uint32_t num_face[]

    const uint32_t num_node_per_edge[]

    const uint32_t max_vertex_per_face[]

    const uint32_t max_node_per_face[]

    const uint32_t * const max_edge_per_face

    const uint32_t * const edge_vertex_order[]

    const uint32_t * const edge_node_order[]

    const uint32_t * const face_vertex_order[]

    const uint32_t * const face_node_order[]

    const uint32_t * const face_edge_order[]

    const uint32_t * const vertex_per_face[]

    const uint32_t * const node_per_face[]

    const uint32_t * const * const edge_per_face

#-----------------------------------------------------------------------------------------------------------------------

cdef extern from "connect/connect.hpp" nogil:

    void connect_boundary_face_to_vertex(
        char* element_name,
        uint32_t num_element,
        uint32_t num_boundary_quad_face,
        uint32_t num_boundary_tri_face,
        const uint32_t * element_to_vertex,
        const int64_t * neighbor,
        uint32_t * boundary_face_to_vertex)

    void connect_element_to_edge(
        uint32_t num_element,
        uint32_t num_edge_per_element,
        uint32_t num_vertex_per_element,
        const uint32_t * edge_vertex_order,
        const uint32_t * element_to_vertex,
        const uint32_t * vertex_to_element_begin,
        const uint32_t * vertex_to_element,
        int64_t * element_to_edge,
        uint32_t * num_internal_edge,
        uint32_t * num_external_edge)

    int64_t connect_element_to_face(
        uint32_t num_element,
        uint32_t num_face_per_element,
        int64_t * neighbor)

    void connect_element_neighbors(
        int topology,
        const uint32_t num_elements,
        const uint32_t num_vertices,
        const uint32_t * element_to_vertex,
        uint32_t * neighbor_elements,
        int8_t * neighbor_faces)

    uint32_t count_vertex_to_element(
        uint32_t num_vertex,
        uint32_t num_vertex_per_element,
        uint32_t num_element,
        const cutil.uint32_ptr element_to_vertex,
        cutil.uint32_ptr vertex_to_element_count)

    uint32_t connect_vertex_to_element(
        uint32_t num_vertex,
        uint32_t num_vertex_per_element,
        uint32_t num_element,
        const cutil.uint32_ptr element_to_vertex,
        const cutil.uint32_ptr vertex_to_element_count,
        cutil.uint32_ptr vertex_to_element)

    void connect_vertex_to_element_face(
        int topology,
        const uint32_t num_elements,
        const uint32_t num_vertices,
        const cutil.uint32_ptr element_to_vertex,
        const cutil.int8_ptr neighbor_faces,
        cutil.uint32_ptr vertex_facet_element,
        cutil.int8_ptr vertex_facet_face)


cdef extern from "connect/global_to_local.hpp" nogil:

    uint32_t compute_global_to_local[T](
        size_t num_entry,
        T max_global,
        const T * element_to_vertex_global,
        T * max_global_index,
        T * min_global_index,
        uint32_t * element_to_vertex_local,
        uint32_t * global_to_local)

    void fill_local_to_global[T](
        T max_global_index,
        T min_global_index,
        const uint32_t * global_to_local,
        T * local_to_global)


cdef extern from "connect/neighbor.hpp" nogil:

    uint32_t neighbor_hex(
        uint32_t num_element,
        uint32_t max_elements_per_vertex,
        uint32_t * elem_set,
        const uint32_t * element_to_vertex,
        const uint32_t * vertex_to_element,
        const uint32_t * vertex_to_element_begin,
        int64_t * neighbor)

    uint32_t neighbor_tet(
        uint32_t num_element,
        uint32_t max_elements_per_vertex,
        uint32_t * elem_set,
        const uint32_t * element_to_vertex,
        const uint32_t * vertex_to_element,
        const uint32_t * vertex_to_element_begin,
        int64_t * neighbor)

    void neighbor_wedge(
        uint32_t num_element,
        uint32_t max_elements_per_vertex,
        uint32_t * elem_set,
        const uint32_t * element_to_vertex,
        const uint32_t * vertex_to_element,
        const uint32_t * vertex_to_element_begin,
        int64_t * neighbor,
        uint32_t * num_quad_faces,
        uint32_t * num_tri_faces)

    uint32_t neighbor_quad(
        uint32_t num_element,
        uint32_t max_elements_per_vertex,
        uint32_t * elem_set,
        const uint32_t * element_to_vertex,
        const uint32_t * vertex_to_element,
        const uint32_t * vertex_to_element_begin,
        int64_t * neighbor)

    uint32_t neighbor_tri(
        uint32_t num_element,
        uint32_t max_elements_per_vertex,
        uint32_t * elem_set,
        const uint32_t * element_to_vertex,
        const uint32_t * vertex_to_element,
        const uint32_t * vertex_to_element_begin,
        int64_t * neighbor)

cdef extern from "connect/connect_util.hpp" nogil:

    void create_boundary_faces_hex(
        uint32_t num_element,
        const uint32_t* element_to_vertex,
        const int64_t* neighbor,
        uint32_t* boundary_face_to_vertex)

    void create_boundary_faces_tet(
        uint32_t num_element,
        const uint32_t* element_to_vertex,
        const int64_t* neighbor,
        uint32_t* boundary_face_to_vertex)

    void create_boundary_faces_wedge(
        uint32_t num_element,
        const uint32_t* element_to_vertex,
        const int64_t* neighbor,
        uint32_t* boundary_quad_to_vertex,
        uint32_t* boundary_tri_to_vertex)

