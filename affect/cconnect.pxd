#!/usr/bin/env python

from libc.stdint cimport int64_t
from libcpp cimport bool

cdef extern from "element_topology.hpp":

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

    const int64_t num_spatial_dim[]

    const int64_t num_vertex[]

    const int64_t num_node[]

    const int64_t num_edge[]

    const int64_t num_face[]

    const int64_t num_node_per_edge[]

    const int64_t max_vertex_per_face[]

    const int64_t max_node_per_face[]

    const int64_t * const max_edge_per_face

    const int64_t * const edge_vertex_order[]

    const int64_t * const edge_node_order[]

    const int64_t * const face_vertex_order[]

    const int64_t * const face_node_order[]

    const int64_t * const face_edge_order[]

    const int64_t * const vertex_per_face[]

    const int64_t * const node_per_face[]

    const int64_t * const * const edge_per_face

#-----------------------------------------------------------------------------------------------------------------------

cdef extern from "connect.hpp":

    void to_zero_based(int64_t n, int64_t* array)

    void to_one_based(int64_t n, int64_t* array)

    void connect_boundary_face_to_vertex(
        char* elementName,
        int64_t numElement,
        int64_t numBoundaryQuadFace,
        int64_t numBoundaryTriFace,
        const int64_t * elementToVertex,
        const int64_t * neighbor,
        int64_t * boundaryFaceToVertex)

    void connect_element_to_edge(
        int64_t numElement,
        int64_t numEdgePerElement,
        int64_t numVertexPerElement,
        const int64_t * edgeVertexOrder,
        const int64_t * elementToVertex,
        const int64_t * vertexToElementBegin,
        const int64_t * vertexToElement,
        int64_t * elementToEdge,
        int64_t * numInternalEdge,
        int64_t * numExternalEdge)

    void connect_element_to_element(
        const char* elementName,
        int64_t numElement,
        int64_t maxElementPerVertex,
        const int64_t * elementToVertex,
        const int64_t * vertexToElementBegin,
        const int64_t * vertexToElement,
        int64_t * neighbor,
        int64_t * numBoundaryQuadFace,
        int64_t * numBoundaryTriFace)

    int64_t connect_element_to_face(
        int64_t numElement,
        int64_t numFacePerElement,
        int64_t * neighbor)

    int64_t count_vertex_to_element(
        int64_t numVertex,
        int64_t numVertexPerElement,
        int64_t numElement,
        const int64_t * elementToVertex,
        int64_t * vertexToElementCount)

    int64_t connect_vertex_to_element(
        int64_t numVertex,
        int64_t numVertexPerElement,
        int64_t numElement,
        const int64_t * elementToVertex,
        const int64_t * vertexToElementCount,
        int64_t * vertexToElement)
