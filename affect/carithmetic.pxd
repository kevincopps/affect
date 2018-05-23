#!/usr/bin/env python

from libc.stdint cimport uint32_t

cimport cutil

cdef extern from "arithmetic/average_element_vertex.hpp" nogil:

    void average_element_vertex[T](
        uint32_t num_vertex_per_element,
        size_t num_elements,
        uint32_t num_node_per_element,
        const T * element_to_node,
        uint32_t num_components,
        const cutil.double_ptr node_values,
        cutil.double_ptr element_value)

cdef extern from "arithmetic/sum_element_vertex.hpp" nogil:

    void sum_element_vertex[T](
        uint32_t num_vertex_per_element,
        size_t num_elements,
        uint32_t num_node_per_element,
        const T * element_to_node,
        uint32_t num_components,
        const cutil.double_ptr node_values,
        cutil.double_ptr element_value)

