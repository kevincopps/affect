#ifndef AFFECT_ARITHMETIC_AVERAGE_ELEMENT_VERTEX_HPP
#define AFFECT_ARITHMETIC_AVERAGE_ELEMENT_VERTEX_HPP

#include <util/aligned_array.hpp>

#include "sum_element_vertex.hpp"

template <typename T>
void average_element_vertex(
    uint32_t num_vertex_per_element,
    size_t num_elements,
    uint32_t num_node_per_element,
    const T * element_to_node,
    uint32_t num_components,
    const aligned::double_ptr __restrict node_values,
    aligned::double_ptr __restrict element_value) {

    sum_element_vertex(num_vertex_per_element,
                       num_elements,
                       num_node_per_element,
                       element_to_node,
                       num_components,
                       node_values,
                       element_value);

    double denominator = (double)num_vertex_per_element;
    aligned::vector_scalar_divide(element_value, num_elements * num_components, denominator);
}

#endif // AFFECT_ARITHMETIC_AVERAGE_ELEMENT_VERTEX_HPP