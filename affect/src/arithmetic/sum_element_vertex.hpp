#ifndef AFFECT_ARITHMETIC_SUM_ELEMENT_VERTEX_HPP
#define AFFECT_ARITHMETIC_SUM_ELEMENT_VERTEX_HPP

#include <util/aligned_array.hpp>

template <typename T>
void sum_element_vertex(
    uint32_t num_vertex_per_element,
    size_t num_elements,
    uint32_t num_node_per_element,
    const T * element_to_node,
    uint32_t num_components,
    const aligned::double_ptr __restrict node_values,
    aligned::double_ptr __restrict element_value) {

    aligned::zero(element_value, num_elements * num_components);

    #pragma omp parallel for simd aligned(element_to_node,node_values,element_value:AFFECT_UTIL_ALIGN) schedule(static)
    for (size_t element = 0; element < num_elements; ++element) {

        const size_t element_value_offset = element * num_components;
        const size_t node_offset = element * num_node_per_element;

        for (uint32_t local_vertex = 0; local_vertex < num_vertex_per_element; ++local_vertex) {

            const size_t node_value_offset = element_to_node[node_offset + local_vertex] * num_components;

            for (uint32_t i = 0; i < num_components; ++i) {
                element_value[element_value_offset + i] += node_values[node_value_offset + i];
            }
        }
    }
}

#endif // AFFECT_ARITHMETIC_SUM_ELEMENT_VERTEX_HPP

