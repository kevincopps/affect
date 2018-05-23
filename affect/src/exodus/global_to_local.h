#include <stdint.h>

extern uint32_t compute_global_to_local64(
    size_t num_entry,
    const int64_t * __restrict element_to_vertex_global,
    size_t * __restrict max_global_index,
    size_t * __restrict min_global_index,
    uint32_t * __restrict element_to_vertex_local,
    uint32_t * __restrict global_to_local); // working space of length at least max_global

extern void fill_local_to_global64(
    size_t max_global_index,
    size_t min_global_index,
    const uint32_t * __restrict global_to_local,
    int64_t * __restrict local_to_global);

extern uint32_t compute_global_to_local32(
    size_t num_entry,
    const int32_t * __restrict element_to_vertex_global,
    size_t * __restrict max_global_index,
    size_t * __restrict min_global_index,
    uint32_t * __restrict element_to_vertex_local,
    uint32_t * __restrict global_to_local); // working space of length at least max_global

extern void fill_local_to_global32(
    size_t max_global_index,
    size_t min_global_index,
    const uint32_t * __restrict global_to_local,
    int32_t * __restrict local_to_global);