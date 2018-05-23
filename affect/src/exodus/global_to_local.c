#include <stdint.h>

#include "array_util.h"

#define ROUND_DOWN(x, s) ((x) & ~((s)-1))

/**
 * Convert global to local indices.
 *
 * The global connectivity is input with signed 64 bit integers, where max value is (2^63 − 1),
 * while the local connectivity is assumed to fit inside an
 * unsigned 32 bit integers (2^32 − 1 or 4,294,967,295).
 *
 * Before calling you must zero global_to_local
 *          ex_aligned_zero_uint32(global_to_local, num_global_nodes);
 *
 */
uint32_t compute_global_to_local64(
    size_t num_entry,
    const int64_t * __restrict element_to_vertex_global,
    size_t * __restrict max_global_index,
    size_t * __restrict min_global_index,
    uint32_t * __restrict element_to_vertex_local,
    uint32_t * __restrict global_to_local ) // working space of length at least num_global_nodes
{
    int64_t k;
    uint32_t num_unique = 0;

    // do we need to worry about possibly exceeding 32 bit limit of unique nodes?
    if (num_entry > (size_t)UINT32_MAX) {

        *max_global_index = -1;
        *min_global_index = INT64_MAX;

        for (uint32_t i = 0; i < num_entry; ++i) { // loop through global connectivity and add entries at locations

            k = element_to_vertex_global[i];

            if (global_to_local[k] == 0) {  // we check the stored values of unique (localID + 1)

                // this is a new global node we haven't seen yet
                // change the connectivity to the new local id
                // increment num_unique by one

                // keep track of the maximum global node we have seen
                *max_global_index = k > *max_global_index ? k : *max_global_index;
                *min_global_index = k < *min_global_index ? k : *min_global_index;

                // check for exceeding the unsigned 32 bit limit
                if (num_unique == UINT32_MAX)
                    break;

                // store the new local id _after_ we increment by one so we can distinguish from zero
                global_to_local[k] = ++num_unique;
            }
        }
        if (num_unique == UINT32_MAX) { // we exceeded limit of 32 bit unsigned integers
            num_unique = 0;
            *max_global_index = -1;
            *min_global_index = INT64_MAX;
            return num_unique;
        }
    } else {
        //
        // no need to worry about exceeding 32 bit limit of unique nodes used by element to vertex array
        //

        /* hand unrolled gives a little speedup because omp does not handle this block */
        const size_t STEP_SIZE = 4;
        register size_t i = 0;
        register int64_t k0, k1, k2, k3;
        register int64_t max = -1;
        register int64_t min = INT64_MAX;

        for(; i < ROUND_DOWN(num_entry, STEP_SIZE); i += STEP_SIZE)
        {
            k0 = element_to_vertex_global[i+0];
            k1 = element_to_vertex_global[i+1];
            k2 = element_to_vertex_global[i+2];
            k3 = element_to_vertex_global[i+3];
            if (global_to_local[k0] == 0) {
                global_to_local[k0] = ++num_unique;
                max = k0 > max ? k0 : max;
                min = k0 < min ? k0 : min;
            }
            if (global_to_local[k1] == 0) {
                global_to_local[k1] = ++num_unique;
                max = k1 > max ? k1 : max;
                min = k1 < min ? k1 : min;
            }
            if (global_to_local[k2] == 0) {
                global_to_local[k2] = ++num_unique;
                max = k2 > max ? k2 : max;
                min = k2 < min ? k2 : min;
            }
            if (global_to_local[k3] == 0) {
                global_to_local[k3] = ++num_unique;
                max = k3 > max ? k3 : max;
                min = k3 < min ? k3 : min;
            }
        }
        for(; i < num_entry; i++) {
            k = element_to_vertex_global[i];
            if (global_to_local[k] == 0) {
                global_to_local[k] = ++num_unique;
                max = k > max ? k : max;
                min = k < min ? k : min;
            }
        }
        *max_global_index = max;
        *min_global_index = min;
    }

    // loop through global connectivity and add entries at locations
    #pragma omp parallel for schedule(static, 1024) shared(element_to_vertex_local, num_entry, global_to_local)
    for (size_t i = 0; i < num_entry; ++i) {
        // change the connectivity to the local id
        element_to_vertex_local[i] = global_to_local[element_to_vertex_global[i]] - 1; // remember we incremented by +1
    }

    return num_unique;
}

// in between calls allocate local_to_global(num_unique * sizeof(uint32_t))

void fill_local_to_global64(
    int64_t max_global_index,
    int64_t min_global_index,
    const uint32_t * __restrict global_to_local,
    int64_t * __restrict local_to_global)  // size num_unique entries that were non-zero in global_to_local
{
    #pragma omp parallel for schedule(static, 1024) shared(local_to_global)
    for (int64_t k = min_global_index; k < max_global_index+1; ++k) {
        if (global_to_local[k] > 0) {
            local_to_global[global_to_local[k] - 1] = k; // remember we started with 1
        }
    }
}


// before calling you must zero global_to_local
//      ex_aligned_zero_uint32(global_to_local, num_global_nodes);

uint32_t compute_global_to_local32(
    size_t num_entry,
    const int32_t * __restrict element_to_vertex_global,
    size_t * __restrict max_global_index,
    size_t * __restrict min_global_index,
    uint32_t * __restrict element_to_vertex_local,
    uint32_t * __restrict global_to_local ) // working space of length at least num_global_nodes
{
    int32_t k;
    uint32_t num_unique = 0;

    // do we need to worry about possibly exceeding 32 bit limit of unique nodes?
    if (num_entry > (size_t)UINT32_MAX) {

        *max_global_index = -1;
        *min_global_index = INT32_MAX;

        for (uint32_t i = 0; i < num_entry; ++i) { // loop through global connectivity and add entries at locations

            k = element_to_vertex_global[i];

            if (global_to_local[k] == 0) {  // we check the stored values of unique (localID + 1)

                // this is a new global node we haven't seen yet
                // change the connectivity to the new local id
                // increment num_unique by one

                // keep track of the maximum global node we have seen
                *max_global_index = k > *max_global_index ? k : *max_global_index;
                *min_global_index = k < *min_global_index ? k : *min_global_index;

                // check for exceeding the unsigned 32 bit limit
                if (num_unique == UINT32_MAX)
                    break;

                // store the new local id _after_ we increment by one so we can distinguish from zero
                global_to_local[k] = ++num_unique;
            }
        }
        if (num_unique == UINT32_MAX) { // we exceeded limit of 32 bit unsigned integers
            num_unique = 0;
            *max_global_index = -1;
            *min_global_index = INT32_MAX;
            return num_unique;
        }
    } else {
        //
        // no need to worry about exceeding 32 bit limit of unique nodes used by element to vertex array
        //

        /* hand unrolled gives a little speedup because omp does not handle this block */
        const size_t STEP_SIZE = 4;
        register size_t i = 0;
        register int32_t k0, k1, k2, k3;
        register int32_t max = -1;
        register int32_t min = INT32_MAX;

        for(; i < ROUND_DOWN(num_entry, STEP_SIZE); i += STEP_SIZE)
        {
            k0 = element_to_vertex_global[i+0];
            k1 = element_to_vertex_global[i+1];
            k2 = element_to_vertex_global[i+2];
            k3 = element_to_vertex_global[i+3];
            if (global_to_local[k0] == 0) {
                global_to_local[k0] = ++num_unique;
                max = k0 > max ? k0 : max;
                min = k0 < min ? k0 : min;
            }
            if (global_to_local[k1] == 0) {
                global_to_local[k1] = ++num_unique;
                max = k1 > max ? k1 : max;
                min = k1 < min ? k1 : min;
            }
            if (global_to_local[k2] == 0) {
                global_to_local[k2] = ++num_unique;
                max = k2 > max ? k2 : max;
                min = k2 < min ? k2 : min;
            }
            if (global_to_local[k3] == 0) {
                global_to_local[k3] = ++num_unique;
                max = k3 > max ? k3 : max;
                min = k3 < min ? k3 : min;
            }
        }
        for(; i < num_entry; i++) {
            k = element_to_vertex_global[i];
            if (global_to_local[k] == 0) {
                global_to_local[k] = ++num_unique;
                max = k > max ? k : max;
                min = k < min ? k : min;
            }
        }
        *max_global_index = max;
        *min_global_index = min;
    }

    // loop through global connectivity and add entries at locations
    #pragma omp parallel for schedule(static, 1024) shared(element_to_vertex_local, num_entry, global_to_local)
    for (size_t i = 0; i < num_entry; ++i) {
        // change the connectivity to the local id
        element_to_vertex_local[i] = global_to_local[element_to_vertex_global[i]] - 1; // remember we incremented by +1
    }

    return num_unique;
}

void fill_local_to_global32(
    int32_t max_global_index,
    int32_t min_global_index,
    const uint32_t * __restrict global_to_local,
    int32_t * __restrict local_to_global)  // size num_unique entries that were non-zero in global_to_local
{
    #pragma omp parallel for schedule(static, 1024) shared(local_to_global)
    for (int32_t k = min_global_index; k < max_global_index+1; ++k) {
        if (global_to_local[k] > 0) {
            local_to_global[global_to_local[k] - 1] = k; // remember we started with 1
        }
    }
}
