#ifndef AFFECT_CONNECT_GLOBAL_TO_LOCAL_HPP
#define AFFECT_CONNECT_GLOBAL_TO_LOCAL_HPP

#include <limits>

#include <util/aligned_array.hpp>

//#define TIMERS_ACTIVE
#include <util/timer.hpp>

//#include "print_array.hpp"

#define AFFECT_CONNECT_ROUND_DOWN(x, s) ((x) & ~((s)-1))

/**
 * Given the global connectivity array, compute a local connectivity array and the global to local ID array.
 *
 * At the same time we are reducing the 64 signed global indices to local indices of 32 bit unsigned integers.
 *
 * Global indices with signed 64 bit integers have max value is (2^63 − 1), whereas
 * local indices must fit with unsignned 32 bit integers (2^32 − 1 or 4,294,967,295).
 *
 * The function will return 0 if the unsigned 32 bit limit is exceeded.
 *
 * @param element_to_vertex_global   (IN) array of element-to-vertex global connectivity,
 *                                        length (numElement*numVertexPerElement)
 * @param num_entry                  (IN) size of the connectivity array element_to_vertex_global
 *                                        (numElement*numVertexPerElement)
 * @param max_global                 (IN) upper bound on the globalID entry values, and
 *                                        size of the global_to_local array
 * @param element_to_vertex_local   (OUT) array of element-to-vertex local connectivity,
 *                                        length (num_entry)
 * @param max_global_index          (OUT) the highest node globalID contained in the element_to_vertex_global array
 * @param global_to_local           (OUT) array for each globalID, the corresponding value of the localID + 1
 *                                        length (max_global)
 * @return                          (OUT) number of unique globalID values,
 *                                        the number of non-zero values in global_to_local array
 *                                        return value is zero if number of local IDs exceeds UINT32_MAX, or
 *                                        if num_entry or max_global is zero.
 */
template <typename T>  // where T is int64 or int32
uint32_t compute_global_to_local(
    size_t num_entry,
    T max_global,
    const T * __restrict element_to_vertex_global,
    T * __restrict max_global_index,
    T * __restrict min_global_index,
    uint32_t * __restrict element_to_vertex_local,
    uint32_t * __restrict global_to_local ) // working space of length at least max_global
{
    T k;
    //uint32_t* local;
    uint32_t num_unique = 0;

    // NOTE: this experimental method using a sort could save memory by doing it in place without as much temporary
    // space, but is slower by a factor of 10 even using the threaded sort in Intel tbb::parallel_sort().
    /*
    START_TIMER(compute_global_to_local_sort);
    compute_global_to_local_sort(
        num_entry,
        max_global,
        element_to_vertex_global,
        max_global_index,
        element_to_vertex_local,
        global_to_local);
    END_TIMER(compute_global_to_local_sort);
    exit(1);
    */

    aligned::zero(global_to_local, max_global);

    // count the unique set of nodes used by the elements
    START_TIMER(compute_global_to_local_count_unique);

    // do we need to worry about possibly exceeding 32 bit limit of unique nodes?
    if (num_entry > (size_t)std::numeric_limits<uint32_t>::max()) {

        *max_global_index = -1;
        *min_global_index = std::numeric_limits<T>::max();

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
                if (num_unique == std::numeric_limits<uint32_t>::max())
                    break;

                // store the new local id _after_ we increment by one so we can distinguish from zero
                global_to_local[k] = ++num_unique;
            }
        }
        if (num_unique == std::numeric_limits<uint32_t>::max()) { // we exceeded limit of 32 bit unsigned integers
            num_unique = 0;
            *max_global_index = -1;
            *min_global_index = std::numeric_limits<T>::max();
            return num_unique;
        }
    } else {
        //
        // no need to worry about exceeding 32 bit limit of unique nodes used by element to vertex array
        //

        /* version with aligned min max calculation done separately from the entries
        // get min max of global vertex IDs used by elements
        auto min_max = aligned::min_max(element_to_vertex_global, num_entry);
        *min_global_index = min_max.first;
        *max_global_index = min_max.second;

        for (size_t i = 0; i < num_entry; ++i) { // loop through global connectivity and add entries at locations
            k = element_to_vertex_global[i];
            if (global_to_local[k] == 0) {  // we check the stored values of unique (localID + 1)
                // this is a new global node we haven't seen yet
                // change the connectivity to the new local id
                // increment num_unique by one
                // store the new local id _after_ we increment by one so we can distinguish from zero
                global_to_local[k] = ++num_unique;
            }
        }
        */

        /* serial version performing the min, max and entries at the same time
        *max_global_index = -1;
        *min_global_index = INT64_MAX;
        for (size_t i = 0; i < num_entry; ++i) { // loop through global connectivity and add entries at locations

            k = element_to_vertex_global[i];

            if (global_to_local[k] == 0) {  // we check the stored values of unique (localID + 1)
                // this is a new global node we haven't seen yet
                // change the connectivity to the new local id
                // increment num_unique by one

                // store the new local id _after_ we increment by one so we can distinguish from zero
                global_to_local[k] = ++num_unique;

                // keep track of the maximum/minimum global node we have seen
                *max_global_index = k > *max_global_index ? k : *max_global_index;
                *min_global_index = k < *min_global_index ? k : *min_global_index;
            }
        }
        */

        /* hand unrolled version */
        const size_t STEP_SIZE = 4;
        register size_t i = 0;
        register T k0, k1, k2, k3;
        register T max = -1;
        register T min = std::numeric_limits<T>::max();

        for(; i < AFFECT_CONNECT_ROUND_DOWN(num_entry, STEP_SIZE); i += STEP_SIZE)
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



        /* threaded version with multiple locks, each locking a section of the array
        auto min_max = aligned::min_max(element_to_vertex_global, num_entry);
        *min_global_index = min_max.first;
        *max_global_index = min_max.second;

        const int lock_count = 4;
        omp_lock_t locks[lock_count];
        for (int i = 0; i < lock_count; ++i)
            omp_init_lock(&locks[i]);

        for (size_t i = 0; i < num_entry; ++i) { // loop through global connectivity and add entries at locations
            k = element_to_vertex_global[i];
            omp_set_lock(&locks[k % lock_count]);
            if (global_to_local[k] == 0) {  // we check the stored values of unique (localID + 1)
                // this is a new global node we haven't seen yet
                // change the connectivity to the new local id
                // increment num_unique by one
                // store the new local id _after_ we increment by one so we can distinguish from zero
                global_to_local[k] = ++num_unique;
            }
            omp_unset_lock(&locks[k % lock_count]);
        }

        for (int i = 0; i < lock_count; ++i)
            omp_destroy_lock(&locks[i]);
        */


    }

    END_TIMER(compute_global_to_local_count_unique);

    START_TIMER(compute_global_to_local_set_local);

    // loop through global connectivity and add entries at locations
    #pragma omp parallel for schedule(static, 1024) shared(element_to_vertex_local, num_entry, global_to_local)
    for (size_t i = 0; i < num_entry; ++i) {
        // change the connectivity to the local id
        element_to_vertex_local[i] = global_to_local[element_to_vertex_global[i]] - 1; // remember we incremented by +1
    }

    END_TIMER(compute_global_to_local_set_local);

    return num_unique;
}


/**
 * Given the global_to_local connectivity array, compute the local_to_global array map.
 *
 * @param max_global_index          (IN) the highest node globalID contained in the element_to_vertex_global array
 * @param max_global_index          (IN) the highest node globalID contained in the element_to_vertex_global array
 * @param global_to_local           (IN) array for each globalID, the corresponding value of the localID + 1
 *                                       length (max_global)
 * @param local_to_global          (OUT) array for each localID, the corresponding globalID
 *                                       length is number of unique globalID values,
 *                                       the number of non-zero values in global_to_local array
 */
template <typename T>
void fill_local_to_global(
    T max_global_index,
    T min_global_index,
    const uint32_t * __restrict global_to_local,
    T * __restrict local_to_global)  // size num_unique entries that were non-zero in global_to_local
{
    START_TIMER(fill_local_to_global);
    #pragma omp parallel for schedule(static, 1024) shared(local_to_global)
    for (int64_t k = min_global_index; k < max_global_index+1; ++k) {
        if (global_to_local[k] > 0) {
            local_to_global[global_to_local[k] - 1] = k; // remember we started with 1
        }
    }
    END_TIMER(fill_local_to_global);
}



/*
#include <algorithm>
#include <vector>
#include <tbb/parallel_sort.h>
void compute_global_to_local_sort(
    int64_t num_entry,
    int64_t max_global,
    const int64_t * element_to_vertex_global,
    int64_t * max_global_index,
    uint32_t * element_to_vertex_local,
    uint32_t * global_to_local ) // working space of length at least max_global
{
    //print_1d_array("element_to_vertex_global", element_to_vertex_global, num_entry);

    std::vector<uint32_t> idx(num_entry);
    std::iota(idx.begin(), idx.end(), 0); // initialize original index locations

    // perform a sort_by_key
    // sort the values (indices) based on comparing values of the key (element_to_vertex_global)

    tbb::parallel_sort(idx.begin(), idx.end(),
              [&element_to_vertex_global](size_t i1, size_t i2) {
                return element_to_vertex_global[i1] < element_to_vertex_global[i2];});

    //std::sort(idx.begin(), idx.end(),
    //          [&element_to_vertex_global](size_t i1, size_t i2) {
    //            return element_to_vertex_global[i1] < element_to_vertex_global[i2];});


    //print_1d_array("idx", idx.data(), idx.size());

    // make element to local vertex array
    // loop through range of equal global vertices,
    // increment num_unique when the global vertex changes
    uint32_t num_unique = 0;
    element_to_vertex_local[idx[num_unique]] = num_unique;
    int64_t previous = element_to_vertex_global[idx[0]];
    for (uint32_t k = 1; k < idx.size(); ++k) {
        if (element_to_vertex_global[idx[k]] != previous) {
            ++num_unique;
            previous = element_to_vertex_global[idx[k]];

            // we could also build the local_to_global right here
        }
        element_to_vertex_local[idx[k]] = num_unique;
    }
    //print_1d_array("element_to_vertex_local", element_to_vertex_local, num_entry);

    // no longer need idx

    auto last = std::unique(idx.begin(), idx.end(),
              [&element_to_vertex_global](size_t i1, size_t i2) {
                return element_to_vertex_global[i1] == element_to_vertex_global[i2];});
    idx.erase(last, idx.end());
    //print_1d_array("idx", idx.data(), idx.size());


    // make local to global array
    std::vector<int64_t> local_to_global(idx.size());
    #pragma omp parallel for
    for (size_t k = 0; k < idx.size(); ++k)
        local_to_global[k] = element_to_vertex_global[ idx[k] ];
    //print_1d_array("local_to_global", local_to_global.data(), local_to_global.size());
}
*/

#undef AFFECT_CONNECT_ROUND_DOWN

#endif  // AFFECT_CONNECT_GLOBAL_TO_LOCAL_HPP