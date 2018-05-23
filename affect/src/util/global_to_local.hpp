#ifndef AFFECT_UTIL_GLOBAL_TO_LOCAL_HPP
#define AFFECT_UTIL_GLOBAL_TO_LOCAL_HPP

#include "aligned_array.hpp"

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
uint32_t compute_global_to_local(
    size_t num_entry,
    int64_t max_global,
    const aligned::int64_ptr __restrict element_to_vertex_global,
    int64_t * max_global_index,
    int64_t * min_global_index,
    aligned::uint32_ptr __restrict element_to_vertex_local,
    aligned::uint32_ptr __restrict global_to_local);


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
void fill_local_to_global(
    int64_t max_global_index,
    int64_t min_global_index,
    const aligned::uint32_ptr __restrict global_to_local,
    aligned::int64_ptr __restrict local_to_global);  // size num_unique entries that were non-zero in global_to_local


#endif  // AFFECT_UTIL_GLOBAL_TO_LOCAL_HPP