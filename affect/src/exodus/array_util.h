#ifndef AFFECT_EXODUS_ARRAY_UTIL_H
#define AFFECT_EXODUS_ARRAY_UTIL_H

#include <stddef.h>
#include <stdint.h>


/**
 * Allocate size bytes of uninitialized storage whose alignment is 64 bytes.
 * The size parameter must be an integral multiple of alignment.
 *
 * Use standard library free() to free these pointers.
 *
 * This is used to allocate arrays aligned for use of SIMD instructions and
 * #pragma omp simd
 *
 * On success, returns the pointer to the beginning of newly allocated memory.
 * The returned pointer must be deallocated with free() or realloc().
 * On failure, returns a null pointer.
 *
 * @param size     (IN) number of bytes to allocate. An integral multiple of alignment (64).
 */
extern void* ex_aligned_allocate(size_t size);

/**
 * Copy N elements of the source vector into the destination vector with a destination stride.
 */
extern void ex_aligned_copy_stride(
    double* __restrict source,
    size_t source_length,
    double* __restrict destination,
    int destination_stride);

/**
 * Subtract one from every entry in the given array.
 *
 * Array must be aligned for use of SIMD instructions.
 *
 * For converting a one-based connectivity array (Fortran style)
 * to zero-based (C style), that is,
 * subtract one from every entry in the given array.
 *
 * @param n     (IN) length of array
 * @param array (INOUT) the one-based integer array to convert
 */
extern void ex_aligned_to_zero_based_int32(int32_t* array, size_t n);
extern void ex_aligned_to_zero_based_int64(int64_t* array, size_t n);


/**
 * Add one to every entry in the given array.
 *
 * Array must be aligned for use of SIMD instructions.
 *
 * For converting onvert a zero-based connectivity array (C style)
 * to one-based (Fortran style), that is,
 *
 *
 * @param n     (IN) length of array
 * @param array (INOUT) a zero-based integer array to convert
 */
extern void ex_aligned_to_one_based_int32(int32_t* array, size_t n);
extern void ex_aligned_to_one_based_int64(int64_t* array, size_t n);

/**
 * Set all entries of array to zero.
 */
extern void ex_aligned_zero_uint32(uint32_t* a, size_t n);

#endif // AFFECT_EXODUS_ARRAY_UTIL_H