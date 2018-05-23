#ifndef AFFECT_UTIL_ALIGNED_ARRAY_HPP
#define AFFECT_UTIL_ALIGNED_ARRAY_HPP

#include <stdlib.h> // for posix_memalign
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <limits>
#include <utility>
#include <omp.h>


// byte alignment good enough for Intel Phi (AVX512 intrinsic extensions).
// also a common cache-line size.
#define AFFECT_UTIL_ALIGN 64

namespace aligned {

typedef int64_t * int64_ptr __attribute__((aligned (AFFECT_UTIL_ALIGN)));
typedef uint_least32_t * uint32_ptr __attribute__((aligned (AFFECT_UTIL_ALIGN)));
typedef int8_t * int8_ptr __attribute__((aligned (AFFECT_UTIL_ALIGN)));
typedef double * double_ptr __attribute__((aligned (AFFECT_UTIL_ALIGN)));

/**
 * Allocate size bytes of uninitialized storage whose alignment is default bytes.
 *
 * This is provided for calling from Cython, inside C++ prefer this is provided for calling from Cython.
 *
 * The size parameter must be an integral multiple of alignment.
 *
 * The number of bytes allocated is sizeof(T)*length.
 *
 * On success, returns the pointer to the beginning of newly allocated memory, or else nullptr.
 * The returned pointer must be deallocated with free() or realloc().
 *
 * @param length     (IN) number of element of the array of type T to allocate,
 */
template <typename T>
inline T* allocate(size_t length) {
    void *mem = NULL;
    int error = posix_memalign(&mem, AFFECT_UTIL_ALIGN, sizeof(T) * length);
    if (error != 0)
        return NULL;
    return static_cast<T*>(mem);
}


inline void* allocate_void(size_t num_bytes) {
    void *mem = NULL;
    int error = posix_memalign(&mem, AFFECT_UTIL_ALIGN, num_bytes);
    if (error != 0)
        return NULL;
    return mem;
}


template <typename T>
inline void fill(T* a, size_t n, const T x) {
    #pragma omp parallel for simd aligned(a:AFFECT_UTIL_ALIGN) shared(a,n) schedule(static)
    for (size_t i = 0; i < n; ++i)
        a[i] = x;
}

// generic type specific
template <typename T>
inline void zero(T* a, size_t n) {
    fill<T>(a, n, static_cast<T>(0));
}

template <typename T>
inline void unaligned_zero(T* a, size_t n) {
    #pragma omp parallel for simd shared(a, n) schedule(static)
    for (size_t i = 0; i < n; ++i)
        a[i] = static_cast<T>(0);
}

template <typename T>
inline T min(const T* a, size_t n) {
    T minimum = std::numeric_limits<T>::max();
#ifdef __clang__ // llvm 4.0.1 does not yet support omp parallel for simd
    #pragma omp simd aligned(a:AFFECT_UTIL_ALIGN) reduction(min:minimum)
    for (size_t i = 0; i < n; ++i) {
        minimum = a[i] < minimum ? a[i] : minimum;
    }
#else
    #pragma omp parallel for simd aligned(a:AFFECT_UTIL_ALIGN) shared(a, n) reduction(min:minimum)
    for (size_t i = 0; i < n; ++i) {
        minimum = a[i] < minimum ? a[i] : minimum;
    }
#endif
    return minimum;
}

template <typename T>
inline T max(const T* a, size_t n) {
    T maximum = std::numeric_limits<T>::min();
#ifdef __clang__ // llvm 4.0.1 does not yet support omp parallel for simd
    #pragma omp simd aligned(a:AFFECT_UTIL_ALIGN) reduction(max:maximum)
    for (size_t i = 0; i < n; ++i) {
        maximum = a[i] > maximum ? a[i] : maximum;
    }
#else
    #pragma omp parallel for simd aligned(a:AFFECT_UTIL_ALIGN) shared(a, n) reduction(max:maximum)
    for (size_t i = 0; i < n; ++i) {
        maximum = a[i] > maximum ? a[i] : maximum;
    }
#endif
    return maximum;
}

template <typename T>
inline std::pair<T,T> min_max(const T* a, size_t n) {
    T minimum = std::numeric_limits<T>::max();
    T maximum = std::numeric_limits<T>::min();
#ifdef __APPLE__ // llvm 4.0.1 does not yet support omp parallel for simd
    #pragma omp simd aligned(a:AFFECT_UTIL_ALIGN) reduction(min:minimum) reduction(max:maximum)
    for (size_t i = 0; i < n; ++i) {
        minimum = a[i] < minimum ? a[i] : minimum;
        maximum = a[i] > maximum ? a[i] : maximum;
    }
#else
    #pragma omp parallel for simd aligned(a:AFFECT_UTIL_ALIGN) shared(a, n) reduction(min:minimum) reduction(max:maximum) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        minimum = a[i] < minimum ? a[i] : minimum;
        maximum = a[i] > maximum ? a[i] : maximum;
    }
#endif
    return std::make_pair(minimum, maximum);
}

// where entries of array 'a' are to be used as an index, check values
template <typename T>
inline int is_index_out_of_range(const T* a, size_t n, size_t minimum, size_t
maximum) {
    std::pair<T,T> min_max_pair = min_max(a, n);
    // some extra checking because T may be signed integer type
    if (min_max_pair.first < 0 || min_max_pair.second < 0)
        return 1;
    if (min_max_pair.first < static_cast<T>(minimum))
        return 1;
    if (static_cast<size_t>(min_max_pair.second) > maximum)
        return 1;
    return 0;
}

// take subarrays at indices, where shape of array a is (num_fields, num_components)
template <typename T, typename U>
inline void take(T* fields, size_t num_components, U* indices, size_t num_indices, U shift, T* out) {
    if (shift == 0) {
        #pragma omp parallel for simd aligned(fields,indices,out:AFFECT_UTIL_ALIGN) \
            shared(fields, num_components, indices, num_indices, shift, out) schedule(static)
        for (size_t i = 0; i < num_indices; ++i) {
            T* fields_start = &fields[indices[i] * num_components];
            T* out_start = &out[i * num_components];
            for (size_t j = 0; j < num_components; ++j)
                out_start[j] = fields_start[j];
        }
    }
    else {
        #pragma omp parallel for simd aligned(fields,indices,out:AFFECT_UTIL_ALIGN) \
            shared(fields, num_components, indices, num_indices, shift, out) schedule(static)
        for (size_t i = 0; i < num_indices; ++i) {
            T* fields_start = &fields[(indices[i] - shift) * num_components];
            T* out_start = &out[i * num_components];
            for (size_t j = 0; j < num_components; ++j)
                out_start[j] = fields_start[j];
        }
    }
}


// x[i] /= y
template <typename T>
inline void vector_scalar_divide(T* x, size_t n, T y) {
#ifdef __clang__ // llvm 4.0.1 does not yet support omp parallel for simd
    #pragma omp simd aligned(x:AFFECT_UTIL_ALIGN)
    for (size_t i = 0; i < n; ++i) {
        x[i] /= y;
    }
#else
    #pragma omp parallel for simd aligned(x:AFFECT_UTIL_ALIGN) shared(x, n, y)
    for (size_t i = 0; i < n; ++i) {
        x[i] /= y;
    }
#endif
}

/*

// Specializations below for some scalar types, measurably faster on LLVM Apple platforms.
#ifdef __APPLE__
    #define MEMSET_SPECIAL(TYPE, NBYTES) \
template<> \
inline void zero_array<TYPE>(TYPE * array, size_t n) { \
    static const TYPE pattern##TYPE = 0; \
    memset_pattern##NBYTES((void *)array, (const void *)&pattern##TYPE, n * sizeof(TYPE)); \
}
#else
    #define MEMSET_SPECIAL(TYPE, NBYTES)
#endif


MEMSET_SPECIAL(int32_t, 4);
MEMSET_SPECIAL(uint_least32_t, 4);
MEMSET_SPECIAL(int64_t, 8);
MEMSET_SPECIAL(uint64_t, 8);
MEMSET_SPECIAL(float, 4);
MEMSET_SPECIAL(double, 8);
*/

} // namespace aligned

#endif // AFFECT_UTIL_ALIGNED_ARRAY_HPP