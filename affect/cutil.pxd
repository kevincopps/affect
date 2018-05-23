#!/usr/bin/env python

from libc.stdint cimport int8_t, uint8_t
from libc.stdint cimport uint32_t
from libc.stdint cimport int64_t


cdef extern from "<utility>" namespace "std":
    # noinspection PyPep8Naming
    cdef cppclass pair[U,V]:
        pair()
        void swap (pair& pr)


cdef extern from "util/aligned_array.hpp" namespace "aligned" nogil:

    # aligned pointer type definitions
    ctypedef int64_t * int64_ptr
    ctypedef uint32_t * uint32_ptr
    ctypedef int8_t * int8_ptr
    ctypedef double * double_ptr

    T* allocate[T](T a, size_t length)
    void* allocate_void(size_t num_bytes)
    void fill[T](T* a, size_t n, const T x)
    void zero[T](T* a, size_t n)
    T min[T](const T* a, size_t n)
    T max[T](const T* a, size_t n)
    pair[T,T] min_max[T](const T* a, size_t n)
    int is_index_out_of_range[T](const T* a, size_t n, size_t minimum, size_t maximum)
    void take[T, U](T* fields, size_t num_components, U* indices, size_t num_indices, U shift, T* out)
    void unaligned_zero[T](T* a, size_t n)


cdef extern from "util/cpu.hpp":
    int can_use_intel_knl_features()
    int initialize_num_threads()


cdef extern from "omp.h":
    int omp_get_num_threads()
    void omp_set_num_threads(int)
    void omp_set_dynamic(int)
