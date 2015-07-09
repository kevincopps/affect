#!/usr/bin/env python

from libc.stdint cimport int64_t

cdef extern from "connect.hpp":

    void to_zero_based(int64_t n, int64_t* array)

    void to_one_based(int64_t n, int64_t* array)
