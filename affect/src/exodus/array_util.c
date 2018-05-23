/*
   Activate posix features, for posix_memalign, and maybe struct timespec
*/
#define _POSIX_C_SOURCE 200112L

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "array_util.h"

/*
   64 byte (not bits) alignment good enough for Intel Phi (AVX512 intrinsic extensions).
*/
#define ARRAY_UTIL_ALIGNMENT 64

/*
   Probably requires compiling with -std=gnu99 on Linux gcc-7
*/
/*
#include <time.h>
#define TIMER_TYPE CLOCK_REALTIME

double time_diff(struct timespec* start, struct timespec* end)
{
    struct timespec temp;
    if ((end->tv_nsec - start->tv_nsec) < 0) {
        temp.tv_sec = end->tv_sec - start->tv_sec - 1;
        temp.tv_nsec = 1000000000+end->tv_nsec - start->tv_nsec;
    } else {
        temp.tv_sec = end->tv_sec - start->tv_sec;
        temp.tv_nsec = end->tv_nsec - start->tv_nsec;
    }
    return (double)temp.tv_sec + (double)temp.tv_nsec * 1E-9;
}
*/


void* ex_aligned_allocate(size_t size) {
    void *mem = NULL;
    int error = posix_memalign(&mem, ARRAY_UTIL_ALIGNMENT, size);
    if (error != 0)
        return NULL;
    return mem;
}

void ex_aligned_copy_stride(double* __restrict source, size_t source_length,
                            double* __restrict destination, int destination_stride) {
/*
    int repeat = 5;
    struct timespec time1, time2;
    double avg_rate = 0.0, peak_rate = 0.0;

    for (int k = 0; k < 5; ++k) {
        clock_gettime(TIMER_TYPE, &time1);
        for(int r=0; r<repeat; r++) {
*/
    #pragma omp parallel for simd aligned(source,destination:ARRAY_UTIL_ALIGNMENT) schedule(static) \
            shared(source,destination) firstprivate(destination_stride)
    for (size_t i = 0; i < source_length; ++i) {
        destination[i*destination_stride] = source[i];
    }
/*
        }
        clock_gettime(TIMER_TYPE, &time2);

        double dtime = time_diff(&time1,&time2);
        double rate = 1E-9 * source_length * sizeof(double) * repeat/dtime;
        avg_rate += rate;
        peak_rate = rate > peak_rate ? rate : peak_rate;
        printf("dtime %f, %f GB/s\n", dtime, rate);
    }
    printf("avg_rate %f GB/s  peak_rate %f GB/s\n", avg_rate/5., peak_rate);
*/
}

void ex_aligned_to_zero_based_int32(int32_t * array, size_t n) {
    #pragma omp parallel for simd aligned(array: ARRAY_UTIL_ALIGNMENT) shared(array) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        --array[i];
    }
}

void ex_aligned_to_zero_based_int64(int64_t * array, size_t n) {

/*
    int repeat = 5;
    struct timespec time1, time2;
    double avg_rate = 0.0, peak_rate = 0.0;

    for (int k = 0; k < 5; ++k) {
        clock_gettime(TIMER_TYPE, &time1);
        for(int r=0; r<repeat; r++) {
*/

    #pragma omp parallel for simd aligned(array: ARRAY_UTIL_ALIGNMENT) shared(array) schedule(static)
    for (size_t i = 0; i < n; ++i) {
        --array[i];
    }
/*
        }
        clock_gettime(TIMER_TYPE, &time2);

        double dtime = time_diff(&time1,&time2);
        double rate = 1E-9 * n * sizeof(double) * repeat/dtime;
        avg_rate += rate;
        peak_rate = rate > peak_rate ? rate : peak_rate;
        printf("dtime %f, %f GB/s\n", dtime, rate);
    }
    printf("avg_rate %f GB/s  peak_rate %f GB/s\n", avg_rate/5., peak_rate);
*/

}

void ex_aligned_to_one_based_int32(int32_t* array, size_t n) {
    #pragma omp simd aligned(array: ARRAY_UTIL_ALIGNMENT)
    for (size_t i = 0; i < n; ++i) {
        ++array[i];
    }
}

void ex_aligned_to_one_based_int64(int64_t* array, size_t n) {
    #pragma omp simd aligned(array: ARRAY_UTIL_ALIGNMENT)
    for (size_t i = 0; i < n; ++i) {
        ++array[i];
    }
}

void ex_aligned_zero_uint32(uint32_t* a, size_t n) {
    #pragma omp parallel for simd aligned(a:ARRAY_UTIL_ALIGNMENT) shared(a,n) schedule(static)
    for (size_t i = 0; i < n; ++i)
        a[i] = 0;
}
