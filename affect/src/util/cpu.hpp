#ifndef AFFECT_CPU_HPP
#define AFFECT_CPU_HPP

/*
 * Test CPU and compile environment for ability to use AVX-512 instructions, specifically AVX-512F (foundation),
 * AVX-512CD (conflict detection), AVX-512ER (exponential and reciprocal) and AVX-512PF (prefetch).
 *
 * If we want an application to run everywhere, in order to use these instructions in a program, we need to make
 * sure that the operating system and the processor have support for them when the application is run.
 *
 * Works with at least the Intel compiler, gcc, clang/LLVM and Microsoft compilers.
 */
int can_use_intel_knl_features(void);


#endif // AFFECT_CPU_HPP