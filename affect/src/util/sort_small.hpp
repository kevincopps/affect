/*
 * Fast efficient sort of small length arrays.
 *
 * Built on the idea of a sorting network with simple swaps.
 *
 * A network sort, also known as an oblivious sort, is a set of comparators.
 * Each comparator exchanges its inputs if the inputs are out of order.
 * Since each comparator has no knowledge of what any other comparator has done,
 * the arrangement of the comparators alone is what guarantees the inputs will be sorted.
 *
 * The order of the comparators can improve performance in the sense that they may occur independently in parallel.
 * Modern compilers can take advantage of this during read/write of the memory.
 *
 * Code based on https://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-uint32_t-array
 *
 * Optimal order of sorting networks up to N=16 can be generated here http://pages.ripco.net/~jgamble/nw.html
 */

// this swap is fast as we can get
#define SWAP(x,y) { uint32_t dx = d[x]; uint32_t dy = d[y]; uint32_t tmp = d[x] = dx < dy ? dx : dy; d[y] ^= dx ^ tmp; }
//#define SWAP(x,y) asm("mov %0, %2 ; cmp %1, %0 ; cmovg %1, %0 ; cmovg %2, %1" : "=r" (x), "=r" (y), "=r" (tmp) : "0" (x), "1" (y) : "cc");

static inline void sort2(uint32_t * d) {
    SWAP(0, 1);
}

static inline void sort3(uint32_t * d) {
    SWAP(1, 2);
    SWAP(0, 2);
    SWAP(0, 1);
}

static inline void sort4(uint32_t * d) {
    SWAP(0, 1);
    SWAP(2, 3);
    SWAP(0, 2);
    SWAP(1, 3);
    SWAP(1, 2);
}

static inline void sort5(uint32_t * d) {
    SWAP(0, 1);
    SWAP(3, 4);
    SWAP(2, 4);
    SWAP(2, 3);
    SWAP(0, 3);
    SWAP(0, 2);
    SWAP(1, 4);
    SWAP(1, 3);
    SWAP(1, 2);
}

static inline void sort6(uint32_t * d) {
    SWAP(1, 2);
    SWAP(4, 5);
    SWAP(0, 2);
    SWAP(3, 5);
    SWAP(0, 1);
    SWAP(3, 4);
    SWAP(1, 4);
    SWAP(0, 3);
    SWAP(2, 5);
    SWAP(1, 3);
    SWAP(2, 4);
    SWAP(2, 3);
}

static inline void sort7(uint32_t * d) {
    SWAP(1, 2);
    SWAP(0, 2);
    SWAP(0, 1);
    SWAP(3, 4);
    SWAP(5, 6);
    SWAP(3, 5);
    SWAP(4, 6);
    SWAP(4, 5);
    SWAP(0, 4);
    SWAP(0, 3);
    SWAP(1, 5);
    SWAP(2, 6);
    SWAP(2, 5);
    SWAP(1, 3);
    SWAP(2, 4);
    SWAP(2, 3);
}

static inline void sort8(uint32_t * d) {
    SWAP(0, 1);
    SWAP(2, 3);
    SWAP(0, 2);
    SWAP(1, 3);
    SWAP(1, 2);
    SWAP(4, 5);
    SWAP(6, 7);
    SWAP(4, 6);
    SWAP(5, 7);
    SWAP(5, 6);
    SWAP(0, 4);
    SWAP(1, 5);
    SWAP(1, 4);
    SWAP(2, 6);
    SWAP(3, 7);
    SWAP(3, 6);
    SWAP(2, 4);
    SWAP(3, 5);
    SWAP(3, 4);
}

/*
 * The max_last function shifts the maximum value in the given array to be the last entry. The remaining previous entries
 * are unsorted.
 *
 * This is like running the first pass of a bubble sort.
 */
static inline void max_last3(uint32_t * d) {
    SWAP(0,1); SWAP(1,2);
}

static inline void max_last4(uint32_t * d) {
    SWAP(0,1); SWAP(1,2); SWAP(2,3);
}

#undef SWAP


