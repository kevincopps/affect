#include <cstdint>

void to_zero_based(int64_t n, int64_t* array) {
    for (; n > 0; ++array, --n)
        --(*array);
    //int64_t* last = array + n;
    //for ( ; array != last; ++array) {
    //    --(*array);
    //}
}

