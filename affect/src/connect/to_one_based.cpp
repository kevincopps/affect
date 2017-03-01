#include <cstdint>

void to_one_based(int64_t n, int64_t* array) {
    int64_t* last = array + n;
    for ( ; array != last; ++array) {
        ++(*array);
    }
}

