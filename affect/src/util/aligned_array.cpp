#include <stdlib.h>

#ifdef __APPLE__
#include <xmmintrin.h>
#endif

#include "aligned_array.hpp"

namespace aligned {

/* this is inligned now in the header
void* aligned_allocate(size_t size) {
#ifdef __APPLE__ // LLVM Clang 4.0.1 does not seem to support C11 function yet
    return _mm_malloc(size, AFFECT_UTIL_ALIGN);
#else
       // Use the standard C11 function.
       // Passing a size which is not an integral multiple of alignment or a alignment
       // which is not valid or not supported by the implementation will probably cause
       // this function to fail and return a null pointer
    return aligned_alloc(AFFECT_UTIL_ALIGN, size);
#endif
}
*/

} // namespace aligned

