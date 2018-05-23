#ifndef AFFECT_INTERSECT_SETS_H
#define AFFECT_INTERSECT_SETS_H

//#include <util/intersection_sse.hpp>
//#include <util/intersection_avx.hpp>
//#include <util/intersection_galloping_sse.hpp>
//#include <util/intersection_galloping_avx2.hpp>
//#include <util/intersection_branchless.hpp>
//#include <util/intersection_scalar.hpp>
//#include <util/intersection_v1_avx2.hpp>
//#include <algorithm>


//#include <iostream>
//inline void print_array(const char* name, const uint32_t* array, size_t length) {
//    std::cout << "    " << name << ":";
//    for (size_t i = 0; i < length; ++i)
//        std::cout << " " << array[i];
//    std::cout << std::endl;
//}

/**
 * Find the intersection of two sets of integers, where each set is given by a range of an array of unique integers in
 * ascending order.
 *
 * @param array1            (IN) the first set, address of an existing array of unique integers in ascending order
 * @param begin1            (IN) the starting index of the subarray of array1
 * @param end1              (IN) one past the end of the subarray of array1
 * @param array2            (IN) the second set, address of an existing array of unique integers in ascending order
 * @param begin2            (IN) the starting index of the subarray of array2
 * @param end2              (IN) one past the end of the subarray of array2
 * @param set              (OUT) the set intersection of the two subarrays, must be long enough to hold the values
 *                               of the intersection
 * @return length          (OUT) the length of the set intersection
 */
inline size_t intersect_sets(const uint32_t * array1, size_t begin1, size_t end1,
                             const uint32_t * array2, size_t begin2, size_t end2,
                             uint32_t * set) {

    // Some alternative implementations were tested for performance below.
    // This the fastest implementation for smaller subarrays (typically of size < 60).
    // And the order of if/else is optimal for when (end2 - begin2) < (end1 - begin1).
    //
    // element_to_element process time (100cm_S.g): 0.497132s laptop
    const uint32_t *const initout(set);
    while ( begin1 < end1 && begin2 < end2 )
        if ( array2[begin2] < array1[begin1] )
            begin2++;
        else if ( array1[begin1] < array2[begin2])
            begin1++;
        else {
            *set++ = array1[ begin1++ ];
            begin2++;
        }
    return (set - initout);


    // element_to_element process time (100cm_S.g): 0.514071s laptop, 0.443495s  desktop
    //return intersect_sets_scalar(array1, begin1, end1, array2, begin2, end2, set);

    // element_to_element process time (100cm_S.g): 0.498139s laptop,
    //return intersect_sets_scalar_alt(array1, begin1, end1, array2, begin2, end2, set);


    // element_to_element process time (100cm_S.g):  laptop,
    // return intersect_sets_scalar_goto(array1, begin1, end1, array2, begin2, end2, set);

    // element_to_element process time (100cm_S.g): 0.511036s laptop,
    //return scalar(array1+begin1, end1 - begin1, array2+begin2, end2 - begin2, set);

    // stl (100cm_S.g):  laptop, 0.495211s  desktop
    //uint32_t * end_result = std::set_intersection(array1+begin1, array1+end1, array2+begin2, array2+end2, set);
    //return end_result - set;

    //
    // branchless
    // element_to_element process time (100cm_S.g): 0.958755s laptop, 0.720802s desktop
    //
    //return intersect_sets_scalar_branchless(array2, begin2, end2, array1, begin1, end1, set);


    // intersect_vector_sse
    // element_to_element process time (100cm_S.g):  desktop
//    const size_t length1 = end1 - begin1;
//    const size_t length2 = end2 - begin2;
//    return intersect_vector_sse(array1 + begin1, length1, array2 + begin2, length2, set);


    //
    // galloping_sse
    // element_to_element process time (100cm_S.g): 0.681228s, 0.540757s desktop
//    const size_t length1 = end1 - begin1;
//    const size_t length2 = end2 - begin2;
//    if (length1 <= length2)
//        return v1(array1 + begin1, length1, array2 + begin2, length2, set);
//    else
//        return v1(array2 + begin2, length2, array1 + begin1, length1, set);

    // vector_avx
    //return intersect_vector_avx(&array1[begin1], end1 - begin1, &array2[begin2], end2 - begin2, set);


    // v1_avx2
    //element_to_element process time (100cm_S.g):  0m2.746s laptop   , 0.546569s desktop
//    if (length1 <= length2)
//        return v1_avx2(set1, length1, set2, length2, set);
//    else
//    return v1_avx2(array2 + begin2, end2 - begin2, array1 + begin1, end1 - begin1, set);

    //return match_scalar(array2 + begin2, end2 - begin2, array1 + begin1, end1 - begin1, set);

}

#endif
