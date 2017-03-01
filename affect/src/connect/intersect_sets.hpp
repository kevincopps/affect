#ifndef AFFECT_INTERSECT_SETS_H
#define AFFECT_INTERSECT_SETS_H


/**
  * find intersection of two ordered sets of integers, these sets lie in
  * sections of the given array. Return value is size of resulting set.
  */
inline int64_t intersect_sets( const int64_t * array1, int64_t begin1, int64_t end1,
                           const int64_t * array2, int64_t begin2, int64_t end2, 
                           int64_t * set) {
  int64_t n = 0;
  while ( begin1 < end1 && begin2 < end2 )
    if ( array1[begin1] < array2[begin2])
      begin1++;
    else if ( array2[begin2] < array1[begin1] )
      begin2++;
    else {
      set[n++] = array1[ begin1++ ];
      begin2++;
    }
  return n;
}


#endif
