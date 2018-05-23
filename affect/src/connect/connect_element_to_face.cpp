#include <algorithm>
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <strstream>

using namespace std;


// replace neighbor array with element-to-face connections
// where interior faces are enumerated first, followed by boundary
// faces
int64_t connect_element_to_face(
  int64_t numElement,
  int64_t numFacePerElement,
  int64_t * neighbor)
{
  int64_t numFace = 0, numInternalFace = 0;

  int64_t * m = neighbor;
  int64_t * m_end = NULL;

  // enumerate and assign all element faces in the interior
  for (int64_t elmt = 0; elmt < numElement; ++elmt) {

    m_end = m + numFacePerElement;

    for ( ; m < m_end; ++m) {
    
      if ( elmt > *m && *m > 0 ) {

        // we are on the higher numbered element on an internal face

        // find the smaller numbered element who has this element as a neighbor
        int64_t * n = neighbor + *m * numFacePerElement;
        const int64_t * n_end = n + numFacePerElement;
        while ( n++ < n_end)
          ; // continue
        if (n == n_end) {
          strstream msg;
          msg << "connect_element_to_face: Unable to find matching face "
              << "of neighbor element " << *m << " for element " << elmt << ".";
          throw runtime_error(msg.str());
        }

        *m = numFace++;
        *n = numFace;
      }
    }
  }
  
  numInternalFace = numFace;

  //
  // now enumerate external, boundary, faces
  //
  m_end = neighbor + numElement * numFacePerElement;

  for ( ; neighbor < m_end; ++neighbor)
    if ( *neighbor == -1 )
      *neighbor = numFace++;

  return numInternalFace; // this is also the starting ID of boundary faces
}


