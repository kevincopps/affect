#include <algorithm>
#include <cassert>
#include <iostream>

#include "connect_util.hpp"
#include "intersect_sets.hpp"
#include "element_topology.hpp"
#include "neighbor_util.hpp"

using namespace std;

#define SET_ELEM_VARS()                                   \
  nbr0 = -1, nbr1 = -1, nbr2 = -1, nbr3 = -1;             \
  localVertex = &elementToVertex[elmt * TET4_num_vertex]; \
  vertex0 = *localVertex++,                               \
  vertex1 = *localVertex++,                               \
  vertex2 = *localVertex++,                               \
  vertex3 = *localVertex


uint32_t neighbor_tet(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet, // working space length 2 * max_elements_per_vertex
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElementBegin,  
  const uint32_t * vertexToElement,
  int64_t * neighbor)
{
  int64_t nbr0, nbr1, nbr2, nbr3;
  size_t localFaces, idx, ndx, m, n, elmt;
  uint32_t * nbrSet = &elemSet[maxElementsPerVertex]; // second half of partition of working space
  const uint32_t * localVertex;
  const uint32_t * TET4_num_vertex_per_face = vertex_per_face[TET4];
  const uint32_t * TET4_face_vertex_order = face_vertex_order[TET4];
  const uint32_t TET4_max_vertex_per_face = max_vertex_per_face[TET4];
  const uint32_t TET4_num_face = num_face[TET4];
  const uint32_t TET4_num_vertex = num_vertex[TET4];

  uint32_t vertex0, vertex1, vertex2, vertex3;
  uint32_t nbrFace0, nbrFace1, nbrFace2, nbrFace3;
  uint32_t numBoundaryFaces = 0;

  std::fill(&neighbor[0], &neighbor[numElement * TET4_num_face], -2);

  for (elmt = 0; elmt < numElement; elmt++) {

    localFaces = elmt * TET4_num_face;
    
    bool doFace0 = -2 == neighbor[ localFaces+0 ],
         doFace1 = -2 == neighbor[ localFaces+1 ],
         doFace2 = -2 == neighbor[ localFaces+2 ],
         doFace3 = -2 == neighbor[ localFaces+3 ];

    // edge1 and 3 pairing
    if ((doFace0 && doFace2) || (doFace1 && doFace3)) {
      SET_ELEM_VARS();
      DO_FACE_PAIR(TET4, 0, 2, 0, 3, 1, 2, 0, 3, 1, 0, 2, 3);
      DO_FACE_PAIR(TET4, 1, 3, 1, 2, 3, 0, 1, 3, 2, 0, 1, 2);
    }
    // edge 2 and 4 pairing
    else if ((doFace0 && doFace1) || (doFace2 && doFace3)) {
      SET_ELEM_VARS();
      DO_FACE_PAIR(TET4, 0, 1, 1, 3, 0, 2, 0, 3, 1, 1, 3, 2);
      DO_FACE_PAIR(TET4, 2, 3, 0, 2, 3, 1, 0, 2, 3, 0, 1, 2);
    }
    // edge 0 and 5 pairing
    else {
      SET_ELEM_VARS();
      DO_FACE_PAIR(TET4, 0, 3, 0, 1, 3, 2, 0, 3, 1, 0, 1, 2);
      DO_FACE_PAIR(TET4, 1, 2, 2, 3, 1, 0, 1, 3, 2, 0, 2, 3);
    }

    FILL_NEIGHBOR(TET4, 0);
    FILL_NEIGHBOR(TET4, 1);
    FILL_NEIGHBOR(TET4, 2);
    FILL_NEIGHBOR(TET4, 3);

  } // loop over elements

  // std::cerr << "neighbor_tet: numBoundaryFaces = " << numBoundaryFaces << std::endl;

  return numBoundaryFaces;
}

