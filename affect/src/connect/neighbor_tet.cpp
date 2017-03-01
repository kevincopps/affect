#include <algorithm>
#include <cassert>

#include "connect_util.hpp"
#include "intersect_sets.hpp"
#include "element_topology.hpp"
#include "neighbor_util.hpp"

using namespace std;

#define SET_ELEM_VARS()                                             \
  nbr0 = -1, nbr1 = -1, nbr2 = -1, nbr3 = -1;                       \
  const int64_t * localVertex = elementToVertex + elmt*TET4_num_vertex; \
  vertex0 = *localVertex++,                                         \
  vertex1 = *localVertex++,                                         \
  vertex2 = *localVertex++,                                         \
  vertex3 = *localVertex


int64_t neighbor_tet(
  int64_t numElement,
  int64_t * elemSet,
  const int64_t * elementToVertex,
  const int64_t * vertexToElementBegin,  
  const int64_t * vertexToElement,
  int64_t * neighbor)
{
    const int64_t TET4_max_vertex_per_face = max_vertex_per_face[TET4];
    const int64_t TET4_num_face = num_face[TET4];
    const int64_t TET4_num_vertex = num_vertex[TET4];
    const int64_t * TET4_num_vertex_per_face = vertex_per_face[TET4];
    const int64_t * TET4_face_vertex_order = face_vertex_order[TET4];

  std::fill(&neighbor[0], &neighbor[numElement * TET4_num_face], -2);

  int64_t numBoundaryFaces = 0;
  int64_t nbr0, nbr1, nbr2, nbr3; 
  int64_t vertex0, vertex1, vertex2, vertex3;
  int64_t nbrFace0, nbrFace1, nbrFace2, nbrFace3;
  int64_t idx, ndx, m, n, nbrSet[2];

  for (int64_t elmt = 0; elmt < numElement; elmt++) {

    int64_t localFaces = elmt * TET4_num_face;
    
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

  return numBoundaryFaces;
}

