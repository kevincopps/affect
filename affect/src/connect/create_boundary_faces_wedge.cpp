#include "element_topology.hpp"


void create_boundary_faces_wedge(
  int64_t numElement,
  const int64_t* elementToVertex,
  const int64_t* neighbor,
  int64_t* boundaryQuadToVertex,
  int64_t* boundaryTriToVertex) {

  const int64_t WEDGE6_num_vertex = num_vertex[WEDGE6];
  const int64_t* localVertexLast = elementToVertex + numElement * WEDGE6_num_vertex;

  for ( ; elementToVertex < localVertexLast; elementToVertex += WEDGE6_num_vertex) {

    bool do0 = *neighbor++ < 0;
    bool do1 = *neighbor++ < 0;
    bool do2 = *neighbor++ < 0;
    bool do3 = *neighbor++ < 0;
    bool do4 = *neighbor++ < 0;

    // the first three are quadrilateral faces
    if (do0) {
      *boundaryQuadToVertex++ = elementToVertex[0];
      *boundaryQuadToVertex++ = elementToVertex[1];
      *boundaryQuadToVertex++ = elementToVertex[4];
      *boundaryQuadToVertex++ = elementToVertex[3];
    }
    if (do1) {
      *boundaryQuadToVertex++ = elementToVertex[1];
      *boundaryQuadToVertex++ = elementToVertex[2];
      *boundaryQuadToVertex++ = elementToVertex[5];
      *boundaryQuadToVertex++ = elementToVertex[4];
    }
    if (do2) {
      *boundaryQuadToVertex++ = elementToVertex[0];
      *boundaryQuadToVertex++ = elementToVertex[3];
      *boundaryQuadToVertex++ = elementToVertex[5];
      *boundaryQuadToVertex++ = elementToVertex[2];
    }

    // the last two are triangular faces
    if (do3) {
      *boundaryTriToVertex++ = elementToVertex[0];
      *boundaryTriToVertex++ = elementToVertex[2];
      *boundaryTriToVertex++ = elementToVertex[1];
    }
    if (do4) {
      *boundaryTriToVertex++ = elementToVertex[3];
      *boundaryTriToVertex++ = elementToVertex[4];
      *boundaryTriToVertex++ = elementToVertex[5];
    }
  }
}

