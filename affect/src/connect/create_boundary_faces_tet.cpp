#include "element_topology.hpp"


void create_boundary_faces_tet(
  uint32_t numElement,
  const uint32_t* elementToVertex,
  const int64_t* neighbor,
  uint32_t* boundaryFaceToVertex) {

  const uint32_t TET4_num_vertex = num_vertex[TET4];
  const uint32_t* localVertexLast = elementToVertex + numElement * TET4_num_vertex;

  for ( ; elementToVertex < localVertexLast; elementToVertex += TET4_num_vertex) {

    bool do0 = *neighbor++ < 0;
    bool do1 = *neighbor++ < 0;
    bool do2 = *neighbor++ < 0;
    bool do3 = *neighbor++ < 0;

    if (do0) {
      *boundaryFaceToVertex++ = elementToVertex[0];
      *boundaryFaceToVertex++ = elementToVertex[1];
      *boundaryFaceToVertex++ = elementToVertex[3];
    }
    if (do1) {
      *boundaryFaceToVertex++ = elementToVertex[1];
      *boundaryFaceToVertex++ = elementToVertex[2];
      *boundaryFaceToVertex++ = elementToVertex[3];
    }
    if (do2) {
      *boundaryFaceToVertex++ = elementToVertex[0];
      *boundaryFaceToVertex++ = elementToVertex[3];
      *boundaryFaceToVertex++ = elementToVertex[2];
    }
    if (do3) {
      *boundaryFaceToVertex++ = elementToVertex[0];
      *boundaryFaceToVertex++ = elementToVertex[2];
      *boundaryFaceToVertex++ = elementToVertex[1];
    }
  }
}

