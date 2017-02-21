#include "element_topology.hpp"


void create_boundary_faces_hex(
  int64_t numElement,
  const int64_t* elementToVertex,
  const int64_t* neighbor,
  int64_t* boundaryFaceToVertex) {

  const int64_t HEX8_num_vertex = num_vertex[HEX8];

  const int64_t* localVertexLast = elementToVertex + numElement * HEX8_num_vertex;

  for ( ; elementToVertex < localVertexLast; elementToVertex += HEX8_num_vertex) {

    bool do0 = *neighbor++ < 0;
    bool do1 = *neighbor++ < 0;
    bool do2 = *neighbor++ < 0;
    bool do3 = *neighbor++ < 0;
    bool do4 = *neighbor++ < 0;
    bool do5 = *neighbor++ < 0;

    if (do0) {
      *boundaryFaceToVertex++ = elementToVertex[0];
      *boundaryFaceToVertex++ = elementToVertex[1];
      *boundaryFaceToVertex++ = elementToVertex[5];
      *boundaryFaceToVertex++ = elementToVertex[4];
    }
    if (do1) {
      *boundaryFaceToVertex++ = elementToVertex[1];
      *boundaryFaceToVertex++ = elementToVertex[2];
      *boundaryFaceToVertex++ = elementToVertex[6];
      *boundaryFaceToVertex++ = elementToVertex[5];
    }
    if (do2) {
      *boundaryFaceToVertex++ = elementToVertex[2];
      *boundaryFaceToVertex++ = elementToVertex[3];
      *boundaryFaceToVertex++ = elementToVertex[7];
      *boundaryFaceToVertex++ = elementToVertex[6];
    }
    if (do3) {
      *boundaryFaceToVertex++ = elementToVertex[0];
      *boundaryFaceToVertex++ = elementToVertex[4];
      *boundaryFaceToVertex++ = elementToVertex[7];
      *boundaryFaceToVertex++ = elementToVertex[3];
    }
    if (do4) {
      *boundaryFaceToVertex++ = elementToVertex[0];
      *boundaryFaceToVertex++ = elementToVertex[3];
      *boundaryFaceToVertex++ = elementToVertex[2];
      *boundaryFaceToVertex++ = elementToVertex[1];
    }
    if (do5) {
      *boundaryFaceToVertex++ = elementToVertex[4];
      *boundaryFaceToVertex++ = elementToVertex[5];
      *boundaryFaceToVertex++ = elementToVertex[6];
      *boundaryFaceToVertex++ = elementToVertex[7];
    }
  }
}

