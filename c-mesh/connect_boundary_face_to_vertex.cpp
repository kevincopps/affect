#include <stdexcept>
#include <strstream>
#include <vector>
#include <cstdint>

#ifdef AFFECT_VERBOSE
#include <iomanip>
#include <iostream>
#endif

#include "connect_util.hpp"
#include "element_topology.hpp"
#include "duration.hpp"
#include "connect.hpp"

using namespace std;


void error_boundary_face_to_vertex( 
  char* elementName,
  int64_t numBoundaryQuadFace,
  int64_t numBoundaryTriFace) {
  
  strstream msg;

  msg << "connect_boundary_face_to_vertex: Unable to handle (" <<
            elementName <<
            ") elements with " <<
            numBoundaryQuadFace <<
            " quadrilateral and " <<
            numBoundaryTriFace <<
            " triangular faces." << endl;  
            
  throw invalid_argument(msg.str());
}

void connect_boundary_face_to_vertex(
  char* elementName, 
  int64_t numElement,
  int64_t numBoundaryQuadFace,
  int64_t numBoundaryTriFace,
  const int64_t * elementToVertex,
  const int64_t * neighbor,
  int64_t * boundaryFaceToVertex)
{
  int64_t numFacePerElement = 0, 
      numVertexPerElement = 0,
      numVertexPerFace = 0, 
      numBoundaryFace = 0;

  if ( is_element_name(elementName,aliases[HEX8]) ) {

    if (numBoundaryTriFace != 0 || numBoundaryQuadFace < 0) {  
      error_boundary_face_to_vertex( 
        elementName,
        numBoundaryQuadFace,
        numBoundaryTriFace);
    }

    START_DURATION(3);

    create_boundary_faces_hex(
      numElement, 
      elementToVertex, 
      neighbor, 
      boundaryFaceToVertex);

    STOP_DURATION(3,create_boundary_faces_hex);
  }
  else if ( is_element_name(elementName,aliases[TET4]) ) {

    if (numBoundaryQuadFace != 0  || numBoundaryTriFace < 0) {  
      error_boundary_face_to_vertex( 
        elementName,
        numBoundaryQuadFace,
        numBoundaryTriFace);
    }

    START_DURATION(3);

    create_boundary_faces_tet(
      numElement, 
      elementToVertex, 
      neighbor, 
      boundaryFaceToVertex);

    STOP_DURATION(3,create_boundary_faces_tet);
    
  }
  else if ( is_element_name(elementName,aliases[WEDGE6]) ) {
  
    if (numBoundaryQuadFace == 0  || numBoundaryTriFace == 0) {  
      error_boundary_face_to_vertex( 
        elementName,
        numBoundaryQuadFace,
        numBoundaryTriFace);
    }

    START_DURATION(3);

    create_boundary_faces_wedge(
      numElement,
      elementToVertex,
      &neighbor[0],
      &boundaryFaceToVertex[0],
      &boundaryFaceToVertex[numBoundaryQuadFace*4]);

    STOP_DURATION(3,create_boundary_faces_wedge);
    
  }
  else {
  
    error_boundary_face_to_vertex( 
      elementName,
      numBoundaryQuadFace,
      numBoundaryTriFace);
  }

# ifdef VERBOSE
  int64_t numFacePerElement = 0;
  int64_t numVertexPerFace  = 0;
  if ( is_element_name(elementName,HEX8_aliases) ) {
    numFacePerElement = HEX8_num_face;
    numVertexPerFace  = HEX8_num_vertex_per_face;
  }
  else if ( is_element_name(elementName,TET4_aliases) ) {
    numFacePerElement = TET4_num_face;
    numVertexPerFace  = TET4_num_vertex_per_face;
  }
  else if ( is_element_name(elementName,WEDGE6_aliases) ) {
    numFacePerElement = WEDGE6_num_face;
    numVertexPerFace  = WEDGE6_num_vertex_per_face;
  }
  else
    error_boundary_face_to_vertex( 
      elementName,
      numBoundaryQuadFace,
      numBoundaryTriFace);
       
  for (int64_t iElem = 0; iElem < numElement; ++iElem) {
    int64_t pos = numFacePerElement * iElem;
    for (int64_t iFace = 0; iFace < numFacePerElement; ++iFace) {
      if (neighbor[pos++] == -2) {
        cout << "element " << iElem << " does not have face " << iFace << endl;
        exit(1);
      }
    }
  }

  cout << "numBoundaryFaces = " << numBoundaryFaces << endl;
  for (int64_t iBFace = 0; iBFace < numBoundaryFaces; ++iBFace) {
    int64_t first = numVertexPerFace * iBFace;
    int64_t last = first + numVertexPerFace;
    for (int64_t i = first; i < last; ++i) {
      cout << setw(4) << boundaryFaceToVertex[i];
    }
    cout << endl;
  }
# endif
}

