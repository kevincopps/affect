#include <stdexcept>
#include <strstream>
#include <vector>

#include "connect_util.hpp"
#include "element_topology.hpp"

using namespace std;


void error_element_to_element(
  const char* elementName,
  int64_t numElement,
  int64_t maxElementPerVertex) {
  
  strstream msg;
  msg << "connect_element_to_element: Unable to handle (" <<
          elementName <<
          ") elements with number of elements " <<
          numElement << " and " <<
          maxElementPerVertex <<
          " max elements per vertex." << endl;
  
  throw invalid_argument(msg.str());  
}

void connect_element_to_element(
  const char* elementName,
  int64_t numElement,
  int64_t maxElementPerVertex,
  const int64_t * elementToVertex,
  const int64_t * vertexToElementBegin,
  const int64_t * vertexToElement,
  int64_t * neighbor,
  int64_t * numBoundaryQuadFaces,
  int64_t * numBoundaryTriFaces)
{
  if (numElement < 1) return;
  
  if ( maxElementPerVertex < 1 ) {
    error_element_to_element(
      elementName,
      numElement,
      maxElementPerVertex);
  }

  vector<int64_t> elemSet(maxElementPerVertex);

  if (is_element_name(elementName,aliases[HEX8])) {
    
    *numBoundaryQuadFaces = neighbor_hex(
      numElement, 
      &elemSet[0],
      elementToVertex,
      vertexToElementBegin,      
      vertexToElement,
      neighbor);
    
  }
  else if (is_element_name(elementName,aliases[TET4])) {
    
    *numBoundaryTriFaces = neighbor_tet(
      numElement, 
      &elemSet[0],
      elementToVertex,
      vertexToElementBegin,
      vertexToElement,
      neighbor);
    
  }
  else if (is_element_name(elementName,aliases[WEDGE6])) {
    
    neighbor_wedge(
      numElement, 
      &elemSet[0],
      elementToVertex,
      vertexToElementBegin,
      vertexToElement,
      neighbor,
      numBoundaryQuadFaces,
      numBoundaryTriFaces);
    
  }  
  else {
    error_element_to_element(
      elementName,
      numElement,
      maxElementPerVertex);    
  }
}

