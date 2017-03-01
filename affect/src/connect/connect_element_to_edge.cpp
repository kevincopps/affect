#include <algorithm>
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <cstdint>

#include "intersect_sets.hpp"
#include "connect_util.hpp"


void connect_element_to_edge(
    int64_t numElement,
    int64_t numEdgePerElement,
    int64_t numVertexPerElement,
    const int64_t * edgeVertexOrder,
    const int64_t * elementToVertex,
    const int64_t * vertexToElementBegin,
    const int64_t * vertexToElement,
    int64_t * elementToEdge,
    int64_t * numInternalEdge,
    int64_t * numExternalEdge)
{
    std::fill(elementToEdge, elementToEdge + numElement * numEdgePerElement, -2);
    *numInternalEdge = 0;
    *numExternalEdge = 0;
    
    const int64_t MAX_ELEMENT_PER_EDGE = 32;
    int64_t m, elemSet[MAX_ELEMENT_PER_EDGE];
    
    const int64_t * localVertex = elementToVertex;
    int64_t * localEdge = elementToEdge;
    const int64_t * localEdgeVertexOrder;

  for (int64_t elmt = 0; elmt < numElement; ++elmt, localVertex += numVertexPerElement) {

    localEdgeVertexOrder = edgeVertexOrder;

    for (int64_t iEdge = 0; iEdge < numEdgePerElement; 
         ++iEdge, localEdge++, localEdgeVertexOrder += 2) {

      if ( -2 != *localEdge ) continue;

      int64_t vertex0 = localVertex[ localEdgeVertexOrder[0] ];
      int64_t vertex1 = localVertex[ localEdgeVertexOrder[1] ];

      m = intersect_sets(vertexToElement,
                        vertexToElementBegin[vertex0],
                        vertexToElementBegin[vertex0+1],
                        vertexToElement,
                        vertexToElementBegin[vertex1],
                        vertexToElementBegin[vertex1+1],
                        elemSet);

      if (m == 0) {

        // this edge is not shared with any other element
        // (it is on the boundary in 2D, or on a non-manifold face in 3D)
        *localEdge = -1;
        ++*numExternalEdge;

      }
      else if (m > MAX_ELEMENT_PER_EDGE) {

        // we have written past the end of the elemSet array, abort.
        throw std::range_error("MAX_ELEMENT_PER_EDGE exceeded. Contact a developer.");

      }
      else {

        *localEdge = *numInternalEdge; *numInternalEdge += 1;

        for (int64_t iNeighbor = 0; iNeighbor < m; ++iNeighbor) {

          int64_t nbr = elemSet[iNeighbor];

          if (nbr != elmt) {

            int64_t nbrEdge = which_edge(
                            edgeVertexOrder,
                            numEdgePerElement,
                            &elementToVertex[nbr * numVertexPerElement],
                            vertex0,vertex1);
            assert( nbrEdge < numEdgePerElement );
            assert(elementToEdge[nbr * numEdgePerElement + nbrEdge] == -2);
            elementToEdge[nbr * numEdgePerElement + nbrEdge] = *localEdge;
          }

        }

      }


    } // iEdge
  } // elmt
}

