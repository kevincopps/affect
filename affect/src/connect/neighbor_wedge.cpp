#include <cassert>
#include <vector>

#include "connect_util.hpp"
#include "intersect_sets.hpp"
#include "element_topology.hpp"
#include "neighbor_util.hpp"

using namespace std;


void neighbor_wedge(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet, // working space length 2 * max_elements_per_vertex
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElementBegin,
  const uint32_t * vertexToElement,
  int64_t * neighbor,
  uint32_t * numQuadFaces,
  uint32_t * numTriFaces)
{
  int64_t nbr0, nbr1, nbr2, nbr3, nbr4;
  size_t localFaces, idx, ndx, m, n, elmt;
  uint32_t * nbrSet = &elemSet[maxElementsPerVertex]; // second half of partition of working space
  const uint32_t * localVertex;
  const uint32_t * WEDGE6_num_vertex_per_face = vertex_per_face[WEDGE6];
  const uint32_t * WEDGE6_face_vertex_order = face_vertex_order[WEDGE6];
  const uint32_t WEDGE6_max_vertex_per_face = max_vertex_per_face[WEDGE6];
  const uint32_t WEDGE6_num_face = num_face[WEDGE6];
  const uint32_t WEDGE6_num_vertex = num_vertex[WEDGE6];

  uint32_t vertex0, vertex1, vertex2, vertex3, vertex4, vertex5;
  uint32_t nbrFace0, nbrFace1, nbrFace2, nbrFace3, nbrFace4;
  *numQuadFaces = 0;
  *numTriFaces = 0;

  std::fill(&neighbor[0], &neighbor[numElement * WEDGE6_num_face], -2);

  // cerr << "neighbor_wedge:" << endl;

  for (elmt = 0; elmt < numElement; elmt++) {

    localFaces = elmt * WEDGE6_num_face;
    
    bool doFace0 = -2 == neighbor[ localFaces+0 ],
         doFace1 = -2 == neighbor[ localFaces+1 ],
         doFace2 = -2 == neighbor[ localFaces+2 ],
         doFace3 = -2 == neighbor[ localFaces+3 ],
         doFace4 = -2 == neighbor[ localFaces+4 ];

    // some extra indicators for the three quad faces
    bool didFace0 = ! doFace0,
         didFace1 = ! doFace1,
         didFace2 = ! doFace2;

    nbr0 = -1, nbr1 = -1, nbr2 = -1, nbr3 = -1, nbr4 = -1;
    localVertex = &elementToVertex[elmt * WEDGE6_num_vertex];
    vertex0 = *localVertex++;
    vertex1 = *localVertex++;
    vertex2 = *localVertex++;
    vertex3 = *localVertex++;
    vertex4 = *localVertex++;
    vertex5 = *localVertex;

    // cerr << "    elmt = " << elmt << endl;
    // cerr << "    doFaceX = " << doFace0 << ", "
    //                          << doFace1 << ", "
    //                          << doFace2 << ", "
    //                          << doFace3 << ", "
    //                          << doFace4 << endl;
    // cerr << "    didFaceX = " << didFace0 << ", "
    //                           << didFace1 << ", "
    //                           << didFace2 << endl;
    // cerr << "    vertexX = " << vertex0 << ", "
    //                          << vertex1 << ", "
    //                          << vertex2 << ", "
    //                          << vertex3 << ", "
    //                          << vertex4 << ", "
    //                          << vertex5 << endl;

    // take care of the first triangular face
    if (doFace3) {

      // cerr << "        inside doFace3 for elmt " << elmt << endl;

      // edge 0
      //M_INTERSECT_VERTEX(I, J)
      m = intersect_sets(vertexToElement,
                        vertexToElementBegin[vertex0],
                        vertexToElementBegin[vertex0+1],
                        vertexToElement,
                        vertexToElementBegin[vertex1],
                        vertexToElementBegin[vertex1+1],
                        elemSet);

      // look for neighbor
      n = intersect_sets(elemSet, 0, m,
                        vertexToElement,
                        vertexToElementBegin[vertex2],
                        vertexToElementBegin[vertex2+1],
                        nbrSet);

      if (n > 1) { // there is a neighbor
        nbr3 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
        nbrFace3 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                    &elementToVertex[nbr3 * WEDGE6_num_vertex],
                    vertex0,vertex1,vertex2);
        assert( nbrFace3 < WEDGE6_num_face);
      }

      if (doFace0) {

        didFace0 = true; // remember that we did this already

        // cerr << "        inside doFace0 (didFace0) for elmt " << elmt << endl;

        // look for neighbor
        n = intersect_sets(elemSet, 0, m,
                          vertexToElement,
                          vertexToElementBegin[vertex4],
                          vertexToElementBegin[vertex4+1],
                          nbrSet);

        if (n > 1) { // there is a neighbor
          nbr0 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
          nbrFace0 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                      &elementToVertex[nbr0 * WEDGE6_num_vertex],
                      vertex1,vertex0,vertex3);
          assert( nbrFace0 < WEDGE6_num_face);
        }
      }
    }

    // take care of the second triangular face
    if (doFace4) {

      // cerr << "        inside doFace4 for elmt " << elmt << endl;

      // edge 5
      m = intersect_sets(vertexToElement,
                        vertexToElementBegin[vertex3],
                        vertexToElementBegin[vertex3+1],
                        vertexToElement,
                        vertexToElementBegin[vertex5],
                        vertexToElementBegin[vertex5+1],
                        elemSet);

      // look for neighbor
      n = intersect_sets(elemSet, 0, m,
                        vertexToElement,
                        vertexToElementBegin[vertex4],
                        vertexToElementBegin[vertex4+1],
                        nbrSet);
      
      if (n > 1) { // there is a neighbor
        nbr4 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
        nbrFace4 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                    &elementToVertex[nbr4 * WEDGE6_num_vertex],
                    vertex3,vertex4,vertex5);
        assert(nbrFace4 < WEDGE6_num_face);
      }

      if (doFace2) {

        didFace2 = true; // remember that we did this already

        // cerr << "        inside doFace2 (didFace2) for elmt " << elmt << endl;

        // look for neighbor
        n = intersect_sets(elemSet, 0, m,
                          vertexToElement,
                          vertexToElementBegin[vertex0],
                          vertexToElementBegin[vertex0+1],
                          nbrSet);

        if (n > 1) { // there is a neighbor
          nbr2 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
          nbrFace2 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                      &elementToVertex[nbr2 * WEDGE6_num_vertex],
                      vertex0,vertex2,vertex5);
          assert(nbrFace2 < WEDGE6_num_face);
        }
      }
    }

    // by this point, we have taken care of faces tri 3, 4, and possibly quad 0, 2
    // do quad 0, 2 here if necessary, and if so, quad 1 as side effect

    if (!didFace0) {

      // cerr << "        inside !didFace0 for elmt " << elmt << endl;

      // edge 7
      m = intersect_sets(vertexToElement,
                        vertexToElementBegin[vertex1],
                        vertexToElementBegin[vertex1+1],
                        vertexToElement,
                        vertexToElementBegin[vertex4],
                        vertexToElementBegin[vertex4+1],
                        elemSet);

      // look for neighbor
      n = intersect_sets(elemSet, 0, m,
                        vertexToElement,
                        vertexToElementBegin[vertex0],
                        vertexToElementBegin[vertex0+1],
                        nbrSet);

      if (n > 1) { // there is a neighbor
        nbr0 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
        nbrFace0 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                    &elementToVertex[nbr0 * WEDGE6_num_vertex],
                    vertex1,vertex0,vertex3);
        assert(nbrFace0 < WEDGE6_num_face);
      }

      if (doFace1) {

        didFace1 = true; // remember we did this

        // cerr << "        inside doFace1 (didFace1) for elmt " << elmt << endl;

        // cerr << "        inside doFace1 for elmt " << elmt << endl;
        // cerr << "        elemSet ";
        // for (uint32_t k = 0; k < m-1; k++)
        //     cerr << elemSet[k] << ", ";
        // cerr << elemSet[m-1] << endl;

        // look for neighbor
        n = intersect_sets(elemSet, 0, m,
                          vertexToElement,
                          vertexToElementBegin[vertex2],
                          vertexToElementBegin[vertex2+1],
                          nbrSet);
        if (n > 1) { // there is a neighbor
          nbr1 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
          nbrFace1 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                      &elementToVertex[nbr1 * WEDGE6_num_vertex],
                      vertex2,vertex1,vertex4);
          assert(nbrFace1 < WEDGE6_num_face);
        }
      }
    }
    
    if (!didFace2) {

      // cerr << "        inside doFace2 for elmt " << elmt << endl;

      // edge 8
      m = intersect_sets(vertexToElement,
                        vertexToElementBegin[vertex2],
                        vertexToElementBegin[vertex2+1],
                        vertexToElement,
                        vertexToElementBegin[vertex5],
                        vertexToElementBegin[vertex5+1],
                        elemSet);

      // look for neighbor
      n = intersect_sets(elemSet, 0, m,
                        vertexToElement,
                        vertexToElementBegin[vertex0],
                        vertexToElementBegin[vertex0+1],
                        nbrSet);

      if (n > 1) { // there is a neighbor
        nbr2 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
        nbrFace2 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                    &elementToVertex[nbr2 * WEDGE6_num_vertex],
                    vertex0,vertex2,vertex5);
        assert(nbrFace2 < WEDGE6_num_face);
      }

      if (!didFace1) {

        // cerr << "        inside !doFace1 for elmt " << elmt << endl;

        // look for neighbor
        n = intersect_sets(elemSet, 0, m,
                          vertexToElement,
                          vertexToElementBegin[vertex1],
                          vertexToElementBegin[vertex1+1],
                          nbrSet);

        if (n > 1) { // there is a neighbor
          nbr1 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
          nbrFace1 = which_face_irregular(
                    WEDGE6_face_vertex_order,
                    WEDGE6_num_face,
                    WEDGE6_max_vertex_per_face,
                    WEDGE6_num_vertex_per_face,
                      &elementToVertex[nbr1 * WEDGE6_num_vertex],
                      vertex2,vertex1,vertex4);
          assert(nbrFace1 < WEDGE6_num_face);
        }
      }
    }

    // by this point we have taken care of faces tri 3, 4, and quad 0, 2, and possibly quad 1
    // do quad 1 if we have not already

    if (!didFace1) {

      // cerr << "        inside (!didFace1) for elmt " << elmt << endl;

      // edge 7
      m = intersect_sets(vertexToElement,
                        vertexToElementBegin[vertex1],
                        vertexToElementBegin[vertex1+1],
                        vertexToElement,
                        vertexToElementBegin[vertex4],
                        vertexToElementBegin[vertex4+1],
                        elemSet);

      // look for neighbor
      n = intersect_sets(elemSet, 0, m,
                        vertexToElement,
                        vertexToElementBegin[vertex2],
                        vertexToElementBegin[vertex2+1],
                        nbrSet);
      if (n > 1) { // there is a neighbor
        nbr1 = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];
        nbrFace1 = which_face_irregular(
                  WEDGE6_face_vertex_order,
                  WEDGE6_num_face,
                  WEDGE6_max_vertex_per_face,
                  WEDGE6_num_vertex_per_face,
                  &elementToVertex[nbr1 * WEDGE6_num_vertex],
                  vertex2,vertex1,vertex4);
        assert(nbrFace1 < WEDGE6_num_face);
      }
    }

    FILL_NEIGHBOR_INCREMENT(WEDGE6, *numQuadFaces, 0);
    FILL_NEIGHBOR_INCREMENT(WEDGE6, *numQuadFaces, 1);
    FILL_NEIGHBOR_INCREMENT(WEDGE6, *numQuadFaces, 2);
    FILL_NEIGHBOR_INCREMENT(WEDGE6, *numTriFaces, 3);
    FILL_NEIGHBOR_INCREMENT(WEDGE6, *numTriFaces, 4);

  } // loop over elements

}

