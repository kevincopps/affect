#include <cassert>
#include <omp.h>
#include <stdio.h>
#include <algorithm>

#include "connect_util.hpp"
#include "intersect_sets.hpp"
#include "element_topology.hpp"
#include "neighbor_util.hpp"

using namespace std;

#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

#define SET_ELEM_VARS()                                             \
  nbr0 = -1, nbr1 = -1, nbr2 = -1, nbr3 = -1, nbr4 = -1, nbr5 = -1; \
  localVertex = &elementToVertex[elmt * HEX8_num_vertex];           \
  vertex0 = *localVertex++, vertex1 = *localVertex++,               \
  vertex2 = *localVertex++, vertex3 = *localVertex++,               \
  vertex4 = *localVertex++, vertex5 = *localVertex++,               \
  vertex6 = *localVertex++, vertex7 = *localVertex


uint32_t neighbor_hex(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet, // working space length 2 * max_elements_per_vertex
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElementBegin,  
  const uint32_t * vertexToElement,
  int64_t * neighbor)
{
  int64_t nbr0, nbr1, nbr2, nbr3, nbr4, nbr5;
  size_t localFaces, idx, ndx, m, n, elmt;
  uint32_t * nbrSet = &elemSet[maxElementsPerVertex]; // second half of partition of working space
  const uint32_t * localVertex;
  const uint32_t * HEX8_num_vertex_per_face = vertex_per_face[HEX8];
  const uint32_t * HEX8_face_vertex_order = face_vertex_order[HEX8];
  const uint32_t HEX8_max_vertex_per_face = max_vertex_per_face[HEX8];
  const uint32_t HEX8_num_face = num_face[HEX8];
  const uint32_t HEX8_num_vertex = num_vertex[HEX8];

  uint32_t vertex0, vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7;
  uint32_t nbrFace0 = 0, nbrFace1 = 0, nbrFace2 = 0, nbrFace3 = 0, nbrFace4 = 0, nbrFace5 = 0;
  uint32_t numBoundaryFaces = 0;

  bool doFace0, doFace1, doFace2, doFace3, doFace4, doFace5;

  std::fill(&neighbor[0], &neighbor[numElement * HEX8_num_face], -2);

  for (elmt = 0; elmt < numElement; elmt++) {

    //printf("neighbor_hex %d.\n", elmt);

    localFaces = elmt * HEX8_num_face;

    doFace0 = -2 == neighbor[ localFaces+0 ];
    doFace1 = -2 == neighbor[ localFaces+1 ];
    doFace2 = -2 == neighbor[ localFaces+2 ];
    doFace3 = -2 == neighbor[ localFaces+3 ];
    doFace4 = -2 == neighbor[ localFaces+4 ];
    doFace5 = -2 == neighbor[ localFaces+5 ];

  //unsigned char doEdge0  = doFace0 && doFace4,
  //     doEdge5  = doFace1 && doFace5,
  //     doEdge11 = doFace2 && doFace3;
  ////if ((doEdge0 && doEdge5) || (doEdge0 && doEdge11) || (doEdge5 && doEdge11)) {
  //if (doEdge0 + doEdge5 + doEdge11 > 1) {

      //SET_ELEM_VARS();
      //DO_FACE_PAIR_WITH_VERTEX(0, 4, 0, 1, 5, 3, 1, 0, 4, 0, 1, 2);
      //DO_FACE_PAIR_WITH_VERTEX(1, 5, 5, 6, 1, 7, 1, 5, 6, 4, 7, 6);
      //DO_FACE_PAIR_WITH_VERTEX(2, 3, 3, 7, 6, 0, 3, 2, 6, 0, 3, 7);
      //// 1.692

  //}
  //else {

  //  unsigned char doEdge1  = doFace1 && doFace4,
  //       doEdge4  = doFace0 && doFace5;
  //  //if ((doEdge1 && doEdge4) || (doEdge1 && doEdge11) || (doEdge4 && doEdge11)) {
  //  if (doEdge1 + doEdge4 + doEdge11 > 1) {

      //SET_ELEM_VARS();
      //DO_FACE_PAIR_WITH_VERTEX(1, 4, 1, 2, 5, 0, 2, 1, 5, 0, 1, 2);
      //DO_FACE_PAIR_WITH_VERTEX(0, 5, 4, 5, 1, 7, 1, 0, 4, 4, 7, 6);
      //DO_FACE_PAIR_WITH_VERTEX(2, 3, 3, 7, 6, 0, 3, 2, 6, 0, 3, 7);
      // 1.603, 24.859 (big-32a.g/PIII)

  //  } 
  //  else {
  //  
    //  unsigned char doEdge6  = doFace2 && doFace5,
    //       doEdge8  = doFace0 && doFace3;
    //  if (doFace1 && doFace4 || doFace2 && doFace5 || doFace0 && doFace3) {
  //    if (doEdge1 + doEdge6 + doEdge8 > 1) {


        nbr0 = -1, nbr1 = -1, nbr2 = -1, nbr3 = -1, nbr4 = -1, nbr5 = -1;
        localVertex = &elementToVertex[elmt * HEX8_num_vertex];
        vertex0 = *localVertex++, vertex1 = *localVertex++,
        vertex2 = *localVertex++, vertex3 = *localVertex++,
        vertex4 = *localVertex++, vertex5 = *localVertex++,
        vertex6 = *localVertex++, vertex7 = *localVertex;

        //SET_ELEM_VARS();

        DO_FACE_PAIR(HEX8, 1, 4, 1, 2, 5, 0, 2, 1, 5, 0, 1, 2);
        DO_FACE_PAIR(HEX8, 2, 5, 6, 7, 2, 4, 3, 2, 6, 4, 7, 6);
        DO_FACE_PAIR(HEX8, 0, 3, 0, 4, 1, 3, 1, 0, 4, 0, 3, 7);
        // 1.593 (tube.g/laptop), 24.469 (big-32a.g/PIII)

  //    }
  //    else {

  //      unsigned char doEdge2  = doFace2 && doFace4,
  //           doEdge5  = doFace1 && doFace5;
  //      //if ((doEdge2 && doEdge5) || (doEdge2 && doEdge8) || (doEdge5 && doEdge8)) {
  //      if (doEdge2 + doEdge5 + doEdge8 > 1) {

            //SET_ELEM_VARS();
            //DO_FACE_PAIR_WITH_VERTEX(2, 4, 2, 3, 6, 1, 3, 2, 6, 0, 1, 2);
            //DO_FACE_PAIR_WITH_VERTEX(1, 5, 5, 6, 1, 4, 2, 1, 5, 4, 7, 6);
            //DO_FACE_PAIR_WITH_VERTEX(0, 3, 0, 4, 1, 3, 1, 0, 4, 0, 3, 7);
            //// 1.673

  //      }
  //      else {

  //        unsigned char doEdge7  = doFace3 && doFace5,
  //            doEdge10 = doFace1 && doFace2;
  //        if (doEdge0 + doEdge7 + doEdge10 > 1) {

            //SET_ELEM_VARS();
            //DO_FACE_PAIR_WITH_VERTEX(0, 4, 0, 1, 5, 3, 1, 0, 4, 0, 1, 2);
            //DO_FACE_PAIR_WITH_VERTEX(3, 5, 4, 7, 0, 5, 0, 3, 7, 4, 7, 6);
            //DO_FACE_PAIR_WITH_VERTEX(1, 2, 2, 6, 1, 3, 2, 1, 5, 3, 2, 6);
            //// 1.723

  //        }
  //        else {

  //          unsigned char doEdge9  = doFace0 && doFace1;
  //          if (doEdge2 + doEdge7 + doEdge9 > 1) {

              //SET_ELEM_VARS();
              //DO_FACE_PAIR_WITH_VERTEX(2, 4, 2, 3, 6, 1, 3, 2, 6, 0, 1, 2);
              //DO_FACE_PAIR_WITH_VERTEX(3, 5, 4, 7, 0, 5, 0, 3, 7, 4, 7, 6);
              //DO_FACE_PAIR_WITH_VERTEX(0, 1, 1, 5, 0, 2, 1, 0, 4, 2, 1, 5);
              //// 1.653 (tube.g/laptop), 26.281 (big-32a.g/PIII)

  //          }
  //          else {

  //            unsigned char doEdge3  = doFace3 && doFace4;
  //            if (doEdge3 + doEdge4 + doEdge10 > 1) {

                //SET_ELEM_VARS();
                //DO_FACE_PAIR_WITH_VERTEX(3, 4, 0, 3, 4, 1, 0, 3, 7, 0, 1, 2);
                //DO_FACE_PAIR_WITH_VERTEX(0, 5, 4, 5, 1, 7, 1, 0, 4, 4, 7, 6);
                //DO_FACE_PAIR_WITH_VERTEX(1, 2, 2, 6, 1, 3, 2, 1, 5, 3, 2, 6);
                //// 1.683

  //            }
  //            else {
  //            
  //              //  (doEdge3 + doEdge6 + doEdge9 > 1)
                //SET_ELEM_VARS();
                //DO_FACE_PAIR_WITH_VERTEX(3, 4, 0, 3, 4, 1, 0, 3, 7, 0, 1, 2);
                //DO_FACE_PAIR_WITH_VERTEX(2, 5, 6, 7, 2, 4, 3, 2, 6, 4, 7, 6);
                //DO_FACE_PAIR_WITH_VERTEX(0, 1, 1, 5, 0, 2, 1, 0, 4, 2, 1, 5);
                //// 1.733

  //            }
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  FILL_NEIGHBOR(HEX8, 0);
  FILL_NEIGHBOR(HEX8, 1);
  FILL_NEIGHBOR(HEX8, 2);
  FILL_NEIGHBOR(HEX8, 3);
  FILL_NEIGHBOR(HEX8, 4);
  FILL_NEIGHBOR(HEX8, 5);

  } // loop over elements

  //printf("neighbor_hex: %lu elements, %d numBoundaryFaces\n", elmt, numBoundaryFaces);

  return numBoundaryFaces;

}

#pragma GCC diagnostic pop