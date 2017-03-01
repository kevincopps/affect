#ifndef AFFECT_NEIGHBOR_UTIL_HPP
#define AFFECT_NEIGHBOR_UTIL_HPP

#define M_INTERSECT_VERTEX(I, J)                         \
  m = intersect_sets(vertexToElement,                    \
                    vertexToElementBegin[vertex##I],     \
                    vertexToElementBegin[vertex##I + 1], \
                    vertexToElement,                     \
                    vertexToElementBegin[vertex##J],     \
                    vertexToElementBegin[vertex##J + 1], \
                    elemSet)

#define N_INTERSECT_VERTEX(I)                            \
  n = intersect_sets(elemSet, 0, m,                      \
                    vertexToElement,                     \
                    vertexToElementBegin[vertex##I],     \
                    vertexToElementBegin[vertex##I + 1], \
                    nbrSet)

#define FIND_NEIGHBOR_FACE(E, N, I, J, K)                   \
  if (n > 1) {                                              \
    nbr##N = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];     \
    nbrFace##N = which_face(                                \
                E##_face_vertex_order,                      \
                E##_num_face,                               \
                E##_max_vertex_per_face,                    \
                &elementToVertex[nbr##N * E##_num_vertex],  \
                vertex##I,vertex##J,vertex##K);             \
    assert( nbrFace##N < E##_num_face );                    \
  }

#define FIND_NEIGHBOR_FACE_IRREGULAR(E, N, I, J, K)         \
  if (n > 1) {                                              \
    nbr##N = nbrSet[0] != elmt ? nbrSet[0] : nbrSet[1];     \
    nbrFace##N = which_face(                                \
                E##_face_vertex_order,                      \
                E##_num_face,                               \
                E##_max_vertex_per_face,                    \
                E##_num_vertex_per_face,                    \
                &elementToVertex[nbr##N * E##_num_vertex],  \
                vertex##I,vertex##J,vertex##K);             \
    assert( nbrFace##N < E##_num_face );                    \
  }

#define DO_FACE_PAIR(E, F1, F2, M1, M2, N1, N2, I1, J1, K1, I2, J2, K2) \
  if (doFace##F1 || doFace##F2) {            \
    M_INTERSECT_VERTEX(M1, M2);              \
    if (doFace##F1) {                        \
      N_INTERSECT_VERTEX( N1 );              \
      FIND_NEIGHBOR_FACE(E, F1, I1, J1, K1); \
    }                                        \
    if (doFace##F2) {                        \
      N_INTERSECT_VERTEX( N2 );              \
      FIND_NEIGHBOR_FACE(E, F2, I2, J2, K2); \
    }                                        \
  }

#define FILL_NEIGHBOR_INCREMENT(E, N, I)         \
  if (doFace##I) {                               \
    idx = localFaces + I;                        \
    if (nbr##I > -1) {                           \
      ndx = nbr##I * E##_num_face + nbrFace##I;  \
      neighbor[idx] = nbr##I ;                   \
      neighbor[ndx] = elmt;                      \
	  }                                          \
    else {                                       \
      ++N;                                       \
      neighbor[idx] = -1;                        \
    }                                            \
  }

#define FILL_NEIGHBOR(E, I)                      \
  FILL_NEIGHBOR_INCREMENT(E, numBoundaryFaces, I)

#endif
