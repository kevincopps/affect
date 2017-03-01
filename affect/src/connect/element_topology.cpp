#include "element_topology.hpp"


const char *names[END_TOPOLOGY] = {
  "hex8","hex20","hex27",
  "pyramid5","pyramid13",
  "quad4","quad8","quad9",
  "quadshell4","quadshell8","quadshell9",
  "tet4","tet8","tet10",
  "tri3","tri4","tri6",
  "trishell3","trishell6",
  "wedge6","wedge15"};

const char *aliases[END_TOPOLOGY] = {
  "hex8|hex",
  "hex20"
  "hex27",
  "pyramid5|pyramid",
  "pyramid13",
  "quad4|quad|quadrilateral4|quadrilateral",
  "quad8|quadrilateral8",
  "quad9|quadrilateral9",
  "shell4|shell|quadshell4|quadshell",
  "shell8|quadshell8",
  "shell9|quadshell9",
  "tetra4|tetra|tet4",
  "tetra8|tet8",
  "tetra10|tet10",
  "triangle3|triangle|tri3|tri",
  "triangle4|tri4",
  "triangle6|tri6",
  "trishell3|trishell|triangleshell3|triangleshell",
  "trishell6|triangleshell6",
  "wedge6|wedge",
  "wedge15"};

// h  h  h  p  p  q  q  q  s  s  s  e  e  e  r  r  r  s  s  w  w

const int64_t num_spatial_dim[END_TOPOLOGY] = {
   3, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 3, 3, 3, 3};

const int64_t num_vertex[END_TOPOLOGY] = {
   8, 8, 8, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 6, 6};

const int64_t num_node[END_TOPOLOGY] = {
   8,20,27, 5,13, 4, 8, 9, 4, 8, 9, 4, 8,10, 3, 4, 6, 3, 6, 6,15};
  
const int64_t num_edge[END_TOPOLOGY] = {
  12,12,12, 8, 8, 4, 4, 4, 4, 4, 4, 6, 6, 6, 3, 3, 3, 3, 3, 9, 9};

const int64_t num_face[END_TOPOLOGY] = {
   6, 6, 6, 5, 5, 1, 1, 1, 2, 2, 2, 4, 4, 4, 1, 1, 1, 2, 2, 5, 5};

const int64_t num_node_per_edge[END_TOPOLOGY] = {
   2, 3, 3, 2, 3, 2, 3, 3, 2, 3, 3, 2, 2, 3, 2, 2, 3, 2, 3, 2, 3};

const int64_t max_vertex_per_face[END_TOPOLOGY] = {
   4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4};

const int64_t max_node_per_face[END_TOPOLOGY] = {
   4, 8, 9, 4, 8, 4, 8, 9, 4, 8, 9, 3, 4, 6, 3, 4, 6, 3, 6, 4, 8};

const int64_t* const max_edge_per_face = max_vertex_per_face;

//---------------------------------------------------------------------
// HEX8
//---------------------------------------------------------------------
const int64_t HEX8_edge_node_order[] = {
  0,1, 1,2, 2,3, 3,0,
  4,5, 5,6, 6,7, 7,4,
  0,4, 1,5, 2,6, 3,7
};

const int64_t HEX8_face_node_order[] = {
    0,1,5,4, 1,2,6,5, 2,3,7,6,
    0,4,7,3, 0,3,2,1, 4,5,6,7
};

const int64_t HEX8_face_edge_order[] = {
  0, 9, 4, 8,   1,10, 5, 9,   2,11, 6,10,
  8, 7,11, 3,   3, 2, 1, 0,   4, 5, 6, 7
};

//---------------------------------------------------------------------
// HEX20
//---------------------------------------------------------------------
const int64_t HEX20_edge_node_order[] = {
  0, 1, 8, 1, 2, 9, 2, 3,10, 3, 0,11,
  4, 5,16, 5, 6,17, 6, 7,18, 7, 4,19,
  0, 4,12, 1, 5,13, 2, 6,14, 3, 7,15
};

const int64_t HEX20_face_node_order[] = {
  0,1,5,4,   8,13,16,12,
  1,2,6,5,   9,14,17,13,
  2,3,7,6,  10,15,18,14,
  0,4,7,3,  12,19,15,11,
  0,3,2,1,  11,10, 9, 8,
  4,5,6,7,  16,17,18,19
};

//---------------------------------------------------------------------
// HEX27
//---------------------------------------------------------------------
const int64_t HEX27_edge_node_order[] = {
  0, 1, 8,  1, 2, 9,  2, 3,10,  3, 0,11,
  4, 5,16,  5, 6,17,  6, 7,18,  7, 4,19,
  0, 4,12,  1, 5,13,  2, 6,14,  3, 7,15
};

const int64_t HEX27_face_node_order[] = {
  0,1,5,4,   8,13,16,12,  25,
  1,2,6,5,   9,14,17,13,  24,
  2,3,7,6,  10,15,18,14,  26,
  0,4,7,3,  12,19,15,11,  23,
  0,3,2,1,  11,10, 9, 8,  21,
  4,5,6,7,  16,17,18,19,  22
};

//---------------------------------------------------------------------
// PYRAMID5
//---------------------------------------------------------------------
const int64_t PYRAMID5_edge_node_order[] = {
  0,1,  1,2,  2,3,  3,0,
  0,4,  1,4,  2,4,  3,4
};

const int64_t PYRAMID5_face_node_order[] = {
  0,1,4,-1,  1,2,4,-1,  2,3,4,-1,  0,4,3,-1,
  0,3,2,1
};

const int64_t PYRAMID5_face_edge_order[] = {
  0,5,4,-1,  1,6,5,-1,  2,7,6,-1,  4,7,3,-1,
  3,2,1,0
};

const int64_t PYRAMID5_node_per_face[] = {
  3, 3, 3, 3, 4
};

//---------------------------------------------------------------------
// PYRAMID13
//---------------------------------------------------------------------
const int64_t PYRAMID13_edge_node_order[] = {
  0,1,5,  1,2,6,   2,3,7,   3,0,8,
  0,4,9,  1,4,10,  2,4,11,  3,4,12,
};

const int64_t PYRAMID13_face_node_order[] = {
   0, 1, 4, 5,10, 9,-1,-1,   1, 2, 4, 6,11,10,-1,-1,
   2, 3, 4, 7,12,11,-1,-1,   0, 4, 3, 9,12, 8,-1,-1,
   0, 3, 2, 1, 8, 7, 6, 5
};

const int64_t PYRAMID13_node_per_face[] = {
  6, 6, 6, 6, 8
};

//---------------------------------------------------------------------
// QUAD4
//---------------------------------------------------------------------
const int64_t QUAD4_edge_node_order[] = {
  0,1, 1,2, 2,3, 3,0
};

const int64_t QUAD4_face_node_order[] = {
  0,1,2,3
};

const int64_t* QUAD4_face_edge_order = QUAD4_face_node_order;

//---------------------------------------------------------------------
// QUAD8
//---------------------------------------------------------------------
const int64_t QUAD8_edge_node_order[] = {
  0,1,4,  1,2,5,  2,3,6,  3,0,7
};

const int64_t QUAD8_face_node_order[] = {
  0,1,2,3,4,5,6,7
};

//---------------------------------------------------------------------
// QUAD9
//---------------------------------------------------------------------
const int64_t QUAD9_edge_node_order[] = {
  0,1,4, 1,2,5, 2,3,6, 3,0,7
};

const int64_t QUAD9_face_node_order[] = {
  0,1,2,3,4,5,6,7,8
};

//---------------------------------------------------------------------
// QUADSHELL4
//---------------------------------------------------------------------
const int64_t* QUADSHELL4_edge_node_order = QUAD4_edge_node_order;

const int64_t QUADSHELL4_face_node_order[] = {
  0,1,2,3, 0,3,2,1
};

const int64_t QUADSHELL4_face_edge_order[] = {
  0,1,2,3, 3,2,1,0
};

//---------------------------------------------------------------------
// QUADSHELL8
//---------------------------------------------------------------------

const int64_t QUADSHELL8_face_node_order[] = {
  0,1,2,3,4,5,6,7,
  0,3,2,1,7,6,5,4
};

//---------------------------------------------------------------------
// QUADSHELL9
//---------------------------------------------------------------------

const int64_t QUADSHELL9_face_node_order[] = {
  0,1,2,3,4,5,6,7,8,
  0,3,2,1,7,6,5,4,8
};

//---------------------------------------------------------------------
// TET4
//---------------------------------------------------------------------
const int64_t TET4_edge_node_order[] = {
  0,1, 1,2, 2,0, 0,3, 1,3, 2,3
};

const int64_t TET4_face_node_order[] = {
  0,1,3, 1,2,3, 0,3,2, 0,2,1
};

const int64_t TET4_face_edge_order[] = {
  0,4,3, 1,5,4, 3,5,2, 2,1,0
};

//---------------------------------------------------------------------
// TET8
//---------------------------------------------------------------------
const int64_t TET8_edge_node_order[] = {
  0,1, 1,2, 2,0, 0,3, 1,3, 2,3
};

const int64_t TET8_face_node_order[] = {
  0,1,3,4, 1,2,3,5, 0,3,2,7, 0,2,1,6 
};

//---------------------------------------------------------------------
// TET10
//---------------------------------------------------------------------
const int64_t TET10_edge_node_order[] = {
  0,1,4, 1,2,5, 2,0,6, 0,3,7, 1,3,8, 2,3,9
};

const int64_t TET10_face_node_order[] = {
  0,1,3,4,8,7,  1,2,3,5,9,8,
  0,3,2,7,9,6,  0,2,1,6,5,4
};

//---------------------------------------------------------------------
// TRI3
//---------------------------------------------------------------------
const int64_t TRI3_edge_node_order[] = {
  0,1, 1,2, 2,0
};

const int64_t TRI3_face_node_order[] = {
  0,1,2
};

const int64_t* TRI3_face_edge_order = TRI3_face_node_order;

//---------------------------------------------------------------------
// TRI4
//---------------------------------------------------------------------
const int64_t TRI4_face_node_order[] = {
  0,1,2,3
};

//---------------------------------------------------------------------
// TRI6
//---------------------------------------------------------------------
const int64_t TRI6_edge_node_order[] = {
  0,1,3,  1,2,4,  2,0,5
};

const int64_t TRI6_face_node_order[] = {
  0,1,2,3,4,5
};

//---------------------------------------------------------------------
// TRISHELL3
//---------------------------------------------------------------------
const int64_t* TRISHELL3_edge_node_order = TRI3_edge_node_order;

const int64_t TRISHELL3_face_node_order[] = {
  0,1,2,
  0,2,1
};

const int64_t TRISHELL3_face_edge_order[] = {
  0,1,2, 2,1,0
};

//---------------------------------------------------------------------
// TRISHELL6
//---------------------------------------------------------------------
const int64_t TRISHELL6_face_node_order[] = {
  0,1,2,3,4,5,
  0,2,1,3,5,4
};

//---------------------------------------------------------------------
// WEDGE6
//---------------------------------------------------------------------
const int64_t WEDGE6_edge_node_order[] = {
  0,1,  1,2, 2,0,
  3,4,  4,5, 5,3, 
  0,3,  1,4, 2,5
};

const int64_t WEDGE6_face_node_order[] = {
  0,1,4,3, 1,2,5,4, 0,3,5,2,
  0,2,1,-1, 3,4,5,-1
};

const int64_t WEDGE6_face_edge_order[] = {
  0,7,3,6,  1,8,4,7, 6,5,8,2,
  2,1,0,-1, 3,4,5,-1
};

const int64_t WEDGE6_node_per_face[] = {
  4, 4, 4, 3, 3
};

//---------------------------------------------------------------------
// WEDGE15
//---------------------------------------------------------------------
const int64_t WEDGE15_edge_node_order[] = {
  0,1, 6,   1,2, 7,  2,0, 8,
  3,4,12,   4,5,13,  5,3,14, 
  0,3, 9,   1,4,10,  2,5,11
};

const int64_t WEDGE15_face_node_order[] = {
  0,1,4,3,  6,10,12, 9,
  1,2,5,4,  7,11,13,10,
  0,3,5,2,  9,14,11, 8,
  0,2,1,    8, 7, 6,-1, -1,
  3,4,5,   12,13,14,-1, -1
};

const int64_t WEDGE15_node_per_face[] = {
  8, 8, 8, 6, 6
};

//---------------------------------------------------------------------

const int64_t* const edge_vertex_order[END_TOPOLOGY] = {
  HEX8_edge_node_order,
  HEX8_edge_node_order,
  HEX8_edge_node_order,
  PYRAMID5_edge_node_order,
  PYRAMID5_edge_node_order,
  QUAD4_edge_node_order,
  QUAD4_edge_node_order,
  QUAD4_edge_node_order,
  QUAD4_edge_node_order,
  QUAD4_edge_node_order,
  QUAD4_edge_node_order,
  TET4_edge_node_order,
  TET4_edge_node_order,
  TET4_edge_node_order,
  TRI3_edge_node_order,
  TRI3_edge_node_order,
  TRI3_edge_node_order,
  TRI3_edge_node_order,
  TRI3_edge_node_order,
  WEDGE6_edge_node_order,
  WEDGE6_edge_node_order
};

const int64_t* const edge_node_order[END_TOPOLOGY] = {
  HEX8_edge_node_order,
  HEX20_edge_node_order,
  HEX27_edge_node_order,
  PYRAMID5_edge_node_order,
  PYRAMID13_edge_node_order,
  QUAD4_edge_node_order,
  QUAD8_edge_node_order,
  QUAD9_edge_node_order,
  QUAD4_edge_node_order,
  QUAD8_edge_node_order,
  QUAD9_edge_node_order,
  TET4_edge_node_order,
  TET8_edge_node_order,
  TET10_edge_node_order,
  TRI3_edge_node_order,
  TRI3_edge_node_order,
  TRI6_edge_node_order,
  TRI3_edge_node_order,
  TRI6_edge_node_order,
  WEDGE6_edge_node_order,
  WEDGE15_edge_node_order
};

const int64_t* const face_vertex_order[END_TOPOLOGY] = {
  HEX8_face_node_order,
  HEX8_face_node_order,
  HEX8_face_node_order,
  PYRAMID5_face_node_order,
  PYRAMID5_face_node_order,
  QUAD4_face_node_order,
  QUAD4_face_node_order,
  QUAD4_face_node_order,
  QUADSHELL4_face_node_order,
  QUADSHELL4_face_node_order,
  QUADSHELL4_face_node_order,
  TET4_face_node_order,
  TET4_face_node_order,
  TET4_face_node_order,
  TRI3_face_node_order,
  TRI3_face_node_order,
  TRI3_face_node_order,
  TRISHELL3_face_node_order,
  TRISHELL3_face_node_order,
  WEDGE6_face_node_order,
  WEDGE6_face_node_order
};

const int64_t* const face_node_order[END_TOPOLOGY] = {
  HEX8_face_node_order,
  HEX20_face_node_order,
  HEX27_face_node_order,
  PYRAMID5_face_node_order,
  PYRAMID13_face_node_order,
  QUAD4_face_node_order,
  QUAD8_face_node_order,
  QUAD9_face_node_order,
  QUADSHELL4_face_node_order,
  QUADSHELL8_face_node_order,
  QUADSHELL9_face_node_order,
  TET4_face_node_order,
  TET8_face_node_order,
  TET10_face_node_order,
  TRI3_face_node_order,
  TRI4_face_node_order,
  TRI6_face_node_order,
  TRISHELL3_face_node_order,
  TRISHELL6_face_node_order,
  WEDGE6_face_node_order,
  WEDGE15_face_node_order
};

const int64_t* const face_edge_order[END_TOPOLOGY] = {
  HEX8_face_edge_order,
  HEX8_face_edge_order,
  HEX8_face_edge_order,
  PYRAMID5_face_edge_order,
  PYRAMID5_face_edge_order,
  QUAD4_face_edge_order,
  QUAD4_face_edge_order,
  QUAD4_face_edge_order,
  QUADSHELL4_face_edge_order,
  QUADSHELL4_face_edge_order,
  QUADSHELL4_face_edge_order,
  TET4_face_edge_order,
  TET4_face_edge_order,
  TET4_face_edge_order,
  TRI3_face_edge_order,
  TRI3_face_edge_order,
  TRI3_face_edge_order,
  TRISHELL3_face_edge_order,
  TRISHELL3_face_edge_order,
  WEDGE6_face_edge_order,
  WEDGE6_face_edge_order
};

const int64_t* const vertex_per_face[END_TOPOLOGY] = {
  0,
  0,
  0,
  PYRAMID5_node_per_face,
  PYRAMID5_node_per_face,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  WEDGE6_node_per_face,
  WEDGE6_node_per_face
};

const int64_t* const node_per_face[END_TOPOLOGY] = {
  0,
  0,
  0,
  PYRAMID5_node_per_face,
  PYRAMID13_node_per_face,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  WEDGE6_node_per_face,
  WEDGE15_node_per_face
};

const int64_t* const * const edge_per_face = vertex_per_face;

int get_topology(const char *name) {
  for (int i = 0; i < END_TOPOLOGY; ++i) {
    if (is_element_name(name, aliases[i])) {
      return i;
    }
  }
  return -1;
}