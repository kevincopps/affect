#ifndef AFFECT_ELEMENT_TOPOLOGY_H
#define AFFECT_ELEMENT_TOPOLOGY_H

#include <cstdint>

// All element topology IDs.
// the value END_TOPOLOGY is used for error checking
enum topology_type {
  HEX8,HEX20,HEX27,
  PYRAMID5,PYRAMID13,
  QUAD4,QUAD8,QUAD9,
  QUADSHELL4,QUADSHELL8,QUADSHELL9,
  TET4,TET8,TET10,
  TRI3,TRI4,TRI6,
  TRISHELL3,TRISHELL6,
  WEDGE6,WEDGE15,
  END_TOPOLOGY};

/**
 * Returns true if a character string is found
 * within a second string (without regard to
 * case).
 *
 * @param s1 what to search for
 * @param s2 the field in which to search 
 */
bool is_element_name(const char *s1, const char *s2);

/**
 * Given the name of an element type, or any valid
 * alias, returns a positive integer ID representing the type of
 * topology, or a negative number if the name is not
 * recognized.
 *
 * @param name an element topology, for example,
 *             "HEX8", "hex8", "TET4", "tetra4", etc.
 */
int get_topology(const char *name);

/**
 * Returns the preferred name given a topology.
 */
extern const char *names[];

/**
 * Returns the possible aliases for names given a topology.
 *
 * A string of names seperated by "|".
 */
extern const char *aliases[];

/**
 * Returns the number of spatial dimensions for a
 * given element topology.
 */
extern const int64_t num_spatial_dim[];

/**
 * Returns the number of vertices in a single element
 * of the given topology.
 * Vertices are end points of edges, and all vertices
 * have an associated node.
 */
extern const int64_t num_vertex[];

/**
 * Returns the number nodes in a single element
 * of the given topology,
 * where a node has unique values of coordinates
 * associated with it. Not every node is necessarily
 * associated with a vertex, i.e., there can exist
 * mid-face or mid-edge nodes.
 */
extern const int64_t num_node[];
  
/**
 * Returns the number of edges of a single element
 * of the given topology.
 */
extern const int64_t num_edge[];

/** 
 * Returns the number of faces of a single element
 * of the given topology. If the topology is 1D
 * or 2D, returns zero.
 */
extern const int64_t num_face[];

/**
 * Returns the number of nodes per each edge
 * of a single element of the given topology.
 */
extern const int64_t num_node_per_edge[];

/**
 * Returns the maximum number of vertices per face of
 * an element of the given topology.
 */
extern const int64_t max_vertex_per_face[];

/**
 * Returns the maximum number of nodes per face of
 * an element of the given topology.
 */
extern const int64_t max_node_per_face[];

/** 
 * Returns the maximum number of edges per face of
 * an element with the given topology.
 */
extern const int64_t* const max_edge_per_face;

/**
 * Returns the enumeration of the edge vertices, for all edges
 * in the given element topology.
 *
 * length 2 * num_edge
 */
extern const int64_t* const edge_vertex_order[];

/**
 * Returns the enumeration of the edge nodes, for all edges
 * in the given element topology, where negative entries
 * mean that no node exists at that position.
 *
 * length num_node_per_edge * num_edge
 */
extern const int64_t* const edge_node_order[];

/**
 * Return the enumeration of the face vertices, for all faces
 * in the given element topology.
 *
 * length num_vertex_per_face * num_face
 */
extern const int64_t* const face_vertex_order[];

/**
 * Return the enumeration of the face nodes, for all faces
 * in the given element topology.
 *
 * length num_vertex_per_face * num_face
 */
extern const int64_t* const face_node_order[];

/**
 * Return the enumeration of the face edges, for all 
 * faces in the given element topology.
 */
extern const int64_t* const face_edge_order[];

/**
 * Return the number of vertices per face, for all
 * faces in the given element topology.
 *
 * length num_face
 */
extern const int64_t* const vertex_per_face[];

/**
 * Return the number of nodes per face, for all
 * faces in the given element topology.
 *
 * length num_face
 */
extern const int64_t* const node_per_face[];

/**
 * Return the number of edge per face, for all
 * faces in the given element topology.
 *
 * length num_face
 */
extern const int64_t* const * const edge_per_face;

/**
 *----------------------------------------------------------
 * NOTES on topology arrays:
 *
 * local sizes of vertices per face -
 * entry 0 returns number of vertices for all faces if homogenous
 *         returns -1 if faces have differing topology
 * has length (nface+1)
 *
 * local sizes of edges per face - 
 * entry 0 returns number of edges for all faces if homogenous
 *         returns -1 if faces have differing topology
 * has length (nface+1)
 *----------------------------------------------------------
 */ 


#endif // AFFECT_ELEMENT_TOPOLOGY_H
