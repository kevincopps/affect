#ifndef AFFECT_CONNECT_H
#define AFFECT_CONNECT_H

#include <cstdint>
#include <vector>
#include <iostream>

#include <util/aligned_array.hpp>

/**
 * Fill an array with the boundary-face-to-vertex connectivity, given
 * the number of boundary faces and both the element-to-vertex and 
 * element-to-element connectivity.
 *
 * @param elementName           (IN) type of element
 * @param numElement            (IN) number of elements
 * @param numBoundaryQuadFace   (IN) number of quadrilateral boundary faces,
 *                                   must equal zero if elementName is TET
 * @param numBoundaryTriFace    (IN) number of triangular boundary faces,
 *                                   must equal zero if elementName is HEX
 * @param elementToVertex       (IN) array of element-to-vertex connectivity,
 *                                   length (numElement*numVertexPerElement)
 * @param neighbor              (IN) an array with an entry for each element 
 *                                   side equal to the ID of the element 
 *                                   neighbor or -1 if no neighbor exists,
 *                                   length (numElement*numFacePerElement)
 * @param boundaryFaceToVertex (OUT) array of boundary-face-to-vertex 
 *                                   connectivity with length 
 *                                   (4*numQuadBoundaryFace + 
 *                                    3*numTriBoundaryFace); for wedge
 *                                   elements, the quad face vertices begin
 *                                   at index zero, and the triangular face
 *                                   vertices begin at index numQuadFaces*4.
 *
 * @see connect_element_to_element
 */
void connect_boundary_face_to_vertex(
    char* elementName,
    uint32_t numElement,
    int32_t numBoundaryQuadFace,
    int32_t numBoundaryTriFace,
    const uint32_t * elementToVertex,
    const int64_t * neighbor,
    uint32_t * boundaryFaceToVertex);

/**
 * Fill a vector with the element-to-edge connectivity, and
 * return a count both of the number of edges shared by two 
 * or more elements, and the number of edges used by only a 
 * single element.
 *
 * @param numElement            (IN) number of elements
 * @param numEdgePerElement     (IN) number of edges in a single element
 * @param numVertexPerElement   (IN) number of vertices in a single element
 * @param edgeVertexOrder       (IN) array of local element-to-vertex enumeration
 * @param elementToVertex       (IN) array of element-to-vertex connectivity,
 *                                   length (numElement*numVertexPerElement)
 * @param vertexToElementBegin  (IN) the index of the first entry in 
 *                                   vertexToElement for a vertex,
 *                                   length (numVertex+1)
 * @param vertexToElement       (IN) the vertex-to-element connectivity,
 *                                   length (vertexToElementBegin[numVertex+1]-1)
 * @param elementToEdge        (OUT) the element-to-edge connectivity
 * @param numInternalEdge      (OUT) number of edges shared by two or more elements
 * @param numExternalEdge      (OUT) number of edges used by only one element
 */  
void connect_element_to_edge(
    uint32_t numElement,
    uint32_t numEdgePerElement,
    uint32_t numVertexPerElement,
    const uint32_t * edgeVertexOrder,
    const uint32_t * elementToVertex,
    const uint32_t * vertexToElementBegin,
    const uint32_t * vertexToElement,
    int64_t * elementToEdge,
    uint32_t * numInternalEdge,
    uint32_t * numExternalEdge);

/**
 * For a single block of elements with 3D topology, return an array of 
 * the element-to-element connectivity, given both the element-to-vertex 
 * connectivity, and the vertex-to-element connectivity.
 *
 * @param numElement            (IN) number of elements
 * @param elemSet               (IN) an array used for working space, 
 *                                   with length equal to the maximum number 
 *                                   of elements connected to a single vertex
 * @param elementToVertex       (IN) array of element to vertex 
 *                                   connectivity, with length (8*numElement)
 * @param vertexToElementBegin  (IN) array with length (numVertices+1),
 *                                   where the ith entry is the starting index 
 *                                   into the vertexToElement array for 
 *                                   the ith element
 * @param vertexToElement       (IN) array of vertex-to-element connectivity
 * @param neighbor             (OUT) on return, the array with an entry for 
 *                                   each element side, containing either the 
 *                                   ID of the element neighbor or -1 if no 
 *                                   neighbor exists,
 *                                   length (numFacePerElement*numElement)
 * @param numBoundaryQuadFace  (OUT) the number of quadrilateral element faces
 *                                   without a neighbor
 * @param numBoundaryTriFace   (OUT) the number of triangular element faces
 *                                   without a neighbor
void connect_element_to_element(
    const char* elementName,
    uint32_t numElement,
    uint32_t maxElementPerVertex,
    const uint32_t * elementToVertex,
    const uint32_t * vertexToElementBegin,
    const uint32_t * vertexToElement,
    int64_t * neighbor,
    uint32_t * numBoundaryQuadFace,
    uint32_t * numBoundaryTriFace);
*/
  
/**
 * Replace the array of element-to-element connectivity with the array of
 * element-to-face connectivity (where interior faces will be enumerated
 * first, followed by the faces on the boundary) and return the number
 * of interior faces.
 *
 * @param numElement           (IN) number of elements
 * @param numFacePerElement    (IN) number of faces in a single element
 * @param neighbor          (INOUT) the array with an entry for each element 
 *                                  side--upon entering this function
 *                                  the entries are either the ID of the element 
 *                                  neighbor or -1 if no neighbor exists--
 *                                  upon return, the entries are face IDs,
 *                                  length (numFacePerElement*numElement)
 *
 * @return                    (OUT) number of internal faces (the starting
 *                                  ID of boundary faces--any face ID
 *                                  less than this number is an internal
 *                                  face).
 *
 * @see connect_element_to_element
 */
int64_t connect_element_to_face(
    uint32_t numElement,
    uint32_t numFacePerElement,
    int64_t * neighbor);

/**
 * For a uniform block of elements, construct an array which
 * contains the sums of counts of the number of elements connected to each vertex.
 * The purpose of this array is to be used as a lookup index into the
 * array holding the set of elements connected to each vertex.
 *
 * Note: the length of the input vertexToElementCount array must be one longer
 * than the number of vertices in the block.
 *
 * @param numVertex             (IN) number of vertices
 * @param numVertexPerElem      (IN) number of vertices per element
 * @param numElement            (IN) number of elements
 * @param elementToVertex       (IN) the element to vertex connectivity,
 *                                   length (numElement*numVertexPerElement)
 * @param vertexToElementCount (OUT) on exit, the ith entry is the sum
 *                                   of the number of elements connected
 *                                   to the vertices in the range [0,i-1].
 *                                   The 0th entry is always zero, and the
 *                                   numVertex+1 entry is the sum of
 *                                   all the elements connected to each the
 *                                   vertices in the block.
 *                                   length (numVertex+1)
 *
 * @return maximum number of elements connected to a single vertex in this block
 */

uint32_t count_vertex_to_element(
    uint32_t numVertex,
    uint32_t numVertexPerElement,
    uint32_t numElement,
    const aligned::uint32_ptr __restrict elementToVertex,
    aligned::uint32_ptr __restrict vertexToElementCount);

/**
 * For a uniform block of elements, construct a packed array of the 
 * vertex-to-element connectivity, it uses as input an array vertexToElementCount
 * as returned from the function count_vertex_to_element. One exit, the
 * list of elements connected to vertex i are contained in the
 * vertexToElementCount(j) where j is in the range [vertexToElementCount[i], vertexToElementCount[i+1]]
 * Note: vertexToElement must be previously allocated.
 *
 * @param numVertex             (IN) number of vertices
 * @param numVertexPerElem      (IN) number of vertices per element
 * @param numElement            (IN) number of elements
 * @param elementToVertex       (IN) the element to vertex connectivity,
 *                                   length (numElement*numVertexPerElement)
 * @param vertexToElementCount (IN)  array as returned from function count_vertex_to_element()
 *                                   length (numVertex+1)
 * @param vertexToElement      (OUT) the vertex-to-element connectivity,
 *                                   length (vertexToElementCount[numVertex+1])
 */
void connect_vertex_to_element(
    uint32_t numVertex,
    uint32_t numVertexPerElement,
    uint32_t numElement,
    const aligned::uint32_ptr __restrict elementToVertex,
    aligned::uint32_ptr __restrict vertexToElementCount,
    aligned::uint32_ptr __restrict vertexToElement );

/**
 *   For a uniform block of elements, construct array of neighbor element IDs and array of local neighbor face IDs
 *   for each element face using the sibling half-facet algorithm.*
 *
 *    Input Args:
 *        topology: type of cells in this block
 *        element_to_vertex(|uint32_2d|): array of shape(num_elements, num_vertex_per_element)
 *
 *    Output Args:
 *        neighbor_elements: array of neighbor element ID for each element face
 *        neighbor_faces: array of neighbor local face ID for each element face
 *
 *    Returns:
 *        num_boundary_faces: number of element faces on the boundary
 */
void connect_element_neighbors(
    int topology,
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const uint32_t * element_to_vertex,
    uint32_t * neighbor_elements,
    int8_t * neighbor_faces);

/**
 * For a uniform block of elements, construct two arrays, storing the last element connected to each vertex and
 * a corresponding local face id.
 *
 * @param num_vertices          (IN) number of vertices
 * @param num_elements          (IN) number of elements
 * @param element_to_vertex     (IN) the element to vertex connectivity, length (numElement*numVertexPerElement)
 * @param neighbor_faces        (IN) array of neighbor local face ID for each element face, as returned by
 *                                   connect_element_neighbors function
 * @param vertex_facet_element (OUT) the vertex-to-last-element containing the vertex connectivity
 *                                   length (num_vertices)
 * @param vertex_facet_face    (OUT) one local face ID of the last element containing the vertex
 */
 void connect_vertex_to_element_face(
    int topology,
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const aligned::uint32_ptr __restrict element_to_vertex,
    const aligned::int8_ptr __restrict neighbor_faces,
    aligned::uint32_ptr __restrict vertex_facet_element,
    aligned::int8_ptr __restrict vertex_facet_face);

#endif
