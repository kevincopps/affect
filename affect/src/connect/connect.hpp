#ifndef AFFECT_CONNECT_H
#define AFFECT_CONNECT_H

#include <cstdint>
#include <vector>

/**
 * Convert a one-based connectivity array (Fortran style)
 * to zero-based (C style), that is, 
 * subtract one from every entry in the given array.
 *
 * @param n     (IN) length of array
 * @param array (INOUT) the one-based integer array to convert
 */
void to_zero_based(int64_t n, int64_t* array);

/**
 * Convert a zero-based connectivity array (C style)
 * to one-based (Fortran style), that is, 
 * subtract add one to every entry in the given array.
 *
 * @param n     (IN) length of array
 * @param array (INOUT) a zero-based integer array to convert
 */
void to_one_based(int64_t n, int64_t* array);

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
    int64_t numElement,
    int64_t numBoundaryQuadFace,
    int64_t numBoundaryTriFace,
    const int64_t * elementToVertex,
    const int64_t * neighbor,
    int64_t * boundaryFaceToVertex);

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
    int64_t numElement,
    int64_t numEdgePerElement,
    int64_t numVertexPerElement,
    const int64_t * edgeVertexOrder,
    const int64_t * elementToVertex,
    const int64_t * vertexToElementBegin,
    const int64_t * vertexToElement,
    int64_t * elementToEdge,
    int64_t * numInternalEdge,
    int64_t * numExternalEdge);

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
 */
void connect_element_to_element(
    const char* elementName,
    int64_t numElement,
    int64_t maxElementPerVertex,
    const int64_t * elementToVertex,
    const int64_t * vertexToElementBegin,
    const int64_t * vertexToElement,
    int64_t * neighbor,
    int64_t * numBoundaryQuadFace,
    int64_t * numBoundaryTriFace);
  
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
    int64_t numElement,
    int64_t numFacePerElement,
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

int64_t count_vertex_to_element(
    int64_t numNode,
    int64_t numVertexPerElement,
    int64_t numElement,
    const int64_t * elementToVertex,
    int64_t * vertexToElementCount);

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
    int64_t numVertex,
    int64_t numVertexPerElement,
    int64_t numElement,
    const int64_t * elementToVertex,
    int64_t * vertexToElementCount,
    int64_t * vertexToElement );


#endif
