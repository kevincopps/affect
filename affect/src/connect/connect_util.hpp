#ifndef AFFECT_CONNECT_UTIL_H
#define AFFECT_CONNECT_UTIL_H

#include <cstdint>

/**
 * Given the array of element-to-element neighbors, 
 * (as output from a previous call to neighbor_hex), 
 * construct the array of face-to-vertex connectivity
 * for the faces of elements with no neighbors.
 *
 * @param numElement      (IN) number of elements
 * @param elementToVertex (IN) array of element-to-vertex connectivity,
 *                             length (8*numElement)
 * @param neighbor        (IN) array of element-to-element connectivity
 *                             with entries of -1 indicating no neighbor exists
 *                             length (6*numElement)
 * @param                (OUT) array of face-to-vertex connectivity
 *                             for all the faces where no neighbor exists
 *                             length (numBoundaryFaces) the return value
 *                             from a previous call to neighbor_hex
 */
void create_boundary_faces_hex(
    int64_t numElement,
    const int64_t* elementToVertex,
    const int64_t* neighbor,
    int64_t* boundaryFaceToVertex);

/**
 * Given the array of element-to-element neighbors, 
 * (as output from a previous call to neighbor_tet), 
 * construct the array of face-to-vertex connectivity
 * for the faces of elements with no neighbors.
 *
 * @param numElement      (IN) number of elements
 * @param elementToVertex (IN) array of element-to-vertex connectivity,
 *                             length (4*numElement)
 * @param neighbor        (IN) array of element-to-element connectivity
 *                             with entries of -1 indicating no neighbor exists
 *                             length (4*numElement)
 * @param                (OUT) array of face-to-vertex connectivity
 *                             for all the faces where no neighbor exists
 *                             length (numBoundaryFaces) the return value
 *                             from a previous call to neighbor_tet
 */
void create_boundary_faces_tet(
    int64_t numElement,
    const int64_t* elementToVertex,
    const int64_t* neighbor,
    int64_t* boundaryFaceToVertex);

/**
 * Given the array of element-to-element neighbors, 
 * (as output from a previous call to neighbor_wedge), 
 * construct the array of face-to-vertex connectivity
 * for the faces of elements with no neighbors.
 *
 * @param numElement      (IN) number of elements
 * @param elementToVertex (IN) array of element-to-vertex connectivity,
 *                             length (6*numElement)
 * @param neighbor        (IN) array of element-to-element connectivity
 *                             with entries of -1 indicating no neighbor exists
 *                             length (5*numElement)
 * @param boundaryQuadToVertex (OUT) array of face-to-vertex connectivity
 *                                   for all quadrilateral faces where no 
 *                                   neighbor exists, length (numQuadFaces) the 
 *                                   value returned from a previous call to 
 *                                   neighbor_wedge
 * @param boundaryTriToVertex  (OUT) array of face-to-vertex connectivity
 *                                   for all triangle faces where no 
 *                                   neighbor exists, length (numTriFaces) the 
 *                                   value returned from a previous call to 
 *                                   neighbor_wedge
 */
void create_boundary_faces_wedge(
    int64_t numElement,
    const int64_t* elementToVertex,
    const int64_t* neighbor,
    int64_t* boundaryQuadToVertex,
    int64_t* boundaryTriToVertex);
  
/**
 * Given the element-to-vertex connectivity, and the vertex-to-element 
 * connectivity, return an array of the element-to-element connectivity.
 * On return, each element side has an entry in the neighbor array with
 * the value equal either to the element ID sharing the side, 
 * or -1 for sides without a neighbor.
 *
 * @param numElement           (IN) number of elements
 * @param elemSet              (IN) an array used for working space, 
 *                                  with length equal to the maximum number 
 *                                  of elements connected to a single vertex
 * @param elementToVertex      (IN) array of element to vertex 
 *                                  connectivity, with length (8*numElement)
 * @param vertexToElement      (IN) array of vertex-to-element connectivity
 * @param vertexToElementBegin (IN) array with length (numVertices+1),
 *                                  where the ith entry is the starting index 
 *                                  into the vertexToElement array for 
 *                                  the ith element
 * @param neighbor            (OUT) on return, the array with an entry for 
 *                                  each element side, containing either the 
 *                                  ID of the element neighbor, or -1 no 
 *                                  neighbor exists, length (6*numElement)
 *
 * @return the number of element sides without neighbors
 */
int64_t neighbor_hex(
    int64_t numElement,
    int64_t * elemSet,
    const int64_t * elementToVertex,
    const int64_t * vertexToElementBegin,
    const int64_t * vertexToElement,
    int64_t * neighbor);

/**
 * Construct element-to-element connectivity.
 * @see neighbor_hex
 */
int64_t neighbor_tet(
    int64_t numElement,
    int64_t * elemSet,
    const int64_t * elementToVertex,
    const int64_t * vertexToElementBegin,
    const int64_t * vertexToElement,
    int64_t * neighbor);

/**
 * Construct element-to-element connectivity.
 * @see neighbor_hex
 */
void neighbor_wedge(
    int64_t numElement,
    int64_t * elemSet,
    const int64_t * elementToVertex,
    const int64_t * vertexToElementBegin,
    const int64_t * vertexToElement,
    int64_t * neighbor,
    int64_t * numQuadFaces,
    int64_t * numTriFaces);

/**
 * Determines the local edge number of an element which shares a set
 * of two vertices.
 *
 * @param edgeVertexOrder an array of the local enumeration of the
 *            element edge vertices, length (2*numEdgePerElement)
 * @param numEdgePerElement number of edges in an element
 * @param localVertex array holding the vertex numbers of the element
 * @param vertex0 a vertex number that the element edge must use
 * @param vertex1 a vertex number that the element edge must use
 *
 * @return local edge number which uses vertex0, and vertex1, 
 *         in the range [0,numEdgePerElement-1], or numEdgePerElement
 *         if there is no such edge
 */
inline int64_t which_edge(
    const int64_t * edgeVertexOrder,
    int64_t numEdgePerElement,
    const int64_t * localVertex,
    int64_t vertex0,
    int64_t vertex1) {

    int64_t edge = 0;
    for ( ; edge < numEdgePerElement; ++edge) {

        int64_t edgeVertex0 = localVertex[ *edgeVertexOrder++ ];
        int64_t edgeVertex1 = localVertex[ *edgeVertexOrder++ ];
        if (edgeVertex0 == vertex1 && edgeVertex1 == vertex0) break;
        if (edgeVertex0 == vertex0 && edgeVertex1 == vertex1) break;

    }
    return edge;
}

/**
 * Determines the local face number of an element which shares a set
 * of three vertices.
 *
 * @param faceVertexOrder an array of the local enumeration of the
 *            element face vertices, length (numVertexPerFace*numFacePerElement)
 * @param numFacePerElement number of faces in an element
 * @param numVertexPerFace number of vertices on an element face
 * @param localVertex array holding the vertex numbers of the element
 * @param vertex0 a vertex number that the element face must use
 * @param vertex1 a vertex number that the element face must use
 * @param vertex2 a vertex number that the element face must use
 *
 * @return local face number which uses vertex0, vertex1, and vertex2,
 *         in the range [0,numFacePerElement-1], or numFacePerElement
 *         if there is no such face
 */
inline int which_face(
    const int64_t * faceVertexOrder,
    const int64_t numFacePerElement,
    int64_t numVertexPerFace,
    const int64_t * localVertex,
    int64_t vertex0,
    int64_t vertex1,
    int64_t vertex2) {

    int64_t vertex;
    bool match0, match1, match2;
    const int64_t * faceVertex;

    for (int64_t face = 0; face < numFacePerElement; ++face) {
        match0 = false;
        match1 = false;
        match2 = false;
        faceVertex = faceVertexOrder + face * numVertexPerFace;
        for (int64_t n = 0; n < numVertexPerFace; ++n) {
            vertex = localVertex[faceVertex[n]];
            if ( !match0 )
                match0 = (vertex0 == vertex);
                if (match0) continue;
            if ( !match1 )
                match1 = (vertex1 == vertex);
                if (match1) continue;
            if ( !match2 )
                match2 = (vertex2 == vertex);
        }
        if (match0 && match1 && match2) return face;
    }
    // if we reached here, we could not find a face
    // assert(false);
    return numFacePerElement;
}

/**
 * Determines the local face number of an element which shares a set
 * of three vertices.
 *
 * @param faceVertexOrder an array of the local enumeration of the
 *            element face vertices, length (numVertexPerFace*numFacePerElement)
 * @param numFacePerElement number of faces in an element
 * @param maxVertexPerFace number of vertices on an element face (4 for wedges/pyramids)
 * @param localVertex array holding the vertex numbers of the element,
 *        contains a -1 when there are less vertices on a face than maxVertexPerFace
 * @param vertex0 a vertex number that the element face must use
 * @param vertex1 a vertex number that the element face must use
 * @param vertex2 a vertex number that the element face must use
 *
 * @return local face number which uses vertex0, vertex1, and vertex2,
 *         in the range [0,numFacePerElement-1], or numFacePerElement
 *         if there is no such face
 */
inline int64_t which_face_irregular(
    const int64_t * faceVertexOrder,
    int64_t numFacePerElement,
    int64_t maxFacePerElement,
    const int64_t * numVertexPerFace,
    const int64_t * localVertex,
    int64_t vertex0,
    int64_t vertex1,
    int64_t vertex2) {

    int64_t vertex;
    bool match0, match1, match2;
    int64_t* vertexOrder = const_cast<int64_t*>(faceVertexOrder);
    int64_t* faceVertex;
    int64_t* lastVertex;

    for (int64_t face = 0; face < numFacePerElement;
        ++face, ++numVertexPerFace, vertexOrder += maxFacePerElement) {
        match0 = false;
        match1 = false;
        match2 = false;
        lastVertex = vertexOrder + *numVertexPerFace;
        for (faceVertex = vertexOrder; faceVertex < lastVertex; ++faceVertex) {
            vertex = localVertex[*faceVertex];
            if ( !match0 )
                match0 = (vertex0 == vertex);
                if (match0) continue;
            if ( !match1 )
                match1 = (vertex1 == vertex);
                if (match1) continue;
            if ( !match2 )
                match2 = (vertex2 == vertex);
        }
        if (match0 && match1 && match2) return face;
    }
    // if we reached here, we could not find a face
    // assert(false);
    return numFacePerElement;
}


#endif
