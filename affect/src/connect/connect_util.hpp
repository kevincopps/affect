#ifndef AFFECT_CONNECT_UTIL_H
#define AFFECT_CONNECT_UTIL_H

#include <cstdint>
#include <iostream>

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
    uint32_t numElement,
    const uint32_t* elementToVertex,
    const int64_t* neighbor,
    uint32_t* boundaryFaceToVertex);

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
    uint32_t numElement,
    const uint32_t* elementToVertex,
    const int64_t* neighbor,
    uint32_t* boundaryFaceToVertex);

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
    uint32_t numElement,
    const uint32_t* elementToVertex,
    const int64_t* neighbor,
    uint32_t* boundaryQuadToVertex,
    uint32_t* boundaryTriToVertex);
  
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
inline uint32_t which_edge(
    const uint32_t * edgeVertexOrder,
    uint32_t numEdgePerElement,
    const uint32_t * localVertex,
    uint32_t vertex0,
    uint32_t vertex1) {

    uint32_t edge = 0;
    for ( ; edge < numEdgePerElement; ++edge) {

        uint32_t edgeVertex0 = localVertex[ *edgeVertexOrder++ ];
        uint32_t edgeVertex1 = localVertex[ *edgeVertexOrder++ ];
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
    const uint32_t * faceVertexOrder,
    const uint32_t numFacePerElement,
    uint32_t numVertexPerFace,
    const uint32_t * localVertex,
    uint32_t vertex0,
    uint32_t vertex1,
    uint32_t vertex2) {

    const uint32_t * faceVertex;
    uint32_t vertex;
    bool match0, match1, match2;

    for (uint32_t face = 0; face < numFacePerElement; ++face) {
        match0 = false;
        match1 = false;
        match2 = false;
        faceVertex = faceVertexOrder + face * numVertexPerFace;
        for (uint32_t n = 0; n < numVertexPerFace; ++n) {
            vertex = localVertex[faceVertex[n]];
            if ( !match0 ) {
                match0 = (vertex0 == vertex);
                if (match0) continue;
            }
            if ( !match1 ) {
                match1 = (vertex1 == vertex);
                if (match1) continue;
            }
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
 *        contains a UINT32_MAX when there are less vertices on a face than maxVertexPerFace
 * @param vertex0 a vertex number that the element face must use
 * @param vertex1 a vertex number that the element face must use
 * @param vertex2 a vertex number that the element face must use
 *
 * @return local face number which uses vertex0, vertex1, and vertex2,
 *         in the range [0,numFacePerElement-1], or numFacePerElement
 *         if there is no such face
 */
inline uint32_t which_face_irregular(
    const uint32_t * faceVertexOrder,
    uint32_t numFacePerElement,
    uint32_t maxVertexPerFace,
    const uint32_t * numVertexPerFace,
    const uint32_t * localVertex,
    uint32_t vertex0,
    uint32_t vertex1,
    uint32_t vertex2) {

    uint32_t vertex;
    bool match0, match1, match2;
    uint32_t* vertexOrder = const_cast<uint32_t*>(faceVertexOrder);
    uint32_t* faceVertex;
    uint32_t* lastVertex;

    // std::cerr << "which_face_irregular: " << std::endl;
    // std::cerr << "    numFacePerElement = " << numFacePerElement << std::endl;
    // std::cerr << "    vertex0 = " << vertex0 << std::endl;
    // std::cerr << "    vertex1 = " << vertex1 << std::endl;
    // std::cerr << "    vertex2 = " << vertex2 << std::endl;

    for (uint32_t face = 0; face < numFacePerElement;
        ++face, ++numVertexPerFace, vertexOrder += maxVertexPerFace) {

        // std::cerr << "        face = " << face << std::endl;
        // std::cerr << "        numVertexPerFace = " << *numVertexPerFace << std::endl;

        match0 = false;
        match1 = false;
        match2 = false;
        lastVertex = vertexOrder + *numVertexPerFace;
        for (faceVertex = vertexOrder; faceVertex < lastVertex; ++faceVertex) {
            vertex = localVertex[*faceVertex];

            // std::cerr << "        vertex = " << vertex << std::endl;

            if ( !match0 ) {
                match0 = (vertex0 == vertex);
                // if (match0) std::cerr << "        match0 = " << vertex0 << std::endl;
                if (match0) continue;
            }
            if ( !match1 ) {
                match1 = (vertex1 == vertex);
                // if (match1) std::cerr << "        match1 = " << vertex1 << std::endl;
                if (match1) continue;
            }
            if ( !match2 ) {
                match2 = (vertex2 == vertex);
                // if (match2) std::cerr << "        match2 = " << vertex2 << std::endl;
            }
        }
        if (match0 && match1 && match2) return face;
    }
    // if we reached here, we could not find a face
    // assert(false);
    return numFacePerElement;
}


#endif
