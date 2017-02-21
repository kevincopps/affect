#ifndef AFFECT_NEIGHBOR_H
#define AFFECT_NEIGHBOR_H


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
  const int64_t * vertexToElement,
  const int64_t * vertexToElementBegin,
  int64_t * neighbor);

int64_t neighbor_tet(
  int64_t numElement, 
  int64_t * elemSet,
  const int64_t * elementToVertex,
  const int64_t * vertexToElement,
  const int64_t * vertexToElementBegin,
  int64_t * neighbor);

void neighbor_wedge(
  int64_t numElement, 
  int64_t * elemSet,
  const int64_t * elementToVertex,
  const int64_t * vertexToElement,
  const int64_t * vertexToElementBegin,
  int64_t * neighbor,
  int64_t * numQuadFaces,
  int64_t * numTriFaces);

int64_t neighbor_quad(
  int64_t numElement, 
  int64_t * elemSet,
  const int64_t * elementToVertex,
  const int64_t * vertexToElement,
  const int64_t * vertexToElementBegin,
  int64_t * neighbor);

int64_t neighbor_tri(
  int64_t numElement, 
  int64_t * elemSet,
  const int64_t * elementToVertex,
  const int64_t * vertexToElement,
  const int64_t * vertexToElementBegin,
  int64_t * neighbor);


#endif
