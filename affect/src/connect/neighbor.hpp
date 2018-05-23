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
 * @param maxElementsPerVertex (IN) maximum number of elements connected to a single vertex
 * @param elemSet              (IN) an array used for working space,
 *                                  with length = 2*maxElementsPerVertex
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
uint32_t neighbor_hex(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet,
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElement,
  const uint32_t * vertexToElementBegin,
  int64_t * neighbor);

uint32_t neighbor_tet(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet,
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElement,
  const uint32_t * vertexToElementBegin,
  int64_t * neighbor);

void neighbor_wedge(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet,
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElement,
  const uint32_t * vertexToElementBegin,
  int64_t * neighbor,
  uint32_t * numQuadFaces,
  uint32_t * numTriFaces);

uint32_t neighbor_quad(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet,
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElement,
  const uint32_t * vertexToElementBegin,
  int64_t * neighbor);

uint32_t neighbor_tri(
  size_t numElement,
  uint32_t maxElementsPerVertex,
  uint32_t * elemSet,
  const uint32_t * elementToVertex,
  const uint32_t * vertexToElement,
  const uint32_t * vertexToElementBegin,
  int64_t * neighbor);

#endif
