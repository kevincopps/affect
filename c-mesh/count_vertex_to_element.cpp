#include <algorithm>
#include <vector>

#ifdef AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT
#include <iostream>
#include <iomanip>
#endif

using namespace std;

int64_t count_vertex_to_element(
  int64_t numVertex,
  int64_t numVertexPerElement,
  int64_t numElement,
  const int64_t * elementToVertex,
  int64_t * vertexToElementCount)
{
  for (int64_t i = 0; i < numVertex + 1; ++i)
    vertexToElementCount[i] = 0;

  // Count number of elements sharing each vertex.
  int64_t elementToVertexLength = numElement * numVertexPerElement;
  for (int64_t i = 0; i < elementToVertexLength; ++i)
    vertexToElementCount[ elementToVertex[i] ]++;

  // Now for each vertex, replace the count of elements sharing it
  // with its beginning index to the lists of its elements.
  int64_t maxElementPerVertex = 0;
  int64_t length = 0, elementsPerNode;
  for (int64_t i = 0; i < numVertex; ++i) {
    elementsPerNode = vertexToElementCount[i];
    if (elementsPerNode > maxElementPerVertex)
      maxElementPerVertex = elementsPerNode;
    vertexToElementCount[i] = length;
    length += elementsPerNode;
  }
  vertexToElementCount[numVertex] = length;

  return maxElementPerVertex;
}

