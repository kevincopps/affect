#include <algorithm>
#include <vector>

#ifdef AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT
#include <iostream>
#include <iomanip>
#endif

using namespace std;


void connect_vertex_to_element(
  int64_t numVertex,
  int64_t numVertexPerElement,
  int64_t numElement,
  const int64_t * elementToVertex,
  int64_t * vertexToElementCount,
  int64_t * vertexToElement)
{
  // Iterate over elements and copy the element to the vertices it is connected to.
  // During the loop we will keep incrementing by one our vertexToElementCount array which acts as a
  // pointer to the beginning of the entries for the ith vertex.
  //
  for (int64_t iElem = 0, start = 0, finish = numVertexPerElement; 
       iElem < numElement; ++iElem, start = finish, finish += numVertexPerElement) {
    for (int64_t j = start; j < finish; ++j)
      vertexToElement[vertexToElementCount[elementToVertex[j]]++] = iElem;
  }

  // reverse the incrementing we did in the previous loops
  for (int64_t i = numVertex-1; i > 0; --i) {
    vertexToElementCount[i] = vertexToElementCount[i-1];
  }
  vertexToElementCount[0] = 0;

#ifdef AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT
  cout << "----Debugging connect_vertex_to_element" << endl;
  //for (int64_t i = 0; i <= numVertex; ++i )
  //  cout << "  " << setw(3) << i << setw(8) << vertexToElementCount[i] << endl;
  for (int64_t i = 0; i < numVertex; ++i ) {
    cout << " vertex " << setw(3) << i << " has elem";
    for (int64_t j = vertexToElementCount[i]; j < vertexToElementCount[i+1]; ++j)
      cout << " " << vertexToElement[j];
    cout << endl;
  }
#endif
  return;
}

