#include <algorithm>
#include <vector>
#include <cassert>
#include <cstdint>

#ifdef AFFECT_VERBOSE
#include <iostream>
#include <iomanip>
#endif

using namespace std;


void connect_edge_to_vertex(
  int64_t numVertex, 
  int64_t numVertexPerFace, 
  int64_t numFaces, 
  const int64_t * faceToVertex,
  vector<int64_t> & edgeToVertex,
  int64_t * edgeBegin)    // working space (numVertex+1)
{
  assert( numVertexPerFace > 2);

  std::fill(edgeBegin, edgeBegin + numVertex + 1, 0);

  // for each edge in a face,
  // add 1 to a count associated with its lowest numbered vertex
  for (int64_t iFace = 0; iFace < numFaces; ++iFace) {
  
    const int64_t * localVertex = faceToVertex + iFace * numVertexPerFace;
  
    for (int64_t iEdge = 0; iEdge < numVertexPerFace; ) {
    
      int64_t vertex0 = localVertex[iEdge++] - 1;
      int64_t vertex1 = localVertex[iEdge % numVertexPerFace] - 1;
      
      edgeBegin[ vertex0 < vertex1 ? vertex0 : vertex1 ]++;
    }
  }
  
  // Now for each vertex, replace the count of edges sharing it
  // with its beginning index to the lists of its edges.
  int64_t length = 0, edgesPerNode;
  for (int64_t i = 0; i < numVertex; ++i) {
    edgesPerNode = edgeBegin[i];
    edgeBegin[i] = length;
    length += edgesPerNode;
  }
  edgeBegin[numVertex] = length;

  // now get memory for the array of vertex-to-vertex connectivity
  vector<int64_t> vertexToVertex(length);

  // now copy names of vertexs to the vertexs it is connected to
  // most, if not all, entries are repeated once  
  for (int64_t iFace = 0; iFace < numFaces; ++iFace) {
  
    const int64_t * localVertex = faceToVertex + iFace * numVertexPerFace;
  
    for (int64_t iEdge = 0; iEdge < numVertexPerFace; ) {
  
      int64_t vertex0 = localVertex[iEdge++] - 1;
      int64_t vertex1 = localVertex[iEdge % numVertexPerFace] - 1;

      if (vertex0 < vertex1)
        vertexToVertex[ edgeBegin[vertex0]++ ] = vertex1;
      else
        vertexToVertex[ edgeBegin[vertex1]++ ] = vertex0;
    }  
  }
  
  // reset the begin indices
  for (int64_t i = numVertex-1; i > 0; --i)
    edgeBegin[i] = edgeBegin[i-1];
  edgeBegin[0] = 0;
  
  vector<int64_t> vertexToVertexCount(numVertex,0);

  // for each vertex, sort the list of vertexs connected to it
  // most, if not all entries are repeated once
  for (int64_t iNode = 0; iNode < numVertex; ++iNode) {
  
    int64_t begin = edgeBegin[iNode],
        end = edgeBegin[iNode+1];
        
    if (begin == end) continue;

    std::sort( &vertexToVertex[ begin ], &vertexToVertex[ end ]);

    int64_t* first = &vertexToVertex[ edgeBegin[iNode] ];
    int64_t* last  = &vertexToVertex[ edgeBegin[iNode+1] ];

    // remove repeated entries
    last = std::unique(first, last);
                   
    vertexToVertexCount[iNode] = (int)(last - first);
  }
  
  // now count number of total edges
  int64_t numEdges = 0;
  vector<int64_t>::const_iterator iCount = vertexToVertexCount.begin(),
                           lastCount = vertexToVertexCount.end();
  while (iCount != lastCount)
    numEdges += *iCount++;
    
  // allocate edge-to-vertex connectivity
  edgeToVertex.reserve(2*numEdges);
  
  for (int64_t iNode = 0; iNode < numVertex; ++iNode) {
  
    int64_t first = edgeBegin[iNode];
    int64_t last = first + vertexToVertexCount[iNode];
  
    for (int64_t j = first; j < last; ++j) {
      edgeToVertex.push_back(iNode);
      edgeToVertex.push_back(vertexToVertex[j]);
    }
  }
}

