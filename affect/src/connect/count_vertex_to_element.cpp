#include <util/aligned_array.hpp>

//#define TIMERS_ACTIVE
#include <util/timer.hpp>

uint32_t count_vertex_to_element(
    uint32_t numVertex,
    uint32_t numVertexPerElement,
    uint32_t numElement,
    const aligned::uint32_ptr __restrict elementToVertex,
    aligned::uint32_ptr __restrict vertexToElementCount)
{
START_TIMER(count_vertex_to_element);

    size_t elementToVertexLength = numElement * numVertexPerElement;
    uint32_t maxElementPerVertex;

    aligned::zero(vertexToElementCount, numVertex + 1);

    // Count number of elements sharing each vertex.
    #pragma omp parallel for schedule(static) shared(vertexToElementCount)
    for (size_t k = 0; k < elementToVertexLength; ++k) {
        #pragma omp atomic update
        vertexToElementCount[ elementToVertex[k] ]++;
    }

    // calculate the maximum count of element sharing a vertex
    maxElementPerVertex = aligned::max(vertexToElementCount, numVertex);

    // Now for each vertex, replace the count of elements sharing it
    // with its beginning index to the lists of its elements.

    uint32_t length = 0, elementsPerNode;

    for (uint32_t i = 0; i < numVertex; ++i) {
        /*
        elementsPerNode = vertexToElementCount[i];
        vertexToElementCount[i] = length;
        length += elementsPerNode;
        */
        length += vertexToElementCount[i];
        vertexToElementCount[i] = length - vertexToElementCount[i];
    }
    vertexToElementCount[numVertex] = length;

END_TIMER(count_vertex_to_element);

    return maxElementPerVertex;
}

