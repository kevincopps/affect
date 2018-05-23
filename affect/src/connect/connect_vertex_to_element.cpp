#include <cstdint>

//#define AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT

#ifdef AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT
#include <iostream>
#include <iomanip>
#endif

#include <util/aligned_array.hpp>

//#define TIMERS_ACTIVE
#include <util/timer.hpp>

// NOTE: thread parallel version is not parallel consistent, i.e., the resulting element IDs in vertexToElement
// may end up in a different order on subsequent executions
//
void connect_vertex_to_element(
    uint32_t numVertex,
    uint32_t numVertexPerElement,
    uint32_t numElement,
    const aligned::uint32_ptr __restrict elementToVertex,
    aligned::uint32_ptr __restrict vertexToElementCount,
    aligned::uint32_ptr __restrict vertexToElement)
{
START_TIMER(connect_vertex_to_element);


    // Iterate over elements and copy the element to the vertices it is connected to.
    // During the loop we will keep incrementing by one our vertexToElementCount array which acts as a
    // pointer to the beginning of the entries for the ith vertex.
    //
    #pragma omp parallel for schedule(static, 1024) shared(vertexToElementCount)
    for (uint32_t iElem = 0; iElem < numElement; ++iElem) {

        const size_t offset = numVertexPerElement * iElem;

        for (uint32_t local_vertex = 0; local_vertex < numVertexPerElement; ++local_vertex) {
            size_t count;
            size_t vertex = elementToVertex[offset + local_vertex];
            #pragma omp atomic capture
            {
                count = vertexToElementCount[vertex];
                vertexToElementCount[vertex]++;
            }
            vertexToElement[count] = iElem;
        }
    }

    // reverse the incrementing we did in the previous loops
    for (uint32_t i = numVertex-1; i > 0; --i) {
        vertexToElementCount[i] = vertexToElementCount[i-1];
    }
    vertexToElementCount[0] = 0;

END_TIMER(connect_vertex_to_element);

#ifdef AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT
    std::cout << "----Debugging connect_vertex_to_element" << std::endl;
    //for (uint32_t i = 0; i <= numVertex; ++i )
    //  cout << "  " << setw(3) << i << setw(8) << vertexToElementCount[i] << endl;
    for (uint32_t i = 0; i < numVertex; ++i ) {
        std::cout << " vertex " << std::setw(3) << i << " has elem";
        for (uint32_t j = vertexToElementCount[i]; j < vertexToElementCount[i+1]; ++j)
            std::cout << " " << vertexToElement[j];
        std::cout << std::endl;
    }
#endif
  return;
}

// serial version is parallel thread consistent
void connect_vertex_to_element_consistent(
    uint32_t numVertex,
    uint32_t numVertexPerElement,
    uint32_t numElement,
    const aligned::uint32_ptr __restrict elementToVertex,
    aligned::uint32_ptr __restrict vertexToElementCount,
    aligned::uint32_ptr __restrict vertexToElement)
{
START_TIMER(connect_vertex_to_element_consistent);

    /* serial version */
    // Iterate over elements and copy the element to the vertices it is connected to.
    // During the loop we will keep incrementing by one our vertexToElementCount array which acts as a
    // pointer to the beginning of the entries for the ith vertex.
    //
    for (uint32_t iElem = 0, start = 0, finish = numVertexPerElement;
        iElem < numElement; ++iElem, start = finish, finish += numVertexPerElement) {
        for (uint32_t j = start; j < finish; ++j) {
            vertexToElement[vertexToElementCount[elementToVertex[j]]++] = iElem;
        }
    }

    // reverse the incrementing we did in the previous loops
    for (uint32_t i = numVertex-1; i > 0; --i) {
        vertexToElementCount[i] = vertexToElementCount[i-1];
    }
    vertexToElementCount[0] = 0;

END_TIMER(connect_vertex_to_element_consistent);

#ifdef AFFECT_DEBUG_CONNECT_VERTEX_TO_ELEMENT
    std::cout << "----Debugging connect_vertex_to_element_consistent" << std::endl;
    //for (uint32_t i = 0; i <= numVertex; ++i )
    //  cout << "  " << setw(3) << i << setw(8) << vertexToElementCount[i] << endl;
    for (uint32_t i = 0; i < numVertex; ++i ) {
        std::cout << " vertex " << std::setw(3) << i << " has elem";
        for (uint32_t j = vertexToElementCount[i]; j < vertexToElementCount[i+1]; ++j)
            std::cout << " " << vertexToElement[j];
        std::cout << std::endl;
    }
#endif
  return;
}


