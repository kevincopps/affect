#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstdint>
#include <iostream>

//#define TIMERS_ACTIVE
#include <util/timer.hpp>
#include <util/sort_small.hpp>
#include <util/aligned_array.hpp>

#include "element_topology.hpp"
#include "print_array.hpp"

/*
  Create an array storing the beginning index to the array of connected faces for each vertex.
  Only anchor vertices will store connected faces.

  An anchor vertex is the maximum vertex ID on any element face.

  On return, for each i ∈ [0, num_vertex) the number of faces connected to the ith vertex is

  anchor_to_face_count[i+1] - anchor_to_face_count[i]

  Arguments:
    num_elements            - (IN) the number of elements in this block
    num_vertices            - (IN) the maximum vertex ID used by this block
    num_vertex_per_element  - (IN) the number of vertices per element in the element_to_vertex array
    num_face_per_element    - (IN) the number of faces per element in this element topology
    max_vertex_per_face     - (IN) the highest number of vertices per each element face in this element topology
    face_vertex_order       - (IN) the enumeration of local vertices for each local face in this element topology
    element_to_vertex       - (IN) array of vertices connected to the elements
    anchor_to_face_count    - (IN/OUT) on entry, an array of zeros of length (num_vertices + 1),
                                       on exit the starting index to a contiguous array of connected faces for each
                                       vertex.
  Return:
    num_faces               - the sum of all faces connected to all anchor vertices (which is also equal to
                              anchor_to_face_count[num_vertices])
 */
uint32_t count_incident_faces_on_anchor_vertices(
  const uint32_t num_elements,
  const uint32_t num_vertices,
  const uint32_t num_vertex_per_element,
  const uint32_t num_face_per_element,
  const uint32_t max_vertex_per_face,
  const uint32_t * face_vertex_order,
  const uint32_t * element_to_vertex,
  uint32_t * anchor_to_face_count)
{
    const uint32_t * element_vertex_end = element_to_vertex + num_elements * num_vertex_per_element;
    const uint32_t * face_vertex_order_end = face_vertex_order + num_face_per_element * max_vertex_per_face;

    // count the faces connected to vertex
    // assumes the anchor_to_face_count must be already initialized to zeros

    // for each element (and its vertices)
    #pragma omp parallel for schedule(static)
    for (const uint32_t * element_vertex_start = element_to_vertex;
                          element_vertex_start < element_vertex_end;
                    element_vertex_start += num_vertex_per_element) {

        // for each face connected to element (and its local vertex ordering)
        for (const uint32_t * face_vertex_order_start = face_vertex_order;
                              face_vertex_order_start < face_vertex_order_end; ) {
            uint32_t anchor_vertex = 0;
            for (uint32_t vertex = 0; vertex < max_vertex_per_face; ++vertex) {
                uint32_t face_vertex = element_vertex_start[ *face_vertex_order_start++ ];
                anchor_vertex = anchor_vertex < face_vertex ? face_vertex : anchor_vertex;
            }
            #pragma omp atomic update
            anchor_to_face_count[anchor_vertex]++;
        }
    }

    uint32_t num_half_facets = 0;

    // for each vertex, replace the count of faces sharing it
    // with its beginning index to the lists of its half-facets.
    for (uint32_t vertex = 0; vertex < num_vertices; ++vertex) {
        num_half_facets += anchor_to_face_count[vertex];
        anchor_to_face_count[vertex] = num_half_facets - anchor_to_face_count[vertex];
    }
    anchor_to_face_count[num_vertices] = num_half_facets;

    return num_half_facets;
}

/*
    Reverse the incrementing of the starting position for each anchor vertex once it has been used to increment
    and insert <element, local_face> and adjacent vertices in the correct position.

    On entering, the anchor vertices are the entries that are non-zero, and the value for the ith anchor vertex
    can be taken from the value for the (i-1)th anchor vertex.
 */
void reverse_anchor_to_face_count(
    uint32_t num_vertices,
    uint32_t * anchor_to_face_count) {

    uint32_t anchor_vertex, previous_anchor;

    // first find the largest vertex ID that is an anchor_vertex (the 0 index is a special case)
    for (anchor_vertex = num_vertices - 1; anchor_vertex > 0; --anchor_vertex)
        if (anchor_to_face_count[anchor_vertex] != 0)
            break;

    // start with last anchor_vertex and search backward (the 0 index is a special case)
    while (anchor_vertex > 0) {

        // search backward for the previous non-zero entry
        for (previous_anchor = anchor_vertex - 1; previous_anchor == 0; --previous_anchor) {
            if (anchor_to_face_count[previous_anchor] != 0)
                break;
        }
        if (anchor_to_face_count[previous_anchor] != 0) {
            anchor_to_face_count[anchor_vertex] = anchor_to_face_count[previous_anchor];
            // now we continue searching backward from this previous non-zero entry
            anchor_vertex = previous_anchor;
        }
        else {
            // there was no previous_anchor, we are at the last anchor vertex
            anchor_to_face_count[anchor_vertex] = 0;
            break;  // exit the while loop
        }
    }
    // In the case that the first anchor vertex is also vertex zero, there is no previous_anchor to search for.
    // Setting to zeroth entry to zero takes care of this case, and has no effect otherwise.
    anchor_to_face_count[0] = 0;
}

/*
    Create arrays of half-facets, <element, local_face> pairs, and their adjacent vertex IDs.

    The anchor_to_face_count array gives the starting address of half-facets connected to an anchor vertex.
 */
void create_half_facets_index(
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const uint32_t num_vertex_per_element,
    const uint32_t num_face_per_element,
    const uint32_t max_vertex_per_face,
    const uint32_t * face_vertex_order,
    const uint32_t * element_to_vertex,
    uint32_t * anchor_to_face_count,
    uint32_t * anchor_face_elements,
    uint32_t * anchor_face_vertices,
    int8_t * anchor_local_faces) {

#define LIMIT_MAX_FACE_VERTICES 4 // large enough for HEX and TET

    const uint32_t * vertex_begin;
    const uint32_t * local_face_vertex;
    uint32_t * anchor_adjacent_start;
    uint32_t element, anchor_vertex, half_facet;
    uint32_t face;
    uint32_t anchor_index, i, adjacent_vertex;
    uint32_t num_adjacent_vertices = max_vertex_per_face - 1;

    assert(max_vertex_per_face <= LIMIT_MAX_FACE_VERTICES);

    for (element = 0, vertex_begin = element_to_vertex;
         element < num_elements;
         ++element, vertex_begin += num_vertex_per_element) {

        for (face = 0, local_face_vertex = face_vertex_order; 
             face < num_face_per_element; 
             ++face, local_face_vertex += max_vertex_per_face) {

            // get the anchor_vertex and its index in the local face vertices
            anchor_index = 0;
            anchor_vertex = vertex_begin[local_face_vertex[anchor_index]];
            for (i = 1; i < max_vertex_per_face; ++i) {
                adjacent_vertex = vertex_begin[local_face_vertex[i]]; 
                if (adjacent_vertex > anchor_vertex) {
                    anchor_index = i;
                    anchor_vertex = adjacent_vertex;
                }
            }

            // position in which to place the half-facet
            half_facet = anchor_to_face_count[anchor_vertex];

            // copy the adjacent vertices other than the one at the anchor_index
            anchor_adjacent_start = anchor_face_vertices + half_facet * num_adjacent_vertices;
            for (i = 0; i < max_vertex_per_face; ++i) {
                if (i != anchor_index)
                    *anchor_adjacent_start++ = vertex_begin[ local_face_vertex[i] ];
            }

            // copy the <element, face> pair
            anchor_face_elements[half_facet] = element;
            anchor_local_faces[half_facet] = (int8_t)face;

            // increment the start position for the next half-facet on the anchor vertex
            anchor_to_face_count[anchor_vertex] += 1;
        }
    }

    // undo the incrementing of the anchor to face starting positions we performed in the previous loop
    reverse_anchor_to_face_count(num_vertices, anchor_to_face_count);
}


void create_half_facets_index_tet(
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const uint32_t * face_vertex_order,
    const uint32_t * element_to_vertex,
    uint32_t * anchor_to_face_count,
    uint32_t * anchor_face_elements,
    uint32_t * anchor_face_vertices,
    int8_t * anchor_local_faces) {

    #pragma omp parallel for schedule(static, 2048) shared(face_vertex_order, element_to_vertex, anchor_to_face_count, \
        anchor_face_elements, anchor_face_vertices, anchor_local_faces)
    for (uint32_t element = 0; element < num_elements; ++element) {

        const uint32_t * local_face_vertex = face_vertex_order;
        const uint32_t * vertex_begin = element_to_vertex + 4 * element;

        for (int8_t face = 0; face < 4; ++face) {

            uint32_t local_face_vertex0 = vertex_begin[*local_face_vertex++];
            uint32_t local_face_vertex1 = vertex_begin[*local_face_vertex++];
            uint32_t local_face_vertex2 = vertex_begin[*local_face_vertex++];

            // get the anchor_vertex and its index in the local face vertices
            uint32_t anchor_index = 0;
            uint32_t anchor_vertex = local_face_vertex0;
            if (local_face_vertex1 > anchor_vertex) {
                anchor_index = 1;
                anchor_vertex = local_face_vertex1;
            }
            if (local_face_vertex2 > anchor_vertex) {
                anchor_index = 2;
                anchor_vertex = local_face_vertex2;
            }

            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wunused-value"  // omp atomic capture seems to cause this warning

            // get position in which to place the half-facet
            // and increment the start position for the next half-facet on the anchor vertex
            uint32_t half_facet;
            #pragma omp atomic capture
            {
                half_facet = anchor_to_face_count[anchor_vertex];
                anchor_to_face_count[anchor_vertex]++;
            }

            // copy the adjacent vertices other than the one at the anchor_index
            uint32_t * v = &anchor_face_vertices[half_facet * 2];
            switch ( anchor_index )
            {
                case 0:
                    v[0] = local_face_vertex1;
                    v[1] = local_face_vertex2;
                    break;
                case 1:
                    v[0] = local_face_vertex0;
                    v[1] = local_face_vertex2;
                    break;
                case 2:
                    v[0] = local_face_vertex0;
                    v[1] = local_face_vertex1;
                    break;
                default:
                    ; // no-op
            }

            // copy the <element, face> pair
            anchor_face_elements[half_facet] = element;
            anchor_local_faces[half_facet] = face;

            #pragma GCC diagnostic pop  // "-Wunused-value"
        }
    }

    // undo the incrementing of the anchor to face starting positions we performed in the previous loop
    reverse_anchor_to_face_count(num_vertices, anchor_to_face_count);
}


void create_half_facets_index_hex(
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const uint32_t * face_vertex_order,
    const uint32_t * element_to_vertex,
    uint32_t * anchor_to_face_count,
    uint32_t * anchor_face_elements,
    uint32_t * anchor_face_vertices,
    int8_t * anchor_local_faces) {

    #pragma omp parallel for schedule(static) shared(face_vertex_order, element_to_vertex, anchor_to_face_count, \
        anchor_face_elements, anchor_face_vertices, anchor_local_faces)
    for (uint32_t element = 0; element < num_elements; ++element) {

        const uint32_t *vertex_begin = element_to_vertex + 8 * element;
        const uint32_t *local_face_vertex = face_vertex_order;
        uint32_t face_vertex0, face_vertex1, face_vertex2, face_vertex3;
        uint32_t anchor_index, anchor_vertex, half_facet;

        for (int8_t face = 0; face < 6; ++face) {

            face_vertex0 = vertex_begin[*local_face_vertex++];
            face_vertex1 = vertex_begin[*local_face_vertex++];
            face_vertex2 = vertex_begin[*local_face_vertex++];
            face_vertex3 = vertex_begin[*local_face_vertex++];

            // get the anchor_vertex and its index in the local face vertices
            uint32_t anchor_index = 0;
            uint32_t anchor_vertex = face_vertex0;
            if (face_vertex1 > anchor_vertex) {
                anchor_index = 1;
                anchor_vertex = face_vertex1;
            }
            if (face_vertex2 > anchor_vertex) {
                anchor_index = 2;
                anchor_vertex = face_vertex2;
            }
            if (face_vertex3 > anchor_vertex) {
                anchor_index = 3;
                anchor_vertex = face_vertex3;
            }

            #pragma GCC diagnostic push
            #pragma GCC diagnostic ignored "-Wunused-value"  // omp atomic capture seems to cause this warning

            // get position in which to place the half-facet
            // and increment the start position for the next half-facet on the anchor vertex
            #pragma omp atomic capture
            {
                half_facet = anchor_to_face_count[anchor_vertex];
                anchor_to_face_count[anchor_vertex]++;
            }

            // copy the adjacent vertices other than the one at the anchor_index
            uint32_t * v = &anchor_face_vertices[half_facet * 3];
            switch ( anchor_index )
            {
                case 0:
                    v[0] = face_vertex1;
                    v[1] = face_vertex2;
                    v[2] = face_vertex3;
                    break;
                case 1:
                    v[0] = face_vertex0;
                    v[1] = face_vertex2;
                    v[2] = face_vertex3;
                    break;
                case 2:
                    v[0] = face_vertex0;
                    v[1] = face_vertex1;
                    v[2] = face_vertex3;
                    break;
                case 3:
                    v[0] = face_vertex0;
                    v[1] = face_vertex1;
                    v[2] = face_vertex2;
                    break;
                default:
                    ; // no-op
            }

            // copy the <element, face> pair
            anchor_face_elements[half_facet] = element;
            anchor_local_faces[half_facet] = face;

            #pragma GCC diagnostic pop  // "-Wunused-value"
        }
    }

    // undo the incrementing of the anchor to face starting positions we performed in the previous loop
    reverse_anchor_to_face_count(num_vertices, anchor_to_face_count);
}

/*
 * For half-facet, sort the list of adjacent vertices.
 */
inline void sort_anchor_face_vertices( const uint32_t max_vertex_per_face,
                                const uint32_t num_half_facets,
                                uint32_t * anchor_face_vertices) {
    uint32_t * anchor_face_vertices_end = anchor_face_vertices + (max_vertex_per_face - 1) * num_half_facets;
    switch (max_vertex_per_face) {
        case 3:
            #pragma omp parallel for schedule(static, 1024)
            for (uint32_t * p = anchor_face_vertices; p < anchor_face_vertices_end; p += 2) {
                sort2(p);
            }
            break;
        case 4:
            #pragma omp parallel for schedule(static, 1024)
            for (uint32_t * p = anchor_face_vertices; p < anchor_face_vertices_end; p += 3) {
                sort3(p);
            }
            break;
        default:
            ; // no-op
    }
}

template <int NUM_ADJACENT_VERTICES>
inline bool equal_adjacent_vertices(const uint32_t * anchor_face_vertices,
                                    uint32_t outer_half_facet, uint32_t inner_half_facet);

template<>
inline bool equal_adjacent_vertices<2>(const uint32_t * anchor_face_vertices,
                                       uint32_t outer_half_facet, uint32_t inner_half_facet) {
    const uint32_t * outer_adjacent = &anchor_face_vertices[2 * outer_half_facet];
    const uint32_t * inner_adjacent = &anchor_face_vertices[2 * inner_half_facet];
    if (*outer_adjacent++ != *inner_adjacent++)
        return false;
    if (*outer_adjacent != *inner_adjacent)
        return false;
    return true;
}

template<>
inline bool equal_adjacent_vertices<3>(const uint32_t * anchor_face_vertices,
                                       uint32_t outer_half_facet, uint32_t inner_half_facet) {
    const uint32_t * outer_adjacent = &anchor_face_vertices[3 * outer_half_facet];
    const uint32_t * inner_adjacent = &anchor_face_vertices[3 * inner_half_facet];
    if (*outer_adjacent++ != *inner_adjacent++)
        return false;
    if (*outer_adjacent++ != *inner_adjacent++)
        return false;
    if (*outer_adjacent != *inner_adjacent)
        return false;
    return true;
}


template <int NUM_ADJACENT_VERTICES>
void connect_neighbors_from_anchor_vertices(
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const uint32_t num_face_per_element,
    const uint32_t * face_vertex_order,
    const uint32_t * element_to_vertex,
    const uint32_t * anchor_to_face_count,
    uint32_t * anchor_face_elements,
    uint32_t * anchor_face_vertices,
    int8_t * anchor_local_faces,
    uint32_t * neighbor_elements,
    int8_t * neighbor_faces) {

    // initialize all neighbor faces to be unknown
    aligned::fill(neighbor_faces, num_elements * num_face_per_element, static_cast<int8_t>(-1));

    // determine first anchor vertex ID, the lowest vertex ID that has half-facets
    uint32_t k;
    for (k = 0; k <= num_vertices; ++k)
        if (anchor_to_face_count[k])  // find the first non-zero
            break;

    #pragma omp parallel for schedule(guided, 256)
    for (uint32_t anchor_vertex = k - 1; anchor_vertex < num_vertices; ++anchor_vertex) {

        uint32_t outer_half_facet = anchor_to_face_count[anchor_vertex];

        // skip this candidate anchor vertex, if we do not have any new half-facets since the previous one
        if (anchor_vertex > 0 && anchor_to_face_count[anchor_vertex-1] == outer_half_facet && outer_half_facet > 0)
            continue;

        // find the next anchor vertex, and ending position of the list of half-faces for this anchor_vertex
        uint32_t anchor_vertex_next;
        for (anchor_vertex_next = anchor_vertex + 1; anchor_vertex_next <= num_vertices; ++anchor_vertex_next)
            if (anchor_to_face_count[anchor_vertex_next] > outer_half_facet)
                break;

        uint32_t inner_half_facet_end = anchor_to_face_count[anchor_vertex_next];
        uint32_t outer_half_facet_end = inner_half_facet_end - 1;

        // outer loop
        for ( ; outer_half_facet < outer_half_facet_end; ++outer_half_facet) {

            uint32_t outer_element = anchor_face_elements[outer_half_facet];
            int8_t outer_face    = anchor_local_faces[outer_half_facet];
            uint32_t outer_face_index = num_face_per_element * outer_element + outer_face;

            // skip this half-facet if we have previously (in the inner loop)
            // set the <neighbor, face> of this <element, face>
            // TODO: handle non-manifold half-faces (more than two allowed)
            if (neighbor_faces[outer_face_index] != -1)
                continue;

            // inner loop
            for (uint32_t inner_half_facet = outer_half_facet + 1;
                 inner_half_facet < inner_half_facet_end;
                 ++inner_half_facet) {

                uint32_t inner_element = anchor_face_elements[inner_half_facet];

                // skip if this is a face on the same element
                if (inner_element == outer_element)
                    continue;

                int8_t inner_face = anchor_local_faces[inner_half_facet];
                uint32_t inner_face_index = num_face_per_element * inner_element + inner_face;

                // skip this half-facet if we have previously (in the inner loop)
                // set the <neighbor, face> of this <element, face>
                // TODO: handle non-manifold half-faces (more than two allowed)
                if (neighbor_faces[inner_face_index] != -1)
                    continue;

                // if adjacent vertices match exchange <outer_element, outer_face> with <inner_element, inner_face>
                if (equal_adjacent_vertices<NUM_ADJACENT_VERTICES>(anchor_face_vertices,
                                                                   outer_half_facet,
                                                                   inner_half_facet)) {
                    neighbor_elements[outer_face_index] = inner_element;
                    neighbor_elements[inner_face_index] = outer_element;
                    neighbor_faces[outer_face_index] = inner_face;
                    neighbor_faces[inner_face_index] = outer_face;
                }
            }
        }
    }
}



/*
    Construct array of neighbor element IDs and array of local neighbor face IDs for each element face using the sibling
    half-facet algorithm.

    Input Args:
        topology: type of cells in this block
        element_to_vertex(|uint32_2d|): array of shape(num_elements, num_vertex_per_element)

    Output Args:
        neighbor_elements: array of neighbor element ID for each element face
        neighbor_faces: array of neighbor local face ID for each element face

    Returns:
        num_boundary_faces: number of element faces on the boundary
*/
void connect_element_neighbors(
    int topology,
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const uint32_t * element_to_vertex,
    uint32_t * neighbor_elements,
    int8_t * neighbor_faces) {

    START_TIMER(element_neighbors);

    // get topology dependent information
    const uint32_t * num_vertex_per_face = vertex_per_face[topology];
    const uint32_t * _face_vertex_order = face_vertex_order[topology];
    const uint32_t _max_vertex_per_face = max_vertex_per_face[topology];
    const uint32_t num_vertex_per_element = num_vertex[topology];
    const int8_t num_face_per_element = num_face[topology];

    // TODO: fix assumptions that do not apply for wedge and pyramid
    if (topology != HEX8 && topology != TET4) {
        std::cerr << "element_neighbors: not implemented for topology of elements other than HEX and TET" << std::endl;
        exit(1);
    }

    // print_2d_array("element_to_vertex", element_to_vertex, num_elements, num_vertex_per_element);
    //std::cerr << std::endl << "count_incident_faces_on_anchor_vertices" << std::endl;

    // allocate space for number of faces connected to anchor vertices
    uint32_t * anchor_to_face_count = aligned::allocate<uint32_t>(num_vertices + 1);
    aligned::zero(anchor_to_face_count, num_vertices + 1);

    START_TIMER(count_incident_faces_on_anchor_vertices);

    uint32_t num_half_facets = count_incident_faces_on_anchor_vertices( num_elements,
                                                                        num_vertices,
                                                                        num_vertex_per_element,
                                                                        num_face_per_element,
                                                                        _max_vertex_per_face,
                                                                        _face_vertex_order,
                                                                        element_to_vertex,
                                                                        anchor_to_face_count);

    END_TIMER(count_incident_faces_on_anchor_vertices);

    //std::cout << "num_half_facets " << num_half_facets << std::endl;
    //print_1d_array("anchor_to_face_count", anchor_to_face_count, num_vertices+1);

    // allocate enough space for the half facets of all anchor vertices, and their adjacent vertices
    uint32_t num_adjacent_vertex = _max_vertex_per_face - 1;

    uint32_t * anchor_face_vertices = aligned::allocate<uint32_t>(num_half_facets * num_adjacent_vertex);
    uint32_t * anchor_face_elements = aligned::allocate<uint32_t>(num_half_facets);
    int8_t * anchor_local_faces = aligned::allocate<int8_t>(num_half_facets);

    aligned::zero(anchor_face_elements, num_half_facets);

    START_TIMER(create_half_facets);

    if (topology == HEX8) {
        create_half_facets_index_hex(
            num_elements,
            num_vertices,
            _face_vertex_order,
            element_to_vertex,
            anchor_to_face_count,
            anchor_face_elements,
            anchor_face_vertices,
            anchor_local_faces);
    }
    else if (topology == TET4) {
        create_half_facets_index_tet(
            num_elements,
            num_vertices,
            _face_vertex_order,
            element_to_vertex,
            anchor_to_face_count,
            anchor_face_elements,
            anchor_face_vertices,
            anchor_local_faces);
    }
    else {
        create_half_facets_index( num_elements,
                        num_vertices,
                        num_vertex_per_element,
                        num_face_per_element,
                        _max_vertex_per_face,
                        _face_vertex_order,
                        element_to_vertex,
                        anchor_to_face_count,
                        anchor_face_elements,
                        anchor_face_vertices,
                        anchor_local_faces);
    }

    END_TIMER(create_half_facets);

    //print_1d_array("anchor_local_faces", anchor_local_faces, num_half_facets);
    //print_1d_array("anchor_face_elements", anchor_face_elements, num_half_facets);
    //print_1d_array("anchor_face_vertices", anchor_face_vertices, num_half_facets * num_adjacent_vertex);

    START_TIMER(sort_anchor_face_vertices);

    sort_anchor_face_vertices(_max_vertex_per_face,
                              num_half_facets,
                              anchor_face_vertices);

    END_TIMER(sort_anchor_face_vertices);

    //print_1d_array("anchor_local_faces", anchor_local_faces, num_half_facets);
    //print_1d_array("anchor_face_elements", anchor_face_elements, num_half_facets);
    //print_1d_array("anchor_face_vertices", anchor_face_vertices, num_half_facets * num_adjacent_vertex);

    START_TIMER(connect_neighbors_from_anchor_vertices);

    if (topology == HEX8) {
        connect_neighbors_from_anchor_vertices<3>( num_elements,
                                                   num_vertices,
                                                   num_face_per_element,
                                                   _face_vertex_order,
                                                   element_to_vertex,
                                                   anchor_to_face_count,
                                                   anchor_face_elements,
                                                   anchor_face_vertices,
                                                   anchor_local_faces,
                                                   neighbor_elements,
                                                   neighbor_faces);
    }
    else if (topology == TET4) {
        connect_neighbors_from_anchor_vertices<2>( num_elements,
                                                   num_vertices,
                                                   num_face_per_element,
                                                   _face_vertex_order,
                                                   element_to_vertex,
                                                   anchor_to_face_count,
                                                   anchor_face_elements,
                                                   anchor_face_vertices,
                                                   anchor_local_faces,
                                                   neighbor_elements,
                                                   neighbor_faces);
    }
    else {
        // TODO: fix assumptions that do not apply for wedge and pyramid
        std::cerr << "element_neighbors: not implemented for topology of elements other than HEX and TET" << std::endl;
        exit(1);
    }

    END_TIMER(connect_neighbors_from_anchor_vertices);

    free(anchor_local_faces);
    free(anchor_face_elements);
    free(anchor_face_vertices);
    free(anchor_to_face_count);

    END_TIMER(element_neighbors);
}
