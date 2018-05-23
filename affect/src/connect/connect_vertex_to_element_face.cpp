#include <cstdint>
#include <util/aligned_array.hpp>

#include "element_topology.hpp"

#define TIMERS_ACTIVE
#include <util/timer.hpp>

/*
Construction of mapping from each vertex to one incident half-facet.

Storage of this (and the given half-facets) allows any adjacency query on the mesh in constant time.

Input: elems: element connectivity
sibhfs: cyclic mappings of sibling half-facets
Output: v2hf: vertex to an incident half-facet

for each elements e in elems do
    for each vertex v of e do
        if v2hf(v) == 0 then
            Set v2hf(v) to (e, f) first facet incident on v in e
        end if
    end for

    // Give border facets higher priorities

    for each facet f in e do
        if sibhfs(e,f) == 0 then
            for each vertex of f do
                Set v2hf(v) to ⟨e, f ⟩
            end for
        end if
    end for
end for

*/


void connect_vertex_to_element_face(
    int topology,
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const aligned::uint32_ptr element_to_vertex,
    const aligned::int8_ptr neighbor_faces,
    aligned::uint32_ptr vertex_facet_element,
    aligned::int8_ptr vertex_facet_face) {

    const uint32_t * const _face_vertex_order = face_vertex_order[topology];
    const uint32_t * const _vertex_last_face = vertex_last_face[topology];
    const uint32_t num_face_per_element = num_face[topology];
    const uint32_t num_vertex_per_element = num_vertex[topology];
    const uint32_t _max_vertex_per_face = max_vertex_per_face[topology];

    START_TIMER(connect_vertex_to_element_face);

    aligned::zero(vertex_facet_face, num_vertices);

    const uint32_t * element_vertex = element_to_vertex;  // start of the element vertices we will increment

    // for each element
    for (size_t element = 0; element < num_elements; ++element) {

        // for each vertex connected to element
        for (uint32_t local_element_vertex = 0; local_element_vertex < num_vertex_per_element; ++local_element_vertex) {

            const uint32_t global_vertex = *element_vertex++;

            // Check if we have a local face stored for this global vertex yet.
            // We use the last local element face connected to the vertex so that if it has been set,
            // it will be non-zero. This is not guaranteed after the end of the next phase below
            // once we overwrite with border facets.
            if (vertex_facet_face[global_vertex] == 0) {

                // store the last incident half facet (element, face)
                vertex_facet_element[global_vertex] = element;
                vertex_facet_face[global_vertex] = _vertex_last_face[local_element_vertex];
            }
        }
    }

    // Give border facets higher priority.
    // For all vertices with an incident border facet,
    // overwrite their connected facet with the border facet.

    const size_t num_facets = num_elements * num_face_per_element;

    //#pragma omp parallel for schedule(guided, 256)   // this makes the results non-deterministic
    for (size_t facet = 0; facet < num_facets; ++facet) {

        if (neighbor_faces[facet] == -1) {

            size_t element = facet / num_face_per_element;
            size_t local_face = facet % num_face_per_element;

            const size_t face_vertices_offset = _max_vertex_per_face * local_face;
            const size_t element_vertices_offset = num_vertex_per_element * element;

            for (uint32_t local_face_vertex = 0; local_face_vertex < _max_vertex_per_face; ++local_face_vertex) {

                const uint32_t global_vertex = element_to_vertex[element_vertices_offset + _face_vertex_order[face_vertices_offset + local_face_vertex]];

                // store the last incident half facet (element, face)

                #pragma omp atomic write
                vertex_facet_element[global_vertex] = element;

                #pragma omp atomic write
                vertex_facet_face[global_vertex] = local_face;
            }
        }
    }
    END_TIMER(connect_vertex_to_element_face);
}

/* version with alignment and restrict pointer aliases
void connect_vertex_to_element_face(
    int topology,
    const uint32_t num_elements,
    const uint32_t num_vertices,
    const aligned::uint32_ptr __restrict__ element_to_vertex,
    const aligned::int8_ptr __restrict__ neighbor_faces,
    aligned::uint32_ptr __restrict__ vertex_facet_element,
    aligned::int8_ptr __restrict__ vertex_facet_face) {

    const uint32_t * __restrict__ _face_vertex_order = static_cast<uint32_t*>(__builtin_assume_aligned(static_cast<const void* const>(face_vertex_order[topology]), AFFECT_UTIL_ALIGN));
    const uint32_t * __restrict__ _vertex_last_face = static_cast<uint32_t*>(__builtin_assume_aligned(static_cast<const void* const>(vertex_last_face[topology]), AFFECT_UTIL_ALIGN));
    const uint32_t num_face_per_element = num_face[topology];
    const uint32_t num_vertex_per_element = num_vertex[topology];
    const uint32_t _max_vertex_per_face = max_vertex_per_face[topology];

    START_TIMER(connect_vertex_to_element_face);

    aligned::zero(vertex_facet_face, num_vertices);
    
    // for each element
    for (size_t element = 0; element < num_elements; ++element) {

        const size_t element_vertices_offset = num_vertex_per_element * element;

        // for each vertex connected to element
        for (uint32_t local_element_vertex = 0; local_element_vertex < num_vertex_per_element; ++local_element_vertex) {

            const uint32_t global_vertex = element_to_vertex[element_vertices_offset + local_element_vertex];

            // Check if we have a local face stored for this global vertex yet.
            // We use the last local element face connected to the vertex so that if it has been set,
            // it will be non-zero. This is not guaranteed after the end of the next phase below
            // once we overwrite with border facets.
            if (vertex_facet_face[global_vertex] == 0) {

                // store the last incident half facet (element, face)
                vertex_facet_element[global_vertex] = element;
                vertex_facet_face[global_vertex] = _vertex_last_face[local_element_vertex];
            }
        }
    }

    // Give border facets higher priority.
    // For all vertices with an incident border facet, overwrite their connected facet with the border facet.
    const size_t num_facets = num_elements * num_face_per_element;

    #pragma omp parallel for schedule(guided, 256)
    for (size_t facet = 0; facet < num_facets; ++facet) {

        if (neighbor_faces[facet] == -1) {

            size_t element = facet / num_face_per_element;
            size_t local_face = facet % num_face_per_element;

            const size_t face_vertices_offset = _max_vertex_per_face * local_face;
            const size_t element_vertices_offset = num_vertex_per_element * element;

            for (uint32_t local_face_vertex = 0; local_face_vertex < _max_vertex_per_face; ++local_face_vertex) {

                const uint32_t global_vertex = element_to_vertex[element_vertices_offset + _face_vertex_order[face_vertices_offset + local_face_vertex]];

                // store the last incident half facet (element, face)

                #pragma omp atomic write
                vertex_facet_element[global_vertex] = element;

                #pragma omp atomic write
                vertex_facet_face[global_vertex] = local_face;
            }
        }
    }
    END_TIMER(connect_vertex_to_element_face);
}
*/