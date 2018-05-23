from .. import exodus
from .. import connect
from .. import display
from .. import util
import os
import numpy
import time


def test_write_scene(edb_display):
    e = edb_display
    util.print_function_starting()
    util.print_bold('    {e.path}')
    blocks = e.element_blocks
    max_node = e.globals.num_nodes()
    global_coordinates = e.nodal.coordinates()

    # # find the largest block of TETs
    # max_num_entries = 0
    # max_block_key = -1
    # for key, block in blocks.items():
    #     topology = connect.CellTopology[block._topology_name.upper()]
    #     if topology == connect.CellTopology.TET4:
    #         num_entries = block.num_entries()
    #         if num_entries > max_num_entries:
    #             max_num_entries = num_entries
    #             max_block_key = key

    face_sets = list()
    all_local_coordinates = list()

    elapsed = 0.0

    items = blocks.connectivity_local_all()

    for key, block, local_connectivity in items:

        util.print_blue(f'    {block} with {local_connectivity.num_unique} nodes')

        topology = connect.CellTopology[block.topology_name]
        if topology == connect.CellTopology.TET4:

            local_nodes = local_connectivity.local_nodes

            # get the connectivity for this block and do the work to get the boundary face to vertex connectivity
            max_element_per_vertex, vertex_to_element_begin, vertex_to_element = connect.vertex_to_element(
                local_connectivity.num_unique, local_nodes)

            neighbor, num_boundary_quad_face, num_boundary_tri_face = connect.element_to_element(
                topology, max_element_per_vertex, local_nodes,
                vertex_to_element_begin, vertex_to_element)

            boundary_face_to_vertex = connect.boundary_face_to_vertex(topology,
                                                                      num_boundary_quad_face,
                                                                      num_boundary_tri_face,
                                                                      local_nodes,
                                                                      neighbor)
            # add to our lists
            face_to_vertex = boundary_face_to_vertex.view()
            face_to_vertex.shape = (num_boundary_tri_face, 3)

            start = time.perf_counter()

            # local_coordinates = global_coordinates[local_connectivity.global_nodes, ...]
            # util.print_array_info('alt_local_coordinates', alt_local_coordinates)

            # local_coordinates = numpy.take(global_coordinates, local_connectivity.global_nodes, axis=0)

            local_coordinates = util.take(global_coordinates, local_connectivity.global_nodes)

            # util.print_array_info('alt_local_coordinates', local_coordinates)
            # numpy.testing.assert_array_equal(local_coordinates, alt_local_coordinates)

            elapsed += time.perf_counter() - start

            face_sets.append((face_to_vertex, local_coordinates))
            all_local_coordinates.append(local_coordinates)

    bbox = display.bounding_box(all_local_coordinates)
    center = numpy.mean(bbox, axis=0)
    print(f'        center: {center}')

    # make an X3d scene file and dump the boundary faces into it
    directory = os.path.dirname(os.path.realpath(__file__))
    path_name = os.path.join(directory, 'test_write_scene')

    print(f'    time elapsed = {elapsed}')

    # display.write_scene(path_name, 'test_write_scene', 'test block boundaries displayed in X3D', face_sets)
