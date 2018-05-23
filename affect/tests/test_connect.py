import sys
import numpy
import time

from .. import exodus
from .. import connect
from .. import util


def get_topology(block):
    topology = connect.CellTopology[block.topology_name]
    if topology == connect.CellTopology.PYRAMID5 or topology == connect.CellTopology.WEDGE6:
        print('    skipping {}'.format(topology.name))
        should_skip = True
    else:
        should_skip = False
    return topology, should_skip


def test_topology(edb_small):
    blocks = edb_small.element_blocks
    util.print_function_starting()
    print('{} element blocks'.format(len(blocks)))
    print('IDs = {}'.format([k for k in blocks.keys()]))
    for key, block in blocks.items():
        print('    block {} has {} entries'.format(block, block._num_entries))
        topology = connect.CellTopology[block.topology_name]
        print('        and {} maps to topology {} ({})'.format(block.topology_name, topology.name, topology.value))


def test_convert_to_local_connectivity(edb_small):
    util.print_function_starting()
    blocks = edb_small.element_blocks
    print('\n{} element blocks'.format(len(blocks)))
    max_node = edb_small.globals.num_nodes()
    for key, block in blocks.items():
        util.print_blue('    {} of {}'.format(block, block.topology_name))
        connectivity = block.connectivity()
        # util.print_array_info('connectivity', connectivity)
        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity(connectivity, max_node)
        util.print_array_info('element_to_vertex_local', element_to_vertex_local)
        util.print_array_info('local_to_global_vertex', local_to_global_vertex)


def test_vertex_to_element(edb_huge):
    e = edb_huge
    util.print_function_starting()
    blocks = e.element_blocks
    max_node = e.globals.num_nodes()
    print('\n{} element blocks'.format(len(blocks)))
    for key, block in blocks.items():
        if key != 36:
            continue
        util.print_blue('    {} of {}'.format(block, block.topology_name))
        connectivity = block.connectivity()
        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity(connectivity, max_node)
        max_vertex = len(local_to_global_vertex)
        # print('    {} unique local vertices'.format(max_vertex))
        # util.print_array_info('local_to_global_vertex', local_to_global_vertex)
        max_element_per_vertex, vertex_to_element_begin, vertex_to_element = connect.vertex_to_element(
            max_vertex, element_to_vertex_local)
        print('    max_element_per_vertex {}'.format(max_element_per_vertex))
        for v in range(10):
            i = vertex_to_element_begin[v]
            j = vertex_to_element_begin[v + 1]
            print('    vertex {} has elements {}'.format(v, vertex_to_element[i:j]))


def test_element_to_element(edb_huge):
    e = edb_huge
    util.print_function_starting()
    blocks = e.element_blocks
    max_node = e.globals.num_nodes()
    print('\n{} element blocks'.format(len(blocks)))
    for key, block in blocks.items():
        # if key != 36:
        #    continue;
        topology, should_skip = get_topology(block)
        if should_skip:
            continue
        util.print_blue('    {} of {}'.format(block, block.topology_name))
        # start_time = time.process_time()
        connectivity = block.connectivity()
        # print('    {:f} block.connectivity()'.format(time.process_time() - start_time))
        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity(connectivity, max_node)
        max_vertex = len(local_to_global_vertex)
        max_element_per_vertex, vertex_to_element_begin, vertex_to_element = connect.vertex_to_element(
            max_vertex, element_to_vertex_local)
        neighbor, num_boundary_quad_face, num_boundary_tri_face = connect.element_to_element(
            topology, max_element_per_vertex, element_to_vertex_local,
            vertex_to_element_begin, vertex_to_element)
        assert neighbor is not None
        assert num_boundary_quad_face + num_boundary_tri_face > 0
        # print('    has on boundary {} quads and {} tris'.format(num_boundary_quad_face, num_boundary_tri_face))
        # util.print_array_info('neighbor', neighbor)


def test_element_neighbors(edb_huge):
    e = edb_huge
    util.print_function_starting()
    blocks = e.element_blocks
    num_global_vertices = e.globals.num_nodes()
    print('\n{} element blocks'.format(len(blocks)))
    print('')
    global_to_local_buffer = util.empty_aligned(num_global_vertices, numpy.uint32)

    for key, block in blocks.items():
        # if key != 36:
        #    continue;

        topology, should_skip = get_topology(block)
        if should_skip:
            continue
        util.print_blue('    {} of {}'.format(block, block.topology_name))
        connectivity = block.connectivity()
        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity_buffer(
            connectivity, num_global_vertices, global_to_local_buffer)
        # util.print_array_info('element_to_vertex_local', element_to_vertex_local)
        # util.print_array_info('local_to_global_vertex', local_to_global_vertex)

        num_vertices = len(local_to_global_vertex)
        # print('    num_vertices = {}'.format(num_vertices))

        neighbor_elements, neighbor_faces = connect.element_neighbors(topology,
                                                                      num_vertices,
                                                                      element_to_vertex_local)
        # util.print_array_info('neighbor_elements', neighbor_elements)
        # util.print_array_info('neighbor_faces', neighbor_faces)
        assert neighbor_elements is not None
        assert neighbor_elements.size > 0
        assert neighbor_faces is not None
        assert neighbor_faces.size > 0


def test_connect_vertex_to_element_face(edb_small):
    e = edb_small
    util.print_function_starting()
    blocks = e.element_blocks
    num_global_vertices = e.globals.num_nodes()
    print('\n{} element blocks'.format(len(blocks)))

    global_to_local_buffer = util.empty_aligned(num_global_vertices, numpy.uint32)
    for key, block in blocks.items():
        if key not in (1, 4):
            continue
        topology, should_skip = get_topology(block)
        if should_skip:
            continue
        util.print_blue('{}'.format(block))
        connectivity = block.connectivity()
        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity_buffer(
            connectivity, num_global_vertices, global_to_local_buffer)
        # util.print_array_info('element_to_vertex_local', element_to_vertex_local)
        # util.print_array_info('local_to_global_vertex', local_to_global_vertex)

        num_vertices = len(local_to_global_vertex)
        # print('    num_vertices = {}'.format(num_vertices))

        neighbor_elements, neighbor_faces = connect.element_neighbors(topology, num_vertices, element_to_vertex_local)
        # util.print_array_info('neighbor_elements', neighbor_elements)
        # util.print_array_info('neighbor_faces', neighbor_faces)

        vertex_facet_element, vertex_facet_face = connect.connect_vertex_to_element_face(topology,
                                                                                         num_vertices,
                                                                                         element_to_vertex_local,
                                                                                         neighbor_faces)
        # check equality of the first ten vertices in some blocks edb_small
        if key == 1:
            numpy.testing.assert_array_equal(vertex_facet_element[:10],
                                             numpy.array([0, 9513, 145, 145, 1, 9105, 16559, 6, 21303, 21277]))
        if key == 4:
            numpy.testing.assert_array_equal(vertex_facet_element[:10],
                                             numpy.array([363, 26, 308, 363, 0, 0, 0, 0, 27, 309]))


def test_connect_vertex_to_element_face_performance(edb_huge):
    e = edb_huge
    util.print_function_starting()
    blocks = e.element_blocks
    num_global_vertices = e.globals.num_nodes()
    print('\n{} element blocks'.format(len(blocks)))

    global_to_local_buffer = util.empty_aligned(num_global_vertices, numpy.uint32)
    for key, block in blocks.items():
        if key != 36:
            continue
        topology, should_skip = get_topology(block)
        if should_skip:
            continue
        util.print_blue('{}'.format(block))

        connectivity = block.connectivity()

        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity_buffer(
            connectivity, num_global_vertices, global_to_local_buffer)

        num_vertices = len(local_to_global_vertex)

        neighbor_elements, neighbor_faces = connect.element_neighbors(topology, num_vertices, element_to_vertex_local)

        vertex_facet_element, vertex_facet_face = connect.connect_vertex_to_element_face(topology,
                                                                                         num_vertices,
                                                                                         element_to_vertex_local,
                                                                                         neighbor_faces)
        assert vertex_facet_element is not None
        assert vertex_facet_face is not None


def test_boundary_face_to_vertex(edb_huge):
    # e = edb_large
    # e = edb_display
    e = edb_huge
    # e = edb_huge_one_block
    util.print_function_starting()
    blocks = e.element_blocks
    max_node = e.globals.num_nodes()
    print('{} element blocks'.format(len(blocks)))
    for key, block in blocks.items():
        util.print_blue('{}'.format(block))
        connectivity = block.connectivity()
        topology = connect.CellTopology[block.topology_name]
        if topology == connect.CellTopology.PYRAMID5:
            print('    SKIP topology = {}'.format(topology.name))
            continue
        else:
            print('    topology = {}'.format(topology.name))
        element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity(connectivity, max_node)
        # util.print_array_info('local_to_global_vertex', local_to_global_vertex)
        min_global_node = local_to_global_vertex.min()
        max_global_node = local_to_global_vertex.max()
        print('    min global node = {} max global node= {}'.format(min_global_node, max_global_node))
        # if max_global_node > max_node:
        #     raise RuntimeError('max_global_node > max_node')
        max_vertex = local_to_global_vertex.size
        max_element_per_vertex, vertex_to_element_begin, vertex_to_element = connect.vertex_to_element(
            max_vertex, element_to_vertex_local)
        print('    max_element_per_vertex {}'.format(max_element_per_vertex))

        neighbor, num_boundary_quad_face, num_boundary_tri_face = connect.element_to_element(
            topology, max_element_per_vertex, element_to_vertex_local,
            vertex_to_element_begin, vertex_to_element)
        print('    has on boundary {} quads and {} tris'.format(num_boundary_quad_face, num_boundary_tri_face))
        # util.print_array_info('neighbor', neighbor)
        boundary_face_to_vertex = connect.boundary_face_to_vertex(topology,
                                                                  num_boundary_quad_face,
                                                                  num_boundary_tri_face,
                                                                  element_to_vertex_local,
                                                                  neighbor)
        # util.print_array_info('boundary_face_to_vertex', boundary_face_to_vertex)
        assert boundary_face_to_vertex is not None
