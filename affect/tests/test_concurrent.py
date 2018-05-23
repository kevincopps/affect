from .. import connect
from .. import util
from concurrent import futures

import time
import numpy


def get_boundary_face_to_vertex(key, connectivity, topology, max_node):
    element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity(connectivity, max_node)
    # util.util.print_array_info('local_to_global_vertex', local_to_global_vertex)
    # min_global_node = local_to_global_vertex.min()
    # max_global_node = local_to_global_vertex.max()
    # print('    min global node = {} max global node= {}'.format(min_global_node, max_global_node))
    # if max_global_node > max_node:
    #     raise RuntimeError('max_global_node > max_node')
    max_vertex = local_to_global_vertex.size
    max_element_per_vertex, vertex_to_element_begin, vertex_to_element = connect.vertex_to_element(
        max_vertex, element_to_vertex_local)
    # print('    max_element_per_vertex {}'.format(max_element_per_vertex))
    neighbor, num_boundary_quad_face, num_boundary_tri_face = connect.element_to_element(
        topology, max_element_per_vertex, element_to_vertex_local,
        vertex_to_element_begin, vertex_to_element)
    # print('    has on boundary {} quads and {} tris'.format(num_boundary_quad_face, num_boundary_tri_face))
    # util.util.print_array_info('neighbor', neighbor)
    boundary_face_to_vertex = connect.boundary_face_to_vertex(topology,
                                                              num_boundary_quad_face,
                                                              num_boundary_tri_face,
                                                              element_to_vertex_local,
                                                              neighbor)
    return key, boundary_face_to_vertex


def process_connectivity(key, connectivity, topology, max_node):
    # print('start')
    shape = connectivity.shape
    num_elements = shape[0]
    start_time = time.process_time()
    key, boundary_face_to_vertex = get_boundary_face_to_vertex(key, connectivity, topology, max_node)
    elapsed_time = time.process_time() - start_time
    return key, num_elements, elapsed_time, boundary_face_to_vertex


def process_connectivity_packed(key, data, topology, max_node):
    print('start')
    start_time = time.process_time()
    connectivity = data.unpack()
    shape = connectivity.shape
    num_elements = shape[0]
    key, boundary_face_to_vertex = get_boundary_face_to_vertex(key, connectivity, topology, max_node)
    elapsed_time = time.process_time() - start_time
    return key, num_elements, elapsed_time, boundary_face_to_vertex


def print_boundary(job):
    key, num_elements, elapsed_time, boundary_face_to_vertex = job.result()

    per_second = num_elements / elapsed_time / 1000000.0
    print(f'    block {key}: {num_elements} elements, {boundary_face_to_vertex.size} boundary vertices, '
          f'{elapsed_time:.4f}s, {per_second:.2f}M elements/s')
    # util.util.print_array_info('boundary_face_to_vertex', boundary_face_to_vertex)


def print_boundary_packed(job):
    key, num_elements, elapsed_time, boundary_face_to_vertex = job.result()

    per_second = num_elements / elapsed_time / 1000000.0
    print(f'    block {key}: {num_elements} elements, {boundary_face_to_vertex.size} boundary vertices, '
          f'{elapsed_time:.4f}s, {per_second:.2f}M elements/s')
    # util.util.print_array_info('boundary_face_to_vertex', boundary_face_to_vertex)


def get_element_neighbors(key, connectivity, topology, max_node):
    element_to_vertex_local, local_to_global_vertex = connect.convert_to_local_connectivity(connectivity, max_node)
    num_vertices = len(local_to_global_vertex)
    neighbor_elements, neighbor_faces = connect.element_neighbors(topology, num_vertices, element_to_vertex_local)
    return key, neighbor_elements, neighbor_faces


def process_element_neighbors_packed(key, data, topology, num_global_vertices):
    # print('start')
    start_time = time.process_time()
    connectivity = data.unpack()
    shape = connectivity.shape
    num_elements = shape[0]
    key, neighbor_elements, neighbor_faces = get_element_neighbors(key, connectivity, topology, num_global_vertices)
    elapsed_time = time.process_time() - start_time
    return key, num_elements, elapsed_time, neighbor_elements, neighbor_faces


def print_element_neighbors_packed(job):
    key, num_elements, elapsed_time, neighbor_elements, neighbor_faces = job.result()

    million_per_second = num_elements / elapsed_time / 1000000.0
    print('    block {}: {} elements, {} neighbor_elements, {:.4f}s, {:.2f}M elements/s'.format(
        key, num_elements, neighbor_elements.size, elapsed_time, million_per_second))


def blocks_in_size_order(blocks):

    def block_num_entries(key_value):
        b = key_value[1]  # the block object
        topology_name = b.topology_name
        # divide by factor of typical millions elements / second
        if topology_name.startswith('TET'):
            factor = 0.9
        elif topology_name.startswith('HEX'):
            factor = 1.75
        elif topology_name.startswith('WEDGE'):
            factor = 1.75
        else:
            factor = 1.00
        return b.num_entries() / factor

    sorted_list = sorted(blocks.items(), key=block_num_entries, reverse=True)
    for key, block in sorted_list:
        yield key, block


def get_topology(block):
    topology = connect.CellTopology[block.topology_name]
    if topology == connect.CellTopology.PYRAMID5 or topology == connect.CellTopology.WEDGE6:
        print('    skipping {}'.format(topology.name))
        should_skip = True
    else:
        should_skip = False
    return topology, should_skip


def test_boundary_face_to_vertex(edb_huge):

    # e = edb_large
    # e = edb_display
    e = edb_huge
    util.print_function_starting()

    blocks = e.element_blocks
    max_node = e.globals.num_nodes()
    print('{} element blocks'.format(len(blocks)))

    with futures.ProcessPoolExecutor(max_workers=len(blocks)) as executor:

        num_job_submitted = 0
        for key, block in blocks_in_size_order(blocks):

            topology, should_skip = get_topology(block)
            if should_skip:
                continue
            util.print_blue('{}'.format(block))
            connectivity = block.connectivity()
            job = executor.submit(process_connectivity, key, connectivity, topology, max_node)
            job.add_done_callback(print_boundary)
            # wait for the first few jobs to start before submitting again
            num_job_submitted += 1
            if num_job_submitted <= 4:
                while not job.running():
                    print('sleep')
                    time.sleep(0.01)


def test_boundary_face_to_vertex_packed(edb_huge):

    # e = edb_large
    # e = edb_display
    e = edb_huge
    util.print_function_starting()

    blocks = e.element_blocks
    max_node = e.globals.num_nodes()
    print('{} element blocks'.format(len(blocks)))

    # with futures.ProcessPoolExecutor(max_workers=len(blocks)) as executor:
    with futures.ProcessPoolExecutor(max_workers=8) as executor:

        for key, block in blocks_in_size_order(blocks):

            topology, should_skip = get_topology(block)
            if should_skip:
                continue
            util.print_blue('{}'.format(block))

            connectivity_compressed = block.connectivity(compress=True)

            job = executor.submit(process_connectivity_packed, key, connectivity_compressed, topology, max_node)
            job.add_done_callback(print_boundary_packed)


def test_element_neighbors_packed(edb_huge):

    # e = edb_large
    # e = edb_display
    e = edb_huge

    util.print_function_starting()

    blocks = e.element_blocks
    num_global_vertices = e.globals.num_nodes()
    print('\n{} element blocks'.format(len(blocks)))

    # with futures.ProcessPoolExecutor(max_workers=len(blocks)) as executor:
    with futures.ProcessPoolExecutor(max_workers=2) as executor:

        for key, block in blocks_in_size_order(blocks):

            topology, should_skip = get_topology(block)
            if should_skip:
                continue
            util.print_blue('{}'.format(block))

            # shape, dtype, data = pack_connectivity(block)
            data = block.connectivity(compress=True)

            job = executor.submit(process_element_neighbors_packed, key, data, topology, num_global_vertices)
            job.add_done_callback(print_element_neighbors_packed)

