from .. import connect
from .. import arithmetic
from .. import util

from concurrent import futures
import time
import numpy


def test_average_element_vertex_global(edb_huge):
    util.print_function_starting()
    e = edb_huge
    blocks = e.element_blocks

    util.print_bold(f'\n{len(blocks)} element blocks')
    # compute an array holding the element centroids
    global_node_coordinates = e.nodal.coordinates()
    for key, block in blocks.items():
        connectivity = block.connectivity()
        util.print_blue(f'    {block}')
        topology = connect.CellTopology[block.topology_name]
        print(f'    {topology.name} with {topology.num_vertex} vertices')
        centroids = arithmetic.average_element_node_values(connectivity, global_node_coordinates, topology.num_vertex)
        print(centroids)


def test_average_element_vertex_local(edb_huge):
    util.print_function_starting()
    e = edb_huge
    blocks = e.element_blocks
    util.print_bold(f'\n{len(blocks)} element blocks')
    # compute an array holding the element centroids
    global_node_coordinates = e.nodal.coordinates()

    # fields = e.nodal.fields
    # util.print_bold('\n{} fields'.format(len(fields)))
    # for field in fields.values():
    #     print('    {}'.format(field))

    for key, block, local in blocks.connectivity_local_all():
        util.print_blue(f'    {block} on {local.num_unique} nodes')
        topology = connect.CellTopology[block.topology_name]
        local_coordinates = util.take(global_node_coordinates, local.global_nodes)
        # util.print_array_info('local_coordinates', local_coordinates)
        centroids = arithmetic.average_element_node_values(local.local_nodes, local_coordinates, topology.num_vertex)
        print(centroids)


def process_average_element_vertex_local(key, topology, local, local_coordinates):
    start_time = time.process_time()
    centroids = arithmetic.average_element_node_values(local.local_nodes, local_coordinates, topology.num_vertex)
    elapsed_time = time.process_time() - start_time
    return key, topology, local.num_unique, elapsed_time, centroids


def print_average_element_vertex_local(job):
    key, topology, num_unique, elapsed_time, centroids = job.result()
    per_second = num_unique / elapsed_time / 1000000.0
    print(f'    block {key:2} {topology.name:>8} {num_unique:7} nodes in {elapsed_time:.4f}s {per_second:5.2f}M nodes/s')


def test_average_element_vertex_local_multiprocess(edb_huge):

    util.print_function_starting()
    e = edb_huge

    blocks = e.element_blocks
    print(f'\n{len(blocks)} element blocks')
    global_node_coordinates = e.nodal.coordinates()

    # with futures.ProcessPoolExecutor(max_workers=len(blocks)) as executor:
    with futures.ProcessPoolExecutor(max_workers=8) as executor:
        for key, block, local in blocks.connectivity_local_all():
            topology = connect.CellTopology[block.topology_name]
            local_coordinates = util.take(global_node_coordinates, local.global_nodes)
            job = executor.submit(process_average_element_vertex_local, key, topology, local, local_coordinates)
            job.add_done_callback(print_average_element_vertex_local)

