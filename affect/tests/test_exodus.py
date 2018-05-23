# test_database.py

from .. import exodus as ex
from .. import util
import sys
import numpy as np
import textwrap
import pickle
import pytest


def get_global_connectivity_from_local(local_connectivity):
    # construct the global connectivity by indexing the local node connectivity into the global node mapping
    local_nodes = local_connectivity.local_nodes
    global_nodes = local_connectivity.global_nodes
    global_connectivity = np.empty_like(local_nodes, dtype=global_nodes.dtype)
    global_connectivity[:, :] = global_nodes[local_nodes[:, :]]
    return global_connectivity


def test_library_version():
    util.print_function_starting()
    print('    ExodusII library API version: {}'.format(ex.library_version()))


def test_non_existent_file():
    util.print_function_starting()
    with pytest.raises(ex.FileNotFound):
        ex.Database("nonexistent.e", ex.Mode.READ_ONLY)


def test_str_repr(edb):
    util.print_function_starting()
    print('    ' + str(edb))


def test_coordinates(edb_huge_one_block):
    util.print_function_starting()
    x = edb_huge_one_block.nodal.coordinates()
    print('    type {} with shape {}'.format(type(x), x.shape))
    

def test_bounding_box(edb):
    util.print_function_starting()
    ndim = edb.globals.dimension()
    coordinates = edb.nodal.coordinates()

    # min_x = np.finfo(np.float64).max
    # min_y = np.finfo(np.float64).max
    # min_z = np.finfo(np.float64).max
    # max_x = np.finfo(np.float64).min
    # max_y = np.finfo(np.float64).min
    # max_z = np.finfo(np.float64).min
    # for point in coordinates:
    #     min_x = min(point[0], min_x)
    #     min_y = min(point[1], min_y)
    #     min_z = min(point[2], min_z)
    #     max_x = max(point[0], max_x)
    #     max_y = max(point[1], max_y)
    #     max_z = max(point[2], max_z)
    # print(np.array([(min_x, min_y, min_z), (max_x, max_y, max_z)]))
    bbox = 0
    if ndim == 2:
        min_x, min_y = np.min(coordinates, axis=0)
        max_x, max_y = np.max(coordinates, axis=0)
        bbox = np.array([(min_x, min_y), (max_x, max_y)])
    elif ndim == 3:
        min_x, min_y, min_z = np.min(coordinates, axis=0)
        max_x, max_y, max_z = np.max(coordinates, axis=0)
        bbox = np.array([(min_x, min_y, min_z), (max_x, max_y, max_z)])
    print(bbox)


def test_element_blocks(edb):
    util.print_function_starting()
    print('    {} element blocks'.format(len(edb.element_blocks)))
    print('    IDs:')
    for line in textwrap.wrap('{}'.format([k for k in edb.element_blocks.keys()])):
        print('        ' + line)
    print('    blocks:')
    for key, block in edb.element_blocks.items():
        print('    {}: {}'.format(key, block))


def test_attributes(edb):
    util.print_function_starting()
    for key, block in edb.element_blocks.items():
        num_attributes = block.num_attributes()
        print('    {}: has {} attributes {}'.format(key, num_attributes, block.attribute_names()))
        if num_attributes:
            print('        {}'.format([block.attribute(k) for k in range(num_attributes)]))


def test_element_block_key_error(edb):
    util.print_function_starting()
    bad_id = -9999999
    element_blocks = edb.element_blocks
    with pytest.raises(ex.EntityKeyError):
        block = element_blocks[bad_id]
        print(block)
    print('    no block with id = {}'.format(bad_id))


def test_block_connectivity(edb):
    util.print_function_starting()
    min_entry = np.iinfo(np.int32).max
    max_entry = -1
    for key, block in edb.element_blocks.items():
        util.print_blue('    block {} connectivity = '.format(key))
        connectivity = block.connectivity()
        print('    connectivity.shape = {}'.format(connectivity.shape), flush=True)
        # print(connectivity)
        block_min_node = np.amin(connectivity)
        block_max_node = np.amax(connectivity)
        print('    block {:5d} first node {:7d} last node {:7d}'.format(key, block_min_node, block_max_node),
              flush=True)
        min_entry = min(min_entry, block_min_node)
        max_entry = max(max_entry, block_max_node)
        partial_connectivity = block.connectivity_partial(0, 1)  # all nodes of 0th element
        np.testing.assert_equal(partial_connectivity, connectivity[0:1, :])
    assert min_entry == 0
    num_nodes = edb.globals.num_nodes()
    assert max_entry == num_nodes - 1


def test_block_connectivity_compression(edb):
    util.print_function_starting()
    for key, block in edb.element_blocks.items():
        x = block.connectivity()
        y_buffer = block.connectivity(compress=True)
        print('    compressed connectivity ({}) buffer'.format(type(y_buffer)))
        y = y_buffer.unpack()
        np.testing.assert_array_equal(x, y)
        print('    block {} connectivity ({} bytes) equal to buffer ({} bytes)'.format(key, x.nbytes, len(y_buffer)))
        break


def test_block_connectivity_local(edb):
    e = edb

    util.print_function_starting()
    for key, block in e.element_blocks.items():
        local_connectivity = block.connectivity_local()  # the LocalConnectivity object
        util.print_blue('    block {} unique {}'.format(block, local_connectivity.num_unique))

        # basic checks on the local and global node IDs
        block_min_node = np.amin(local_connectivity.local_nodes)
        block_max_node = np.amax(local_connectivity.local_nodes)
        assert block_min_node == 0
        assert block_max_node == local_connectivity.global_nodes.size - 1
        global_min_node = np.amin(local_connectivity.global_nodes)
        global_max_node = np.amax(local_connectivity.global_nodes)
        assert global_min_node == local_connectivity.min_global
        assert global_max_node == local_connectivity.max_global
        assert global_max_node <= (e.globals.num_nodes() - 1)

        # read the global connectivity from file and check the constructed version is identical
        global_connectivity = get_global_connectivity_from_local(local_connectivity)
        connectivity = block.connectivity()
        np.testing.assert_equal(global_connectivity, connectivity)
        break


def test_block_connectivity_local_compression(edb):
    util.print_function_starting()
    for key, block in edb.element_blocks.items():
        local_connectivity = block.connectivity_local(compress=True)  # the LocalConnectivity object
        print('    block {} compressed connectivity buffer'.format(key))
        global_connectivity = get_global_connectivity_from_local(local_connectivity)
        # read the global connectivity from file and check the constructed version is identical
        connectivity = block.connectivity()
        np.testing.assert_equal(global_connectivity, connectivity)
        break


def test_block_connectivity_local_compression_pickleable(edb):
    util.print_function_starting()
    for key, block in edb.element_blocks.items():
        c = block.connectivity_local(compress=True)  # the LocalConnectivity object
        pickle_string = pickle.dumps(c)
        size_pickle_string = sys.getsizeof(pickle_string)
        print('    sizeof pickle string = {:.0f}k bytes'.format(size_pickle_string / 1024))
        local_connectivity = pickle.loads(pickle_string)
        assert local_connectivity.min_global == c.min_global
        assert local_connectivity.max_global == c.max_global
        assert local_connectivity.num_unique == c.num_unique
        local_nodes = local_connectivity.local_nodes
        global_nodes = local_connectivity.global_nodes
        size_arrays = local_nodes.nbytes + global_nodes.nbytes
        print('    sizeof arrays = {:.0f}k bytes'.format(size_arrays / 1024))
        print('    size ratio {:5.3f}'.format(size_pickle_string / size_arrays))
        np.testing.assert_equal(local_nodes, c.local_nodes)
        np.testing.assert_equal(global_nodes, c.global_nodes)
        break


def test_block_connectivity_local_all(edb):
    e = edb
    util.print_function_starting()
    items = e.element_blocks.connectivity_local_all()
    print('   iterator is type = {}'.format(type(items)), flush=True)
    for key, block, local_connectivity in items:
        util.print_blue('    block {} unique {}'.format(block, local_connectivity.num_unique))
        # check reconstructed global connectivity is identical to that read from the database
        global_connectivity = get_global_connectivity_from_local(local_connectivity)
        connectivity = block.connectivity()
        np.testing.assert_equal(global_connectivity, connectivity)


def test_side_sets(edb):
    util.print_function_starting()
    for key, side_set in edb.side_sets.items():
        elements, sides = side_set.entries()
        assert isinstance(elements, np.ndarray)
        assert isinstance(sides, np.ndarray)
        num_entries = side_set.num_entries()
        assert elements.shape == (num_entries,)  # a side is an element and a local side index
        assert sides.shape == (num_entries,)
        print('    side set {:3d} with {:5d} sides'.format(key, num_entries))


def test_field_array():
    util.print_function_starting()
    pos_info = ex.Field('position', ('x', 'y', 'z'))
    print('    pos_info:', pos_info)
    arr = np.zeros(3, dtype=np.double)
    print('    arr:', arr)
    obj = ex.FieldArray(arr, pos_info)
    assert str(type(obj)) == '<class \'affect.exodus.FieldArray\'>'
    assert obj.info.name == 'position'
    assert len(obj.info.components) == 3
    assert obj.shape == (3,)
    assert obj.dtype == np.double
    v = obj[1:]
    assert str(type(v)) == '<class \'affect.exodus.FieldArray\'>'
    assert len(v.info.components) == 3


def assert_field_variable_length(entity_dict):
    util.print_function_starting()
    entity_fields = entity_dict.fields
    assert entity_dict.num_variables() == sum(len(f.components) for f in entity_fields.values())


def print_entity_fields_info(entity_dict):
    util.print_function_starting()
    util.print_blue('    {} fields:'.format(entity_dict.__class__.__name__))
    assert_field_variable_length(entity_dict)
    for v in entity_dict.fields.values():
        print('    ' + str(v))


def test_fields(edb):
    util.print_function_starting()
    print_entity_fields_info(edb.globals)
    print_entity_fields_info(edb.nodal)
    print_entity_fields_info(edb.element_blocks)


def test_node_variable_at_all_times(edb):
    util.print_function_starting()
    num_variables = edb.nodal.num_variables()
    names = edb.nodal.variable_names()
    num_time_steps = edb.globals.num_times()
    num_nodes = edb.globals.num_nodes()
    print('    num_times = {}'.format(num_time_steps))
    for i in range(num_variables):
        node = num_nodes - 1  # the last node
        values = edb.nodal.variable_at_times(i, node, 0, -1)
        assert isinstance(names, list)
        assert isinstance(values, np.ndarray)
        print('    variable {} at last node ({}) is {} of {} with shape {}:'.format(
            names[i], node, type(values), values.dtype, values.shape))
        # print(values)


def test_node_field_at_all_times(edb):
    util.print_function_starting()
    fields = edb.nodal.fields
    num_nodes = edb.globals.num_nodes()
    node = num_nodes - 1
    for f in fields.values():
        util.print_blue('field {}:'.format(f))
        node_field_at_times = edb.nodal.field_at_times(f, node, 0, -1)
        print('    for node_{} is with shape {}:'.format(node, node_field_at_times.shape))
        print(node_field_at_times)


def test_times(edb):
    util.print_function_starting()
    num_times = edb.globals.num_times()
    print('    num_times = {}'.format(num_times))
    times = edb.globals.times()
    print('{}'.format(times))


def test_info_records(edb):
    util.print_function_starting()
    records = edb.info_records
    if len(records):
        print('   length: {}'.format(len(records)))
        #for record in records:
        #    print(f'    {record}')


def test_qa_records(edb):
    util.print_function_starting()
    records = edb.qa_records
    if len(records):
        for record in records:
            print('  {}'.format(record))


def test_summary(edb):
    util.print_function_starting()
    print(edb.summary())


def test_global(edb):
    util.print_function_starting()
    for k, v in edb.globals.__dict__.items():
        if k.startswith('__'):
            pass
        else:
            print('    {}: {}'.format(k, v))


def test_node_map_key_error(edb):
    util.print_function_starting()
    bad_id = -9999999
    node_maps = edb.node_maps
    with pytest.raises(ex.EntityKeyError):
        node_map = node_maps[bad_id]
        print(node_map)
    print('    no map with id = {}'.format(bad_id))


def test_face_set_key_error(edb):
    util.print_function_starting()
    bad_id = -9999999
    face_sets = edb.face_sets
    with pytest.raises(ex.EntityKeyError):
        face_set = face_sets[bad_id]
        print(face_set)
    print('    no set with id = {}'.format(bad_id))


def test_coordinates_local(edb):
    e = edb
    util.print_function_starting()

    block_key = 1
    mean_a = 0.0

    # read range of local coordinates for only one block
    block = e.element_blocks[block_key]
    local_connectivity = block.connectivity_local()  # the LocalConnectivity object
    util.print_blue('    block {} {}'.format(block_key, block.topology_name))
    local_coordinates = e.nodal.coordinates_local(local_connectivity)
    mean_b = np.mean(local_coordinates, axis=0)
    print('    mean local coordinates', mean_b)

    # read all global coordinates once and use numpy.take in order to index
    global_coordinates = e.nodal.coordinates()
    for key, block in e.element_blocks.items():
        if key != block_key:
            continue
        local_connectivity = block.connectivity_local()  # the LocalConnectivity object
        util.print_blue('    block {} {}'.format(key, block.topology_name))
        local_coordinates = global_coordinates.take(local_connectivity.global_nodes, axis=0)
        mean_a = np.mean(local_coordinates, axis=0)
        print('    mean local coordinates', mean_a)
    np.testing.assert_array_equal(mean_a, mean_b)
