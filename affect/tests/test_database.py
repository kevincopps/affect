# test_database.py

from .. import exodus as ex
import numpy as np
import textwrap
import pytest


class ConsoleCode:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARK_CYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


def print_bold(*args):
    print(ConsoleCode.BOLD + ' '.join(map(str, args)) + ConsoleCode.END)
    

def print_blue(*args):
    print(ConsoleCode.BLUE + ' '.join(map(str, args)) + ConsoleCode.END)


def test_library_version():
    print_bold('\nlibrary_version:')
    print('    ExodusII library API version: {}'.format(ex.library_version()))


def test_non_existent_file():
    with pytest.raises(ex.FileNotFound):
        ex.Database("nonexistent.e", ex.Mode.READ_ONLY)


def test_str_repr(edb):
    print_bold('\nstr_repr:')
    print('    ' + str(edb))


def test_coordinates(edb):
    print_bold('\ncoordinates')
    x = edb.nodal.coordinates()
    print('    type {} with shape {}'.format(type(x), x.shape))
    

def test_bounding_box(edb):
    print_bold('\nbounding_box')
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
    print_bold('\nelement_blocks:')
    print('    {} element blocks'.format(len(edb.element_blocks)))
    print('    IDs:')
    for line in textwrap.wrap('{}'.format([k for k in edb.element_blocks.keys()])):
        print('        ' + line)
    print('    blocks:')
    for key, block in edb.element_blocks.items():
        print('    {}: {}'.format(key, block))


def test_attributes(edb):
    print_bold('\nblock attributes:')
    for key, block in edb.element_blocks.items():
        num_attributes = block.num_attributes()
        print('    {}: has {} attributes {}'.format(key, num_attributes, block.attribute_names()))
        if num_attributes:
            print('        {}'.format([block.attribute(k) for k in range(num_attributes)]))


def test_element_block_key_error(edb):
    print_bold('\nelement_block_key_error:')
    bad_id = -9999999
    element_blocks = edb.element_blocks
    with pytest.raises(ex.EntityKeyError):
        block = element_blocks[bad_id]
        print(block)
    print('    no block with id = {}'.format(bad_id))


def test_block_connectivity(edb):
    print_bold('\nblock connectivity:')
    min_entry = np.iinfo(np.int32).max
    max_entry = -1
    for key, block in edb.element_blocks.items():
        connectivity = block.connectivity()
        # print_blue('    block {} connectivity = '.format(key))
        # print(type(connectivity))
        # print('connectivity.shape = {}'.format(connectivity.shape))
        # print(connectivity)
        block_min_node = np.amin(connectivity)
        block_max_node = np.amax(connectivity)
        print('    block {:5d} first node {:7d} last node {:7d}'.format(key, block_min_node, block_max_node))
        min_entry = min(min_entry, block_min_node)
        max_entry = max(max_entry, block_max_node)
        partial_connectivity = block.partial_connectivity(0, 1)  # all nodes of 0th element
        np.testing.assert_equal(partial_connectivity, connectivity[0:1, :])
        print('    partial_connectivity(0,1) == connectivity[0:1, :]')
    assert min_entry == 0
    num_nodes = edb.globals.num_nodes()
    assert max_entry == num_nodes - 1


def test_side_sets(edb):
    print_bold('\nside sets:')
    for key, side_set in edb.side_sets.items():
        sides = side_set.entries()
        assert isinstance(sides, np.ndarray)
        num_entries = side_set.num_entries()
        assert sides.shape == (num_entries, 2)  # a side is an element and a local side index
        print('    side set {:3d} with {:5d} sides'.format(key, num_entries))


def test_field_array():
    print_bold('\nfield_array:')
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
    entity_fields = entity_dict.fields
    assert entity_dict.num_variables() == sum(len(f.components) for f in entity_fields.values())


def print_entity_fields_info(entity_dict):
    print_blue('    {} fields:'.format(entity_dict.__class__.__name__))
    assert_field_variable_length(entity_dict)
    for v in entity_dict.fields.values():
        print('    ' + str(v))


def test_fields(edb):
    print_bold('\ntest_fields:')
    print_entity_fields_info(edb.globals)
    print_entity_fields_info(edb.nodal)
    print_entity_fields_info(edb.element_blocks)


def test_node_variable_at_all_times(edb):
    print_bold('\nnode_variable_at_all_times:')
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
    print_bold('\nnode_field_at_all_times:')
    fields = edb.nodal.fields
    num_nodes = edb.globals.num_nodes()
    node = num_nodes - 1
    for f in fields.values():
        print_blue('field {}:'.format(f))
        node_field_at_times = edb.nodal.field_at_times(f, node, 0, -1)
        print('    for node_{} is with shape {}:'.format(node, node_field_at_times.shape))
        print(node_field_at_times)


def test_times(edb):
    print_bold('\ntimes:')
    num_times = edb.globals.num_times()
    print('    num_times = {}'.format(num_times))
    times = edb.globals.times()
    print('{}'.format(times))


def test_info_records(edb):
    records = edb.info_records
    if len(records):
        print_bold('\ninfo_records:')
        print('   length: {}'.format(len(records)))


def test_qa_records(edb):
    records = edb.qa_records
    if len(records):
        print_bold('\nqa_records: ({})'.format(len(records)))
    for record in records:
        print('  {}'.format(record))


def test_summary(edb):
    print_bold('\nsummary:')
    print(edb.summary())


def test_global(edb):
    print_bold('\nglobal:')
    for k, v in edb.globals.__dict__.items():
        if k.startswith('__'):
            pass
        else:
            print('    {}: {}'.format(k, v))


def test_node_map_key_error(edb):
    print_bold('\nnode_map_key_error:')
    bad_id = -9999999
    node_maps = edb.node_maps
    with pytest.raises(ex.EntityKeyError):
        node_map = node_maps[bad_id]
        print(node_map)
    print('    no map with id = {}'.format(bad_id))


def test_face_set_key_error(edb):
    print_bold('\nface_set_key_error:')
    bad_id = -9999999
    face_sets = edb.face_sets
    with pytest.raises(ex.EntityKeyError):
        face_set = face_sets[bad_id]
        print(face_set)
    print('    no set with id = {}'.format(bad_id))
