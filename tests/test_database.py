import sys
import numpy as np
sys.path.append("../affect")
import exodus
from nose.tools import raises, assert_equal

class TestDatabase():


    e = None  # database for most of the test functions

    def setUp(self):

        exodus.debug_messages(exodus.VERBOSE|exodus.DEBUG)

        # read a mesh we will use for many of the tests
        base = "../../meshes/"
        #file = "cube_1M_elem.e"
        #file = "contact_puzzle.e"
        #file = "thermal/francis-w76-ISLloc1ht.e"
        #file = "thermal/thermal_final.e"
        #file = "impact_stage/impact-stage-history.e"
        #file = "lapjoint_hex/lapjoint_hex.e"
        file = "large25m/b6112_unstr_out.e"

        path = base + file
        self.e = exodus.Database(path)
        print str(type(self.e))

    def tearDown(self):
        pass

    def test_library_version(self):
        print 'Exodus library API version: {}'.format(exodus.library_version())

    @raises(IOError)
    def test_non_existent(self):
        y = exodus.Database("nonexistent.e", exodus.OPEN_READ_ONLY)

    def test_str_repr(self):
        print str(self.e)

    def test_get_coord(self):
        coords = self.e.nodal.coordinates()
        print 'nodal.coordinates is type {} with shape {}'.format(type(coords), coords.shape)
        #print coords

    def test_element_blocks(self):
        print '{} element blocks'.format(len(self.e.element_blocks))
        print 'IDs = {}'.format(self.e.element_blocks.keys())
        for key, block in self.e.element_blocks.iteritems():
            print str(block)

    def test_connectivity(self):
        min_entry = np.iinfo(np.int32).max
        max_entry = -1
        for key, block in self.e.element_blocks.iteritems():
            connectivity = block.connectivity()
            #print 'block {} connectivity = '.format(key)
            #print type(connectivity)
            #print 'connectivity.shape = {}'.format(connectivity.shape)
            #print connectivity
            min_entry = min(min_entry, np.amin(connectivity))
            max_entry = max(max_entry, np.amax(connectivity))
        assert_equal(min_entry, 0, 'minimum node entry != 0 in connectivity')
        num_nodes = self.e.globals.num_nodes()
        assert_equal(max_entry, num_nodes-1, 'maximum node entry != num_nodes-1 in connectivity')

    def test_get_side_sets(self):
        for key, side_set in self.e.side_sets.iteritems():
            sides = side_set.entries()
            print 'side_set {} with {} sides is type {} with shape {}'.format(
                key, side_set.num_entries(), type(sides), sides.shape)
            #print sides

    def test_field_array(self):
        pos_info = exodus.Field('position', ('x', 'y', 'z'))
        print 'pos_info:', pos_info
        arr = np.zeros(3, dtype=np.double)
        print 'arr:', arr
        obj = exodus.FieldArray(arr, pos_info)
        assert_equal(str(type(obj)), '<class \'exodus.FieldArray\'>')
        assert_equal(obj.info.name, 'position')
        assert_equal(len(obj.info.components), 3)
        assert_equal(obj.shape, (3,))
        assert_equal(obj.dtype, np.double)
        v = obj[1:]
        assert_equal(str(type(v)), '<class \'exodus.FieldArray\'>')
        assert_equal(len(v.info.components), 3)

    def assert_field_variable_length(self, entity_dict):
        entity_fields = entity_dict.fields
        assert_equal(entity_dict.num_variables(), sum(len(f.components) for f in entity_fields.itervalues()))

    def print_entity_fields_info(self, entity_dict):
        print '\n{} fields:'.format(entity_dict.__class__.__name__)
        self.assert_field_variable_length(entity_dict)
        for v in entity_dict.fields.itervalues():
            print v

    def test_fields(self):
        self.print_entity_fields_info(self.e.globals)
        self.print_entity_fields_info(self.e.nodal)
        self.print_entity_fields_info(self.e.element_blocks)

    def test_get_node_variable_at_all_times(self):
        num_variables = self.e.nodal.num_variables()
        names = self.e.nodal.variable_names()
        num_time_steps = self.e.globals.num_times()
        num_nodes = self.e.globals.num_nodes()
        print 'num_times = {}'.format(num_time_steps)
        for i in range(num_variables):
            node = num_nodes - 1 # the last node
            values = self.e.nodal.variable_at_times(i, node, 0, -1)
            print 'variable {} for node_{} is type {} with shape {}:'.format(names[i], node,
                                                                            type(values), values.shape)
            #print values

    def test_node_field_at_all_times(self):
        fields = self.e.nodal.fields
        num_nodes = self.e.globals.num_nodes()
        node = num_nodes - 1
        for f in fields.itervalues():
            print 'field {}:'.format(f)
            node_field_at_times = self.e.nodal.field_at_times(f, node, 0, -1)
            print '    for node_{} is with shape {}:'.format(node, node_field_at_times.shape)
            print node_field_at_times

    def test_num_times(self):
        num_times = self.e.globals.num_times()
        print "num_times = {}".format(num_times)

    def test_all_times(self):
        times = self.e.globals.times()
        print 'times:'
        print times

    def test_get_info_records(self):
        records = self.e.info_records
        if len(records):
            print 'len(info_records): {}'.format(len(records))
        #for record in records:
        #    print '  {}'.format(record)

    def test_get_qa_records(self):
        records = self.e.qa_records
        if len(records):
            print 'len(qa_records): {}'.format(len(records))
            print 'qa records:'
        for record in records:
            print '  {}'.format(record)

    def test_get_summary(self):
        print self.e.summary()