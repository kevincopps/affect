import sys
import os
import numpy as np
import textwrap
import pytest
from .. import exodus as ex
from .. import connect


class TestDisplay():


    e = None  # database for most of the test functions

    def setUp(self):

        exodus.debug_messages(exodus.VERBOSE|exodus.DEBUG)

        # read a mesh we will use for many of the tests
        #base = "../../meshes/"
        #file = "cube_1M_elem.e"
        #file = "contact_puzzle.e"
        #file = "thermal/francis-w76-ISLloc1ht.e"
        #file = "thermal/thermal_final.e"
        #file = "impact_stage/impact-stage-history.e"
        #file = "lapjoint_hex/lapjoint_hex.e"
        #file = "large25m/b6112_unstr_out.e"

        base = './'
        file = 'one_hex.e'

        path = base + file
        self.e = exodus.Database(path)
        print(str(type(self.e)))

    def tearDown(self):
        pass

    def test_vertex_to_element(self):
        blocks = self.e.element_blocks
        print '{} element blocks'.format(len(blocks))
        print 'IDs = {}'.format(blocks.keys())
        num_nodes = self.e.globals.num_nodes()
        print 'num_nodes = {}'.format(num_nodes)
        for key, block in blocks.iteritems():
            print str(block)
            print str(block._topology_name)
            connectivity = block.connectivity()
            print 'block {} connectivity = '.format(key)
            print type(connectivity)
            print 'connectivity.shape = {}'.format(connectivity.shape)
            print connectivity
            print
            topology = connect.get_topology(block._topology_name)
            num_elements = block._num_entries
            vertex_begin = connect.vertex_to_element(num_nodes, topology, num_elements, connectivity)
            print type(vertex_begin)
            print vertex_begin

