import numpy as np
import timeit
from .. import exodus as ex


BASE_PATH = '../../../meshes/'
# read a mesh we will use for many of the tests
# FILE = "cube_1M_elem.e"
# FILE = "contact_puzzle.e"
FILE = 'thermal/francis-w76-ISLloc1ht.e'
# FILE = "thermal/thermal_final.e"
# FILE = "impact_stage/impact-stage-history.e"
# FILE = "lapjoint_hex/lapjoint_hex.e"
# FILE = "large25m/b6112_unstr_out.e"
FILE_PATH = BASE_PATH + FILE


def side_set_entries():
    e = ex.Database(FILE_PATH)
    side_sets = e.side_sets
    # keys = side_sets.keys()
    side_set = side_sets[4]
    entries = side_set.entries()


def block_connectivity():
    e = ex.Database(FILE_PATH)
    block_id = 17
    block = e.element_blocks[block_id]
    # num_elem = block.num_entries()
    # print 'block_{} with {}'.format(block_id, num_elem)
    connectivity = block.connectivity
    

def bounding_box():
    e = ex.Database(FILE_PATH)
    ndim = e.globals.dimension()
    coordinates = e.nodal.coordinates()

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

    if ndim == 2:
        min_x, min_y = np.min(coordinates, axis=0)
        max_x, max_y = np.max(coordinates, axis=0)
        bbox = np.array([(min_x, min_y), (max_x, max_y)])
    elif ndim == 3:
        min_x, min_y, min_z = np.min(coordinates, axis=0)
        max_x, max_y, max_z = np.max(coordinates, axis=0)
        bbox = np.array([(min_x, min_y, min_z), (max_x, max_y, max_z)])


if __name__ == '__main__':

    t = timeit.timeit("block_connectivity()", setup="from __main__ import block_connectivity",number=1000)
    print('block_connectivity() time = {}'.format(t))