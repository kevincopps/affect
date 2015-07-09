import sys
sys.path.append("../affect")
import exodus

def side_set_entries():
    e = exodus.Database("francis-w76-ISLloc1ht.e")
    side_sets = e.side_sets
    #keys = side_sets.keys()
    side_set = side_sets[4]
    entries = side_set.entries()

def block_connectivity():
    e = exodus.Database("francis-w76-ISLloc1ht.e")
    block_id = 17
    block = e.element_blocks[block_id]
    #num_elem = block.num_entries()
    #print 'block_{} with {}'.format(block_id, num_elem)
    connectivity = block.connectivity

if __name__ == '__main__':
    import timeit
    t = timeit.timeit("block_connectivity()", setup="from __main__ import block_connectivity",number=1000)
    print 'block_connectivity() time = {}'.format(t)