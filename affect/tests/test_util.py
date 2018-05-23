from .. import util
import pickle
import numpy as np


def test_byte_align():
    a = np.empty((2, 3, 2), dtype=int)
    b = util.byte_align(a)
    np.testing.assert_equal(a, b)
    assert util.is_byte_aligned(b)


def test_empty_aligned():
    shape = (2, 3, 2)
    d = util.empty_aligned(shape, dtype=np.double)
    assert d.size == 12
    assert d.shape == shape
    assert d.dtype == np.double
    assert util.is_byte_aligned(d)
    e = util.empty_aligned(shape, dtype='int32', n=128)
    assert e.size == 12
    assert e.shape == shape
    assert e.dtype == np.int32
    assert util.is_byte_aligned(e, 128)


def test_zeros_aligned():
    shape = (2, 3)
    x = util.zeros_aligned(shape, dtype=np.uint32, n=32)
    assert x.size == 6
    assert x.shape == shape
    assert x.dtype == np.uint32
    assert util.is_byte_aligned(x, 32)
    assert not np.any(x)


def test_compress():
    x = np.array([[1, 2.0], [3.0, 4]])
    b = util.compress(x)
    y = util.decompress(b, x.shape, x.dtype)
    assert util.is_byte_aligned(y)
    np.testing.assert_array_equal(x, y)


def test_compressed_array_from_array():
    x = np.array([[1, 2.0], [3.0, 4]])
    b = util.CompressedArray(x)
    pickle_string = pickle.dumps(b)  # also test pickling/unpickling
    c = pickle.loads(pickle_string)
    y = c.unpack()
    np.testing.assert_array_equal(x, y)


def test_compressed_array_from_pointer():
    x = np.array([[1, 2.0], [3.0, 4]])
    # normally the following constructor is only called from within Cython
    b = util.CompressedArray(shape=(2, 2), dtype=x.dtype, pointer=x.__array_interface__['data'][0])
    pickle_string = pickle.dumps(b)
    c = pickle.loads(pickle_string)
    y = c.unpack()
    np.testing.assert_array_equal(x, y)


def test_compressed_array_shape():
    x = np.array([[1, 2.0], [3.0, 4], [5, 6.0]])
    b = util.CompressedArray(x)
    assert x.shape == b.shape


def test_compressed_array_dtype():
    x = np.array([[1, 2.0], [3.0, 4], [5, 6.0]])
    b = util.CompressedArray(x)
    assert x.dtype == b.dtype


def test_compressed_array_size():
    x = np.array([[1, 2.0], [3.0, 4], [5, 6.0]])
    b = util.CompressedArray(x)
    assert x.size == b.size


def test_compressed_array_itemsize():
    x = np.array([[1, 2.0], [3.0, 4], [5, 6.0]])
    b = util.CompressedArray(x)
    assert x.itemsize == b.itemsize


def test_compressed_array_nbytes():
    x = np.array([[1, 2.0], [3.0, 4], [5, 6.0]])
    b = util.CompressedArray(x)
    assert x.nbytes == b.nbytes


def test_compressed_array_len():
    x = np.random.random((128, 128))
    b = util.CompressedArray(x)
    print('\n    16384 (float64) compression ratio = {:6.3f}'.format(len(b) / x.nbytes))
    assert len(b) < x.nbytes  # rough guess of compression ratio of 0.25
    y = np.random.randint(16384, size=(128, 128))
    c = util.CompressedArray(y)
    print('    16384 (int) compression ratio = {:6.3f}'.format(len(c) / y.nbytes))
    assert len(c) < y.nbytes  # rough guess of compression ratio of 0.25


def test_take_indices():
    x = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], [10.0, 11.0, 12.0]])
    util.print_array_info('x', x)
    i = np.array([1, 2])
    util.print_array_info('i', i)
    y = util.take(x, i)
    util.print_array_info('y', y)
    z = x[i, :]  # fancy indexing
    util.print_array_info('z', z)
    np.testing.assert_array_equal(y, z)


