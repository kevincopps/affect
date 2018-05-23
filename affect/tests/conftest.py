# configuration and fixtures for pytest tests

from .. import exodus as ex
import os
import pytest
import numpy

CURRENT_PATH = os.path.dirname(__file__)
MESH_PATH = os.path.abspath(os.path.join(CURRENT_PATH, '../../../meshes/'))  # large mesh test files
DATA_PATH = os.path.abspath(os.path.join(CURRENT_PATH, './data/'))  # large mesh test files


def open_database(file_name: str, directory=MESH_PATH) -> ex.Database:
    """
    Open database file that exists in the MESH_PATH and print message to stdout.

    Turn on verbose debug messages.

    Returns:
        exodus - the database instance.
    """
    path = os.path.join(directory, file_name)
    ex.debug_messages(ex.Messages.VERBOSE | ex.Messages.DEBUG)  # turn on capturing debug messages
    e = ex.Database(path)
    print(f'\nOpening {file_name} {str(type(e))}\n')
    return e


def release_database(file_name: str, exodus: ex.Database):
    """
    Release the database file resource and set debug messages to previous default.

    Args:
        file_name: name of the file opened previously
        exodus: database instance
    """
    print(f'\nClosing {file_name}')
    exodus.close()
    ex.debug_messages(ex.Messages.DEFAULT)  # reset debug messages off


@pytest.fixture(scope="session")
def edb():
    """
    Return a single open Exodus database from a collection of various Exodus database files on disk.
    Returns: exodus
    """

    # read a mesh we will use for many of the tests
    # file_name = 'cube_1M_elem.e'
    # file_name = 'contact_puzzle.e'
    file_name = 'thermal/francis-w76-ISLloc1ht.e'
    # file_name = 'thermal/thermal_final.e'
    # file_name = 'impact_stage/impact-stage-history.e'
    # file_name = 'lapjoint_hex/lapjoint_hex.e'
    # file_name = 'large25m/b6112_unstr_out.e'
    # file_name = 'my_mock_aff_1_scale_13.g'

    e = open_database(file_name)
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)


@pytest.fixture(scope="session")
def edb_display():
    """
    Return a single open Exodus database from a collection of various Exodus database files on disk.

    Returns: exodus
    """

    # read a mesh we will use for many of the tests
    # file_name = 'cube_1M_elem.e'
    file_name = 'contact_puzzle.e'
    # file_name = 'thermal/francis-w76-ISLloc1ht.e'
    # file_name = 'thermal/thermal_final.e'
    # file_name = 'impact_stage/impact-stage-history.e'
    # file_name = 'lapjoint_hex/lapjoint_hex.e'
    # file_name = 'large25m/b6112_unstr_out.e'
    e = open_database(file_name)
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)


@pytest.fixture(scope="session")
def edb_frf_database():
    """
    Return a single open Exodus database from a collection of various Exodus database files on disk.

    Returns:
        exodus
    """
    file_name = 'p1f.exo'
    subdirectory = 'SRS-FRF-example/model'
    base = os.path.abspath(os.path.join(MESH_PATH, subdirectory))
    e = open_database(file_name, base, )
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)

@pytest.fixture(scope="session")
def frf_data_dict():
    """

    Returns:
        loaded - dict of ndarray of same length dtype=numpy.float64_t, with keys times, force_z, and acceleration_z
    """
    path = os.path.abspath(os.path.join(DATA_PATH, 'frf_test_data.npz'))
    loaded = numpy.load(path)  # load compressed arrays from .npz file
    yield loaded

@pytest.fixture(scope="session")
def edb_small():
    """
    Return a single open Exodus database from a collection of various Exodus database files on disk.
    
    Returns:
        exodus
    """
    # file_name = 'one_hex.exo'
    # file_name = 'test.exo'
    file_name = 'tet_hex.exo'
    # file_name = 'tet_pyramid.exo'
    # file_name = 'tet_wedge.exo'

    e = open_database(file_name, DATA_PATH, )
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)


@pytest.fixture(scope="session")
def edb_medium():
    # file_name = 'my_mock_aff_1_scale_13.g'
    # file_name = 'cylinder.g'
    # file_name = 'lapjoint_hex/lapjoint_hex.e'
    file_name = 'sphere_void_layers4.g'
    e = open_database(file_name)
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)


@pytest.fixture(scope="session")
def edb_large():
    """
    Return a single open Exodus database with 1M elements or more from a collection of Exodus database files on disk.

    Returns:
        exodus
    """

    file_name = 'cube_1M_elem.e'

    e = open_database(file_name)
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)


@pytest.fixture(scope="session")
def edb_huge():
    """
    Return a single open Exodus database with 25M elements or more from a collection of Exodus database files on disk.

    Returns:
        exodus
    """
    file_name = 'large25m/b6112_unstr_out.e'  # has pyramids
    e = open_database(file_name)
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)


@pytest.fixture(scope="session")
def edb_huge_one_block():
    """
    Return a single open Exodus database with 25M elements or more from a collection of Exodus database files on disk.

    Returns:
        exodus
    """
    # file_name = 'large25m/b6112_unstr_out.e' # has pyramids
    # file_name = 'LNG_mossT23.g' # seems corrupt in name dimension
    # file_name = 'dowding-fc-model.e' # 21913 nodes lots of field data
    file_name = '100cm_S.g'

    e = open_database(file_name)
    yield e  # everything after yield serves as tear down and reset
    release_database(file_name, e)
