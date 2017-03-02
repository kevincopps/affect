# conftest.py
# configuration and fixtures for pytest tests

import os
import textwrap
import pytest
from .. import exodus as ex


@pytest.fixture(scope="session")
def edb():
    """
    Return a single open Exodus database from a collection of various Exodus database files on disk.
    :return: exodus
    """

    dir_name = os.path.dirname(__file__)
    goal_dir = os.path.join(dir_name, '../../../meshes/')
    base_path = os.path.abspath(goal_dir)

    # read a mesh we will use for many of the tests
    # file_name = 'cube_1M_elem.e'
    # file_name = 'contact_puzzle.e'
    # file_name = 'thermal/francis-w76-ISLloc1ht.e'
    # file_name = 'thermal/thermal_final.e'
    file_name = 'impact_stage/impact-stage-history.e'
    # file_name = 'lapjoint_hex/lapjoint_hex.e'
    # file_name = 'large25m/b6112_unstr_out.e'

    path = os.path.join(base_path, file_name)
    ex.debug_messages(ex.VERBOSE | ex.DEBUG)  # turn on capturing debug messages
    e = ex.Database(path)
    print('\nOpening {} {}\n'.format(file_name, str(type(e))))
    # use yield instead of return
    yield e
    # everything after yield serves as tear down and reset
    print('\nClosing {}\n'.format(file_name))
    e.close()
    ex.debug_messages(ex.DEFAULT)  # reset debug messages off



