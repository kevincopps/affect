#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import getpass
import glob
import os
import numpy as np

# platform specific header and library directories
from sys import platform as _platform
if _platform == "linux" or _platform == "linux2":
    # linux
    pass
elif _platform == "darwin":
    os.environ["CC"] = "clang"
    os.environ["CXX"] = "clang"
    user = getpass.getuser()
    other_base = "/Users/" + user + "/Library/Enthought/Canopy_64bit/User"
    other_include = other_base + "/include"
    other_library = other_base + "/lib"
elif _platform == "win32":
    # Windows...
    pass

source_files = ["exodus.pyx"]
source_files += glob.glob('../c-mesh/*.cpp')
mesh_include = '../c-mesh'

extensions = [
    Extension("exodus",
              include_dirs = [mesh_include,other_include,np.get_include()],
              sources=source_files,
              libraries=["exoIIv2c"],
              library_dirs = [other_library],
              extra_compile_args=["-I../../exodus/exodus/cbind/include/","-stdlib=libc++","-mmacosx-version-min=10.7"],
              extra_link_args=["-L../../exodus/exodus/cbind/"],
              language="c++",
             ),
]
setup(
    ext_modules = cythonize(extensions),
)

# setup(
#     cmdclass = {'build_ext': build_ext},
#     ext_modules = [
#         Extension("exodus", ["exodus.pyx"],
#                   include_dirs = [mesh_include,other_include,np.get_include()],
#                   libraries=["exoIIv2c"],
#                   library_dirs = [other_library],
#                   extra_compile_args=["-I../../exodus/exodus/cbind/include/"],
#                   extra_link_args=["-L../../exodus/exodus/cbind/"]
#                   #define_macros=[("NPY_NO_DEPRECATED_API",None)]
#                   #define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_8_API_VERSION')]
#                   )
#     ]
# )