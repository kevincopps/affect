#!/usr/bin/env python

from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext
import pkg_resources
import distutils.ccompiler
import distutils.sysconfig
import glob
import os
import sys
import multiprocessing.pool

long_description = """
Affect is a library for processing computer simulation data on unstructured grids.
"""

class BuildExtensions(build_ext):
    """
    Subclass setuptools build_ext command. Does the following
    1) it makes sure numpy is available
    2) it injects numpy's core/include directory in the include_dirs parameter of all extensions
    3) it runs the original build_ext command
    """

    def run(self):
        # According to
        # https://pip.pypa.io/en/stable/reference/pip_install.html#installation-order
        # at this point we can be sure pip has already installed numpy
        numpy_includes = pkg_resources.resource_filename('numpy', 'core/include')
        for ext in self.extensions:
            if hasattr(ext, 'include_dirs'):
                if numpy_includes not in ext.include_dirs:
                    ext.include_dirs.append(numpy_includes)
        build_ext.run(self)

def get_netcdf_include():
    # setup and distutils compiler does not find the base include path in anaconda installs
    include = 'netcdf.h'
    start = os.path.abspath(os.path.join(distutils.sysconfig.get_python_inc(), '..'))
    for root, dirs, files in os.walk(start):
        for file in files:
            if file == include:
                return root
    print('Could not find {}'.format(include))
    sys.exit(1)


NUMBER_PARALLEL_COMPILES = 8

# platform specific header and library directories
other_library = ''
other_link_args = ''
other_compile_args = ''
exodus_compile_args = ''

connect_source_files = ['affect/connect.pyx']
connect_source_files += glob.glob('affect/src/connect/*.cpp')
connect_include = ['affect/src/connect']

exodus_source_files = ['affect/exodus.pyx']
exodus_source_files += glob.glob('affect/src/exodus/*.c')
exodus_include = ['affect/src/exodus']

python_base = sys.prefix
_platform = sys.platform

if _platform == 'linux' or _platform == 'linux2':
    os.environ["CC"] = 'gcc'
    os.environ["CXX"] = 'g++'
    other_library = []  #['/usr/local/opt/llvm/lib']  # location of libiomp5 (however, it may be in anaconda)
    other_link_args = []  #['-mmacosx-version-min=10.11']
    exodus_include.append(get_netcdf_include())
    other_compile_args = ['-fopenmp']
    exodus_compile_args = ['-Dexodus_EXPORTS', '-Wno-sign-compare', '-Wno-uninitialized']
    connect_compile_args = ['-std=c++0x', '-Wno-deprecated', '-Wno-unused-variable',
                            '-Wno-uninitialized']

elif _platform == 'darwin':
    # Prerequisites:
    #   brew install llvm (until Apple "clang -fopenmp" will be supported)
    #   export PATH="/usr/local/opt/llvm/bin:${PATH}"
    #   may also have to "export DYLD_LIBRARY_PATH=/usr/local/lib"?
    os.environ["CC"] = 'clang'
    os.environ["CXX"] = 'clang'
    other_library = ['/usr/local/opt/llvm/lib']  # location of libiomp5 (however, it may be in anaconda)
    other_link_args = ['-mmacosx-version-min=10.11']
    other_include = []
    other_compile_args = ['-stdlib=libc++', '-mmacosx-version-min=10.11', '-fopenmp']
    exodus_compile_args = ['-Dexodus_EXPORTS', '-Wno-unused-function', '-Wno-sometimes-uninitialized',
                           '-Wno-unreachable-code', '-Wno-sign-compare']
    connect_compile_args = ['-Wno-unused-function', '-Wno-unneeded-internal-declaration', '-Wno-unused-variable']
elif _platform == 'win32':
    os.environ["CC"] = 'gcc'
    os.environ["CXX"] = 'gcc'


extensions = [
    Extension('affect.exodus',
              sources=exodus_source_files,
              include_dirs=exodus_include,
              libraries=['iomp5', 'netcdf'],
              library_dirs=other_library,
              extra_compile_args=other_compile_args+exodus_compile_args,
              extra_link_args=other_link_args,
              language='c',
              ),
    Extension('affect.connect',
              include_dirs=connect_include,
              sources=connect_source_files,
              libraries=['iomp5'],
              library_dirs=other_library,
              extra_compile_args=other_compile_args+connect_compile_args,
              extra_link_args=other_link_args,
              language='c++',
              ),
]


def parallel_c_compile(self, sources, output_dir=None, macros=None, include_dirs=None, debug=0, extra_preargs=None,
                       extra_postargs=None, depends=None):
    """
    Enable parallel compiles of C/C++ during development. A monkey-patch of the distutils.ccompiler.
    """
    # those lines are copied from distutils.ccompiler.CCompiler directly
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(output_dir, macros, include_dirs, sources,
                                                                          depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    def _single_compile(obj):
        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(NUMBER_PARALLEL_COMPILES).imap(_single_compile, objects))
    return objects

# Monkey patch to allow parallel compile
distutils.ccompiler.CCompiler.compile = parallel_c_compile

install_requires = ['cython >= 0.24.1', 'hdf4 >= 4.2.12', 'hdf5 >= 1.8.17', 'netcdf4 >= 1.2.4', 'numexpr >= 2.6.1',
                    'numpy >= 1.11.1']

setup(
    name='affect',
    description='Affect - Processing Computational Simulations',
    long_description=long_description,
    license='MIT',
    version='0.1',
    author='Kevin Copps',
    author_email='kdcopps@sandia.gov',
    maintainer='Kevin Copps',
    maintainer_email='kdcopps@sandia.gov',
    url='https://github.com/kdcopps/affect',
    packages=['affect'],
    classifiers=['Programming Language :: Python :: 3', ],
    setup_requires=['setuptools>=18.0', 'numpy', 'cython', 'pytest-runner'],
    install_requires=install_requires,
    tests_require=['pytest'],
    zip_safe=False,
    cmdclass={'build_ext': BuildExtensions},
    ext_modules=extensions
)
