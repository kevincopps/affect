#!/usr/bin/env python

from setuptools import setup
from Cython.Distutils import Extension
from Cython.Distutils.build_ext import build_ext
# from setuptools.extension import Extension
# from setuptools.command.build_ext import build_ext
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


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace('-Wstrict-prototypes', '')


class BuildExtensions(build_ext):
    """
    Subclass setuptools build_ext command. Does the following
    1) it makes sure numpy is available
    2) it injects numpy's core/include directory in the include_dirs parameter of all extensions
    3) it runs the original build_ext command
    """

    def run(self):
        # According to https://pip.pypa.io/en/stable/reference/pip_install.html#installation-order
        # at this point we can be sure pip has already installed numpy
        from numpy import get_include
        numpy_includes = get_include()
        # Alternatively
        # numpy_includes = pkg_resources.resource_filename('numpy', 'core/include')
        for ext in self.extensions:
            if hasattr(ext, 'include_dirs'):
                if numpy_includes not in ext.include_dirs:
                    ext.include_dirs.append(numpy_includes)
        build_ext.run(self)


def get_netcdf_include():
    # setup and distutils compiler does not find the base include path in anaconda installs
    include = 'netcdf.h'
    python_include = distutils.sysconfig.get_python_inc()
    start = os.path.abspath(os.path.join(python_include, '..'))
    for root, dirs, files in os.walk(start):
        for file in files:
            if file == include:
                return root
    raise RuntimeError('Error: setup.py could not find {}'.format(include))

def get_cpu_options(platform):
    # return a clang/gcc compile flag with the highest level of sse, avx instructions supported on the compile CPU
    flag = ['']
    if platform == 'darwin':
        return ['-march=native']
    elif platform.startswith('linux'):
        return ['-march=native']
    else:
        raise Exception(f'Unknown platform ({platform}')
    return flag


# platform specific header and library directories
other_include = ''
other_library = ''
other_link_args = ''
other_compile_args = ''

util_compile_args = ''
exodus_compile_args = ''
connect_compile_args = ''
arithmetic_compile_args = ''

util_source_files = ['affect/util.pyx']
util_source_files += glob.glob('affect/src/util/*.cpp')
util_include = ['affect/src/util', 'affect/src']

arithmetic_source_files = ['affect/arithmetic.pyx']
arithmetic_source_files += glob.glob('affect/src/arithmetic/*.cpp')
arithmetic_include = ['affect/src/arithmetic', 'affect/src']

connect_source_files = ['affect/connect.pyx']
connect_source_files += glob.glob('affect/src/connect/*.cpp')
connect_include = ['affect/src/connect', 'affect/src']

dynamics_source_files = ['affect/dynamics.pyx']

exodus_source_files = ['affect/exodus.pyx']
exodus_source_files += glob.glob('affect/src/exodus/*.c')
exodus_include = ['affect/src/exodus']

python_base = sys.prefix
_platform = sys.platform

if _platform == 'linux' or _platform == 'linux2':
    os.environ['CC'] = 'gcc'
    os.environ['CXX'] = 'g++'
    other_include = []
    other_library = []  # ['/usr/local/opt/llvm/lib']  # location of libiomp5 (however, it may be in anaconda)
    other_link_args = ['-fopenmp']  # ['-mmacosx-version-min=10.11']
    exodus_include.append(get_netcdf_include())

    other_compile_args = ['-fopenmp']
    other_compile_args += get_cpu_options(_platform)
    exodus_compile_args = ['-Dexodus_EXPORTS', '-std=c11', '-Wno-sign-compare', '-Wno-uninitialized']
    connect_compile_args = ['-std=c++14', '-Wno-deprecated', '-Wno-unused-variable', '-Wno-uninitialized']
    arithmetic_compile_args = ['-std=c++14', '-Wno-deprecated', '-Wno-unused-variable', '-Wno-uninitialized']
    util_compile_args = ['-std=c++14', '-Wno-unused-function']

elif _platform == 'darwin':

    other_library = ['/usr/local/opt/llvm/lib']  # location of libiomp5 (however, it may be in anaconda)
    other_link_args = ['-mmacosx-version-min=10.11']
    other_include = []

    # -fslp-vectorize-aggressive
    # is a relatively minor performance improvement for long compile/link times
    # other_compile_args = ['-stdlib=libc++', '-mmacosx-version-min=10.11', '-fopenmp', '-march=corei7-avx', 
    #                       '-mavx2', '-fslp-vectorize-aggressive']
    other_compile_args = ['-mmacosx-version-min=10.11', '-fopenmp']
    other_compile_args += get_cpu_options(_platform)
    exodus_compile_args = ['-Dexodus_EXPORTS', '-Wno-unused-function', '-Wno-sometimes-uninitialized',
                           '-Wno-unreachable-code', '-Wno-sign-compare']
    connect_compile_args = ['-std=c++14', '-Wno-unused-function', '-Wno-unused-variable', '-Wno-deprecated']
    arithmetic_compile_args = ['-std=c++14', '-Wno-unused-function', '-Wno-unused-variable', '-Wno-deprecated']
    util_compile_args = ['-std=c++14', '-Wno-unused-function']

    os.environ['CC'] = 'gcc-7'
    os.environ['CXX'] = 'g++-7'

    # Clang Prerequisites:
    #   install the Xcode and command line tools:
    #       xcode-select --install
    #   brew install llvm (until Apple "clang -fopenmp" is supported)
    #   export PATH="/usr/local/opt/llvm/bin:${PATH}"
    #   in the past we have had to also "export DYLD_LIBRARY_PATH=/usr/local/lib"?
    # To use the bundled libc++ please add the following LDFLAGS:
    #   LDFLAGS = "-L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib"
    # For compilers to find this software you may need to set:
    #   LDFLAGS: -L/usr/local/opt/llvm/lib
    #   CPPFLAGS: -I/usr/local/opt/llvm/include
    #
    # What target flags does a machine support? Try the following:
    #    'echo | clang -E - -march=native -###'
    #
    # os.environ["CC"] = 'clang'
    # os.environ["CXX"] = 'clang'
    # # OpenMP with offloading to Nvidia GPU
    # # -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda --cuda-path=$CUDA_HOME
    # other_compile_args += '-stdlib=libc++'
    # connect_compile_args += ['-Wno-unneeded-internal-declaration', '-Wshorten-64-to-32',
    #                          '-Rpass-analysis=loop-vectorize']
    # util_compile_args += ['-Wshorten-64-to-32']


elif _platform == 'win32':
    os.environ['CC'] = 'gcc'
    os.environ['CXX'] = 'gcc'


extensions = [
    Extension('affect.exodus',
              sources=exodus_source_files,
              include_dirs=other_include + exodus_include,
              libraries=['iomp5', 'netcdf'],
              library_dirs=other_library,
              extra_compile_args=other_compile_args + exodus_compile_args,
              extra_link_args=other_link_args,
              language='c',
              ),
    Extension('affect.connect',
              include_dirs=other_include + connect_include,
              sources=connect_source_files,
              libraries=['iomp5'],
              library_dirs=other_library,
              extra_compile_args=other_compile_args + connect_compile_args,
              extra_link_args=other_link_args,
              language='c++',
              # undef_macros=['NDEBUG']  # enable assert(), etc.
              ),
    Extension('affect.arithmetic',
              include_dirs=other_include + arithmetic_include,
              sources=arithmetic_source_files,
              libraries=['iomp5'],  # 'tbb'],
              library_dirs=other_library,
              extra_compile_args=other_compile_args + arithmetic_compile_args,
              extra_link_args=other_link_args,
              language='c++',
              # undef_macros=['NDEBUG']  # enable assert(), etc.
              ),
    Extension('affect.dynamics',
              sources=dynamics_source_files,
              library_dirs=other_library,
              extra_compile_args=other_compile_args,
              extra_link_args=other_link_args,
              language='c',
              ),
    Extension('affect.util',
              include_dirs=other_include + util_include,
              sources=util_source_files,
              library_dirs=other_library,
              extra_compile_args=other_compile_args+util_compile_args,
              extra_link_args=other_link_args,
              language='c++',
              # undef_macros=['NDEBUG']  # enable assert(), etc.
              ),
]

#
# global cython directives
#
for e in extensions:
    e.cython_directives = {'embedsignature': True}

    # build with support coverage or profiling using
    # "python setup.py --linetrace build_ext -i"
    if '--linetrace' in sys.argv:
        e.define_macros = [('CYTHON_TRACE_NOGIL', '1')]
        e.cython_directives = {'linetrace': True, 'binding': True, 'profile': True, 'embedsignature': True}

if '--linetrace' in sys.argv:
    sys.argv.remove('--linetrace')

NUMBER_PARALLEL_COMPILES = 16

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


# monkey patch compile to enable parallel compile
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
    setup_requires=['setuptools>=18.0', 'numpy', 'cython', 'netcdf4', 'pytest-runner'],
    install_requires=install_requires,
    tests_require=['pytest'],
    zip_safe=False,
    cmdclass={'build_ext': BuildExtensions},  # use our own subclass
    ext_modules=extensions
)

# alternatively
# cmdclass={'build_ext': Cython.Build.build_ext},
