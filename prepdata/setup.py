from distutils.core import setup
from Cython.Build import cythonize
import numpy

# python setup.py build_ext --inplace

setup(
    ext_modules  = cythonize("prepdata_utils.pyx"),
    include_dirs = [numpy.get_include()]
)
