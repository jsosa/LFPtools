from setuptools import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

binaries = ['bin/lfp-fixelevs',
            'bin/lfp-getbankelevs',
            'bin/lfp-getbedelevs',
            'bin/lfp-getdepths',
            'bin/lfp-getslopes',
            'bin/lfp-getwidths',
            'bin/lfp-prepdata',
            'bin/lfp-rasterresample',
            'bin/lfp-split',
            'bin/lfp-getinflows',
            ]

ext_modules = [
    Extension("lfptools.prepdata_utils", ["lfptools/prepdata_utils.pyx"],
              include_dirs=[numpy.get_include()],
              )
]

setup(
    name='lfptools',
    version='0.1',
    description='A software suite to easily & quickly prepare large scale LISFLOOD-FP models using freely available data',
    url='http://github.com/jsosa/LFPtools',
    author='Jeison Sosa',
    author_email='sosa.jeison@gmail.com',
    license='MIT',
    packages=['lfptools'],
    zip_safe=False,
    scripts=binaries,
    ext_modules=cythonize(ext_modules),
)
