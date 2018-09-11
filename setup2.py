from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("chim_search2.pyx")
)
