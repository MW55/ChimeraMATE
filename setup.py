from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize('chimeramate_cy.pyx'))
