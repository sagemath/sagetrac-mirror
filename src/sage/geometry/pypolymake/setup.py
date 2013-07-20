# Run as:
#    python setup.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import os
try:
    SAGE_LOCAL = os.environ['SAGE_LOCAL']
    SAGE_ROOT = os.environ['SAGE_ROOT']
except KeyError:
    print("pypolymake installation failed!")
    print("SAGE_LOCAL and SAGE_ROOT should be defined.")
    print("Perhaps you are not inside a Sage shell.")

import sys, os.path as path
sys.path = [SAGE_ROOT+'/devel/sage'] + sys.path

ext_modules = cythonize("*.pyx", cpp=True, cplus=True, language='c++',
        include_path=[SAGE_ROOT+'/devel/sage'])
polymake = ext_modules[0]

#polymake = Extension('polymake',
#        ['polymake.pyx'],
#        language='c++',
#        include_dirs=[SAGE_LOCAL+'/include'],
#        libraries=['polymake', 'gmp'],
#        cmd_class = {'build_ext':build_ext}
#        )

setup(
  name = 'pypolymake',
  #cmdclass = {'build_ext': build_ext},
  ext_modules = [polymake],
        include_dirs=[path.join(SAGE_LOCAL,'include'),
            path.join(SAGE_LOCAL,'include','csage')]
  )
