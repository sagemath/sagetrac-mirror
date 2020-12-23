from Cython.Build.Dependencies import default_create_extension
import numpy

def create_extension(template, kwds):
    # Add numpy and source folder to the include search path used by the compiler
    # This is a workaround for https://github.com/cython/cython/issues/1480

    include_dirs = kwds.get('include_dirs', []) + [numpy.get_include(), 'src', '.']
    kwds['include_dirs'] = include_dirs
    return default_create_extension(template, kwds)
