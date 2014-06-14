"""
Context manager for chdir
"""

import contextlib
import os
 
@contextlib.contextmanager
def chdir(dirname=None):
    cwd = os.getcwd()
    try:
        if dirname is not None:
            os.chdir(dirname)
        yield
    finally:
        os.chdir(cwd)
