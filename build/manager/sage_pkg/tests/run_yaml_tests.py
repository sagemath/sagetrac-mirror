"""
Run the PyYaml tests
"""

import sys
import os

# The Package Manager root directory
PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PY2 = bool(sys.version_info[0] == 2)

def run_PyYaml_tests(PyYaml_dir):
    PyYaml_dir = os.path.abspath(PyYaml_dir)
    if PY2:
        compat = os.path.join(PKG_DIR, 'compat', 'python-2')
        tests = os.path.join(PyYaml_dir, 'tests', 'lib')
    else:
        compat = os.path.join(PKG_DIR, 'compat', 'python-3')
        tests = os.path.join(PyYaml_dir, 'tests', 'lib3')
    sys.path = [compat, tests] + list(sys.path)
    cwd = os.getcwd()
    try:
        os.chdir(PyYaml_dir)
        __import__('test_all')
        sys.modules['test_all'].main()
    finally:
        os.chdir(cwd)

