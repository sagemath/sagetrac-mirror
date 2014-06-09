"""
Adjust sys.path
"""


import os
import sys

APP_DIR = os.path.dirname(__file__)

def append_compat(path):
    path = os.path.join(APP_DIR, 'compat', path)
    path = os.path.abspath(path)
    if path not in sys.path:
        sys.path.append(path)


if sys.version_info[0] == 2:
    if sys.version_info[1] == 6:
        append_compat('python-2.6')
    append_compat('python-2')
elif sys.version_info[0] == 3:
    append_compat('python-3')



