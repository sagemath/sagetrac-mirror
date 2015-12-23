"""
Walk a file tree and read the beginning of each file
"""

import os

from sage.tests.local.test_base import TestException, FileMagicException
from sage.tests.local.binary_audit import ELFBinaryFile


FileClasses = (
    ELFBinaryFile,
)


def recursive_check(root_path):
    """
    Recursively iterate through subdirectories and apply checks

    INPUT:

    - ``root_path`` -- string. The root directory

    OUTPUT:

    Boolean. Whether all tests were successfull. Diagnostics about
    failed tests is printed to stdout.

    EXAMPLES::

        sage: from sage.tests.local.recursive_check import recursive_check
        sage: recursive_check(SAGE_LOCAL)   # long time
        True
    """
    result = True
    for path, dirs, files in os.walk(root_path):
        for filename in files:
            fqn = os.path.join(path, filename)
            if os.path.islink(fqn):
                continue
            if not os.path.isfile(fqn):
                continue
            with open(fqn, 'rb') as f:
                head = f.read(512)
            for cls in FileClasses:
                try:
                    cls(fqn, head).perform_tests()
                except FileMagicException:
                    pass
                except TestException as exc:
                    print('File {0}: {1}'.format(fqn, exc))
                    result = False
    return result

