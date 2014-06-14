"""
Utility functions
"""

import os
import binascii

from sage_pkg.config import config



def sha_to_binary(sha_ascii):
    assert len(sha_ascii) == 40
    return binascii.unhexlify(sha_ascii)


def sha_to_ascii(sha_binary):
    assert len(sha_binary) == 20
    return binascii.exlify(sha_binary)


def full_split_repo_path(path):
    """
    Utility function to split path into list

    INPUT:

    - ``path`` -- string. Path relative to the repo root or absolute
      path contained in the repo root.

    EXAMPLES::

        >>> import os
        >>> from sage_pkg.git.util import full_split_repo_path
        >>> full_split_repo_path('a/b/c/d')
        ['a', 'b', 'c', 'd']
        >>> full_split_repo_path('a/b')
        ['a', 'b']
        >>> full_split_repo_path('a')
        ['a']
        >>> full_split_repo_path('a/')
        ['a']
        >>> full_split_repo_path(os.path.join(config.path.root, 'a', 'b'))
        ['a', 'b']
    """
    if os.path.isabs(path):
        path = os.path.abspath(path)
        from sage_pkg.config import config
        if path.startswith(config.path.root):
            path = path[len(config.path.root):].lstrip(os.path.sep)
        else:
            raise ValueError('path must be relative or absolute and start with the package path')
    path = path.rstrip(os.path.sep)
    dirs = []
    while len(path) > 0:
        path, directory = os.path.split(path)
        dirs = [directory] + dirs
    return dirs
