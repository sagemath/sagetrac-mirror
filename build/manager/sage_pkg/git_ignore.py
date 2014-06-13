"""
Interface to dot-gitignore
"""


import os
import glob


_gitignore = None

def get_gitignore():
    """
    Return the content of the gitignore file

    EXAMPLES::

        >>> from sage_pkg.git_ignore import get_gitignore
        >>> get_gitignore()   # doctest: +ELLIPSIS
        ['*.pyc', ...
    """
    global _gitignore
    if _gitignore is not None:
        return _gitignore
    from sage_pkg.config import config
    filename = os.path.join(config.path.root, '.gitignore')
    try:
        with open(filename, 'r') as f:
            gitignore = f.read()
    except OSError:
        gitignore = ''
    _gitignore = []
    for line in gitignore.splitlines():
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        _gitignore.append(line)
    return _gitignore


        
    
def is_ignored(name):
    """
    Return whether ``name`` is matched by ``.gitignore``
    
    INPUT:

    - ``name`` -- string.

    OUTPUT:

    Boolean.

    EXAMPLES::

        >>> from sage_pkg.git_ignore import is_ignored
        >>> is_ignored('file.py')
        False
        >>> is_ignored('dir/file.py')
        False
        >>> is_ignored('dir/file.pyc')
        True
        >>> is_ignored('dir/file.py~')
        True
    """
    for pattern in get_gitignore():
        if glob.fnmatch.fnmatch(name, pattern):
            return True
    return False
