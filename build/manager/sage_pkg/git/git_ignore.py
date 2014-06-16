"""
Interface to dot-gitignore
"""


import os
import glob
from sage_pkg.config import config
from sage_pkg.logger import logger
from sage_pkg.git.util import full_split_repo_path
from .gitignore_pattern import GitIgnorePattern


_gitignore = dict()

def get_gitignore(path_components):
    r"""
    Return the content of the gitignore file

    INPUT:

    - ``path_component`` -- tuple of strings. Nested directory names
      relative to the repostiory root.

    EXAMPLES::

        >>> from sage_pkg.git.git_ignore import get_gitignore
        >>> get_gitignore(())   # doctest: +ELLIPSIS
        (gitignore:...
        >>> get_gitignore(('src', 'c_lib'))   # doctest: +ELLIPSIS
        (gitignore:...
    """
    global _gitignore
    assert isinstance(path_components, tuple)
    try:
        return _gitignore[path_components]
    except KeyError:
        pass
    from sage_pkg.config import config
    path = [config.path.root] + list(path_components) + ['.gitignore']
    filename = os.path.join(*path)
    try:
        with open(filename, 'r') as f:
            gitignore = f.read()
    except (OSError, IOError):
        gitignore = ''
    result = []
    root = os.path.join(config.path.root, *path_components)
    for line in gitignore.splitlines():
        pattern = GitIgnorePattern(line, root=root)
        if pattern:
            result.append(pattern)
    result = tuple(result)
    _gitignore[path_components] = result
    return result


def is_ignored(name):
    """
    Return whether ``name`` is matched by ``.gitignore``
    
    INPUT:

    - ``name`` -- string. File name. Either relative to the
      repository root or absolute and inside the repo.

    OUTPUT:

    Boolean.

    EXAMPLES::

        >>> from sage_pkg.git.git_ignore import is_ignored
        >>> is_ignored('file.py')
        False
        >>> is_ignored('dir/file.py')
        False
        >>> is_ignored('dir/file.pyc')
        True
        >>> is_ignored('dir/file.py~')
        True
        >>> is_ignored('src/c_lib/libcsage.so')
        True
        >>> is_ignored(os.path.join(config.path.root, 'src/c_lib/libcsage.so'))
        True
    """
    parts = full_split_repo_path(name)
    fullname = os.path.join(config.path.root, *parts)    # includes directory name relative to root
    path_components = ()
    for part in parts:
        for pattern in get_gitignore(path_components):
            if pattern.matches(fullname):
                return True
        path_components += (part,)
    logger.info('file is not ignored: %s', name)
    return False
