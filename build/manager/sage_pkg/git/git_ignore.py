"""
Interface to dot-gitignore
"""


import os
import glob
from sage_pkg.git.util import full_split_repo_path



_gitignore = dict()

def get_gitignore(path_components):
    """
    Return the content of the gitignore file

    INPUT:

    - ``path_component`` -- tuple of strings. Nested directory names
      relative to the repostiory root.

    EXAMPLES::

        >>> from sage_pkg.git.git_ignore import get_gitignore
        >>> get_gitignore([])   # doctest: +ELLIPSIS
        ('*.pyc', ...
        >>> get_gitignore(['src', 'c_lib'])   # doctest: +ELLIPSIS
        ('*.os', '*.so', '*.dylib', '.sconsign.dblite')
    """
    global _gitignore
    path_components = tuple(path_components)
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
    for line in gitignore.splitlines():
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith('#'):
            continue
        result.append(line)
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
    filename = parts[-1]               # just the file name
    path_components = []
    fullname = os.path.join(*parts)    # includes directory name relative to root
    for part in parts:
        for pattern in get_gitignore(path_components):
            if pattern.startswith(os.path.sep):
                if glob.fnmatch.fnmatch(fullname, pattern):
                    return True
            else:
                if glob.fnmatch.fnmatch(filename, pattern):
                    return True
        path_components.append(part)
    return False
