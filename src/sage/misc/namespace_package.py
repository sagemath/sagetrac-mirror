"""
Utility functions for namespace packages in Sage
"""
from importlib import import_module
import os, glob

def install_doc(package, doc):
    """
    Install the docstring ``doc`` to the package.

    TESTS:

        sage: from sage.misc.namespace_package import install_doc
        sage: install_doc('sage', 'hello')
        sage: from inspect import getdoc
        sage: getdoc(sage)
        'hello'
    """
    pkg = import_module(package)
    pkg.__doc__ = doc         # enable sage.package?
    pkg.getdoc = lambda: doc  # enable help(sage.package)


def is_package_or_sage_namespace_package_dir(path):
    """
    Return whether ``path`` is a directory that contains a Python package.

    Ordinary Python packages are recognized by the presence of `__init__.py`.

    Implicit namespace packages (PEP 420) are only recognized if they
    follow the conventions of the Sage library, i.e., the directory contains
    a file ``all.py`` or a file matching the pattern ``all__*.py``
    such as ``all__sagemath_categories.py``.

    EXAMPLES:

    :mod:`sage.cpython` is an ordinary package::

        sage: from sage.misc.namespace_package import is_package_or_sage_namespace_package_dir
        sage: d = os.path.dirname(sage.cpython.__file__); d
        '.../sage/cpython'
        sage: is_package_or_sage_namespace_package_dir(d)
        True

    :mod:`sage` is designated to become an implicit namespace package::

        sage: d = os.path.dirname(sage.env.__file__); d
        '.../sage'
        sage: is_package_or_sage_namespace_package_dir(d)
        True

    Not a package::

        sage: d = os.path.join(os.path.dirname(sage.symbolic.__file__), 'ginac'); d
        '.../sage/symbolic/ginac'
        sage: is_package_or_sage_namespace_package_dir(d)
        False
    """
    if os.path.exists(os.path.join(path, '__init__.py')):  # ordinary package
        return True
    if os.path.exists(os.path.join(path, 'all.py')):       # complete namespace package
        return True
    for _ in glob.iglob(os.path.join(path, 'all__*.py')):
        return True                                        # partial namespace package
    return False
