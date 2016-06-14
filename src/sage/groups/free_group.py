r"""
TESTS::

    sage: from sage.groups.free_group import FreeGroup
    sage: F = FreeGroup(('a','b'))
    doctest:...: DeprecationWarning: 
    Importing FreeGroup from here is deprecated. If you need to use it,
    please import it directly from sage.groups.free_groups.free_group
    See http://trac.sagemath.org/20154 for details.
    sage: F
    Free Group on generators {a, b}
"""
from sage.misc.lazy_import import lazy_import
lazy_import('sage.groups.free_groups.free_group', '*', deprecation=20154)
