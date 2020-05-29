r"""
Catalog of Path Tableaux

The ``path_tableaux`` object may be used to access examples of various algebras
currently implemented in Sage. Using tab-completion on this object is an
easy way to discover and quickly create the path tableaux that are available
(as listed here).

Let ``<tab>`` indicate pressing the tab key.  So begin by typing
``algebras.<tab>`` to the see the currently implemented named path tableaux.

- :class:`algebras.ArikiKoike
"""

from sage.combinat.path_tableaux.dyck_path import DyckPaths, DyckPath

from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.path_tableaux.dyck_path', 'DyckPaths', 'DyckPath')