r"""
Library of Interesting Groups

Type ``groups.matrix.<tab>`` to access examples
of groups implemented as permutation groups.
"""
from __future__ import absolute_import

# groups imported here will be available
# via  groups.matrix.<tab>
#
# Do not use this file for code
#
# If you import a new group, then add an
# entry to the list in the module-level
# docstring of groups/groups_catalog.py

from .all import GL, SL, Sp, SU, GU, SO, GO
from .all import QuaternionMatrixGroupGF3 as QuaternionGF3

