"""
Type ``groups.presentation.<tab>`` to access examples
of groups implemented as finite presentations (quotients of
free groups).
"""
from __future__ import absolute_import

# groups imported here will be available
# via  groups.presentation.<tab>
#
# Do not use this file for code
#
# If you import a new group, then add an
# entry to the list in the module-level
# docstring of groups/groups_catalog.py

from .finitely_presented_named import DihedralPresentation as Dihedral
from .finitely_presented_named import CyclicPresentation as Cyclic
from .finitely_presented_named import DiCyclicPresentation as DiCyclic
from .finitely_presented_named import FinitelyGeneratedAbelianPresentation as FGAbelian
from .finitely_presented_named import KleinFourPresentation as KleinFour
from .finitely_presented_named import SymmetricPresentation as Symmetric
from .finitely_presented_named import QuaternionPresentation as Quaternion
from .finitely_presented_named import AlternatingPresentation as Alternating
from .finitely_presented_named import BinaryDihedralPresentation as BinaryDihedral


