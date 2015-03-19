from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class CIS(UniqueRepresentation, Parent):

    def __init__(self):
        Parent.__init__(self, category=CycleIndexSeries())
        self._sym_ = SymmetricFunctions(QQ)