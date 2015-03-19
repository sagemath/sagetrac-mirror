from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.partition import Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class CISSets(UniqueRepresentation, Parent):

    def __init__(self):
        Parent.__init__(self, category=CycleIndexSeries())

    def Frobenius_characteristic(self, n):
        h = SymmetricFunctions(QQ).h()
        return h.monomial(Partition([n]))

    def _repr_(self):
        return "Z_E"