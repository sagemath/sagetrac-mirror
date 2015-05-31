# -*- coding: utf-8 -*-
"""
The species `C`, of *cycles*.

References
----------

 _[BBL] Combinatorial species and tree-like structures,
 François Bergeron, Gilbert Labelle and Pierre Leroux,
 1998, Cambridge University Press

"""
# *******************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# *******************************************************************************
from itertools import permutations
from sage.combinat.species2 import SpeciesDesign
from sage.combinat.species2.cycle_index_series.cycles import CyclesCIS
from sage.combinat.structures import Structure
from sage.misc.misc import uniq
from sage.sets.set import Set
from sage.structure.list_clone import ClonableArray


class Cycle(Structure, ClonableArray):

    def __init__(self, parent=None, cyc=None):

        m = cyc[:].index(min(cyc))
        cyc = tuple(cyc[m:] + cyc[:m])
        ClonableArray.__init__(self, parent, cyc)

        Structure.__init__(self, parent)

    def check(self):
        assert(len(uniq(self)) == len(self)), "all labels should be distinct"


class Cycles(SpeciesDesign):

    def grading(self, cyc):
        return len(cyc)

    def transport(self, sigma):
        return lambda cyc: self._element_constructor_(map(sigma, cyc))

    def _repr_(self):
        return "Cycles"
    
    def cycle_index_series(self):
        return CyclesCIS()

    Element = Cycle

    class Structures(SpeciesDesign.Structures):

        def __iter__(self):
            U = self.finite_set()

            if U.cardinality() == 0:
                #yield self._element_constructor_(())
                return

            mi = min(U)
            S = U.difference(Set([mi]))
            for sig in permutations(S):
                yield self._element_constructor_((mi,) + sig)