# -*- coding: utf-8 -*-
"""
Restriction of cycle index series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998
"""
#*****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.rings.infinity import Infinity
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class RestrictionCIS(UniqueRepresentation, Parent):
    """
    Restriction of cycle index series
    """

    def __init__(self, ZF, min=0, max=Infinity):
        # we suppose the simplification have been done by the species...
        Parent.__init__(self, category=CycleIndexSeries())
        self._ZF_ = ZF
        self._min_ = min
        self._max_ = max

    def Frobenius_characteristic(self, n):
        if n < self._min_ or n > self._max_:
            return CycleIndexSeries().zero().Frobenius_characteristic(n)
        return self._ZF_.Frobenius_characteristic(n)