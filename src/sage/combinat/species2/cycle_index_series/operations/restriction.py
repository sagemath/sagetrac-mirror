# -*- coding: utf-8 -*-
"""
Restriction of cycle index series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.species2.cycle_index_series import CIS
from sage.rings.infinity import Infinity


class RestrictionCIS(CIS):
    """
    Restriction of cycle index series

    TESTS::

        sage: ZP = Permutations().cycle_index_series().restriction(min=3)
        sage: TestSuite(ZP).run()

    """

    def __init__(self, ZF, min=0, max=Infinity):
        # we suppose the simplification have been done by the species...
        CIS.__init__(self)
        self._ZF_ = ZF
        self._min_ = min
        self._max_ = max

    def Frobenius_characteristic(self, n):
        if n < self._min_ or n > self._max_:
            return CycleIndexSeries().zero().Frobenius_characteristic(n)
        return self._ZF_.Frobenius_characteristic(n)

    def generating_series(self):
        return self._ZF_.generating_series().restricted(min=self._min_, max=self._max_)

    def type_generating_series(self):
        return self._ZF_.type_generating_series().restricted(min=self._min_, max=self._max_)

    def _repr_(self):
        if self._min_ != 0 and self._max_ != Infinity:
            return repr(self._ZF_) + "[" + repr(self._min_) + ", " + repr(self._max_) + "]"
        elif self._min_ != 0:
            return repr(self._ZF_) + "(≥" + repr(self._min_) + ")"
        else:
            return repr(self._ZF_) + "(≤" + repr(self._max_) + ")"