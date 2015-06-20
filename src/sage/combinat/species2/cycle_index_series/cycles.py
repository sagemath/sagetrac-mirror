# -*- coding: utf-8 -*-
"""
Cycle index series of the species of cyclic permutations

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
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2.cycle_index_series.operations.partitional_composite import ZeroCIS
from sage.combinat.species2.formal_power_series.cycles import CyclesEPS, CyclesOPS
from sage.rings.arith import euler_phi


class CyclesCIS(CIS):
    """

    Cycle index series of the species of cycles

    """

    def __init__(self):
        CIS.__init__(self)
        self._gen = self._gen_()
        self._gen.send(None)

    def _gen_(self):

        p = self._sym_.p()


        self._index_ = [ZeroCIS(p, 0)]
        self._last_index_ = 0

        while True:
            k = (yield)
            while self._last_index_ < k:
                self._last_index_ += 1
                self._index_.append(EulerCIS(self._last_index_))

    def Frobenius_characteristic(self, n):
        """
        MATH::

            Z_C(p_1, p_2, \cdots) = \sum_{k\geqslant 1}
            \frac{\phi(k)}{k}\log\left(
                \frac{1}{1 - p_k}
            \right) = \sum_{k, m \geqslant 1} \frac{\phi(k)p_k^m}{km}\,.

        """
        self._gen.send(n)
        return sum(self._index_[k].Frobenius_characteristic(n)
                   for k in range(n+1))

    def _repr_(self):
        """
        TEST::

        """
        return "ZC"

    def generating_series(self):
        return CyclesEPS()

    def type_generating_series(self):
        return CyclesOPS()

class EulerCIS(CIS):

    def __init__(self, k):
        CIS.__init__(self)
        self._k_ = k
        self._k_euler = euler_phi(k)

    def Frobenius_characteristic(self, n):
        p = self._sym_.p()
        if n % self._k_ == 0:
            m = n / self._k_
            return p[self._k_]**m * self._k_euler / (self._k_ * m)
        return p.zero()