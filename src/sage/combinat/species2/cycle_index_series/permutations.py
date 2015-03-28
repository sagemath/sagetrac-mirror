# -*- coding: utf-8 -*-
"""
Cycle index series of the species of permutations

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
from sage.combinat.partition import Partitions
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2.formal_power_series.permutations import PermutationsEPS, PermutationsOPS


class PermutationsCIS(CIS):
    """
    Cycle index series of the species of permutations

    TESTS::

        sage: from sage.combinat.species2.cycle_index_series.permutations import PermutationsCIS
        sage: TestSuite(PermutationsCIS()).run()

    """
    def Frobenius_characteristic(self, n):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.permutations import PermutationsCIS
            sage: ZP = PermutationsCIS()
            sage: list(ZP.generating_series().coefficients(6))
            [1, 1, 2, 6, 24, 120, 720]
            sage: list(ZP.type_generating_series().coefficients(19))
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490]

        """
        p = self._sym_.p()
        return p.sum_of_monomials(Partitions(n))

    def _repr_(self):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.permutations import PermutationsCIS
            sage: PermutationsCIS()
            ZP

        """
        return "ZP"

    def generating_series(self):
        return PermutationsEPS()

    def type_generating_series(self):
        return PermutationsOPS()