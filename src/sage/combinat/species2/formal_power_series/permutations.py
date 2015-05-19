# -*- coding: utf-8 -*-
"""
Formal power series of the species of permutations

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
from sage.combinat.partition import cached_number_of_partitions
from sage.combinat.species2.formal_power_series import FPS
from sage.rings.arith import factorial


class PermutationsEPS(FPS):
    """
    Formal power series of the structures of species of permutations

    TESTS::

        sage: from sage.combinat.species2.cycle_index_series.permutations import PermutationsCIS
        sage: TestSuite(PermutationsCIS()).run()

    """
    
    def coefficient(self, n):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.permutations import PermutationsCIS
            sage: ZP = PermutationsCIS()
            sage: list(ZP.generating_series().coefficients(6))
            [1, 1, 2, 6, 24, 120, 720]
            sage: list(ZP.type_generating_series().coefficients(19))
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490]

        """
        return factorial(n)

    def _repr_(self):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.permutations import PermutationsCIS
            sage: PermutationsFPS()
            ZP

        """
        return "P"

    def _valuation_(self):
        """
        The valuation of the permutations generating series is `0`.
        """
        return 0

class PermutationsOPS(FPS):
    """
    Formal power series of the isomorphism type of species of permutations
    """

    def coefficient(self, n):
        """
        TESTS::


        """
        return cached_number_of_partitions(n)