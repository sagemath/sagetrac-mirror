# -*- coding: utf-8 -*-
"""
Cycle index series of the species of compositions

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
from sage.rings.arith import multinomial


class CompositionsCIS(CIS):
    """
    Cycle index series of the species of compositions

    TESTS::

        sage: from sage.combinat.species2.cycle_index_series.compositions import CompositionsCIS
        sage: TestSuite(CompositionsCIS()).run()

    """

    def Frobenius_characteristic(self, n):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.compositions import CompositionsCIS
            sage: ZC = CompositionsCIS()
            sage: list(ZC.generating_series().coefficients(6))
            [1, 1, 3, 13, 75, 541, 4683]
            sage: list(ZC.type_generating_series().coefficients(6))
            [1, 1, 2, 4, 8, 16, 32]

        """
        h = self._sym_.h()
        return sum(multinomial(lambd.to_exp()) * h(lambd)
                   for lambd in Partitions(n))

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.species2.cycle_index_series.compositions import CompositionsCIS
            sage: CompositionsCIS()
            ZC

        """
        return "ZC"
