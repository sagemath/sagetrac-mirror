# -*- coding: utf-8 -*-
"""
Formal power series of the species of cycles

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
from sage.categories.formal_power_series import ExponentialPowerSeries, OrdinaryPowerSeries
from sage.combinat.species2.formal_power_series import FPS
from sage.rings.arith import factorial
from sage.rings.integer import Integer


class CyclesEPS(FPS):
    """
    Formal power series of the structures of species of cycles

    TESTS::

        sage: from sage.combinat.species2.cycle_index_series.cycles import CyclesEPS
        sage: TestSuite(PermutationsEPS()).run()

    """

    def __init__(self):
        FPS.__init__(self, category=ExponentialPowerSeries())
    
    def coefficient(self, n):
        """
        TEST::



        """
        if n == 0: return Integer(0)
        else:      return factorial(n-1)

    def _repr_(self):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.cycles import CyclesEPS
            sage: CyclesEPS()
            C

        """
        return "C"

    def valuation(self):
        """
        The valuation of the cycles generating series is `1`.
        """
        return Integer(1)

class CyclesOPS(FPS):
    """
    Formal power series of the isomorphism type of species of permutations
    """

    def __init__(self):
        FPS.__init__(self, category=OrdinaryPowerSeries())

    def coefficient(self, n):
        """
        TESTS::


        """
        if n == 0: return Integer(0)
        else:      return Integer(1)

    def valuation(self):
        """
        The valuation of the cycles generating series is `1`.
        """
        return 1