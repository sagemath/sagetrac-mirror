# -*- coding: utf-8 -*-
"""
Default implementation of the cycle index series of a species

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
from collections import defaultdict
from itertools import imap
from sage.combinat.partition import Partitions
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2 import partition_to_permutation
from sage.rings.integer import Integer
from sage.sets.set import Set


class genericCIS(CIS):
    """
    Default implementation of a cycle index series

    Let `F` be a species. This class implements the basic definition of `Z_F`:

    MATH::

        Z_F(p_1, p_2, \cdots) = \sum_{n \geqslant 0} \frac{1}{n!} \left(
                    \sum_{\sigma \in S_n} \mathtt{fix}\, F[\sigma]\, p_1^{\sigma_1} p_2^{\sigma_2} \cdots
                \right) \,.

    (section 1.2, _[BBL])

    """

    def __init__(self, F):
        """
        Default implementation of a cycle index series

        :param F: a species

        """
        CIS.__init__(self)
        self._F_ = F

    def _repr_(self):
        s = repr(self._F_)
        return "Z"+(s if len(s) < 5 else "{" + s[:3] + "..}")

    def Frobenius_characteristic(self, n):
        """
        Compute `ch(Z_F, n)`

        MATH::

            ch(Z_F, n) = \frac{1}{n!} \left(
                    \sum_{\sigma \in S_n} \mathtt{fix}\, F[\sigma]\, p_1^{\sigma_1} p_2^{\sigma_2} \cdots
                \right) \,.

        """
        p = self._sym_.p()
        fix = defaultdict(Integer)
        for struct in self._F_.structures(Set(range(1, n+1))):
            for pi in Partitions(n):
                sigma = partition_to_permutation(pi)
                if struct == self._F_.transport(sigma)(struct):
                    fix[pi] += 1
        return p.sum_of_terms(imap(lambda (pi, c): (pi, c / pi.aut()), fix.iteritems()))
