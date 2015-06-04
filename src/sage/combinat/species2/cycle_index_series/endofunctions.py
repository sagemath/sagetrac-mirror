# -*- coding: utf-8 -*-
"""
Cycle index series of the species of endofunctions

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998

.. [CL89] Calcul combinatoire du nombre d'endofunctions et
  d'arborescences laissées fixes par une permutation,
  Ivan Constantineau and Jacques Labelle
  Annales des sciences mathématiques du Québec, 1989

.. [Lab86] Some new computational methods in the theory of species
  Gilbert Labelle
  Springer, 1986

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from itertools import imap
from sage.combinat.partition import Partitions
from sage.combinat.species2.cycle_index_series import CIS
from sage.misc.misc_c import prod
from sage.rings.arith import divisors
from sage.rings.integer import Integer


class EndofunctionsCIS(CIS):
    """
    Cycle index series of the species of endofunctions

    MATH::

        Z_{\mathtt{End}}(p_1,p_2,p_3, \cdots) =
            \sum_{n \geqslant 0} \frac{1}{n!}\sum_{\sigma \in \S_n}
                \mathtt{fix}\, \mathtt{End}[\sigma] p_1^{\sigma^{(1)}}p_2^{\sigma^{(2)}}\cdots p_n^{\sigma^{(n)}}\,,

    where `\sigma^{(i)}` is a the number of cycles of length `k` in `\sigma`; and

    MATH::

        \mathtt{fix}\, \mathtt{End}[\sigma(\lambda)] =
            \prod_{k=1}^{n} \left(
                \sum_{d | k} d\sigma^{(d)}
            \right)^{\sigma^{(k)}}\,,

    (see _[Lab86] and _[CL89]).

    TESTS::

        sage: from sage.combinat.species2.cycle_index_series.endofunctions import EndofunctionsCIS
        sage: TestSuite(EndofunctionsCIS()).run()

    """
    def Frobenius_characteristic(self, n):
        """

        MATH::

            \mathtt{ch}(Z_{\mathtt{End}}; n) = \sum_{n \geqslant 0} \frac{1}{n!}\sum_{\sigma \in \S_n}
                \mathtt{fix}\, \mathtt{End}[\sigma] p_1^{\sigma^{(1)}}p_2^{\sigma^{(2)}}\cdots p_n^{\sigma^{(n)}}\,,

        where

        MATH::

            \mathtt{fix}\, \mathtt{End}[\sigma(\lambda)] =
                \prod_{k=1}^{n} \left(
                    \sum_{d | k} d\sigma^{(d)}
                \right)^{\sigma^{(k)}}\,,

        (see _[Lab86] and _[CL89]).

        TEST::

            sage: from sage.combinat.species2.cycle_index_series.endofunctions import EndofunctionsCIS
            sage: ZEnd = EndofunctionsCIS()
            sage: ZEnd.Frobenius_characteristic(4)
            32/3*p[1, 1, 1, 1] + 4*p[2, 1, 1] + 2*p[2, 2] + 4/3*p[3, 1] + p[4]

            sage: list(ZEnd.egs().coefficients(9)) # indirect doctest
            [1, 1, 4, 27, 256, 3125, 46656, 823543, 16777216, 387420489]
            sage: list(ZEnd.ogs().coefficients(9)) # indirect doctest
            [1, 1, 3, 7, 19, 47, 130, 343, 951, 2615]

        """
        def fix(pi):
            piexp = pi.to_exp()

            accpr = Integer(1)
            k = 1
            for nk in piexp:
                if nk == 0:
                    k += 1
                    continue

                accpr *= sum(d * piexp[d-1] for d in divisors(k)) ** nk
                k += 1

            return accpr


        p = self._sym_.p()

        return p.sum_of_terms(imap(lambda pi: (pi, fix(pi) / pi.aut()), Partitions(n)))

    def _repr_(self):
        """
        TEST::

            sage: from sage.combinat.species2.cycle_index_series.endofunctions import EndofunctionsCIS
            sage: EndofunctionsCIS()
            ZEnd

        """
        return "ZEnd"
