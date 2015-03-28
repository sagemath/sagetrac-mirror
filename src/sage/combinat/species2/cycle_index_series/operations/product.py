# -*- coding: utf-8 -*-
"""
Product of cycle index series

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
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.combinat.species2.formal_power_series.operations.product import OrdinaryProd
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.misc_c import prod
from sage.rings.integer import Integer


class Prod(CIS):
    """
    The product of cycle index series.

    MATH::

        Z_{F \cdot G}(p_1, p_2, p_3, \cdots) = Z_F(p_1, p_2, p_3, \cdots) Z_G(p_1, p_2, p_3, \cdots)\,.

    Properties:

     - Neutral element: `Z_F \cdot Z_1 = Z_1 \cdot Z_F = Z_F`,
     - Absorbing element: `Z_F \cdot Z_0 = Z_0 \cdot Z_F = Z_0`,
     - Distributive: `Z_F \cdot (Z_G + Z_G) = Z_F \cdot Z_G + Z_F \cdot Z_G`,
     - Associative: `Z_F \cdot (Z_G \cdot Z_H) = (Z_F \cdot Z_G) \cdot Z_H`,
     - Commutative: `Z_F \cdot Z_G = Z_G \cdot Z_F`.

    EXAMPLE::

        sage: ZP = Permutations().cycle_index_series()
        sage: ZS = SetPartitions().cycle_index_series()
        sage: ZP * ZS
        Z{Set..}⋅ZP
        sage: sum(ZP.Frobenius_characteristic(k) * ZS.Frobenius_characteristic(3 - k) for k in range(4))
        23/6*p[1, 1, 1] + 9/2*p[2, 1] + 5/3*p[3]
        sage: (ZP * ZS).Frobenius_characteristic(3)
        23/6*p[1, 1, 1] + 9/2*p[2, 1] + 5/3*p[3]

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: Z1 = CycleIndexSeries().one()
        sage: ZP = Permutations().cycle_index_series()
        sage: ZP * Z1 == ZP
        True
        sage: Z1 * ZP == ZP
        True

        sage: Z0 = CycleIndexSeries().zero()
        sage: ZP * Z0 == Z0
        True
        sage: Z0 * ZP == Z0
        True

        sage: ZC = SetPartitions().cycle_index_series()
        sage: ZX = CycleIndexSeries().singletons()
        sage: (ZP + ZC) * ZX == ZP * ZX + ZC * ZX
        True

        sage: (ZP * ZC) * ZX == ZP * (ZC * ZX)
        True

        sage: ZP * ZC == ZC * ZP
        True

        sage: TestSuite(ZP * ZC).run()

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):

        # commutativity
        dic_cis = defaultdict(Integer)

        for i, (ZF, nf) in enumerate(args):
            # neutral element
            if ZF == CycleIndexSeries().one():
                continue
            # absorbing element
            elif ZF == CycleIndexSeries().zero():
                if nf != 0:
                    return CycleIndexSeries().zero()
                # else continue
            # distributivity ZF ⋅ (ZG + ZH) ⋅ ZR |--> ZF ⋅ ZG ⋅ ZR + ZF ⋅ ZH ⋅ ZR
            elif isinstance(ZF, Add):
                return Add(*tuple((Prod(*(tuple(dic_cis.items()) + ((ZG, ng),) + args[i+1:])), 1)
                                  for (ZG, ng) in ZF._dic_cis_.iteritems()))
            # associativity ZF ⋅ (ZG ⋅ ZH) |--> ZF ⋅ ZG ⋅ ZH
            elif isinstance(ZF, Prod):
                for (ZG, ng) in ZF._dic_cis_.iteritems():
                    dic_cis[ZG] += ng
            # otherwise
            else:
                dic_cis[ZF] += nf

        # # cleanning ##
        for (ZF, nf) in list(dic_cis.items()):
            if nf == 0:
                del dic_cis[ZF]
            if nf < 0:
                raise NotImplementedError("Virtual species I suppose... please implement...")

        if len(dic_cis.keys()) == 0:
            return CycleIndexSeries().one()
        elif len(dic_cis.keys()) == 1 and dic_cis.values()[0] == 1:
            return dic_cis.keys()[0]
        else:
            return super(Prod, cls).__classcall__(cls, tuple(dic_cis.items()))

    def __init__(self, dic_cis):
        CIS.__init__(self)
        self._dic_cis_ = dict(dic_cis)

    def _repr_(self):
        return "⋅".join(imap(lambda (ZF, nf): repr(ZF) + ("^%d" % nf if nf != 1 else ""), self._dic_cis_.iteritems()))

    @cached_method
    def Frobenius_characteristic(self, n):
        """
        MATH::

            [n](Z_F \cdot Z_G) = \sum_{i + j = n} [i]Z_F \times [j]Z_G

        :param n: a non-negative integer
        """

        def rec_prod(cis, n):
            ZF = cis[0]
            if len(cis) == 1:
                return ZF.Frobenius_characteristic(n)

            acc = h.zero()
            for k in range(ZF._valuation_(), n+1):
                chZFk = ZF.Frobenius_characteristic(k)
                if chZFk != h.zero():
                    acc += chZFk * rec_prod(cis[1:], n-k)
            return acc

        cis = reduce(lambda e, f: e+f, tuple((ZF,)*nf for ZF, nf in self._dic_cis_.iteritems()), ())
        h = cis[0].Frobenius_characteristic(0).parent()

        return rec_prod(cis, n)

    @cached_method
    def generating_series(self):
        from sage.combinat.species2.formal_power_series.operations.product import ExponentialProd
        return ExponentialProd(*map(lambda (ZF, nf): (ZF.generating_series(), nf), self._dic_cis_.iteritems()))

    @cached_method
    def type_generating_series(self):
        from sage.combinat.species2.formal_power_series.operations.product import OrdinaryProd
        return OrdinaryProd(*map(lambda (ZF, nf): (ZF.type_generating_series(), nf), self._dic_cis_.iteritems()))

    def is_pointing_of(self, ZG):
        """
        Test if `Z_F` (*self*) is the pointing of `Z_G`.

        In other terms, this method tests if

        MATH::

            Z_F = Z_X \cdot Z_G'\,.

        :param ZG: a cycle index series

        """
        # TODO: this method is not consistent... should be implement in other operator class...
        from sage.combinat.species2.cycle_index_series.operations.derivative import Derivative
        Zsing = CycleIndexSeries().singletons()
        if len(self._dic_cis_) == 2 and Zsing in self._dic_cis_.keys() and self._dic_cis_[Zsing] == 1:
            ZGd = Derivative(ZG)
            return ZGd in self._dic_cis_.keys() and self._dic_cis_[ZGd] == 1
        return False
