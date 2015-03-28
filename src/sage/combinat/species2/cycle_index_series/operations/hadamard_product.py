# -*- coding: utf-8 -*-
"""
Hadamard product of cycle index series

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
from sage.categories.formal_power_series import ExponentialPowerSeries
from sage.combinat.partition import Partitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.misc import exists
from sage.rings.integer import Integer


class HadamardProduct(CIS):
    """
    The *Hadamard product* of cycles index series

    MATH::

        (Z_F \times Z_G)(p_1, p_2, p_3, \cdots) = Z_F(p_1, p_2, p_3, \cdots) \times Z_G(p_1, p_2, p_3, \cdots)\,.

    Properties:

     - Commutative: `Z_F \times Z_G = Z_G \times Z_F`,
     - Absorbing element: ` Z_0 \times Z_F = Z_F \times Z_0 = Z_0`,
     - Neutral element: `Z_E \times Z_F = Z_F \times Z_E = Z_F`,
     - Distributive: `(Z_F + Z_G) \times Z_H = Z_F \times Z_H + Z_G \times Z_H`.

    EXAMPLE::

        sage: ZP = Permutations().cycle_index_series()
        sage: ZS = SetPartitions().cycle_index_series()
        sage: ZPS = ZP.hadamard_product(ZS); ZPS
        Z{Set..} x ZP
        sage: ZP.Frobenius_characteristic(3)
        p[1, 1, 1] + p[2, 1] + p[3]
        sage: ZS.Frobenius_characteristic(3)
        5/6*p[1, 1, 1] + 3/2*p[2, 1] + 2/3*p[3]
        sage: ZPS.Frobenius_characteristic(3)
        5*p[1, 1, 1] + 3*p[2, 1] + 2*p[3]

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: CIS = CycleIndexSeries()
        sage: Z0, ZE, ZX = CIS.zero(), CIS.sets(), CIS.singletons()
        sage: ZP = Permutations().cycle_index_series()
        sage: ZS = SetPartitions().cycle_index_series()

        sage: ZS.hadamard_product(ZP) == ZP.hadamard_product(ZS)
        True

        sage: Z0.hadamard_product(ZP)
        Z0

        sage: ZP.hadamard_product(ZE)
        ZP

        sage: (ZP + ZS).hadamard_product(ZX) == ZP.hadamard_product(ZX) + ZS.hadamard_product(ZX)
        True

        sage: TestSuite(ZP.hadamard_product(ZS)).run()

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):

        # commutativity
        dic_cis = defaultdict(Integer)

        for i, (ZF, nf) in enumerate(args):
            # neutral element
            if ZF == CycleIndexSeries().sets():
                continue
            # absorbing element
            elif ZF == CycleIndexSeries().zero():
                if nf != 0:
                    return CycleIndexSeries().zero()
                # else continue
            # distributive
            elif isinstance(ZF, Add):
                return Add(*tuple((HadamardProduct(*(tuple(dic_cis.items()) + ((ZG, ng),) + args[i+1:])), 1)
                                  for (ZG, ng) in ZF._dic_cis_.iteritems()))
            # associative
            elif isinstance(ZF, HadamardProduct):
                for (ZG, ng) in ZF._dic_cis_.iteritems():
                    dic_cis[ZG] += ng
            # otherwise
            else:
                dic_cis[ZF] += nf

        # # cleaning ##
        for (ZF, nf) in list(dic_cis.items()):
            if nf == 0:
                del dic_cis[ZF]
            if nf < 0:
                raise NotImplementedError("Virtual species I suppose... please implement...")

        if len(dic_cis.keys()) == 0:
            return CycleIndexSeries().sets()
        elif len(dic_cis.keys()) == 1 and dic_cis.values()[0] == 1:
            return dic_cis.keys()[0]
        else:
            return super(HadamardProduct, cls).__classcall__(cls, tuple(dic_cis.items()))

    def __init__(self, dic_cis):
        CIS.__init__(self)
        self._dic_cis_ = dict(dic_cis)
        self._p_ = SymmetricFunctions(self._dic_cis_.keys()[0].Frobenius_characteristic(0).base_ring()).p()

    def Frobenius_characteristic(self, n):
        """
        TESTS::

            sage: from sage.categories.cycle_index_series import CycleIndexSeries
            sage: CIS = CycleIndexSeries()
            sage: ZX = CIS.singletons()
            sage: ZXX = ZX.Hadamard_product(ZX); ZXX
            Z_X^{x2}
            sage: [ZXX.Frobenius_characteristic(n) for n in range(10)]
            [0, p[1], 0, 0, 0, 0, 0, 0, 0, 0]

        """
        if self._valuation_() > n:
            return self._sym_.zero()
        p = self._p_
        list_ch = [(p(ZF.Frobenius_characteristic(n)), nf) for (ZF, nf) in self._dic_cis_.iteritems()]

        if exists(list_ch, lambda (chZF, nf): chZF == p.zero())[0]:
            return p.zero()

        acc = p.zero()
        for pi in Partitions(n):

            accn = 0
            tmp = Integer(1)
            zpi = pi.aut()

            for (chZF, nf) in list_ch:
                coefPi = chZF.coefficient(pi)
                if coefPi == 0:
                    tmp = Integer(0)
                    break
                tmp *= (coefPi * zpi) ** nf
                accn += nf

            acc += tmp / zpi * p.monomial(pi)

        return acc

    @cached_method
    def generating_series(self):
        from sage.combinat.species2.formal_power_series.operations.hadamard_product import HadamardProduct
        return HadamardProduct(map(lambda (ZF, nf): (ZF.generating_series(), nf), self._dic_cis_))

    def _repr_(self):
        return " x ".join(imap(lambda (ZF, nf): repr(ZF) + ("^{x%d}" % nf if nf != 1 else ""),
                               self._dic_cis_.iteritems()))