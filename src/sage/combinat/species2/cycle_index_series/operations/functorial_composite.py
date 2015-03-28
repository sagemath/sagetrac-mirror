# -*- coding: utf-8 -*-
"""
Functorial composite of cycle index series

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
from sage.combinat.partition import Partition, Partitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.combinat.species2.cycle_index_series.operations.hadamard_product import HadamardProduct
from sage.combinat.species2.cycle_index_series.operations.product import Prod
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.arith import lcm, divisors, moebius
from sage.rings.integer import Integer


class FunctorialComposite(CIS):
    """
    Functorial composite of cycle index series.

    MATH::

        (Z_F \Box Z_G)(p_1, p_2, p_3, \cdots) = ??

    Properties:

     - Associative: `Z_F \Box (Z_G \Box Z_H) = (Z_F \Box Z_G) \Box Z_H`,
     - Neutral element: `Z_E^\bullet (= Z_X Z_E')`,
     - Distributive on the right by the operation of hadamard product:
     `(Z_F \times Z_G) \Box Z_H = (Z_F \Box Z_H) \times (Z_G \Box Z_H)`,
     - Distributive on the right by the sum of cycle index series:
     `(Z_F + Z_G) \Bod Z_H = Z_F \Box Z_H + Z_G \Box Z_H`.


    The following test verifies that the property

    MATH::

        F \times F = ((X + X^2) \cdot E) \Box F

    is satisfied.

    EXAMPLES::

        sage: ZP = Permutations().cycle_index_series()
        sage: ZS = SetPartitions().cycle_index_series()
        sage: ZPS = ZP.functorial_composite(ZS); ZPS
        ZP▢ Z{Set..}
        sage: ZP[ZS]
        ZP▢ Z{Set..}
        sage: ZPS.Frobenius_characteristic(3)
        20*p[1, 1, 1] + 3*p[2, 1] + 2/3*p[3]

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZP = Permutations().cycle_index_series()
        sage: ZPP = ZP.Hadamard_product(ZP)
        sage: CIS = CycleIndexSeries()
        sage: ZE = CIS.sets()
        sage: ZX = CIS.singletons()
        sage: ZX2 = ZX*ZX
        sage: ZT = ((ZX + ZX2) * ZE).functorial_composite(ZP)
        sage: for n in range(4): assert(ZT.Frobenius_characteristic(n) == ZPP.Frobenius_characteristic(n))

        sage: ZP.functorial_composite(ZC).functorial_composite(ZE) == ZP.functorial_composite(ZC.functorial_composite(ZE))
        True

        sage: ZP.functorial_composite(ZE.pointing()) == ZP
        True
        sage: ZE.pointing().functorial_composite(ZP) == ZP
        True

        sage: ZP.Hadamard_product(ZC).functorial_composite(ZE) == (ZP.functorial_composite(ZE)).Hadamard_product(ZC.functorial_composite(ZE))
        True

        sage: TestSuite(ZP.functorial_composite(ZC)).run()

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, ZF, ZG):
        # associativity
        if isinstance(ZG, FunctorialComposite):
            return FunctorialComposite(FunctorialComposite(ZF, ZG._ZF_), ZG._ZG_)
        # neutral element : E• = X⋅E'
        if isinstance(ZF, Prod) and ZF.is_pointing_of(CycleIndexSeries().sets()):
            return ZG
        if isinstance(ZG, Prod) and ZG.is_pointing_of(CycleIndexSeries().sets()):
            return ZF
        # distributive right x
        if isinstance(ZF, HadamardProduct):
            return HadamardProduct(*tuple((FunctorialComposite(ZH, ZG), nh)
                                          for (ZH, nh) in ZF._dic_cis_.iteritems()))
        # distributive right +
        if isinstance(ZF, Add):
            return Add(*tuple((FunctorialComposite(ZH, ZG), nh)
                              for (ZH, nh) in ZF._dic_cis_.iteritems()))
        # otherwise
        return super(FunctorialComposite, cls).__classcall__(cls, ZF, ZG)

    def __init__(self, ZF, ZG):
        assert(ZG._valuation_() > 0)
        CIS.__init__(self)
        self._ZF_, self._ZG_ = ZF, ZG
        self._p_ = SymmetricFunctions(ZF.Frobenius_characteristic(0).base_ring()).p()

    def cycle_type_of_Gsigma(self, sigma):
        """
        This method compute the cycle `[(G[\sigma])_1, (G[\sigma])_2, (G[\sigma])_3, \cdots]` of the permutation
        `G[\sigma]` where:

        MATH::

            (G[\sigma])_k = \frac{1}{k} \sum_{d | k} \mu\left( \frac{k}{d} \right) \mathtt{fix}\, G[\sigma^d]\,.

        (section 2.2 - Proposition 3, _[BBL])

        NOTE: there `\sigma` is a partition.

        TESTS::

            sage: ZP = SetPartitions().cis()
            sage: from sage.combinat.species2.cycle_index_series.compositions import CompositionsCIS
            sage: ZC = CompositionsCIS()
            sage: ZCP = ZC.functorial_composite(ZP)
            sage: [ZCP.cycle_type_of_Gsigma(pi) for pi in Partitions(3)]
            [[3, 1, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]

        """
        def fix(ZG, sigma, kd):
            tau = sigma.power(kd)
            return p(ZG.Frobenius_characteristic(sigma.size())).coefficient(tau) * tau.aut()

        p = self._p_
        ZG = self._ZG_

        res = []
        for k in range(1, min(ZG.generating_series().coefficient(sigma.size()), lcm(sigma[:])) + 1):
            acc = sum(moebius(k/d) * fix(ZG, sigma, d)
                      for d in divisors(k)
                     ) / Integer(k)
            res.extend([k]*int(acc))
        res.reverse()
        return Partition(res)

    def _repr_(self):
        return repr(self._ZF_) + "▢ " + repr(self._ZG_)

    def generating_series(self):
        """
        TESTS::

            sage: Sp = Species()
            sage: E = Sp.sets()
            sage: G = (E*E).functorial_composite(E.restricted(min=2, max=2)*E) # simple graphs
            sage: G.generating_series()
        """
        return self._ZF_.generating_series().functorial_composite(self._ZG_.generating_series())

    def Frobenius_characteristic(self, n):
        """

        TESTS::


            sage: m = SymmetricFunctions(QQ).m()
            sage: Sp = Species()
            sage: E = Sp.sets()
            sage: G = (E*E).functorial_composite(E.restricted(min=2, max=2)*E) # simple graphs
            sage: for n in range(6): m(G.cis().Frobenius_characteristic(n))
            m[]
            m[1]
            2*m[1, 1] + 2*m[2]
            8*m[1, 1, 1] + 6*m[2, 1] + 4*m[3]
            64*m[1, 1, 1, 1] + 40*m[2, 1, 1] + 28*m[2, 2] + 20*m[3, 1] + 11*m[4]
            1024*m[1, 1, 1, 1, 1] + 576*m[2, 1, 1, 1] + 336*m[2, 2, 1] + 240*m[3, 1, 1] + 148*m[3, 2] + 90*m[4, 1] + 34*m[5]

        """
        p = self._p_
        res = p.zero()
        ZFm = p(self._ZF_.Frobenius_characteristic(self.cycle_type_of_Gsigma(Partitions(n).first()).size()))
        for sigma in Partitions(n):
            tau = self.cycle_type_of_Gsigma(sigma)
            fixFtau = ZFm.coefficient(tau) * tau.aut()
            res += fixFtau / sigma.aut() * p(sigma)
        return res