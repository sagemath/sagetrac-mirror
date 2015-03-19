# -*- coding: utf-8 -*-
"""
Functorial composite of cycle index series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998
"""
#*****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.partition import Partition, Partitions
from sage.combinat.permutation import Permutation
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2 import partition_to_permutation
from sage.combinat.species2.cycle_index_series.operations.hadamard_product import HadamardProduct
from sage.combinat.species2.cycle_index_series.operations.product import Prod
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.arith import lcm, divisors, moebius
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class FunctorialComposite(UniqueRepresentation, Parent):
    """
    Functorial composite of cycle index series.

    MATH::

        (Z_F \Box Z_G)(p_1, p_2, p_3, \cdots) = ??

    Properties:

     - Associative: `Z_F \Box (Z_G \Box Z_H) = (Z_F \Box Z_G) \Box Z_H`,
     - Neutral element: `Z_E^\bullet (= Z_X Z_E')`,
     - Distributive on the right by the operation of hadamard product:
     `(Z_F \times Z_G) \Box Z_H = (Z_F \Box Z_H) \times (Z_G \Box Z_H)`.


    The following test verifies that the property

    MATH::

        F \times F = ((X + X^2) \cdot E) \Box F

    is satisfied.

    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: from sage.combinat.species2.cycle_index_series.cis_compositions import CISCompo
        sage: ZC = CISCompo()
        sage: ZCC = ZC.Hadamard_product(ZC)
        sage: CIS = CycleIndexSeries()
        sage: ZE = CIS.sets()
        sage: ZX = CIS.singletons()
        sage: ZX2 = ZX*ZX
        sage: ZT = ((ZX + ZX2) * ZE).functorial_composition(ZC)
        sage: for n in range(4): assert(ZT.Frobenius_characteristic(n) == ZCC.Frobenius_characteristic(n))

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, ZF, ZG, **opts):
        # associativity
        if isinstance(ZG, FunctorialComposite):
            return FunctorialComposite(FunctorialComposite(ZF, ZG._ZF_), ZG._ZG_)
        # neutral element : E• = X⋅E'
        if isinstance(ZF, Prod) and ZF.is_pointing_of(CycleIndexSeries().sets()):
            return ZG
        if isinstance(ZG, Prod) and ZG.is_pointing_of(CycleIndexSeries().sets()):
            return ZF
        # distributive right
        if isinstance(ZF, HadamardProduct):
            return HadamardProduct(*tuple((FunctorialComposite(ZH, ZG), nh)
                                          for (ZH, nh) in ZF._dic_cis_.iteritems()))
        # otherwise
        return super(FunctorialComposite, cls).__classcall__(cls, ZF, ZG)

    def __init__(self, ZF, ZG):
        Parent.__init__(self, category=CycleIndexSeries())
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
        """

        def pow(sigma, k):
            if k == 1: return sigma
            else: return pow(Permutation([sigma(i) for i in range(1, len(sigma)+1)]), k-1)

        def fix(ZG, sigma, d):
            tau = Partition(map(len, pow(partition_to_permutation(sigma), k).to_cycles()))
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


    def Frobenius_characteristic(self, n):
        p = self._p_
        res = p.zero()
        for sigma in Partitions(n):
            tau = self.cycle_type_of_Gsigma(sigma)
            fixFtau = p(self._ZF_.Frobenius_characteristic(tau.size())).coefficient(tau) * tau.aut()
            res += fixFtau / sigma.aut() * p(sigma)
        return res