# -*- coding: utf-8 -*-
"""
(Partitional) composite of cycle index series

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
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.cycle_index_series import CIS
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.combinat.species2.cycle_index_series.operations.product import Prod
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.misc_c import prod
from sage.rings.integer import Integer


class Composite(CIS):
    """
    (Partitional) composite/plethystic substitution of cycle index series

    MATH::

        Z_{F \circ G}(p_1, p_2, p_3, \cdots) = Z_F(Z_G(p_1, p_2, p_3, \cdots), Z_G(p_2, p_4, p_6, \cdots), \cdots)\,.

    Properties:

     - Associative: `Z_F(Z_G(Z_H)) = Z_F(Z_G)(Z_H)`,
     - Distributive: `(Z_F \cdot Z_G)(Z_H) = Z_F(Z_H) \cdot Z_G(Z_H)`,
     - Distributive: `(Z_F + Z_G)(Z_H) = Z_F(Z_H) + Z_G(Z_H)`,
     - Neutral element: `Z_F(Z_X) = Z_X(Z_F) = Z_F`.

    EXAMPLE::

        sage: ZP = Permutations().cycle_index_series()
        sage: ZS = SetPartitions().cycle_index_series().restriction(min=1)
        sage: ZPS = ZP(ZS); ZPS
        ZP(Z{Set..}(≥1))
        sage: ZP.composite(ZS)
        ZP(Z{Set..}(≥1))
        sage: ZPS.Frobenius_characteristic(3)
        23/6*p[1, 1, 1] + 9/2*p[2, 1] + 5/3*p[3]

    TESTS::

        sage: ZP = Permutations().cycle_index_series()
        sage: ZS = SetPartitions().cycle_index_series().restriction(min=1)
        sage: from sage.combinat.species2.cycle_index_series.compositions import CompositionsCIS
        sage: ZC = CompositionsCIS().restriction(min=1)
        sage: ZP(ZS(ZC)) == ZP(ZS)(ZC)
        True

        sage: ZS = SetPartitions().cycle_index_series()
        sage: (ZP * ZS)(ZC)
        ZP(ZC(≥1))⋅Z{Set..}(ZC(≥1))

        sage: (ZP + ZS)(ZC)
        ZP(ZC(≥1)) + Z{Set..}(ZC(≥1))

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZX = CycleIndexSeries().singletons()
        sage: ZP(ZX)
        ZP
        sage: ZX(ZP)
        ZP

        sage: TestSuite(ZX(ZP)).run()

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):
        ZF, ZG = args
        # neutral element
        if ZF == CycleIndexSeries().singletons():
            return ZG
        if ZG == CycleIndexSeries().singletons():
            return ZF
        # associative: ZF(ZG(ZH)) |--> ZF(ZG)(ZH)
        if isinstance(ZG, Composite):
            ZG, ZH = ZG._ZF_, ZG._ZG_
            return Composite(Composite(ZF, ZG), ZH)
        # distributive 1
        if isinstance(ZF, Prod):
            return Prod(*[(Composite(ZH, ZG), nh) for (ZH, nh) in ZF._dic_cis_.iteritems()])
        # distributive 2
        if isinstance(ZF, Add):
            return Add(*[(Composite(ZH, ZG), nh) for (ZH, nh) in ZF._dic_cis_.iteritems()])
        # # otherwise ##
        return super(Composite, cls).__classcall__(cls, ZF, ZG)

    def __init__(self, ZF, ZG):
        assert(ZG._valuation_() > 0), "The cycle index series should satisfy: `[0]ZG = 0`."
        CIS.__init__(self)
        self._ZF_ = ZF
        self._ZG_ = ZG

        self._gen = self._gen_()
        self._gen.send(None)

    def _gen_(self):
        ZF = self._ZF_
        ZG = self._ZG_
        p = SymmetricFunctions(ZF.Frobenius_characteristic(0).base_ring()).p()

        self._index_ = [ZeroCIS(p, ZF.Frobenius_characteristic(0))] * (self._valuation_()+1)
        self._last_index_ = self._valuation_()
        
        while True:
            n = (yield)
            while self._last_index_ < n:
                self._last_index_ += 1
                self._index_.append(CycleIndexSeries().zero())

                for pi, coef in p(ZF.Frobenius_characteristic(self._last_index_)):
                    if coef == p.zero():
                        continue
                    acc = ZeroCIS(p, coef)
                    acc *= prod(StretchCIS(ZG, part) for part in pi)
                    self._index_[-1] += acc

    @cached_method
    def Frobenius_characteristic(self, n):
        self._gen.send(n)
        return sum(self._index_[k].Frobenius_characteristic(n) for k in range(n+1))

    def generating_series(self):
        return self._ZF_.generating_series().substitution(self._ZG_.generating_series())

    # TODO: Implement `\tilde{F\circ G}(t) = Z_{F}(\tilde{G}(t), \tilde{G}(t^2), \cdots)`
    # more efficient or not? I suppose

    def _repr_(self):
        return repr(self._ZF_) + "(" + repr(self._ZG_) + ")"


class StretchCIS(CIS):

    def __init__(self, ZF, k):
        CIS.__init__(self)
        self._ZF_ = ZF
        self._k_ = Integer(k)
        ###
        self._p_ = SymmetricFunctions(ZF.Frobenius_characteristic(0).base_ring()).p()

    @cached_method
    def Frobenius_characteristic(self, n):
        k = self._k_
        p = self._p_
        ZF = self._ZF_

        return p[k].plethysm(ZF.Frobenius_characteristic(int(n/k))) \
            if k.divides(n) else p.zero()


class ZeroCIS(CIS):

    def __init__(self, p, coeff):
        CIS.__init__(self)
        self._coeff_ = coeff
        self._p_ = p

    def Frobenius_characteristic(self, n):
        return self._p_(self._coeff_) \
            if n == 0 else self._p_.zero()
