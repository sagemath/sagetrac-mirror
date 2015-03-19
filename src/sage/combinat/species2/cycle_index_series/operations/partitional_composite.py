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
from __builtin__ import staticmethod
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.combinat.species2.cycle_index_series.operations.product import Prod
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class Composite(UniqueRepresentation, Parent):
    """
    (Partitional) composite of cycle index series

    Properties:

     - Associative: `Z_F(Z_G(Z_H)) = Z_F(Z_G)(Z_H)`,
     - Distributive: `(Z_F \cdot Z_G)(Z_H) = Z_F(Z_H) \cdot Z_G(Z_H)`,
     - Distributive: `(Z_F + Z_G)(Z_H) = Z_F(Z_H) + Z_G(Z_H)`,
     - Neutral element: `Z_F(Z_X) = Z_X(Z_F) = Z_F`.

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        ZF, ZG = args
        # neutral element
        if ZF == CycleIndexSeries().singleton():
            return ZG
        if ZG == CycleIndexSeries().singleton():
            return ZF
        # associative: ZF(ZG(ZH)) |--> ZF(ZG)(ZH)
        if isinstance(ZG, Composite):
            ZG = ZG._F_
            ZH = ZG._G_
            return Composite(Composite(ZF, ZG), ZH)
        # distributive 1
        if isinstance(ZF, Prod):
            return Prod((Composite(ZH, ZG), nh) for (ZH, nh) in ZF._dic_cis_.iteritems())
        # distributive 2
        if isinstance(ZF, Add):
            return Add((Composite(ZH, ZG), nh) for (ZH, nh) in ZF._dic_cis_.iteritems())
        ## otherwise ##
        return super(Composite, cls).__classcall__(cls, ZF, ZG)

    def __init__(self, ZF, ZG):
        ZG0 = ZG.Frobenius_characteristic(0)
        assert(ZG0 == ZG0.parent().zero()), "The cycle index series should satisfy: `[0]Z_G = 0`."
        Parent.__init__(self, category=CycleIndexSeries())
        self._ZF_ = ZF
        self._ZG_ = ZG

        self._gen = self._gen_()
        self._gen.send(None)

    def _gen_(self):
        ZF = self._ZF_
        ZG = self._ZG_
        p = SymmetricFunctions(ZF.Frobenius_characteristic(0).base_ring()).p()

        self._index_ = [ZeroCIS(p, ZF.Frobenius_characteristic(0))]
        self._last_index_ = 0
        
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

    def _repr_(self):
        return repr(self._ZF_) + "(" + repr(self._ZG_) + ")"

class StretchCIS(UniqueRepresentation, Parent):

    def __init__(self, ZF, k):
        Parent.__init__(self, category=CycleIndexSeries())
        self._ZF_ = ZF
        self._k_ = Integer(k)
        ###
        self._p_ = SymmetricFunctions(ZF.Frobenius_characteristic(0).base_ring()).p()


    @cached_method
    def Frobenius_characteristic(self, n):
        k = self._k_
        p = self._p_
        ZF = self._ZF_

        return p[k].plethysm(ZF.Frobenius_characteristic(n/k)) \
               if k.divides(n) else p.zero()

class ZeroCIS(UniqueRepresentation, Parent):

    def __init__(self, p, coeff):
        Parent.__init__(self, category=CycleIndexSeries())
        self._coeff_ = coeff
        self._p_ = p

    def Frobenius_characteristic(self, n):
        return self._p_(self._coeff_) \
               if n == 0 else self._p_.zero()
