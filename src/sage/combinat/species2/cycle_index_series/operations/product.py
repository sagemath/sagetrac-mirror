# -*- coding: utf-8 -*-
"""
Product of cycle index series

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
from collections import defaultdict
from itertools import imap
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class Prod(UniqueRepresentation, Parent):
    """
    The product of cycle index series.

    Properties:

     - Neutral element: `Z_F \cdot Z_1 = Z_1 \cdot Z_F = Z_F`,
     - Absorbing element: `Z_F \cdot Z_0 = Z_0 \cdot Z_F = Z_0`,
     - Distributive: `Z_F \cdot (Z_G + Z_G) = Z_F \cdot Z_G + Z_F \cdot Z_G`,
     - Associative: `Z_F \cdot (Z_G \cdot Z_H) = (Z_F \cdot Z_G) \cdot Z_H`,
     - Commutative: `Z_F \cdot Z_G = Z_G \cdot Z_F`.

    TESTS::

        sage: from sage.combinat.species2.cycle_index_series.cis_sets import CISSets
        sage: ZE = CISSets()
        sage: ZEE = ZE * ZE
        sage: [ZEE.generating_series().coefficient(n) for n in range(10)]
        [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZX = CycleIndexSeries().singleton()
        sage: ZX * ZE
        Z_E⋅Z_X
        sage: ZX * ZE == ZE * ZX
        True
        sage: ZE * ZX * ZE == ZE * ZE * ZX
        True

    """
    # TODO test

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **options):

        # commutativity
        dic_cis = defaultdict(Integer)

        for i, (ZF, nf) in enumerate(args):
            # neutral element
            if ZF == CycleIndexSeries().one():
                continue
            # absorbing element
            elif ZF == CycleIndexSeries().zero():
                if nf != 0: return CycleIndexSeries().zero()
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

        ## cleanning ##
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
        Parent.__init__(self, category=CycleIndexSeries())
        self._dic_cis_ = dict(dic_cis)

    def _repr_(self):
        return "⋅".join(imap(lambda (ZF, nf): repr(ZF) + ("^%d"%nf if nf != 1 else ""), self._dic_cis_.iteritems()))

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
            for k in range(n+1):
                chZFk = ZF.Frobenius_characteristic(k)
                if chZFk != h.zero():
                    acc += chZFk * rec_prod(cis[1:], n-k)
            return acc

        cis = reduce(lambda e,f: e+f, tuple((ZF,)*nf for ZF, nf in self._dic_cis_.iteritems()), ())
        h = cis[0].Frobenius_characteristic(0).parent()

        return rec_prod(cis, n)

    def is_pointing_of(self, ZG):
        """
        Test if `Z_F` (*self*) is the pointing of `Z_G`.

        In other terms, this method tests if

        MATH::

            Z_F = Z_X \cdot Z_G'\,.

        :param ZG: a cycle index series

        """
        from sage.combinat.species2.cycle_index_series.operations.derivative import Derivative
        Zsing = CycleIndexSeries().singletons()
        if len(self._dic_cis_) == 2 and Zsing in self._dic_cis_.keys() and self._dic_cis_[Zsing] == 1:
            ZGd = Derivative(ZG)
            return ZGd in self._dic_cis_.keys() and self._dic_cis_[ZGd] == 1
        return False
