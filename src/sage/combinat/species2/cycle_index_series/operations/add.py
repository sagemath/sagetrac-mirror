# -*- coding: utf-8 -*-
"""
Sum of cycle index series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998

"""
# ******************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from collections import defaultdict
from itertools import imap
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.categories.formal_power_series import OrdinaryPowerSeries, ExponentialPowerSeries
from sage.combinat.species2.cycle_index_series import CIS
from sage.misc.cachefunc import cached_method
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer import Integer


class Add(CIS):
    """
    The sum of cycle index series.

    Properties:

     - Neutral element: `Z_F + Z_0 = Z_0 + Z_F = Z_F`,
     - Associative: `Z_F + (Z_G + Z_H) = (Z_F + Z_G) + Z_H`,
     - Commutative: `Z_F + Z_G = Z_G + Z_F`.

    EXAMPLE:

        sage: ZP = Permutations().cycle_index_series(); ZP
        ZP
        sage: ZC = SetPartitions().cycle_index_series(); ZC
        Z{Set..}
        sage: ZP + ZC   # random order
        Z{Set..} + ZP


    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZE = CycleIndexSeries().sets()
        sage: ZX = CycleIndexSeries().singletons()
        sage: ZX + ZE
        ZE + ZX
        sage: ZX + ZE == ZE + ZX
        True

        sage: ZZ = CycleIndexSeries().zero()
        sage: ZE + ZZ == ZE
        True

        sage: ZP = Permutations().cycle_index_series()
        sage: ZC = SetPartitions().cycle_index_series()
        sage: ZE + (ZP + ZC) == (ZE + ZP) + ZC
        True

        sage: TestSuite(ZP + ZC).run()

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args):
        # args = ((Z_F, nf), (Z_G, ng), ...) means nf.Z_F + ng.Z_G + ...
        dic_cis = defaultdict(Integer)
        for (ZF, nf) in args:
            # Neutral element:
            if ZF == CycleIndexSeries().zero():
                continue
            # Associativity
            elif isinstance(ZF, Add):
                for (ZG, ng) in ZF._dic_cis_.iteritems():
                    dic_cis[ZG] += ng * nf
            else:
                dic_cis[ZF] += nf

        # ### simplify (virtual species) ####
        for (ZF, nf) in list(dic_cis.items()):
            if nf == 0:
                del dic_cis[ZF]
        # ###################################

        if len(dic_cis.keys()) == 0:
            return CycleIndexSeries().zero()
        elif len(dic_cis.keys()) == 1 and dic_cis[dic_cis.keys()[0]] == 1:
            return dic_cis.keys()[0]
        else:
            return super(Add, cls).__classcall__(cls, tuple(dic_cis.items()))

    def __init__(self, dic_cis):
        CIS.__init__(self)
        self._dic_cis_ = dict(dic_cis)

    def _repr_(self):
        return " + ".join(imap(lambda (ZF, nf): (repr(nf) + "⋅" if nf != 1 else "") +
                                                repr(ZF), self._dic_cis_.iteritems()))

    def Frobenius_characteristic(self, n):
        """
        MATH::

            [n](Z_F + Z_G) = [n]Z_F + [n]Z_G\,.

        :param n: an non-negative integer

        """
        if self._valuation_() > n:
            return self._sym_.zero()
        return sum(nf * ZF.Frobenius_characteristic(n)
                   for (ZF, nf) in self._dic_cis_.iteritems())

    @cached_method
    def generating_series(self):
        from sage.combinat.species2.formal_power_series.operations.add import Add
        return Add(*map(lambda (ZF, nf): (ZF.generating_series(), nf), self._dic_cis_.iteritems()),
                   category=ExponentialPowerSeries())

    @cached_method
    def type_generating_series(self):
        from sage.combinat.species2.formal_power_series.operations.add import Add
        return Add(*map(lambda (ZF, nf): (ZF.type_generating_series(), nf), self._dic_cis_.iteritems()),
                   category=OrdinaryPowerSeries())