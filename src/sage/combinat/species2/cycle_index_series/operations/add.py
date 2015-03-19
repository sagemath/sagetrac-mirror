# -*- coding: utf-8 -*-
"""
Sum of cycle index series

Reference
---------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux
  Cambridge University Press, 1998
"""
#******************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from collections import defaultdict
from itertools import imap
from sage.categories.cycle_index_series import CycleIndexSeries
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class Add(UniqueRepresentation, Parent):
    """
    The sum of cycle index series.

    Properties:

     - Neutral element: `Z_F + Z_0 = Z_0 + Z_F = Z_F`,
     - Associative: `Z_F + (Z_G + Z_H) = (Z_F + Z_G) + Z_H`,
     - Commutative: `Z_F + Z_G = Z_G + Z_F`.


    TESTS::

        sage: from sage.categories.cycle_index_series import CycleIndexSeries
        sage: ZE = CycleIndexSeries().sets()
        sage: ZX = CycleIndexSeries().singleton()
        sage: ZX + ZE
        Z_E + Z_X
        sage: ZX + ZE == ZE + ZX
        True

        sage: ZZ =CycleIndexSeries().zero()
        sage: ZE + ZZ == ZE
    """
    # TODO test

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
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

        #### simplify (virtual species) ####
        for (ZF, nf) in list(dic_cis.items()):
            if nf == 0:
                del dic_cis[ZF]
        ####################################

        if len(dic_cis.keys()) == 0:
            return CycleIndexSeries.zero()
        elif len(dic_cis.keys()) == 1 and dic_cis[dic_cis.keys()[0]] == 1:
            return dic_cis.keys()[0]
        else:
            return super(Add, cls).__classcall__(cls, tuple(dic_cis.items()))

    def __init__(self, dic_cis):
        Parent.__init__(self, category=CycleIndexSeries())
        self._dic_cis_ = dict(dic_cis)

    def _repr_(self):
        return " + ".join(imap(lambda (ZF, nf): (repr(nf) + "⋅" if nf != 1 else "") + repr(ZF), self._dic_cis_.iteritems()))

    def Frobenius_characteristic(self, n):
        """
        MATH::

            [n](Z_F + Z_G) = [n]Z_F + [n]Z_G\,.

        :param n: an non-negative integer

        """
        return sum(nf * ZF.Frobenius_characteristic(n)
                   for (ZF, nf) in self._dic_cis_.iteritems())

