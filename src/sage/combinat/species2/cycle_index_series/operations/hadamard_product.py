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
from sage.combinat.partition import Partitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.misc import exists
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class HadamardProduct(UniqueRepresentation, Parent):
    """
    The *Hadamard product* of cycles index series

    MATH::

        (Z_F \times Z_G)(p_1, p_2, p_3, \cdots) = Z_F(p_1, p_2, p_3, \cdots) \times Z_G(p_1, p_2, p_3, \cdots)\,.

    Properties:

     - Commutative: `Z_F \times Z_G = Z_G \times Z_F`,
     - Absorbing element: ` Z_0 \times Z_F = Z_F \times Z_0 = Z_0`,
     - Neutral element: `Z_E \times Z_F = Z_F \times Z_E = Z_F`,
     - Distributive: `(Z_F + Z_G) \times Z_H = Z_F \times Z_H + Z_G \times Z_H`.

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):

        # commutativity
        dic_cis = defaultdict(Integer)

        for i, (ZF, nf) in enumerate(args):
            # neutral element
            if ZF == CycleIndexSeries().sets():
                continue
            # absorbing element
            elif ZF == CycleIndexSeries().zero():
                if nf != 0: return CycleIndexSeries().zero()
                # else continue
            # distributive
            elif isinstance(ZF, Add):
                return Add(*tuple(HadamardProduct(*(tuple(dic_cis.items()) + ((ZG, ng),) + args[i+1:]))
                                  for (ZG, ng) in ZF._dic_cis_.iteritems()))
            # associative
            elif isinstance(ZF, HadamardProduct):
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
            return CycleIndexSeries().sets()
        elif len(dic_cis.keys()) == 1 and dic_cis.values()[0] == 1:
            return dic_cis.keys()[0]
        else:
            return super(HadamardProduct, cls).__classcall__(cls, tuple(dic_cis.items()))

    def __init__(self, dic_cis):
        Parent.__init__(self, category=CycleIndexSeries())
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

    def _repr_(self):
        return " x ".join(imap(lambda (ZF, nf): repr(ZF) + ("^{x%d}"%nf if nf != 1 else ""),
                               self._dic_cis_.iteritems()))