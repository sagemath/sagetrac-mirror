# -*- coding: utf-8 -*-
r"""
Derivative of cycle index series

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
from sage.combinat.partition import Partition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.species2.cycle_index_series.operations.partitional_composite import Composite
from sage.combinat.species2.cycle_index_series.operations.add import Add
from sage.combinat.species2.cycle_index_series.operations.product import Prod
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.integer import Integer
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


class Derivative(UniqueRepresentation, Parent):
    """
    The *derivative* of `Z_F`.

    MATH::

        Z_F'(p_1, p_2, p_3, \cdots) = \left(\frac{\partial}{\partial p_1} Z_F\right)(p_1, p_2, p_3, \cdots)\,.

    Properties:

     - Distributive: `(Z_F \cdot Z_G)' = Z_F'\cdot Z_G + Z_F \cdot Z_G'`,
     - Distributive: `(Z_F + Z_G)' = Z_F' + Z_G'`,
     - Chain rule: `(Z_F(Z_G))' = Z_F'(Z_G) \cdot Z_G'`,
     - "Neutral" element: `Z_E' = Z_E`??
    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, ZF, **opts):
        # Z_E
        if ZF == CycleIndexSeries().sets():
            return ZF
        # distributive: (ZF ⋅ ZG)' --> ZF ⋅ ZG' + ZF' ⋅ ZG
        if isinstance(ZF, Prod):
            return Add(
                *[(Prod(*([(ZH, nh) if ZH != ZG else (ZH, nh-1)
                        for (ZH, nh) in ZF._dic_cis_.iteritems()] +
                          [(Derivative(ZG), Integer(1))])), ng)
                for (ZG, ng) in ZF._dic_cis_.iteritems()]
            )
        # distributive: (ZF + ZG)' --> ZF' + ZG'
        if isinstance(ZF, Add):
            return Add(*[(Derivative(ZG), ng) for (ZG, ng) in ZF._dic_cis_.iteritems()])
        # chain rule: (ZF(ZG))' --> ZF'(ZG) ⋅ ZG'
        if isinstance(ZF, Composite):
            # ZF = ZG(ZH)
            ZG, ZH = ZF._ZF_, ZF._ZG_
            return Prod((Composite(Derivative(ZG), ZH), Integer(1)), (Derivative(ZG), Integer(1)))
        # otherwise
        return super(Derivative, cls).__classcall__(cls, ZF)

    def __init__(self, ZF):
        self._ZF_ = ZF
        self._p_ = SymmetricFunctions(ZF.Frobenius_characteristic(0).base_ring()).p()

    def _repr_(self):
        return repr(self._ZF_) + "'"

    def Frobenius_characteristic(self, n):
        p = self._p_
        ch = p(self._ZF_.Frobenius_characteristic(n+1))
        return ch.map_item(lambda pi, c: (Partition(pi[:-1]), Integer(pi.to_exp()[0])*c))
        # TODO ...

