"""
Benkart-Kang-Kashiwara crystals for the general-linear Lie superalgebra
"""

#*****************************************************************************
#       Copyright (C) 2017 Franco Saliola <saliola@gmail.com>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.regular_supercrystals import RegularSuperCrystals
from sage.combinat.partition import _Partitions
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.letters import CrystalOfBKKLetters
from sage.combinat.crystals.tensor_product import CrystalOfWords
from sage.combinat.crystals.tensor_product_element import CrystalOfBKKTableauxElement


class CrystalOfBKKTableaux(CrystalOfWords):
    """
    Crystal of tableaux for type `A(m,n)`.

    EXAMPLES::

        sage: from sage.combinat.crystals.bkk_crystals import CrystalOfBKKTableaux
        sage: T = CrystalOfBKKTableaux(['A', [1,1]], [2,1])
    """
    @staticmethod
    def __classcall_private__(cls, ct, shape):
        ct = CartanType(ct)
        shape = _Partitions(shape)
        if len(shape) > ct.m + 1 and shape[ct.m] > ct.n + 1:
            raise ValueError("invalid hook shape")
        return super(CrystalOfBKKTableaux, cls).__classcall__(cls, ct, shape)

    def __init__(self, ct, shape):
        r"""
        EXAMPLES::

            sage: from sage.combinat.crystals.bkk_crystals import CrystalOfBKKTableaux
            sage: T = CrystalOfBKKTableaux(['A', [1,1]], [2,1])
            sage: T
            Crystal of BKK tableaux of skew shape [2, 1] of gl(2|2)

            sage: TestSuite(T).run()
        """
        self._shape = shape
        self._cartan_type = ct
        m = ct.m + 1
        n = ct.n + 1
        C = CrystalOfBKKLetters(ct)
        tr = shape.conjugate()
        mg = []
        for i,col_len in enumerate(tr):
            for j in range(col_len - m):
                mg.append(C(i+1))
            for j in range(max(0, m - col_len), m):
                mg.append(C(-j-1))
        mg = list(reversed(mg))
        Parent.__init__(self, category=RegularSuperCrystals())
        self.module_generators = (self.element_class(self, mg),)

    def _repr_(self):
        m = self._cartan_type.m + 1
        n = self._cartan_type.n + 1
        return "Crystal of BKK tableaux of skew shape {} of gl({}|{})".format(self.shape(), m, n)

    def shape(self):
        r"""
        EXAMPLES::

            sage: from bkk_tableaux import BKKTableaux
            sage: B = BKKTableaux([1,1,1], 2, 2)
            sage: B.shape()
            [1, 1, 1] / []
            sage: B = BKKTableaux(([3,1,1], [2,1]), 2, 2)
            sage: B.shape()
            [3, 1, 1] / [2, 1]
        """
        return self._shape

    class Element(CrystalOfBKKTableauxElement):
        pass

