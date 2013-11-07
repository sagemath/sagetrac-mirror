r"""
Crystals of Gelfand-Tsetlin patterns

AUTHORS:

- Travis Scrimshaw: Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from copy import deepcopy
from sage.combinat.root_system.cartan_type import CartanType
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.gelfand_tsetlin_patterns import GelfandTsetlinPattern

class CrystalOfGelfandTsetlinPatternsElement(GelfandTsetlinPattern):
    r"""
    A Gelfand-Tsetlin pattern in a crystal.
    """
    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        Let `T` be a Gelfand-Tsetlin pattern, we have `e_i(T)` as increasing
        position `j` of the row of length `i` by `1` where `j` is the leftmost
        position such that `a_{kj} - a_{k+1,j} \leq a_{k-1,j-1} - a_{k,j-1}`
        where `k = n - i + 1`.

        EXAMPLES::

            sage: C = CrystalOfGelfandTsetlinPatterns(2, [0,1,2])
            sage: mg = C.module_generators[0]
            sage: mg.e(1)
            sage: mg.f(1).e(1) == mg
            True
            sage: C([[0,1,2],[1,1],[1]]).e(2)
            [[0, 1, 2], [1, 2], [1]]
        """
        gt = deepcopy(list(self))
        if i == 1:
            if gt[-1][0] == gt[-2][1]:
                return None
            gt[-1][0] += 1
            return self.__class__(self.parent(), gt)
        row = gt[-i]
        if row[0] > gt[-i-1][0]:
            row[0] += 1
            return self.__class__(self.parent(), gt)
        for j in range(1, len(row)):
            if row[j] - gt[-i+1][j-1] <= gt[-i-1][j] - row[j-1]:
                row[j] += 1
                return self.__class__(self.parent(), gt)
        return None

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        Let `T` be a Gelfand-Tsetlin pattern, we have `f_i(T)` as decreasing
        position `j` of the row of length `i` by `1` where `j` is the rightmost
        position such that `a_{kj} - a_{k+1,j} > \sum_{m=\ell}^{j-1} a_{k-1,m}
        - a_{k,m}` where `k = n - i + 1` and `\ell` is the first value such
        that `a_{k\ell} - a_{k+1,\ell} > 0`.

        EXAMPLES::

            sage: C = CrystalOfGelfandTsetlinPatterns(2, [0,1,2])
            sage: mg = C.module_generators[0]
            sage: mg.f(1)
            [[0, 1, 2], [1, 2], [1]]
            sage: mg.f(1).e(1) == mg
            True
            sage: mg.f_string([1,2,2])
            [[0, 1, 2], [0, 1], [1]]
        """
        gt = deepcopy(list(self))
        if i == 1:
            if gt[-1][0] == gt[-2][1]:
                return None
            gt[-1][0] -= 1
            return self.__class__(self.parent(), gt)
        row = gt[-i]
        for j in reversed(range(1, len(row))):
            if row[j] - gt[-i+1][j-1] > gt[-i-1][j] - row[j-1]:
                row[j] -= 1
                return self.__class__(self.parent(), gt)
        if row[0] > gt[-i-1][0]:
            row[0] -= 1
            return self.__class__(self.parent(), gt)
        return None

#    def epsilon(self, i):
#        r"""
#        Return the value of `\varepsilon_i` of ``self``.
#
#        EXAMPLES::
#
#            sage: C = CrystalOfGelfandTsetlinPatterns(2, [0,1,2])
#            sage: mg = C.module_generators[0]
#            sage: mg.epsilon(1)
#            0
#            sage: mg.epsilon(2)
#            0
#            sage: mg.f(1).epsilon(1)
#            1
#            sage: mg.f(1).epsilon(2)
#            0
#        """
#        if i == 1:
#            return self[-2][1] - self[-1][0]
#        row = self[-i]
#        diff = self[-i-1][-1] - row[-1]
#        for j in range(len(row)-1):
#            diff += min(self[-i-1][j+1], self[-i+1][j]) - row[j]
#        return diff
#
#    def phi(self, i):
#        r"""
#        Return the value of `\varphi_i` of ``self``.
#
#        EXAMPLES::
#
#            sage: C = CrystalOfGelfandTsetlinPatterns(2, [0,1,2])
#            sage: mg = C.module_generators[0]
#            sage: mg.phi(1)
#            1
#            sage: mg.phi(2)
#            1
#            sage: mg.f(1).phi(1)
#            0
#            sage: mg.f(1).phi(2)
#            2
#        """
#        if i == 1:
#            return self[-1][0] - self[-2][0]
#        row = self[-i]
#        diff = row[0] - self[-i-1][0]
#        for j in range(1, len(row)):
#            diff += row[j] - max(self[-i-1][j], self[-i+1][j-1])
#        return diff

class CrystalOfGelfandTsetlinPatterns(Parent, UniqueRepresentation):
    r"""
    The crystal of Gelfand-Tsetlin patterns which is equivalent to the
    tableaux model under the bijection.

    INPUT:

    - ``n`` -- type `A_n^{(1)}`
    - ``top_row`` -- The top row which is equivalent to the dominant weight

    EXAMPLES::
    """
    @staticmethod
    def __classcall_private__(cls, n, top_row):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: C = CrystalOfGelfandTsetlinPatterns(3, [1,1,3])
            sage: C2 = CrystalOfGelfandTsetlinPatterns(int(3), (0,1,1,3))
            sage: C is C2
            True
        """
        top_row = list(top_row)
        if len(top_row) > n+1:
            raise ValueError("The top row can have at most %s entries"%(n+1))
        for i in range(len(top_row)-1):
            if top_row[i] > top_row[i+1]:
                raise ValueError("The top row must be non-decreasing")
        while len(top_row) <= n:
            top_row.insert(0, 0)
        return super(CrystalOfGelfandTsetlinPatterns, cls).__classcall__(cls,n,tuple(top_row))

    def __init__(self, n, top_row):
        r"""
        EXAMPLES::

            sage: C = CrystalOfGelfandTsetlinPatterns(2, [1,1,3])
            sage: TestSuite(C).run()
        """
        self._cartan_type = CartanType(['A',n])
        self._top_row = top_row
        Parent.__init__(self, category=ClassicalCrystals())
        gt = [list(top_row[i:]) for i in range(len(top_row))]
        self.module_generators = (self.element_class(self, gt),)

    def _element_constructor_(self, data):
        r"""
        Construct an element of ``self`` from ``data``.

        INPUT:

        - ``data`` -- A Gelfand-Tsetlin pattern

        EXAMPLES::

            sage: C = CrystalOfGelfandTsetlinPatterns(2, [1, 2])
            sage: gt = C([[0,1,2],[1,1],[1]]); gt
            [[0, 1, 2], [1, 1], [1]]
            sage: gt.parent()
            Crystal of Gelfand-Tsetlin patterns of type ['A', 2] with top row (0, 1, 2)
        """
        return self.element_class(self, data)

    Element = CrystalOfGelfandTsetlinPatternsElement

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: CrystalOfGelfandTsetlinPatterns(2, [1, 2])
            Crystal of Gelfand-Tsetlin patterns of type ['A', 2] with top row (0, 1, 2)
        """
        return "Crystal of Gelfand-Tsetlin patterns of type %s with top row %s"%(self._cartan_type, self._top_row)

