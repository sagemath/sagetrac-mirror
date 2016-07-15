"""
Shifted Tableaux

AUTHORS:

- Zachary Hamaker, Tobias Johnson (2016): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2016 Zachary Hamaker <your email>
#                          Tobias Johnson <your email>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.functions.other import factorial
from sage.misc.misc_c import prod
from sage.misc.prandom import randrange
from sage.structure.list_clone import ClonableArray
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.permutation import Permutation
from sage.combinat.posets.posets import Poset
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

class ShiftedTableau(ClonableArray):
    """
    A shifted tableau.

    EXAMPLES::

        sage: T = ShiftedTableaux([4,2])
        sage: T([[1,2,4,5],[3,6]])[1]
        [3, 6]
        sage: len(T([[1,2,4,5],[3,6]]))
        2
        sage: T([[1,2,4,5],[3,6]])
        [[1, 2, 4, 5], [3, 6]]
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, t):
        """
        This ensures that a shifted tableau is only ever constructed
        as an ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: data = [[1,2,4,5],[3,6]]
            sage: t = ShiftedTableau(data)
            sage: T = ShiftedTableaux([4,2])
            sage: t == T(data)
            True
        """
        if isinstance(t, cls):
            return t
        shape = [len(row) for row in t]
        return ShiftedTableaux(shape)(t)

    def check(self):
        """
        Check that ``self`` is a valid shifted tableaux.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: t = T([[1,2,4,5],[3,6]])
            sage: t.check()
        """
        for i,row in enumerate(self):
            if i > 0:
                if not all(val > self[i-1][j+1] for j,val in enumerate(row)):
                    raise ValueError("non-increasing column")
            if not all(row[j] < row[j+1] for j in range(len(row)-1)):
                raise ValueError("non-increasing row")

    def _latex_(self):
        r"""
        Return LaTex code for ``self``.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: latex(T([[1,2,4,5],[3,6]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{2}&\lr{4}&\lr{5}\\\cline{1-4}
            &\lr{3}&\lr{6}\\\cline{2-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        L = [[None]*i + row for (i,row) in enumerate(self)]
        return tex_from_array(L)

    def pp(self):
        """
        Print out a nice version of ``self``.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: T([[1,2,4,5],[3,6]]).pp()
            1  2  4  5
               3  6
        """
        for row_num, row in enumerate(self):
            print "   "*row_num + "".join("{!r:3}".format(n) for n in row)

    # Perhaps call this number_of_yang_baxter_configurations (or positions)?
    def yang_baxter_moves(self):
        """
        Count the number of Yang-Baxter-type configurations in ``self``.

        A Yang-Baxter-type configuration is a subtableau of the form::

             i  i+1
                i+2

        along the superdiagonal.

        EXAMPLES::

            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.yang_baxter_moves()
            1
        """
        return sum(1 if pos else 0 for pos in self.yang_baxter_positions())

    def yang_baxter_positions(self):
        """
        Return a vector of ``True`` and ``False`` values giving
        the locations of Yang-Baxter-type configurations in ``self``.

        .. SEEALSO::

            :meth:`yang_baxter_moves`

        EXAMPLES::

            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.yang_baxter_positions()
            [False, True, False]
        """
        return [self[r+1][0] == self[r][0] + 2
                for r in range(len(self) - 1)]

    def position(self, value):
        """
        Return the coordinates of the cell of ``self`` labelled value.

        EXAMPLES::

            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.position(6)
            (1, 1)
            sage: t.position(12)
            (3, 0)
            sage: t.position(22)
            Traceback (most recent call last):
            ...
            ValueError: 22 is not in shifted tableau
        """
        for i in range(len(self)):
            if value in self[i]:
                return (i,self[i].index(value))
        raise ValueError("{} is not in shifted tableau".format(value))

    def shape(self):
        """
        Return the shape of ``self``.

        EXAMPLES::

            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.shape()
            (6, 4, 3, 1)
        """
        return self.parent()._shape

class ShiftedTableaux(UniqueRepresentation, Parent):
    """
    Set of shifted tableaux of a fixed shape.
    """
    @staticmethod
    def __classcall_private__(cls, shape):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: T1 = ShiftedTableaux([5,2,1])
            sage: T2 = ShiftedTableaux((5,2,1))
            sage: T1 is T2
            True
        """
        shape = tuple(shape)
        if not all(shape[i] > shape[i+1] > 0 for i in range(len(shape)-1)):
            raise ValueError("invalid shape")
        return super(ShiftedTableaux, cls).__classcall__(cls, shape)

    def __init__(self, shape):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: T = ShiftedTableaux([5,2,1])
            sage: TestSuite(T).run()
        """
        self._shape = shape
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T[0:4]
            [[[1, 2, 3, 4, 5, 6], [7, 8, 9, 10], [11, 12, 13], [14]],
             [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10], [11, 12, 14], [13]],
             [[1, 2, 3, 4, 5, 6], [7, 8, 9, 11], [10, 12, 13], [14]],
             [[1, 2, 3, 4, 5, 6], [7, 8, 9, 11], [10, 12, 14], [13]]]
        """
        P = shape_to_poset(self._shape)
        L = P.linear_extensions()
        for l in L:
            t = linear_extension_to_tableau(l,self._shape)
            yield self.element_class(self,t)

    def cardinality(self):
        """
        Return the number of shifted tableaux with the shape of ``shape``.

        This uses the shifted hook-length formula originally due to Thrall in 1952.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: T.cardinality()
            5
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T.cardinality()
            1716
        """
        n = sum(self._shape)
        denom = prod(len(hook(self._shape, (i, j))) + 1
                     for i, row_len in enumerate(self._shape)
                     for j in range(row_len))
        return factorial(n) // denom

    def random_element(self):
        """
        Return a random shifted tableau with the shape of ``self``.

        This uses a variant of the Greene-Nijenhuis-Wilf hook-walk
        algorithm in [Sagan80]_.

        REFERENCES:

        .. [Sagan80] Sagan 1980.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: T.random_element()  # random
            [[1, 2, 4, 5], [3, 6]]
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T.random_element()  # random
            [[1, 2, 3, 5, 7, 10], [4, 6, 9, 12], [8, 11, 14], [13]]
        """
        tableau = [ [0] * row_length for row_length in self._shape ]
        N = sum(self._shape)
        current_shape = list(self._shape)
        for next_number in reversed(range(1, N+1)):
            cell = random_cell(current_shape)
            hook_coords = hook(current_shape, cell)
            while hook_coords:
                cell = hook_coords[ randrange(len(hook_coords)) ]
                hook_coords = hook(current_shape, cell)
            # done with hook walk, so place next_number in tableau
            tableau[cell[0]][cell[1]] = next_number
            current_shape[cell[0]] -= 1

        return self.element_class(self,tableau)

    Element = ShiftedTableau


#####################################################################
## Helper functions

def hook(shape, pos):
    """
    Determines the coordinates of the cells in the hook corresponding to
    the cell indexed by the tuple ``pos`` in the shifted shape ``shape``.

    EXAMPLES::

        sage: from sage.combinat.shifted_tableaux import hook
        sage: hook([6,4,3,1], (1,1))
        [(1, 2), (1, 3), (2, 0), (3, 0)]
    """
    # arm first:
    hook_coordinates = [(pos[0], col) for col in range(pos[1]+1, shape[pos[0]])]

    # now the leg:
    (r, c) = pos
    r, c = r+1, c-1
    nrows = len(shape)
    while c >= 0 and r < nrows and c < shape[r]:
        hook_coordinates.append( (r, c) )
        r, c = r+1, c-1
    if c < 0 and r < len(shape):
        hook_coordinates.extend( (r, j) for j in range(shape[r]) )
    return hook_coordinates

def random_cell(shape):
    """
    Return a uniformly random cell from the given shifted tableau shape.

    EXAMPLES:

        sage: from sage.combinat.shifted_tableaux import random_cell
        sage: [random_cell([6,4,3,1]) for i in range(6)]  # random
        [(2, 2), (1, 1), (1, 2), (1, 1), (2, 1), (3, 0)]
    """
    rnd = randrange(0, sum(shape))
    row = 0
    cells_so_far = 0
    while cells_so_far + shape[row] <= rnd:
        cells_so_far += shape[row]
        row += 1
    return (row, rnd - cells_so_far)

def shape_to_poset(shape):
    r"""
    Convert the shifted shape ``shape`` into a poset with elements `1,2,...,n`.

    EXAMPLES::

        sage: from sage.combinat.shifted_tableaux import shape_to_poset
        sage: shape_to_poset([4,2]).cover_relations()
        [[1, 2], [2, 3], [2, 5], [3, 4], [3, 6], [5, 6]]
    """
    elts = range(1, sum(shape)+1)
    rels = [(i, i+1) for i in range(1, shape[0])]
    tot = shape[0]
    prev = shape[0]
    for row_len in shape[1:]:
        rels += [(i, i+1) for i in range(tot+1, tot+row_len)]
        rels += [(i-prev+1, i) for i in range(tot+1, tot+1+row_len)]
        tot += row_len
        prev = row_len
    return Poset((elts, rels))

def linear_extension_to_tableau(linear_extension, shape):
    """
    Convert a linear extension for the poset ``shape_to_poset(shape)``
    into a shifted tableau of shape ``shape``.

    EXAMPLES::

        sage: from sage.combinat.shifted_tableaux import shape_to_poset
        sage: from sage.combinat.shifted_tableaux import linear_extension_to_tableau
        sage: P = shape_to_poset([4,2])
        sage: [linear_extension_to_tableau(L,[4,2]) for L in P.linear_extensions()]
        [[[1, 2, 3, 4], [5, 6]],
         [[1, 2, 3, 5], [4, 6]],
         [[1, 2, 3, 6], [4, 5]],
         [[1, 2, 4, 5], [3, 6]],
         [[1, 2, 4, 6], [3, 5]]]
    """
    L = list(Permutation(linear_extension).inverse()) # tableau in a list, rows not separate
    A = [sum(shape[:i]) for i in range(len(shape)+1)] # determine indices where new rows start
    return [L[A[i]:A[i+1]] for i in range(len(A)-1)]

