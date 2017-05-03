"""
Shifted Tableaux

Shifted tableaux arise in the projective representation theory of the symmetric
groups and other related contexts [HH92]_. The shifted tableaux are indexed by
**strict partitions**, or partitions with distinct parts. Unlike 
:class:`~sage.combinat.tableau.Tableau`, which are drawn as left justified arrays, shifted tableaux are drawn
with row `r+1` shifted one position to the right of row `r`. For example, the
standard shifted tableaux of shape `(4,1)` are::

    1 2 3 4    1 2 3 5    1 2 4 5
      5          4          3   

EXAMPLES::

    sage: t = ShiftedTableau([[1,2,3],[4,5]]); t
    [[1, 2, 3], [4, 5]]
    sage: t.pp()
    1  2  3
       4  5
    sage: t(0,0)
    1
    sage: t(1,0)
    4
    sage: t.shape()
    [3, 2]
    sage: t.size()
    5
    sage: t.entries()
    (1, 2, 3, 4, 5)
    sage: t.parent()
    Shifted tableaux of shape [3, 2]
    sage: t.category()
    Category of elements of Shifted tableaux of shape [3, 2]
    sage: ShiftedTableaux([4,1])
    Shifted tableaux of shape [4, 1]
    sage: ShiftedTableaux([4,1])[:]
    [[[1, 2, 3, 4], [5]], [[1, 2, 3, 5], [4]], [[1, 2, 4, 5], [3]]]

AUTHORS:

- Zachary Hamaker, Tobias Johnson (2016): Initial version

- Travis Scrimshaw (2016): Put into category framework

- Andrew Mathas (2016): Added extended functionality and parent classes

Element classes:

* :class:`ShiftedTableau`

Factory classes:

* :class:`ShiftedTableaux`

Parent classes:

* :class:`ShiftedTableaux_all`
* :class:`ShiftedTableaux_size`
* :class:`ShiftedTableaux_shape`
"""

#*****************************************************************************
#       Copyright (C) 2016 Zachary Hamaker <zachary.hamaker at gmail.com>
#                          Tobias Johnson <tobias.johnson at nyu.edu>
#                          Travis Scrimshaw <tcscrims at gmail.com>
#                          Andrew Mathas <andrew.mathas atsydney.edu.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import

from .partition import Partition, StrictPartitions
from .permutation import Permutation
from .posets.posets import Poset
from .tableau import Tableaux

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.functions.other import factorial
from sage.misc.cachefunc import cached_method
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.misc_c import prod
from sage.misc.prandom import randrange
from sage.rings.integer import Integer
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class ShiftedTableau(ClonableArray):
    """
    A shifted tableau.

    EXAMPLES::

        sage: T = ShiftedTableaux([4,2])
        sage: T([[1,2,4,5],[3,6]])[1]
        (3, 6)
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

    def __init__(self, parent, t):
        r"""
        Initialize a shifted tableau.

        TESTS::

            sage: s = ShiftedTableau([[1,2],[3]])
            sage: t = ShiftedTableaux([2,1])([[1,2],[3]])
            sage: s==t
            True
            sage: t.parent()
            Shifted tableaux of shape [2, 1]
            sage: s.parent()
            Shifted tableaux of shape [2, 1]
            sage: r = ShiftedTableaux([2,1])(s); r.parent()
            Shifted tableaux of shape [2, 1]
            sage: s is t # identical shifted tableaux are distinct objects
            False

        A shifted tableau is shallowly immutable. The entries
        themselves may be mutable objects, though in that case the
        resulting ShiftedTableau should be unhashable.

            sage: T = ShiftedTableau([[1,2],[3]])
            sage: t0 = T[0]
            sage: t0[1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
            sage: T[0][1] = 5
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        if isinstance(t, ShiftedTableau):
            # Since we are (supposed to be) immutable, we can share the underlying data
            ClonableArray.__init__(self, parent, t)
            return

        # Normalize t to be a list of tuples.
        t = [tuple(_) for _ in t]

        ClonableArray.__init__(self, parent, t)
        # This dispatches the input verification to the :meth:`check`
        # method.

    def check(self):
        """
        Check that ``self`` is a valid shifted tableaux.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: t = T([[1,2,4,5],[3,6]])
            sage: t.check()
        """
        if [len(_) for _ in self] not in StrictPartitions():
            raise ValueError('shape must be a strict partition')
        entries = sorted(sum((list(_) for _ in self), []))
        if entries != range(1, len(entries)+1):
            raise ValueError('must contain the numbers 1,2,...,<size>')
        for i,row in enumerate(self):
            if i > 0:
                if not all(val > self[i-1][j+1] for j,val in enumerate(row)):
                    raise ValueError('non-increasing column')
            if not all(row[j] < row[j+1] for j in range(len(row)-1)):
                raise ValueError('non-increasing row')

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: t = ShiftedTableau([[1,2,3],[4,5]])
            sage: ShiftedTableaux.options.display="list"
            sage: t
            [[1, 2, 3], [4, 5]]
            sage: ShiftedTableaux.options.display="array"
            sage: t
              1  2  3
              4  5
            sage: ShiftedTableaux.options.display="compact"; t
            1,2,3/4,5
            sage: ShiftedTableaux.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        """
        Return a string representation of ``self`` as a list.

        EXAMPLES::

            sage: T = ShiftedTableau([[1,2,3],[4,5]])
            sage: T._repr_list()
            '[[1, 2, 3], [4, 5]]'
        """
        return repr([list(_) for _ in self])

    def _repr_compact(self):
        """
        Return a compact string representation of ``self``.

        EXAMPLES::

            sage: ShiftedTableau([[1,2,3],[4,5]])._repr_compact()
            '1,2,3/4,5'
            sage: ShiftedTableau([])._repr_compact()
            '-'
        """
        if not self:
            return '-'
        return '/'.join(','.join('%s'%r for r in row) for row in self)

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: t = ShiftedTableau([[1,2,3,6],[4,5]])
            sage: print(t._repr_diagram())
              1  2  3  6
                 4  5
            sage: ShiftedTableaux.options(convention="french")
            sage: print(t._repr_diagram())
                 4  5
              1  2  3  6
            sage: ShiftedTableaux.options._reset()

        TESTS:

            sage: ShiftedTableau([[1,2,3,4,5,6,7,8,9,10,11],[12,13,14,15,16],[17,18,19]]).pp()
            1  2  3  4  5  6  7  8  9 10 11
              12 13 14 15 16
                 17 18 19
        """
        if not self:
            return "  -"

        # Get the widths of the columns
        str_tab = [[str(data) for data in row] for row in self]
        col_widths = [2]*len(str_tab[0])
        for r, row in enumerate(str_tab):
            for i,e in enumerate(row):
                col_widths[r+i] = max(col_widths[r+i], len(e))

        if self.parent().options.convention == "French":
            str_tab.reverse()
            return '\n'.join('{}{:>{}}'.format(' '*(len(str_tab)-row_index), '', sum(col_widths[row_index:]))
                             + ' '.join('{entry:>{width}}'.format(entry=e,width=col_widths[i])
                                        for i,e in enumerate(row))
                             for row_index, row in enumerate(str_tab))
        else:
            return '\n'.join('{}{:>{}}'.format(' '*(row_index+1), '', sum(col_widths[:row_index]))
                             + ' '.join('{entry:>{width}}'.format(entry=e,width=col_widths[i])
                                        for i,e in enumerate(row))
                             for row_index, row in enumerate(str_tab))

    def pp(self):
        """
        Print out a nice version of ``self``.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: T([[1,2,4,5],[3,6]]).pp()
            1  2  4  5
               3  6
        """
        print(self._repr_diagram())

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
            sage: T.parent().convention='French'
            sage: latex(T([[1,2,4,5],[3,6]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{2}&\lr{4}&\lr{5}\\\cline{1-4}
            &\lr{3}&\lr{6}\\\cline{2-3}
            \end{array}$}
            }
        """
        return self.parent().options._dispatch(self, '_latex_', 'latex')

    _latex_list = _repr_list

    def _latex_diagram(self):
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
        L = [[None]*i + list(row) for (i,row) in enumerate(self)]
        return tex_from_array(L)

    def __call__(self, *cell):
        r"""
        Function call of ``self``.

        INPUT:

        - ``cell`` -- a pair of integers, tuple, or list specifying a cell in
          the tableau

        OUTPUT:

        - the value in the corresponding cell

        EXAMPLES::

            sage: t = ShiftedTableau([[1,2,3],[4,5]])
            sage: t(1,0)
            4
            sage: t((1,0))
            4
            sage: t(3,3)
            Traceback (most recent call last):
            ...
            IndexError: the cell (3,3) is not contained in [[1, 2, 3], [4, 5]]
        """
        try:
            i,j = cell
        except ValueError:
            i,j = cell[0]

        try:
            return self[i][j]
        except IndexError:
            raise IndexError("the cell (%d,%d) is not contained in %s"%(i, j, repr(self)))

    def __eq__(self, other):
        r"""
        Check whether ``self`` is equal to ``other``.

        .. TODO::

            This overwrites the equality check of
            :class:`~sage.structure.list_clone.ClonableArray` in order to
            circumvent the coercion framework.  Eventually this should be
            solved more elegantly, for example along the lines of what was
            done for `k`-tableaux.

            For now, two elements are equal if their underlying defining
            lists compare equal.

        INPUT:

        ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        Boolean.

        TESTS::

            sage: t = ShiftedTableau([[1,2]])
            sage: t == 0
            False
            sage: t == ShiftedTableaux([2])([[1,2]])
            True
        """
        if isinstance(other, ShiftedTableau):
            return list(self) == list(other)
        else:
            return list(self) == other

    def __ne__(self, other):
        r"""
        Check whether ``self`` is unequal to ``other``.

        See the documentation of :meth:`__eq__`.

        INPUT:

        ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        Boolean.

        TESTS::

            sage: ShiftedTableau([[1,2],[3]]) !=[]
            True
        """
        if isinstance(other, ShiftedTableau):
            return list(self) != list(other)
        else:
            return list(self) != other

    def shape(self):
        """
        Return the shape of the shifted tableau ``self``.

        EXAMPLES::

            sage: ShiftedTableau([[1, 2, 3, 5, 8], [4, 6], [7]]).shape()
            [5, 2, 1]
        """
        return self.parent().shape()

    def size(self):
        """
        Return the size of the shape of the shifted tableau ``self``.

        EXAMPLES::

            sage: ShiftedTableau([[1, 2, 4, 7], [3, 5], [6]]).size()
            7
        """
        return self.parent().size()

    def entries(self):
        """
        Return the tuple of all entries of ``self``, in the order obtained
        by reading across the rows from top to bottom (in English
        notation).

        EXAMPLES::

            sage: ShiftedTableau([[1, 2, 4, 6], [3, 5],[7]]).entries()
            (1, 2, 4, 6, 3, 5, 7)
        """
        return sum(self, ())

    def entry(self, cell):
        """
        Returns the entry of cell ``cell`` in the tableau ``self``.

        Here, ``cell`` should be given as a tuple `(i,j)` of zero-based
        coordinates (so the northwesternmost cell in English notation
        is `(0,0)`).

        EXAMPLES::

            sage: ShiftedTableau([[1,2,4,5],[3,6]]).entry( (0,0) )
            1
            sage: ShiftedTableau([[1,2,4,6],[3,5]]).entry( (0,1) )
            2
        """
        i,j = cell
        return self[i][j]

    # Perhaps call this number_of_yang_baxter_configurations (or positions)?
    def yang_baxter_moves(self):
        """
        Count the number of Yang-Baxter-type configurations in ``self``.

        A Yang-Baxter-type configuration is a subtableau of the form::

             i  i+1
                i+2

        along the superdiagonal.

        EXAMPLES::

            sage: t = ShiftedTableau([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
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
        Return the coordinates of the cell of ``self`` labeled value.

        EXAMPLES::

            sage: t = ShiftedTableau([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
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


class ShiftedTableaux(UniqueRepresentation, Parent):
    r"""
    A factory for the various classes of shifted standard tableaux.

    INPUT:

    - either a non-negative integer (possibly specified with the keyword
      ``n``) or a partition

    OUTPUT:

    - with no argument, the class of all standard tableaux

    - with a non-negative integer argument, ``n``, the class of all standard
      tableaux of size ``n``

    - with a partition argument, the class of all standard tableaux of that
      shape

    A standard tableau is a semistandard tableaux which contains each of the
    entries from 1 to ``n`` exactly once.

    All classes of standard tableaux are iterable.

    EXAMPLES::

        sage: ST = ShiftedTableaux(4); ST
        Shifted tableaux of size 4
        sage: ST.first()
        [[1, 2, 3, 4]]
        sage: ST.last()
        [[1, 2, 4], [3]]
        sage: ST.cardinality()
        3
        sage: ST.list()
        [[[1, 2, 3, 4]], [[1, 2, 3], [4]], [[1, 2, 4], [3]]]

    .. SEEALSO::

        - :class:`ShiftedTableau`
        - :class:`StandardTableau`

    TESTS::

        sage: ShiftedTableaux()([])
        []
        sage: ST = ShiftedTableaux([3,2]); ST
        Shifted tableaux of shape [3, 2]
        sage: ST.first()
        [[1, 2, 3], [4, 5]]
        sage: ST.last()
        [[1, 2, 4], [3, 5]]
        sage: ST.cardinality()
        2
        sage: ST.list()
        [[[1, 2, 3], [4, 5]], [[1, 2, 4], [3, 5]]]
    """
    # use Tableaux options
    options = Tableaux.options

    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`ShiftedTableaux` for
        more information.

        TESTS::

            sage: ShiftedTableaux()
            Shifted tableaux
            sage: ShiftedTableaux(3)
            Shifted tableaux of size 3
            sage: ShiftedTableaux([2,1])
            Shifted tableaux of shape [2, 1]
            sage: ShiftedTableaux(0)
            Shifted tableaux of size 0

            sage: ShiftedTableaux(-1)
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
            sage: ShiftedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: the argument must be a non-negative integer or a partition
        """
        if args:
            size = args[0]
        elif 'size' in kwargs:
            size = kwargs[n]
        else:
            size = None

        if size is None:
            return ShiftedTableaux_all()

        elif size in StrictPartitions():
            return ShiftedTableaux_shape(Partition(size))

        if not isinstance(size, (int, Integer)) or size < 0:
            raise ValueError("the argument must be a non-negative integer or a partition")

        return ShiftedTableaux_size(size)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: [[1,1],[2,3]] in ShiftedTableaux()
            False
            sage: [[1,2],[3,4]] in ShiftedTableaux()
            False
            sage: [[1,3],[2]]   in ShiftedTableaux()
            False
            sage: [[1,5],[2]]   in ShiftedTableaux()
            False
            sage: [] in ShiftedTableaux()
            True

        Check that :trac:`14145` is fixed::

            sage: 1 in ShiftedTableaux()
            False
        """
        if isinstance(x, ShiftedTableau) or x == []:
            return True
        elif not isinstance(x, list):
            return False

        # have a list, so check if it is shifted standard
        entries = sorted(sum((list(_) for _ in x), []))
        if entries != range(1, len(entries)+1):
            return False  # must contain 1,2...,n

        if [len(_) for _ in x] not in StrictPartitions():
            return False # must have strict partition shape

        for i,row in enumerate(x):
            if i > 0:
                if not all(val > x[i-1][j+1] for j,val in enumerate(row)):
                    return False  # increasing down columns
            if not all(row[j] < row[j+1] for j in range(len(row)-1)):
                return False  # increasing along rows

        return True

    _is_a = __contains__

    def an_element(self):
        r"""
        Return a particular shifted tableaux in the class.

        TESTS::

            sage: ShiftedTableaux().an_element()
            []
            sage: ShiftedTableaux(4).an_element()
            [[1, 2, 3, 4]]
            sage: ShiftedTableaux([3,1]).an_element()
            [[1, 2, 3], [4]]
        """
        return self[0]

class ShiftedTableaux_all(ShiftedTableaux, DisjointUnionEnumeratedSets):
    """
    The class of all shifted tableaux.
    """
    def __init__(self):
        r"""
        Initializes the class of all shifted tableaux.

        TESTS::

            sage: TestSuite( ShiftedTableaux() ).run()

        Check that containment is correctly inherited::

            sage: [[1,2],[3]] in ShiftedTableaux()
            True
            sage: [[1,3],[2]] in ShiftedTableaux()
            False
            sage: [[1,1],[2]] in ShiftedTableaux()
            False
            sage: [[1,2],[3,4]] in ShiftedTableaux()
            False
        """
        DisjointUnionEnumeratedSets.__init__(
                self,
                Family(NonNegativeIntegers(), ShiftedTableaux_size),
                facade=True, keepkey=False
        )

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedTableaux()    # indirect doctest
            Shifted tableaux
        """
        return "Shifted tableaux"

    def __getitem__(self, r):
        r"""
        Return the ``r``th item of ``self``.

        The default implementation of ``__getitem__`` for infinite
        enumerated sets does not allow for finite slices so we override it.

        EXAMPLES::

            sage: ShiftedTableaux()[10]
            [[1, 2, 3, 5], [4]]
            sage: ShiftedTableaux()[10:16]  
            [[[1, 2, 3, 5], [4]],
             [[1, 2, 4, 5], [3]],
             [[1, 2, 3], [4, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 2, 3, 4, 5, 6]],
             [[1, 2, 3, 4, 5], [6]]]
            sage: ShiftedTableaux()[:5]
            [[], [[1]], [[1, 2]], [[1, 2, 3]], [[1, 2], [3]]]

        TESTS::

            sage: ShiftedTableaux()[5:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set

            sage: ShiftedTableaux()[:]
            Traceback (most recent call last):
            ...
            ValueError: infinite set
        """
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start = 0 if r.start is None else r.start
            stop = r.stop
            if stop is None and not self.is_finite():
                raise ValueError('infinite set')
        else:
            raise ValueError('r must be an integer or a slice')
        count = 0
        tabs = []
        for t in self:
            if count == stop:
                break
            if count >= start:
                tabs.append(t)
            count += 1

        # this is to cope with empty slices endpoints like [:6] or [:]
        if count == stop or stop is None:
            return tabs
        raise IndexError('value out of range')


class ShiftedTableaux_size(ShiftedTableaux, DisjointUnionEnumeratedSets):
    """
    Shifted tableaux of fixed size `n`.

    .. WARNING::

        Input is not checked; please use :class:`ShiftedTableaux` to ensure
        the options are properly parsed.
    """
    def __init__(self, size):
        r"""
        Initializes the class of all shifted tableaux of size ``n``.

        TESTS::

            sage: TestSuite( ShiftedTableaux(0) ).run()
            sage: TestSuite( ShiftedTableaux(3) ).run()
        """
        self._size = Integer(size)
        from .partition import StrictPartitions_size
        DisjointUnionEnumeratedSets.__init__(self,
                family=Family(StrictPartitions_size(self._size), ShiftedTableaux_shape),
                facade=True, keepkey=False
        )

    def _is_a(self, x):
        """
        TESTS::

            sage: ST4 = ShiftedTableaux(4)
            sage: all([st in ST4 for st in ST4])
            True
            sage: ST5 = ShiftedTableaux(5)
            sage: filter(lambda x: x in ST4, ST5)
            []
        """
        if isinstance(x, ShiftedTableau):
            return sum(map(len, x)) == self._size

        return ShiftedTableaux.__contains__(self, x) and sum(map(len, x)) == self._size

    __contains__ = _is_a

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedTableaux(3)    # indirect doctest
            Shifted tableaux of size 3
        """
        return "Shifted tableaux of size %s" % self._size

    def size(self):
        """
        Return the size of the shifted tableaux ``self``.

        EXAMPLES::

            sage: ShiftedTableaux(6).size()
            6
        """
        return self._size

    def random_element(self):
        r"""
        Return a random ``ShiftedTableau``.

        EXAMPLES::

            sage: ShiftedTableaux(5).random_element() # random
            [[1, 4, 5], [2], [3]]
            sage: ShiftedTableaux(0).random_element()
            []
            sage: ShiftedTableaux(1).random_element()
            [[1]]

        TESTS::

            sage: all([ShiftedTableaux(10).random_element() in ShiftedTableaux(10)
            ....:      for i in range(20)])
            True
        """
        return ShiftedTableaux_shape(StrictPartitions(self._size).random_element()).random_element()


class ShiftedTableaux_shape(ShiftedTableaux):
    """
    Shifted tableaux of a fixed shape.
    """
    Element = ShiftedTableau

    def __init__(self, shape):
        r"""
        Initializes the class of semistandard tableaux of shape ``p`` and no
        maximum entry.

        TESTS::

            sage: TestSuite( ShiftedTableaux([3,2,1]) ).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._shape = shape

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedTableaux([3,2,1])    # indirect doctest
            Shifted tableaux of shape [3, 2, 1]
        """
        return "Shifted tableaux of shape {}".format(self._shape)

    def __contains__(self, x):
        """
        TESTS::

            sage: ST421 = ShiftedTableaux([4,2,1])
            sage: all([st in ST421 for st in ST421])
            True
            sage: ST42 = ShiftedTableaux([4,2])
            sage: filter(lambda x: x in ST42, ST421)
            []
        """
        if isinstance(x, ShiftedTableau):
            return [len(row) for row in x] == self._shape

        return (ShiftedTableaux.__contains__(self, x)
                and [len(row) for row in x] == self._shape)

    def _element_constructor_(self, t):
        r"""
        Constructs an object from ``t`` as an element of ``self``, if
        possible.

        INPUT:

        - ``t`` -- data which can be interpreted as a tableau

        OUTPUT:

        - the corresponding tableau object

        TESTS::

            sage: ShiftedTableaux([3])([[1,2,3]]).parent() is ShiftedTableaux([3])
            True
            sage: ShiftedTableaux([3])([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not an element of Shifted tableaux of shape [3]
        """
        if not t in self:
            raise ValueError("{} is not an element of {}".format(t, self))

        return self.element_class(self, t)

    def shape(self):
        """
        Return the shape of the shifted tableaux ``self``.

        EXAMPLES::

            sage: ShiftedTableaux([6,4,3,1]).shape()
            [6, 4, 3, 1]
        """
        return self._shape

    @cached_method
    def size(self):
        """
        Return the shape of the shifted tableaux ``self``.

        EXAMPLES::

            sage: ShiftedTableaux([6,4,3,1]).size()
            14
        """
        return self._shape.size()

    def cardinality(self):
        """
        Return the number of shifted tableaux with the shape of ``shape``.

        This uses the shifted hook-length formula originally due to
        Thrall in 1952.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: T.cardinality()
            5
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T.cardinality()
            1716
        """
        n = sum(self._shape)
        denom = prod(len(shifted_hook_cells(self._shape, (i, j) )) + 1
                     for i, row_len in enumerate(self._shape)
                     for j in range(row_len))
        return factorial(n) // denom

    def random_element(self):
        """
        Return a random shifted tableau with the shape of ``self``.

        This uses a variant of the Greene-Nijenhuis-Wilf hook-walk
        algorithm in [Sag1980]_.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: T.random_element()  # random
            [[1, 2, 4, 5], [3, 6]]
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T.random_element()  # random
            [[1, 2, 3, 5, 7, 10], [4, 6, 9, 12], [8, 11, 14], [13]]
        """
        tableau = [[0] * row_length for row_length in self._shape]
        N = sum(self._shape)
        current_shape = list(self._shape)
        for next_number in reversed(range(1, N+1)):
            cell = random_shifted_cell(current_shape)
            hook_coords = shifted_hook_cells(current_shape, cell)
            while hook_coords:
                cell = hook_coords[ randrange(len(hook_coords)) ]
                hook_coords = shifted_hook_cells(current_shape, cell)
            # done with hook walk, so place next_number in tableau
            tableau[cell[0]][cell[1]] = next_number
            current_shape[cell[0]] -= 1

        return self.element_class(self, tableau)

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
        if self._shape == []:
            yield ShiftedTableau([])
            return

        for l in self.strict_partition_poset().linear_extensions():
            yield self.element_class(self, self.linear_extension_to_shifted_tableau(l))

    def strict_partition_poset(self):
        r"""
        Convert the shifted shape ``shape`` into a poset with
        elements `1, 2, \ldots, n`.

        EXAMPLES::

            sage: ShiftedTableaux([4,2]).strict_partition_poset().cover_relations()
            [[1, 2], [2, 3], [2, 5], [3, 4], [3, 6], [5, 6]]
        """
        elts = range(1, sum(self._shape)+1)
        rels = [(i, i+1) for i in range(1, self._shape[0])]
        tot = self._shape[0]
        prev = self._shape[0]
        for row_len in self._shape[1:]:
            rels += [(i, i+1) for i in range(tot+1, tot+row_len)]
            rels += [(i-prev+1, i) for i in range(tot+1, tot+1+row_len)]
            tot += row_len
            prev = row_len
        return Poset((elts, rels))

    def linear_extension_to_shifted_tableau(self, linear_extension):
        """
        Convert a linear extension for the poset
        :meth:`strict_partition_poset()` into a shifted tableau
        of shape ``shape``.

        EXAMPLES::

            sage: T = ShiftedTableaux([4,2])
            sage: P = T.strict_partition_poset()
            sage: [T.linear_extension_to_shifted_tableau(L) for L in P.linear_extensions()]
            [[[1, 2, 3, 4], [5, 6]],
             [[1, 2, 3, 5], [4, 6]],
             [[1, 2, 3, 6], [4, 5]],
             [[1, 2, 4, 5], [3, 6]],
             [[1, 2, 4, 6], [3, 5]]]
        """
        # tableau in a list, rows not separate
        L = list(Permutation(linear_extension).inverse())
        # determine indices where new rows start
        A = [sum(self._shape[:i]) for i in range(len(self._shape)+1)]
        return [L[A[i]:A[i+1]] for i in range(len(A)-1)]

#####################################################################
## Helper functions

def shifted_hook_cells(shape, pos):
    """
    Determines the coordinates of the cells in the hook corresponding to
    the cell indexed by the tuple ``pos``.

    EXAMPLES::

        sage: from sage.combinat.tableau_shifted import shifted_hook_cells
        sage: shifted_hook_cells([6,4,3,1], (1,1))
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

def random_shifted_cell(shape):
    """
    Return a uniformly random cell in shifted tableaux in ``self``.

    EXAMPLES::

        sage: from sage.combinat.tableau_shifted import random_shifted_cell
        sage: [random_shifted_cell([6,4,3,1]) for i in range(6)]  # random
        [(2, 2), (1, 1), (1, 2), (1, 1), (2, 1), (3, 0)]
    """
    rnd = randrange(0, sum(shape))
    row = 0
    cells_so_far = 0
    while cells_so_far + shape[row] <= rnd:
        cells_so_far += shape[row]
        row += 1
    return (row, rnd - cells_so_far)

