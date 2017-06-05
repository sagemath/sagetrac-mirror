r"""
Set Partitions

AUTHORS:

- Mike Hansen

- MuPAD-Combinat developers (for algorithms and design inspiration).

- Travis Scrimshaw (2013-02-28): Removed ``CombinatorialClass`` and added
  entry point through :class:`SetPartition`.

This module defines a class for immutable partitioning of a set. For
mutable version see :func:`DisjointSet`.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
#*****************************************************************************
from __future__ import print_function, absolute_import, division
from six.moves import range
from six import add_metaclass

from sage.sets.set import Set, is_Set

import itertools

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.rings.infinity import infinity
from sage.rings.integer import Integer

from sage.combinat.misc import IterableFunctionCall
from sage.combinat.combinatorial_map import combinatorial_map
import sage.combinat.subset as subset
from sage.combinat.partition import Partition, Partitions
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.combinat import bell_number, stirling_number2
from sage.combinat.permutation import Permutation
from functools import reduce


@add_metaclass(InheritComparisonClasscallMetaclass)
class SetPartition(ClonableArray):
    """
    A partition of a set.

    A set partition `p` of a set `S` is a partition of `S` into subsets
    called parts and represented as a set of sets. By extension, a set
    partition of a nonnegative integer `n` is the set partition of the
    integers from 1 to `n`. The number of set partitions of `n` is called
    the `n`-th Bell number.

    There is a natural integer partition associated with a set partition,
    namely the nonincreasing sequence of sizes of all its parts.

    There is a classical lattice associated with all set partitions of
    `n`. The infimum of two set partitions is the set partition obtained
    by intersecting all the parts of both set partitions. The supremum
    is obtained by transitive closure of the relation `i` related to `j`
    if and only if they are in the same part in at least one of the set
    partitions.

    We will use terminology from partitions, in particular the *length* of
    a set partition `A = \{A_1, \ldots, A_k\}` is the number of parts of `A`
    and is denoted by `|A| := k`. The *size* of `A` is the cardinality of `S`.
    We will also sometimes use the notation `[n] := \{1, 2, \ldots, n\}`.

    EXAMPLES:

    There are 5 set partitions of the set `\{1,2,3\}`::

        sage: SetPartitions(3).cardinality()
        5

    Here is the list of them::

        sage: SetPartitions(3).list()
        [{{1, 2, 3}},
         {{1}, {2, 3}},
         {{1, 3}, {2}},
         {{1, 2}, {3}},
         {{1}, {2}, {3}}]

    There are 6 set partitions of `\{1,2,3,4\}` whose underlying partition is
    `[2, 1, 1]`::

        sage: SetPartitions(4, [2,1,1]).list()
        [{{1}, {2}, {3, 4}},
         {{1}, {2, 4}, {3}},
         {{1}, {2, 3}, {4}},
         {{1, 4}, {2}, {3}},
         {{1, 3}, {2}, {4}},
         {{1, 2}, {3}, {4}}]

    Since :trac:`14140`, we can create a set partition directly by
    :class:`SetPartition`, which creates the base set by taking the
    union of the parts passed in::

        sage: s = SetPartition([[1,3],[2,4]]); s
        {{1, 3}, {2, 4}}
        sage: s.parent()
        Set partitions
    """
    @staticmethod
    def __classcall_private__(cls, parts, check=True):
        """
        Create a set partition from ``parts`` with the appropriate parent.

        EXAMPLES::

            sage: s = SetPartition([[1,3],[2,4]]); s
            {{1, 3}, {2, 4}}
            sage: s.parent()
            Set partitions
        """
        P = SetPartitions()
        return P.element_class(P, parts)

    def __init__(self, parent, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1,3],[2,4]])
            sage: TestSuite(s).run()
            sage: SetPartition([])
            {}
        """
        self._latex_options = {}
        ClonableArray.__init__(self, parent, sorted(map(Set, s), key=min))

    def check(self):
        """
        Check that we are a valid ordered set partition.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.check()
        """
        assert self in self.parent()

    def __hash__(self):
        """
        Return the hash of ``self``.

        The parent is not included as part of the hash.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = SetPartition([[1], [2,3], [4]])
            sage: B = P([[1], [2,3], [4]])
            sage: hash(A) == hash(B)
            True
        """
        return sum(hash(x) for x in self)

    def __eq__(self, y):
        """
        Check equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = SetPartition([[1], [2,3], [4]])
            sage: B = P([[1], [2,3], [4]])
            sage: A == B
            True
            sage: C = P([[2, 3], [1], [4]])
            sage: A == C
            True
            sage: D = P([[1], [2, 4], [3]])
            sage: A == D
            False
        """
        if not isinstance(y, SetPartition):
            return False
        return list(self) == list(y)

    def __lt__(self, y):
        """
        Check that ``self`` is less than ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: A < B
            True
            sage: C = P([[1,2,4], [3]])
            sage: B < C
            True
            sage: B < B
            False
            sage: D = P([[1,4], [2], [3]])
            sage: E = P([[1,4], [2,3]])
            sage: D < E
            True
            sage: F = P([[1,2,4], [3]])
            sage: E < C
            False
            sage: A < E
            True
            sage: A < C
            True
        """
        if not isinstance(y, SetPartition):
            return False
        return [sorted(_) for _ in self] < [sorted(_) for _ in y]

    def __gt__(self, y):
        """
        Check that ``self`` is greater than ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: B > A
            True
            sage: A > B
            False
        """
        if not isinstance(y, SetPartition):
            return False
        return [sorted(_) for _ in self] > [sorted(_) for _ in y]

    def __le__(self, y):
        """
        Check that ``self`` is less than or equals ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: A <= B
            True
            sage: A <= A
            True
        """
        return self == y or self < y

    def __ge__(self, y):
        """
        Check that ``self`` is greater than or equals ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: B >= A
            True
            sage: B >= B
            True
        """
        return self == y or self > y

    def _cmp_(self, y):
        """
        Return the result of ``cmp``.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: A < B
            True
            sage: cmp(A, B)
            -1
            sage: cmp(B, A)
            1
        """
        if self < y:
            return -1
        if self > y:
            return 1
        return 0

    __cmp__ = _cmp_

    def __mul__(self, other):
        r"""
        The product of the set partitions ``self`` and ``other``.

        The product of two set partitions `B` and `C` is defined as the
        set partition whose parts are the nonempty intersections between
        each part of `B` and each part of `C`. This product is also
        the infimum of `B` and `C` in the classical set partition
        lattice (that is, the coarsest set partition which is finer than
        each of `B` and `C`). Consequently, ``inf`` acts as an alias for
        this method.

        .. SEEALSO::

            :meth:`sup`

        EXAMPLES::

            sage: x = SetPartition([ [1,2], [3,5,4] ])
            sage: y = SetPartition(( (3,1,2), (5,4) ))
            sage: x * y
            {{1, 2}, {3}, {4, 5}}

            sage: S = SetPartitions(4)
            sage: sp1 = S([[2,3,4], [1]])
            sage: sp2 = S([[1,3], [2,4]])
            sage: s = S([[2,4], [3], [1]])
            sage: sp1.inf(sp2) == s
            True

        TESTS:

        Here is a different implementation of the ``__mul__`` method
        (one that was formerly used for the ``inf`` method, before it
        was realized that the methods do the same thing)::

            sage: def mul2(s, t):
            ....:     temp = [ss.intersection(ts) for ss in s for ts in t]
            ....:     temp = filter(lambda x: x != Set([]), temp)
            ....:     return s.__class__(s.parent(), temp)

        Let us check that this gives the same as ``__mul__`` on set
        partitions of `\{1, 2, 3, 4\}`::

            sage: all( all( mul2(s, t) == s * t for s in SetPartitions(4) )
            ....:      for t in SetPartitions(4) )
            True
        """
        new_composition = []
        for B in self:
           for C in other:
                BintC = B.intersection(C)
                if BintC:
                    new_composition.append(BintC)
        return SetPartition(new_composition)

    inf = __mul__

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: S([[1,3],[2,4]])
            {{1, 3}, {2, 4}}
        """
        return '{' + ', '.join(('{' + repr(sorted(x))[1:-1] + '}' for x in self)) + '}'

    def set_latex_options(self, **kwargs):
        r"""
        Set the latex options for use in the ``_latex_`` function

        - ``tikz_scale`` -- (default: 1) scale for use with tikz package

        - ``plot`` -- (default: ``None``) ``None`` returns the set notation,
          ``linear`` returns a linear plot, ``cyclic`` returns a cyclic
          plot

        - ``color`` -- (default: ``'black'``) the arc colors

        - ``fill`` -- (default: ``False``) if ``True`` then fills ``color``,
          else you can pass in a color to alter the fill color -
          *only works with cyclic plot*

        - ``show_labels`` -- (default: ``True``) if ``True`` shows labels -
          *only works with plots*

        - ``radius`` -- (default: ``"1cm"``) radius of circle for cyclic
          plot - *only works with cyclic plot*

        - ``angle`` -- (default: 0) angle for linear plot

        EXAMPLES::

            sage: SP = SetPartition([[1,6], [3,5,4]])
            sage: SP.set_latex_options(tikz_scale=2,plot='linear',fill=True,color='blue',angle=45)
            sage: SP.set_latex_options(plot='cyclic')
            sage: SP.latex_options()
            {'angle': 45,
             'color': 'blue',
             'fill': True,
             'plot': 'cyclic',
             'radius': '1cm',
             'show_labels': True,
             'tikz_scale': 2}

        """
        valid_args = ['tikz_scale', 'plot', 'color', 'fill', 'show_labels',
                      'radius', 'angle']

        for key in kwargs:
            if key not in valid_args:
                raise ValueError("unknown keyword argument: %s"%key)
            if key == 'plot':
                if not (kwargs['plot'] == 'cyclic'
                        or kwargs['plot'] == 'linear'
                        or kwargs['plot'] is None):
                    raise ValueError("plot must be None, 'cyclic', or 'linear'")

        for opt in kwargs:
            self._latex_options[opt] = kwargs[opt]

    def latex_options(self):
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary. The default values are set using the global options.

        Options can be found in :meth:`set_latex_options`

        EXAMPLES::

            sage: SP = SetPartition([[1,6], [3,5,4]]); SP.latex_options()
            {'angle': 0,
             'color': 'black',
             'fill': False,
             'plot': None,
             'radius': '1cm',
             'show_labels': True,
             'tikz_scale': 1}
        """
        opts = self._latex_options.copy()
        if "tikz_scale" not in opts:
            opts["tikz_scale"] = 1
        if "plot" not in opts:
            opts["plot"] = None
        if "color" not in opts:
            opts['color'] = 'black'
        if "fill" not in opts:
            opts["fill"] = False
        if "show_labels" not in opts:
            opts['show_labels'] = True
        if "radius" not in opts:
            opts['radius'] = "1cm"
        if "angle" not in opts:
            opts['angle'] = 0
        return opts

    def _latex_(self):
        r"""
        Return a `\LaTeX` string representation of ``self``.

        EXAMPLES::

            sage: x = SetPartition([[1,2], [3,5,4]])
            sage: latex(x)
            \{\{1, 2\}, \{3, 4, 5\}\}

            sage: x.set_latex_options(plot='linear', angle=25, color='red')
            sage: latex(x)
            \begin{tikzpicture}[scale=1]
            \node[below=.05cm] at (0,0) {$1$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (0) at (0,0) {};
            \node[below=.05cm] at (1,0) {$2$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (1) at (1,0) {};
            \node[below=.05cm] at (2,0) {$3$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (2) at (2,0) {};
            \node[below=.05cm] at (3,0) {$4$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (3) at (3,0) {};
            \node[below=.05cm] at (4,0) {$5$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (4) at (4,0) {};
            \draw[color=red] (1) to [out=115,in=65] (0);
            \draw[color=red] (3) to [out=115,in=65] (2);
            \draw[color=red] (4) to [out=115,in=65] (3);
            \end{tikzpicture}

            sage: p = SetPartition([['a','c'],['b',1],[20]])
            sage: p.set_latex_options(plot='cyclic', color='blue', fill=True, tikz_scale=2)
            sage: latex(p)
            \begin{tikzpicture}[scale=2]
            \draw (0,0) circle [radius=1cm];
            \node[label=90:1] (0) at (90:1cm) {};
            \node[label=18:20] (1) at (18:1cm) {};
            \node[label=-54:a] (2) at (-54:1cm) {};
            \node[label=-126:b] (3) at (-126:1cm) {};
            \node[label=-198:c] (4) at (-198:1cm) {};
            \draw[-,thick,color=blue,fill=blue,fill opacity=0.1] (2.center) -- (4.center) -- cycle;
            \draw[-,thick,color=blue,fill=blue,fill opacity=0.1] (0.center) -- (3.center) -- cycle;
            \draw[-,thick,color=blue,fill=blue,fill opacity=0.1] (1.center) -- cycle;
            \fill[color=black] (0) circle (1.5pt);
            \fill[color=black] (1) circle (1.5pt);
            \fill[color=black] (2) circle (1.5pt);
            \fill[color=black] (3) circle (1.5pt);
            \fill[color=black] (4) circle (1.5pt);
            \end{tikzpicture}
        """
        latex_options = self.latex_options()
        if latex_options["plot"] is None:
            return repr(self).replace("{",r"\{").replace("}",r"\}")

        from sage.misc.latex import latex
        latex.add_package_to_preamble_if_available("tikz")
        res = "\\begin{{tikzpicture}}[scale={}]\n".format(latex_options['tikz_scale'])

        cardinality = self.base_set_cardinality()
        from sage.rings.integer_ring import ZZ
        if all(x in ZZ for x in self.base_set()):
            sort_key = ZZ
        else:
            sort_key = str
        base_set = sorted(self.base_set(), key=sort_key)
        color = latex_options['color']

        # If we want cyclic plots
        if latex_options['plot'] == 'cyclic':
            degrees = 360 // cardinality
            radius = latex_options['radius']

            res += "\\draw (0,0) circle [radius={}];\n".format(radius)

            # Add nodes
            for k,i in enumerate(base_set):
                location = (cardinality - k) * degrees - 270
                if latex_options['show_labels']:
                    res += "\\node[label={}:{}]".format(location, i)
                else:
                    res += "\\node"
                res += " ({}) at ({}:{}) {{}};\n".format(k, location, radius)

            # Setup partitions
            for partition in sorted(self, key=str):
                res += "\\draw[-,thick,color="+color
                if latex_options['fill'] is not False:
                    if isinstance(latex_options['fill'], str):
                        res += ",fill=" + latex_options['fill']
                    else:
                        res += ",fill={},fill opacity=0.1".format(color)
                res += "] "
                res += " -- ".join("({}.center)".format(base_set.index(j))
                                   for j in sorted(partition, key=sort_key))
                res += " -- cycle;\n"

            # Draw the circles on top
            for k in range(len(base_set)):
                res += "\\fill[color=black] ({}) circle (1.5pt);\n".format(k)

        # If we want line plots
        elif latex_options['plot'] == 'linear':
            angle = latex_options['angle']
            # setup line
            for k,i in enumerate(base_set):
                if latex_options['show_labels']:
                    res += "\\node[below=.05cm] at ({},0) {{${}$}};\n".format(k, i)
                res += "\\node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] "
                res += "({k}) at ({k},0) {{}};\n".format(k=k)

            # setup arcs
            for partition in sorted(self, key=str):
                p = sorted(partition, key=sort_key)
                if len(p) <= 1:
                    continue
                for k in range(1, len(p)):
                    res += "\\draw[color={}] ({})".format(color, base_set.index(p[k]))
                    res += " to [out={},in={}] ".format(90+angle, 90-angle)
                    res += "({});\n".format(base_set.index(p[k-1]))
        else:
            raise ValueError("plot must be None, 'cyclic', or 'linear'")

        res += "\\end{tikzpicture}"
        return res

    cardinality = ClonableArray.__len__

    def sup(self, t):
        """
        Return the supremum of ``self`` and ``t`` in the classical set
        partition lattice.

        The supremum of two set partitions `B` and `C` is obtained as the
        transitive closure of the relation which relates `i` to `j` if
        and only if `i` and `j` are in the same part in at least
        one of the set partitions `B` and `C`.

        .. SEEALSO::

            :meth:`__mul__`

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: sp1 = S([[2,3,4], [1]])
            sage: sp2 = S([[1,3], [2,4]])
            sage: s = S([[1,2,3,4]])
            sage: sp1.sup(sp2) == s
            True
        """
        res = Set(list(self))
        for p in t:
            inters = Set([x for x in list(res) if x.intersection(p) != Set([])])
            res = res.difference(inters).union(_set_union(inters))
        return SetPartition(res)

    def pipe(self, other):
        r"""
        Return the pipe of the set partitions ``self`` and ``other``.

        The pipe of two set partitions is defined as follows:

        For any integer `k` and any subset `I` of `\ZZ`, let `I + k`
        denote the subset of `\ZZ` obtained by adding `k` to every
        element of `k`.

        If `B` and `C` are set partitions of `[n]` and `[m]`,
        respectively, then the pipe of `B` and `C` is defined as the
        set partition

        .. MATH::

            \{ B_1, B_2, \ldots, B_b,
            C_1 + n, C_2 + n, \ldots, C_c + n \}

        of `[n+m]`, where `B = \{ B_1, B_2, \ldots, B_b \}` and
        `C = \{ C_1, C_2, \ldots, C_c \}`. This pipe is denoted by
        `B | C`.

        EXAMPLES::

            sage: SetPartition([[1,3],[2,4]]).pipe(SetPartition([[1,3],[2]]))
            {{1, 3}, {2, 4}, {5, 7}, {6}}
            sage: SetPartition([]).pipe(SetPartition([[1,2],[3,5],[4]]))
            {{1, 2}, {3, 5}, {4}}
            sage: SetPartition([[1,2],[3,5],[4]]).pipe(SetPartition([]))
            {{1, 2}, {3, 5}, {4}}
            sage: SetPartition([[1,2],[3]]).pipe(SetPartition([[1]]))
            {{1, 2}, {3}, {4}}
        """
        # Note: GIGO if self and other are not standard.
        parts = list(self)
        n = self.base_set_cardinality()
        for newpart in other:
            raised_newpart = Set([i + n for i in newpart])
            parts.append(raised_newpart)
        return SetPartition(parts)

    @combinatorial_map(name='shape')
    def shape(self):
        r"""
        Return the integer partition whose parts are the sizes of the sets
        in ``self``.

        EXAMPLES::

            sage: S = SetPartitions(5)
            sage: x = S([[1,2], [3,5,4]])
            sage: x.shape()
            [3, 2]
            sage: y = S([[2], [3,1], [5,4]])
            sage: y.shape()
            [2, 2, 1]
        """
        return Partition(sorted(map(len, self), reverse=True))

    # we define aliases for shape()
    shape_partition = shape
    to_partition = shape

    @combinatorial_map(name='to permutation')
    def to_permutation(self):
        r"""
        Convert a set partition of `\{1,...,n\}` to a permutation by considering
        the blocks of the partition as cycles.

        The cycles are such that the number of excedences is maximised, that is,
        each cycle is of the form `(a_1,a_2, ...,a_k)` with `a_1<a_2<...<a_k`.

        EXAMPLES::

            sage: s = SetPartition([[1,3],[2,4]])
            sage: s.to_permutation()
            [3, 4, 1, 2]

        """
        return Permutation(tuple( map(tuple, self.standard_form()) ))

    def standard_form(self):
        r"""
        Return ``self`` as a list of lists.

        When the ground set is totally ordered, the elements of each
        block are listed in increasing order.

        This is not related to standard set partitions (which simply
        means set partitions of `[n] = \{ 1, 2, \ldots , n \}` for some
        integer `n`) or standardization (:meth:`standardization`).

        EXAMPLES::

            sage: [x.standard_form() for x in SetPartitions(4, [2,2])]
            [[[1, 2], [3, 4]], [[1, 3], [2, 4]], [[1, 4], [2, 3]]]

        TESTS::

            sage: SetPartition([(1, 9, 8), (2, 3, 4, 5, 6, 7)]).standard_form()
            [[1, 8, 9], [2, 3, 4, 5, 6, 7]]
        """
        return [sorted(_) for _ in self]

    def apply_permutation(self, p):
        r"""
        Apply ``p`` to the underlying set of ``self``.

        INPUT:

        - ``p`` -- a permutation

        EXAMPLES::

            sage: x = SetPartition([[1,2], [3,5,4]])
            sage: p = Permutation([2,1,4,5,3])
            sage: x.apply_permutation(p)
            {{1, 2}, {3, 4, 5}}
            sage: q = Permutation([3,2,1,5,4])
            sage: x.apply_permutation(q)
            {{1, 4, 5}, {2, 3}}
        """
        return self.__class__(self.parent(), [Set(map(p, B)) for B in self])

    def is_noncrossing(self):
        r"""
        Check if ``self`` is noncrossing.

        EXAMPLES::

            sage: x = SetPartition([[1,2],[3,4]])
            sage: x.is_noncrossing()
            True
            sage: x = SetPartition([[1,3],[2,4]])
            sage: x.is_noncrossing()
            False

        AUTHOR: Florent Hivert
        """
        l = list(self)
        mins = [min(_) for _ in l]
        maxs = [max(_) for _ in l]

        for i in range(1, len(l)):
            for j in range(i):
                poss = [mins[i], maxs[i], mins[j], maxs[j]]
                possort = sorted(poss)
                cont = [possort.index(mins[i]), possort.index(maxs[i])]
                if cont == [0,2] or cont == [1,3]:
                    return False
                if (cont == [0,3] and
                    any(mins[j] < x < maxs[j] for x in l[i])):
                    return False
                if (cont == [1,2] and
                    any(mins[i] < x < maxs[i] for x in l[j])):
                    return False
        return True

    def is_atomic(self):
        """
        Return if ``self`` is an atomic set partition.

        A (standard) set partition `A` can be split if there exist `j < i`
        such that `\max(A_j) < \min(A_i)` where `A` is ordered by minimal
        elements. This means we can write `A = B | C` for some nonempty set
        partitions `B` and `C`. We call a set partition *atomic* if it
        cannot be split and is nonempty. Here, the pipe symbol
        `|` is as defined in method :meth:`pipe`.

        EXAMPLES::

            sage: SetPartition([[1,3], [2]]).is_atomic()
            True
            sage: SetPartition([[1,3], [2], [4]]).is_atomic()
            False
            sage: SetPartition([[1], [2,4], [3]]).is_atomic()
            False
            sage: SetPartition([[1,2,3,4]]).is_atomic()
            True
            sage: SetPartition([[1, 4], [2], [3]]).is_atomic()
            True
            sage: SetPartition([]).is_atomic()
            False
        """
        if len(self) == 0:
            return False
        maximum_so_far = max(self[0])
        for S in self[1:]:
            if maximum_so_far < min(S):
                return False
            maximum_so_far = max(maximum_so_far, max(S))
        return True

    def base_set(self):
        """
        Return the base set of ``self``, which is the union of all parts
        of ``self``.

        EXAMPLES::

            sage: SetPartition([[1], [2,3], [4]]).base_set()
            {1, 2, 3, 4}
            sage: SetPartition([[1,2,3,4]]).base_set()
            {1, 2, 3, 4}
            sage: SetPartition([]).base_set()
            {}
        """
        return reduce(lambda x,y: x.union(y), self, Set([]))

    def base_set_cardinality(self):
        """
        Return the cardinality of the base set of ``self``, which is the sum
        of the sizes of the parts of ``self``.

        This is also known as the *size* (sometimes the *weight*) of
        a set partition.

        EXAMPLES::

            sage: SetPartition([[1], [2,3], [4]]).base_set_cardinality()
            4
            sage: SetPartition([[1,2,3,4]]).base_set_cardinality()
            4
        """
        return sum(len(x) for x in self)

    size = base_set_cardinality

    def standardization(self):
        """
        Return the standardization of ``self``.

        Given a set partition `A = \{A_1, \ldots, A_n\}` of an ordered
        set `S`, the standardization of `A` is the set partition of
        `\{1, 2, \ldots, |S|\}` obtained by replacing the elements of
        the parts of `A` by the integers `1, 2, \ldots, |S|` in such
        a way that their relative order is preserved (i. e., the
        smallest element in the whole set partition is replaced by
        `1`, the next-smallest by `2`, and so on).

        EXAMPLES::

            sage: SetPartition([[4], [1, 3]]).standardization()
            {{1, 2}, {3}}
            sage: SetPartition([[4], [6, 3]]).standardization()
            {{1, 3}, {2}}
            sage: SetPartition([]).standardization()
            {}
        """
        if len(self) == 0:
            return self
        temp = [list(_) for _ in self]
        mins = [min(p) for p in temp]
        over_max = max([max(p) for p in temp]) + 1
        ret = [[] for i in range(len(temp))]
        cur = 1
        while min(mins) != over_max:
            m = min(mins)
            i = mins.index(m)
            ret[i].append(cur)
            cur += 1
            temp[i].pop(temp[i].index(m))
            if len(temp[i]) != 0:
                mins[i] = min(temp[i])
            else:
                mins[i] = over_max
        return SetPartition(ret)

    def restriction(self, I):
        """
        Return the restriction of ``self`` to a subset ``I``
        (which is given as a set or list or any other iterable).

        EXAMPLES::

            sage: A = SetPartition([[1], [2,3]])
            sage: A.restriction([1,2])
            {{1}, {2}}
            sage: A.restriction([2,3])
            {{2, 3}}
            sage: A.restriction([])
            {}
            sage: A.restriction([4])
            {}
        """
        ret = []
        for part in self:
            newpart = [i for i in part if i in I]
            if len(newpart) != 0:
                ret.append(newpart)
        return SetPartition(ret)

    def ordered_set_partition_action(self, s):
        r"""
        Return the action of an ordered set partition ``s`` on ``self``.

        Let `A = \{A_1, A_2, \ldots, A_k\}` be a set partition of some
        set `S` and `s` be an ordered set partition (i.e., set composition)
        of a subset of `[k]`. Let `A^{\downarrow}` denote the standardization
        of `A`, and `A_{\{ i_1, i_2, \ldots, i_m \}}` denote the sub-partition
        `\{A_{i_1}, A_{i_2}, \ldots, A_{i_m}\}` for any subset
        `\{i_1, \ldots, i_m\}` of `\{1, \ldots, k\}`. We define the set
        partition `s(A)` by

        .. MATH::

            s(A) = A_{s_1}^{\downarrow} | A_{s_2}^{\downarrow} | \cdots
            | A_{s_q}^{\downarrow}.

        where `s = (s_1, s_2, \ldots, s_q)`. Here, the pipe symbol
        `|` is as defined in method :meth:`pipe`.

        This is `s[A]` in section 2.3 in [LM2011]_.

        INPUT:

        - ``s`` -- an ordered set partition with base set a subset
          of `\{1, \ldots, k\}`

        EXAMPLES::

            sage: A = SetPartition([[1], [2,4], [3]])
            sage: s = OrderedSetPartition([[1,3], [2]])
            sage: A.ordered_set_partition_action(s)
            {{1}, {2}, {3, 4}}
            sage: s = OrderedSetPartition([[2,3], [1]])
            sage: A.ordered_set_partition_action(s)
            {{1, 3}, {2}, {4}}

        We create Figure 1 in [LM2011]_ (we note that there is a typo in the
        lower-left corner of the table in the published version of the
        paper, whereas the arXiv version gives the correct partition)::

            sage: A = SetPartition([[1,3], [2,9], [4,5,8], [7]])
            sage: B = SetPartition([[1,3], [2,8], [4,5,6], [7]])
            sage: C = SetPartition([[1,5], [2,8], [3,4,6], [7]])
            sage: s = OrderedSetPartition([[1,3], [2]])
            sage: t = OrderedSetPartition([[2], [3,4]])
            sage: u = OrderedSetPartition([[1], [2,3,4]])
            sage: A.ordered_set_partition_action(s)
            {{1, 2}, {3, 4, 5}, {6, 7}}
            sage: A.ordered_set_partition_action(t)
            {{1, 2}, {3, 4, 6}, {5}}
            sage: A.ordered_set_partition_action(u)
            {{1, 2}, {3, 8}, {4, 5, 7}, {6}}
            sage: B.ordered_set_partition_action(s)
            {{1, 2}, {3, 4, 5}, {6, 7}}
            sage: B.ordered_set_partition_action(t)
            {{1, 2}, {3, 4, 5}, {6}}
            sage: B.ordered_set_partition_action(u)
            {{1, 2}, {3, 8}, {4, 5, 6}, {7}}
            sage: C.ordered_set_partition_action(s)
            {{1, 4}, {2, 3, 5}, {6, 7}}
            sage: C.ordered_set_partition_action(t)
            {{1, 2}, {3, 4, 5}, {6}}
            sage: C.ordered_set_partition_action(u)
            {{1, 2}, {3, 8}, {4, 5, 6}, {7}}

        REFERENCES:

        .. [LM2011] \A. Lauve, M. Mastnak. *The primitives and antipode in
           the Hopf algebra of symmetric functions in noncommuting variables*.
           Advances in Applied Mathematics. **47** (2011). 536-544.
           :arxiv:`1006.0367v3` :doi:`10.1016/j.aam.2011.01.002`.
        """
        cur = 1
        ret = []
        for part in s:
            sub_parts = [list(self[i-1]) for i in part] # -1 for indexing
            # Standardizing sub_parts (the cur variable not being reset
            # to 1 gives us the offset we want):
            mins = [min(_) for _ in sub_parts]
            over_max = max(map(max, sub_parts)) + 1
            temp = [[] for i in range(len(part))]
            while min(mins) != over_max:
                m = min(mins)
                i = mins.index(m)
                temp[i].append(cur)
                cur += 1
                sub_parts[i].pop(sub_parts[i].index(m))
                if len(sub_parts[i]) != 0:
                    mins[i] = min(sub_parts[i])
                else:
                    mins[i] = over_max
            ret += temp
        return SetPartition(ret)

    def refinements(self):
        """
        Return a list of refinements of ``self``.

        .. SEEALSO::

            :meth:`coarsenings`

        EXAMPLES::

            sage: SetPartition([[1,3],[2,4]]).refinements()
            [{{1, 3}, {2, 4}},
             {{1, 3}, {2}, {4}},
             {{1}, {2, 4}, {3}},
             {{1}, {2}, {3}, {4}}]
            sage: SetPartition([[1],[2,4],[3]]).refinements()
            [{{1}, {2, 4}, {3}}, {{1}, {2}, {3}, {4}}]
            sage: SetPartition([]).refinements()
            [{}]
        """
        L = [SetPartitions(part) for part in self]
        return [SetPartition(sum(map(list, x), [])) for x in itertools.product(*L)]

    def coarsenings(self):
        """
        Return a list of coarsenings of ``self``.

        .. SEEALSO::

            :meth:`refinements`

        EXAMPLES::

            sage: SetPartition([[1,3],[2,4]]).coarsenings()
            [{{1, 2, 3, 4}}, {{1, 3}, {2, 4}}]
            sage: SetPartition([[1],[2,4],[3]]).coarsenings()
            [{{1, 2, 3, 4}},
             {{1}, {2, 3, 4}},
             {{1, 3}, {2, 4}},
             {{1, 2, 4}, {3}},
             {{1}, {2, 4}, {3}}]
            sage: SetPartition([]).coarsenings()
            [{}]
        """
        SP = SetPartitions(len(self))
        def union(s):
            # Return the partition obtained by combining, for every
            # part of s, those parts of self which are indexed by
            # the elements of this part of s into a single part.
            ret = []
            for part in s:
                cur = Set([])
                for i in part:
                    cur = cur.union(self[i-1]) # -1 for indexing
                ret.append(cur)
            return ret
        return [SetPartition(union(s)) for s in SP]

    def strict_coarsenings(self):
        r"""
        Return all strict coarsenings of ``self``.

        Strict coarsening is the binary relation on set partitions
        defined as the transitive-and-reflexive closure of the
        relation `\prec` defined as follows: For two set partitions
        `A` and `B`, we have `A \prec B` if there exist parts
        `A_i, A_j` of `A` such that `\max(A_i) < \min(A_j)` and
        `B = A \setminus \{A_i, A_j\} \cup \{ A_i \cup A_j \}`.

        EXAMPLES::

            sage: A = SetPartition([[1],[2,3],[4]])
            sage: A.strict_coarsenings()
            [{{1}, {2, 3}, {4}}, {{1, 2, 3}, {4}}, {{1, 4}, {2, 3}},
             {{1}, {2, 3, 4}}, {{1, 2, 3, 4}}]
            sage: SetPartition([[1],[2,4],[3]]).strict_coarsenings()
            [{{1}, {2, 4}, {3}}, {{1, 2, 4}, {3}}, {{1, 3}, {2, 4}}]
            sage: SetPartition([]).strict_coarsenings()
            [{}]
        """
        # This is more or less generic code for computing a
        # transitive-and-reflexive closure by depth-first search.
        todo = [self]
        visited = set([self])
        ret = [self]
        while todo:
            A = todo.pop()
            for i, part in enumerate(A):
                for j, other in enumerate(A[i+1:]):
                    if max(part) < min(other):
                        next = A[:i]
                        next.append(part.union(other))
                        next += A[i+1:i+1+j] + A[i+j+2:]
                        next = SetPartition(next)
                        if next not in visited:
                            todo.append(next)
                            visited.add(next)
                            ret.append(next)
        return ret

    def arcs(self):
        r"""
        Return ``self`` as a list of arcs.

        Consider a set partition `X = \{X_1, X_2, \ldots, X_n\}`,
        then the arcs are the pairs `i - j` such that `i,j \in X_k`
        for some `k`.

        EXAMPLES::

            sage: A = SetPartition([[1],[2,3],[4]])
            sage: A.arcs()
            [(2, 3)]
            sage: B = SetPartition([[1,3,6,7],[2,5],[4]])
            sage: B.arcs()
            [(1, 3), (3, 6), (6, 7), (2, 5)]
        """
        arcs = []
        for p in self:
            p = sorted(p)
            for i in range(len(p)-1):
                arcs.append((p[i], p[i+1]))
        return arcs

    def plot(self, angle=None, color='black', base_set_dict=None):
        r"""
        Return a plot of ``self``.

        INPUT:

        - ``angle`` -- (default: `\pi/4`) the angle at which the arcs take off
          (if angle is negative, the arcs are drawn below the horizontal line)

        - ``color`` -- (default: ``'black'``) color of the arcs

        - ``base_set_dict`` -- (optional) dictionary with keys elements
          of :meth:`base_set()` and values as integer or float

        EXAMPLES::

            sage: p = SetPartition([[1,10,11],[2,3,7],[4,5,6],[8,9]])
            sage: p.plot()
            Graphics object consisting of 29 graphics primitives

        .. PLOT::

            p = SetPartition([[1,10,11],[2,3,7],[4,5,6],[8,9]])
            sphinx_plot(p.plot())

        ::

            sage: p = SetPartition([[1,3,4],[2,5]])
            sage: print(p.plot().description())
            Point set defined by 1 point(s):    [(0.0, 0.0)]
            Point set defined by 1 point(s):    [(1.0, 0.0)]
            Point set defined by 1 point(s):    [(2.0, 0.0)]
            Point set defined by 1 point(s):    [(3.0, 0.0)]
            Point set defined by 1 point(s):    [(4.0, 0.0)]
            Text '1' at the point (0.0,-0.1)
            Text '2' at the point (1.0,-0.1)
            Text '3' at the point (2.0,-0.1)
            Text '4' at the point (3.0,-0.1)
            Text '5' at the point (4.0,-0.1)
            Arc with center (1.0,-1.0) radii (1.41421356237,1.41421356237)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
            Arc with center (2.5,-0.5) radii (0.707106781187,0.707106781187)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
            Arc with center (2.5,-1.5) radii (2.12132034356,2.12132034356)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
            sage: p = SetPartition([['a','c'],['b','d'],['e']])
            sage: print(p.plot().description())
            Point set defined by 1 point(s):  [(0.0, 0.0)]
            Point set defined by 1 point(s):    [(1.0, 0.0)]
            Point set defined by 1 point(s):    [(2.0, 0.0)]
            Point set defined by 1 point(s):    [(3.0, 0.0)]
            Point set defined by 1 point(s):    [(4.0, 0.0)]
            Text 'a' at the point (0.0,-0.1)
            Text 'b' at the point (1.0,-0.1)
            Text 'c' at the point (2.0,-0.1)
            Text 'd' at the point (3.0,-0.1)
            Text 'e' at the point (4.0,-0.1)
            Arc with center (1.0,-1.0) radii (1.41421356237,1.41421356237)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
            Arc with center (2.0,-1.0) radii (1.41421356237,1.41421356237)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
            sage: p = SetPartition([['a','c'],['b','d'],['e']])
            sage: print(p.plot(base_set_dict={'a':0,'b':1,'c':2,'d':-2.3,'e':5.4}).description())
            Point set defined by 1 point(s):    [(-2.3, 0.0)]
            Point set defined by 1 point(s):    [(0.0, 0.0)]
            Point set defined by 1 point(s):    [(1.0, 0.0)]
            Point set defined by 1 point(s):    [(2.0, 0.0)]
            Point set defined by 1 point(s):    [(5.4, 0.0)]
            Text 'a' at the point (0.0,-0.1)
            Text 'b' at the point (1.0,-0.1)
            Text 'c' at the point (2.0,-0.1)
            Text 'd' at the point (-2.3,-0.1)
            Text 'e' at the point (5.4,-0.1)
            Arc with center (-0.65,-1.65) radii (2.33345237792,2.33345237792)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
            Arc with center (1.0,-1.0) radii (1.41421356237,1.41421356237)
             angle 0.0 inside the sector (0.785398163397,2.35619449019)
        """
        from sage.plot.graphics import Graphics
        from sage.plot.point import point
        from sage.plot.text import text
        from sage.plot.arc import arc
        from sage.functions.other import sqrt
        from sage.symbolic.constants import pi
        from sage.functions.trig import tan, sin
        from sage.functions.generalized import sgn

        diag = Graphics()
        sorted_vertices_list = list(self.base_set())
        sorted_vertices_list.sort()

        if angle is None:
            angle = pi / 4

        if base_set_dict is not None:
            vertices_dict = base_set_dict
        else:
            vertices_dict = {val: pos for pos,val in enumerate(sorted_vertices_list)}

        for elt in vertices_dict:
            pos = vertices_dict[elt]
            diag += point((pos,0), size=30, color=color)
            diag += text(elt, (pos, -sgn(angle)*0.1), color=color)
            # TODO: change 0.1 to something proportional to the height of the picture

        for (k,j) in self.arcs():
            pos_k,pos_j = float(vertices_dict[k]),float(vertices_dict[j])
            center = ((pos_k+pos_j) / 2, -abs(pos_j-pos_k) / (2*tan(angle)))
            r1 = abs((pos_j-pos_k) / (2*sin(angle)))
            sector = (sgn(angle) * (pi/2 - angle), sgn(angle) * (pi/2 + angle))
            diag += arc(center=center, r1=r1, sector=sector, color=color)

        diag.axes(False)
        return diag

class SetPartitions(UniqueRepresentation, Parent):
    r"""
    An (unordered) partition of a set `S` is a set of pairwise
    disjoint nonempty subsets with union `S`, and is represented
    by a sorted list of such subsets.

    ``SetPartitions(s)`` returns the class of all set partitions of the set
    ``s``, which can be given as a set or a string; if a string, each
    character is considered an element.

    ``SetPartitions(n)``, where ``n`` is an integer, returns the class of
    all set partitions of the set `\{1, 2, \ldots, n\}`.

    You may specify a second argument `k`. If `k` is an integer,
    :class:`SetPartitions` returns the class of set partitions into `k` parts;
    if it is an integer partition, :class:`SetPartitions` returns the class of
    set partitions whose block sizes correspond to that integer partition.

    The Bell number `B_n`, named in honor of Eric Temple Bell,
    is the number of different partitions of a set with `n` elements.

    EXAMPLES::

        sage: S = [1,2,3,4]
        sage: SetPartitions(S,2)
        Set partitions of {1, 2, 3, 4} with 2 parts
        sage: SetPartitions([1,2,3,4], [3,1]).list()
        [{{1}, {2, 3, 4}}, {{1, 3, 4}, {2}}, {{1, 2, 4}, {3}}, {{1, 2, 3}, {4}}]
        sage: SetPartitions(7, [3,3,1]).cardinality()
        70

    In strings, repeated letters are not considered distinct as of
    :trac:`14140`::

        sage: SetPartitions('abcde').cardinality()
        52
        sage: SetPartitions('aabcd').cardinality()
        15

    REFERENCES:

    - :wikipedia:`Partition_of_a_set`
    """
    @staticmethod
    def __classcall_private__(cls, s=None, part=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: T = SetPartitions([1,2,3,4])
            sage: S is T
            True
        """
        if s is None:
            return SetPartitions_all()
        if isinstance(s, (int, Integer)):
            s = frozenset(range(1, s+1))
        else:
            try:
                if s.cardinality() == infinity:
                    raise ValueError("The set must be finite")
            except AttributeError:
                pass
            s = frozenset(s)

        if part is not None:
            if isinstance(part, (int, Integer)):
                if len(s) < part:
                    raise ValueError("part must be <= len(set)")
                else:
                    return SetPartitions_setn(s, part)
            else:
                if part not in Partitions(len(s)):
                    raise ValueError("part must be a partition of %s"%len(s))
                else:
                    return SetPartitions_setparts(s, Partition(part))
        else:
            return SetPartitions_set(s)

    def __contains__(self, x):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: SA = SetPartitions()
            sage: all(sp in SA for sp in S)
            True
            sage: Set([Set([1,2]),Set([3,7])]) in SA
            True
            sage: Set([Set([1,2]),Set([2,3])]) in SA
            False
            sage: Set([]) in SA
            True
        """
        # x must be a set
        if not (isinstance(x, (SetPartition, set, frozenset)) or is_Set(x)):
            return False

        # Check that all parts are disjoint
        base_set = reduce( lambda x,y: x.union(y), map(Set, x), Set([]) )
        if len(base_set) != sum(map(len, x)):
            return False

        # Check to make sure each element of x is a set
        for s in x:
            if not (isinstance(s, (set, frozenset)) or is_Set(s)):
                return False

        return True

    def _element_constructor_(self, s):
        """
        Construct an element of ``self`` from ``s``.

        INPUT:

        - ``s`` -- A set of sets

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: elt = S([[1,3],[2,4]]); elt
            {{1, 3}, {2, 4}}
            sage: P = SetPartitions()
            sage: P(elt).parent() is P
            True
            sage: S = SetPartitions([])
            sage: S([])
            {}
        """
        if isinstance(s, SetPartition):
            if s.parent() is self:
                return s
            if isinstance(s.parent(), SetPartitions):
                return self.element_class(self, list(s))
            raise ValueError("cannot convert %s into an element of %s"%(s, self))
        return self.element_class(self, s)

    Element = SetPartition

    def _iterator_part(self, part):
        """
        Return an iterator for the set partitions with block sizes
        corresponding to the partition ``part``.

        INPUT:

        -  ``part`` -- a :class:`Partition` object

        EXAMPLES::

            sage: S = SetPartitions(3)
            sage: it = S._iterator_part(Partition([1,1,1]))
            sage: list(sorted(map(list, next(it))))
            [[1], [2], [3]]
            sage: S21 = SetPartitions(3,Partition([2,1]))
            sage: len(list(S._iterator_part(Partition([2,1])))) == S21.cardinality()
            True
        """
        nonzero = []
        expo = [0] + part.to_exp()

        for i in range(len(expo)):
            if expo[i] != 0:
                nonzero.append([i, expo[i]])

        taillesblocs = [(x[0])*(x[1]) for x in nonzero]

        blocs = OrderedSetPartitions(self._set, taillesblocs)

        for b in blocs:
            lb = [IterableFunctionCall(_listbloc, nonzero[i][0], nonzero[i][1], b[i]) for i in range(len(nonzero))]
            for x in itertools.product(*lb):
                yield _union(x)

    def is_less_than(self, s, t):
        r"""
        Check if `s < t` in the refinement ordering on set partitions.

        This means that `s` is a refinement of `t` and satisfies
        `s \neq t`.

        A set partition `s` is said to be a refinement of a set
        partition `t` of the same set if and only if each part of
        `s` is a subset of a part of `t`.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1,3],[2,4]])
            sage: t = S([[1],[2],[3],[4]])
            sage: S.is_less_than(t, s)
            True
            sage: S.is_less_than(s, t)
            False
            sage: S.is_less_than(s, s)
            False
        """
        if hasattr(s.parent(), "_set"):
            S = s.parent()._set
        else:
            S = s.base_set()
        if hasattr(t.parent(), "_set"):
            T = t.parent()._set
        else:
            T = t.base_set()
        if S != T:
            raise ValueError("cannot compare partitions of different sets")

        if s == t:
            return False

        for p in s:
            x = p[0]
            for t_ in t:
                if x in t_:
                    break
            for p_ in p:
                if p_ not in t_:
                    return False
        return True

    lt = is_less_than

    def is_strict_refinement(self, s, t):
        r"""
        Return ``True`` if ``s`` is a strict refinement of ``t`` and
        satisfies `s \neq t`.

        A set partition `s` is said to be a strict refinement of a set
        partition `t` of the same set if and only if one can obtain
        `t` from `s` by repeatedly combining pairs of parts whose
        convex hulls don't intersect (i. e., whenever we are combining
        two parts, the maximum of each of them should be smaller than
        the minimum of the other).

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1],[2],[3],[4]])
            sage: t = S([[1,3],[2,4]])
            sage: u = S([[1,2,3,4]])
            sage: S.is_strict_refinement(s, t)
            True
            sage: S.is_strict_refinement(t, u)
            False
            sage: A = SetPartition([[1,3],[2,4]])
            sage: B = SetPartition([[1,2,3,4]])
            sage: S.is_strict_refinement(s, A)
            True
            sage: S.is_strict_refinement(t, B)
            False
        """
        if hasattr(s.parent(), "_set"):
            S = s.parent()._set
        else:
            S = frozenset(s.base_set())
        if hasattr(t.parent(), "_set"):
            T = t.parent()._set
        else:
            T = frozenset(t.base_set())
        if S != T:
            raise ValueError("cannot compare partitions of different sets")

        if s == t:
            return False

        for p in t:
            L = [x for x in list(s) if x.issubset(p)]
            if sum(len(x) for x in L) != len(p) \
                    or any(max(L[i]) > min(L[i+1]) for i in range(len(L)-1)):
                return False
        return True

class SetPartitions_all(SetPartitions):
    r"""
    All set partitions.
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions()
            sage: TestSuite(S).run()
        """
        SetPartitions.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SetPartitions()
            Set partitions
        """
        return "Set partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = SetPartitions().__iter__()
            sage: [next(it) for x in range(10)]
            [{}, {{1}}, {{1, 2}}, {{1}, {2}}, {{1, 2, 3}}, {{1}, {2, 3}},
             {{1, 3}, {2}}, {{1, 2}, {3}}, {{1}, {2}, {3}}, {{1, 2, 3, 4}}]
        """
        n = 0
        while True:
            for x in SetPartitions_set(frozenset(range(1, n+1))):
                yield self.element_class(self, list(x))
            n += 1

class SetPartitions_set(SetPartitions):
    """
    Set partitions of a fixed set `S`.
    """
    @staticmethod
    def __classcall_private__(cls, s):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SetPartitions(set([2,1,4]))
            sage: S2 = SetPartitions([4,1,2])
            sage: S3 = SetPartitions((1,2,4))
            sage: S1 is S2, S1 is S3
            (True, True)
        """
        return super(SetPartitions_set, cls).__classcall__(cls, frozenset(s))

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions(3)
            sage: TestSuite(S).run()
            sage: SetPartitions(0).list()
            [{}]
            sage: SetPartitions([]).list()
            [{}]
        """
        self._set = s
        SetPartitions.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions([1,2,3])
            Set partitions of {1, 2, 3}
        """
        return "Set partitions of %s"%(Set(self._set))

    def __contains__(self, x):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: all(sp in S for sp in S)
            True
            sage: SetPartition([[1,3],[2,4]]) in SetPartitions(3)
            False
            sage: SetPartition([[1,3],[2,4]]) in SetPartitions(4, [3,1])
            False
            sage: SetPartition([[2],[1,3,4]]) in SetPartitions(4, [3,1])
            True
        """
        # Must pass the general check
        if not SetPartitions.__contains__(self, x):
            return False

        # Make sure that the number of elements match up
        if sum(map(len, x)) != len(self._set):
            return False

        # Make sure that the union of all the sets is the original set
        if reduce(lambda u, s: u.union(Set(s)), x, Set([])) != Set(self._set):
            return False

        return True

    def cardinality(self):
        """
        Return the number of set partitions of the set `S`.

        The cardinality is given by the `n`-th Bell number where `n` is the
        number of elements in the set `S`.

        EXAMPLES::

            sage: SetPartitions([1,2,3,4]).cardinality()
            15
            sage: SetPartitions(3).cardinality()
            5
            sage: SetPartitions(3,2).cardinality()
            3
            sage: SetPartitions([]).cardinality()
            1
        """
        return bell_number(len(self._set))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SetPartitions(3).list()
            [{{1, 2, 3}}, {{1}, {2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{1}, {2}, {3}}]
        """
        for p in Partitions(len(self._set)):
            for sp in self._iterator_part(p):
                yield self.element_class(self, sp)

    def base_set(self):
        """
        Return the base set of ``self``.

        EXAMPLES::

            sage: SetPartitions(3).base_set()
            {1, 2, 3}
        """
        return Set(self._set)

    def base_set_cardinality(self):
        """
        Return the cardinality of the base set of ``self``.

        EXAMPLES::

            sage: SetPartitions(3).base_set_cardinality()
            3
        """
        return len(self._set)

class SetPartitions_setparts(SetPartitions_set):
    r"""
    Class of all set partitions with fixed partition sizes corresponding to
    an integer partition `\lambda`.
    """
    @staticmethod
    def __classcall_private__(cls, s, parts):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SetPartitions(4, [2,2])
            sage: T = SetPartitions([1,2,3,4], Partition([2,2]))
            sage: S is T
            True
        """
        if isinstance(s, (int, Integer)):
            s = list(range(1, s + 1))
        return super(SetPartitions_setparts, cls).__classcall__(cls, frozenset(s), Partition(parts))

    def __init__(self, s, parts):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: TestSuite(S).run()
        """
        SetPartitions_set.__init__(self, s)
        self.parts = parts

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions(4, [2,2])
            Set partitions of {1, 2, 3, 4} with sizes in [2, 2]
        """
        return "Set partitions of %s with sizes in %s"%(Set(self._set), self.parts)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        This algorithm counts for each block of the partition the
        number of ways to fill it using values from the set.  Then,
        for each distinct value `v` of block size, we divide the result by
        the number of ways to arrange the blocks of size `v` in the
        set partition.

        For example, if we want to count the number of set partitions
        of size 13 having [3,3,3,2,2] as underlying partition we
        compute the number of ways to fill each block of the
        partition, which is `\binom{13}{3} \binom{10}{3} \binom{7}{3}
        \binom{4}{2}\binom{2}{2}` and as we have three blocks of size
        `3` and two blocks of size `2`, we divide the result by
        `3!2!` which gives us `600600`.

        EXAMPLES::

            sage: SetPartitions(3, [2,1]).cardinality()
            3
            sage: SetPartitions(13, Partition([3,3,3,2,2])).cardinality()
            600600

        TESTS::

            sage: all((len(SetPartitions(size, part)) == SetPartitions(size, part).cardinality() for size in range(8) for part in Partitions(size)))
            True
            sage: sum((SetPartitions(13, p).cardinality() for p in Partitions(13))) == SetPartitions(13).cardinality()
            True
        """
        from sage.misc.misc_c import prod

        remaining_subset_size = Integer(len(self._set))
        cardinal = Integer(1)
        for subset_size in self.parts:
            cardinal *= remaining_subset_size.binomial(subset_size)
            remaining_subset_size -= subset_size

        repetitions = (Integer(rep).factorial()
                       for rep in self.parts.to_exp_dict().values()
                       if rep != 1)
        cardinal /= prod(repetitions)
        return Integer(cardinal)

    def __iter__(self):
        """
        An iterator for all the set partitions of the given set with
        the given sizes.

        EXAMPLES::

            sage: SetPartitions(3, [2,1]).list()
            [{{1}, {2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}]
        """
        for sp in self._iterator_part(self.parts):
            yield self.element_class(self, sp)

    def __contains__(self, x):
        """
        Check containment.

        TESTS::

            sage: S = SetPartitions(4, [3,1])
            sage: Set([Set([1,2,3]), Set([4])]) in S
            True
            sage: Set([Set([1,3]), Set([2,4])]) in S
            False
            sage: Set([Set([1,2,3,4])]) in S
            False
        """
        if not SetPartitions_set.__contains__(self, x):
            return False
        return sorted(map(len, x)) == list(reversed(self.parts))

class SetPartitions_setn(SetPartitions_set):
    @staticmethod
    def __classcall_private__(cls, s, n):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SetPartitions(set([2,1,4]), 2)
            sage: S2 = SetPartitions([4,1,2], 2)
            sage: S3 = SetPartitions((1,2,4), 2)
            sage: S1 is S2, S1 is S3
            (True, True)
        """
        return super(SetPartitions_setn, cls).__classcall__(cls, frozenset(s), n)

    def __init__(self, s, n):
        """
        TESTS::

            sage: S = SetPartitions(5, 3)
            sage: TestSuite(S).run()
        """
        self.n = n
        SetPartitions_set.__init__(self, s)

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions(5, 3)
            Set partitions of {1, 2, 3, 4, 5} with 3 parts
        """
        return "Set partitions of %s with %s parts"%(Set(self._set), self.n)

    def cardinality(self):
        """
        The Stirling number of the second kind is the number of partitions
        of a set of size `n` into `k` blocks.

        EXAMPLES::

            sage: SetPartitions(5, 3).cardinality()
            25
            sage: stirling_number2(5,3)
            25
        """
        return stirling_number2(len(self._set), self.n)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SetPartitions(3).list()
            [{{1, 2, 3}}, {{1}, {2, 3}}, {{1, 3}, {2}}, {{1, 2}, {3}}, {{1}, {2}, {3}}]
        """
        for p in Partitions(len(self._set), length=self.n):
            for sp in self._iterator_part(p):
                yield self.element_class(self, sp)

    def __contains__(self, x):
        """
        Check containment.

        TESTS::

            sage: S = SetPartitions(4, 2)
            sage: Set([Set([1,2,3]), Set([4])]) in S
            True
            sage: Set([Set([1,3]), Set([2,4])]) in S
            True
            sage: Set([Set([1,2,3,4])]) in S
            False
        """
        if not SetPartitions_set.__contains__(self, x):
            return False
        return len(x) == self.n

def _listbloc(n, nbrepets, listint=None):
    r"""
    Decompose a set of `n \times n` ``brepets`` integers (the list
    ``listint``) in ``nbrepets`` parts.

    It is used in the algorithm to generate all set partitions.

    .. WARNING::

        Internal function that is not to be called by the user.

    EXAMPLES::

        sage: list(sage.combinat.set_partition._listbloc(2,1))
        [{{1, 2}}]
        sage: l = [Set([Set([3, 4]), Set([1, 2])]), Set([Set([2, 4]), Set([1, 3])]), Set([Set([2, 3]), Set([1, 4])])]
        sage: list(sage.combinat.set_partition._listbloc(2,2,[1,2,3,4])) == l
        True
    """
    if isinstance(listint, (int, Integer)) or listint is None:
        listint = Set(range(1,n+1))

    if nbrepets == 1:
        yield Set([listint])
        return

    l = sorted(listint)
    smallest = Set(l[:1])
    new_listint = Set(l[1:])

    f = lambda u, v: u.union(_set_union([smallest,v]))

    for ssens in subset.Subsets(new_listint, n-1):
        for z in _listbloc(n, nbrepets-1, new_listint-ssens):
            yield f(z,ssens)

def _union(s):
    """
    TESTS::

        sage: s = Set([ Set([1,2]), Set([3,4]) ])
        sage: sage.combinat.set_partition._union(s)
        {1, 2, 3, 4}
    """
    result = Set([])
    for ss in s:
        result = result.union(ss)
    return result

def _set_union(s):
    """
    TESTS::

        sage: s = Set([ Set([1,2]), Set([3,4]) ])
        sage: sage.combinat.set_partition._set_union(s)
        {{1, 2, 3, 4}}
    """
    result = Set([])
    for ss in s:
        result = result.union(ss)
    return Set([result])

def cyclic_permutations_of_set_partition(set_part):
    """
    Returns all combinations of cyclic permutations of each cell of the
    set partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.set_partition import cyclic_permutations_of_set_partition
        sage: cyclic_permutations_of_set_partition([[1,2,3,4],[5,6,7]])
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]
    """
    return list(cyclic_permutations_of_set_partition_iterator(set_part))

def cyclic_permutations_of_set_partition_iterator(set_part):
    """
    Iterates over all combinations of cyclic permutations of each cell
    of the set partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.set_partition import cyclic_permutations_of_set_partition_iterator
        sage: list(cyclic_permutations_of_set_partition_iterator([[1,2,3,4],[5,6,7]]))
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]
    """
    from sage.combinat.permutation import CyclicPermutations
    if len(set_part) == 1:
        for i in CyclicPermutations(set_part[0]):
            yield [i]
    else:
        for right in cyclic_permutations_of_set_partition_iterator(set_part[1:]):
            for perm in CyclicPermutations(set_part[0]):
                yield [perm] + right
