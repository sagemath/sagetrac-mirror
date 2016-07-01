r"""
The Directed Convex Polyominoes
===============================

The goal of this module is to give some tools to manipulate the
directed convex polyominoes.
"""
#******************************************************************************
#  Copyright (C) 2014 Adrien Boussicault (boussica@labri.fr),
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy
)
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.global_options import GlobalOptions
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_union_enumerated_sets \
    import DisjointUnionEnumeratedSets
from sage.rings.integer import Integer
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from copy import deepcopy
from sage.matrix.constructor import matrix
from sage.combinat.combinat import catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.partition import Partition
from sage.combinat.partition import Partitions
from sage.combinat.parallelogram_polyomino import ParallelogramPolyomino
from sage.combinat.parallelogram_polyomino import _drawing_tool
from sage.functions.other import binomial

default_tikz_options = dict(
    scale=1, line_size=1, point_size=3.5,
    color_line='black', color_point='black',
    color_bounce_0='red', color_bounce_1='blue',
    translation=[0, 0], rotation=0,
    mirror=None
)

DirectedConvexPolyominoesOptions = GlobalOptions(
    name='directed convex polyominoes',
    doc=r"""
    """,
    end_doc=r"""
    """,
    tikz_options=dict(
        default=default_tikz_options,
        description='the tikz options',
        checker=lambda x: Set(x.keys()).issubset(
            Set(
                [
                    'scale', 'line_size', 'point_size',
                    'color_line', 'color_point', 'translation', 'mirror',
                    'rotation', 'color_bounce_0', 'color_bounce_1'
                ]
            )
        )
    ),
    drawing_components=dict(
        default=dict(diagram=True),
        description='Different tree-like tableaux components to draw',
        checker=lambda x: Set(x.keys()).issubset(
            Set(['diagram', 'bounce_0', 'bounce_1'])
        )
    ),
    display=dict(
        default="list",
        values=dict(
            list='displayed as list',
            drawing='as a drawing'
        )
    ),
    latex=dict(
        default="drawing",
        values=dict(
            list='displayed as list',
            drawing='as a drawing'
        )
    )
)


def _maximal_cut(pp):
    if pp.size() == 0:
        return [0, 0]
    h = pp.heights()[-1] - 1
    w = pp.widths()[-1] - 1
    return [h, w]


def _maximal_partition_cut(pp):
    [h, w] = _maximal_cut(pp)
    return [w for i in range(h)]


def _include(pp1, pp2):
    if len(pp1) > len(pp2):
        return False
    for i in range(len(pp1)):
        if pp1[i] > pp2[i]:
            return False
    return True


class DirectedConvexPolyomino(ClonableList):
    r"""
    The class of directed convex polyominoes.
    """
    __metaclass__ = InheritComparisonClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        """
        return cls._auto_parent._element_constructor_(*args, **opts)

    @lazy_class_attribute
    def _auto_parent(cls):
        r"""
        """
        return DirectedConvexPolyominoes()

    def check(self):
        r"""
        This method raises an error if the internal data of the class does not
        represent a parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
        """
        pp = self.parallelogram_polyomino()
        cut = self.cut()
        maximal_cut = _maximal_partition_cut(pp)
        if not _include(cut, maximal_cut):
            raise ValueError(
                "The cut must be included in %s." % (str(maximal_cut))
            )

    def __hash__(self):
        r"""
        Return the hash code of the parallelogram polyomino

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: hash(dcp) == hash((pp, cut))
            True

            sage: DCPS = DirectedConvexPolyominoes(7)
            sage: D = { DCPS[0]: True, DCPS[1]: True }
            sage: D[DCPS[0]] = False
            sage: import pprint
            sage: pp = pprint.PrettyPrinter()
            sage: pp.pprint(D)
            {[[[0, 0, 0, 0, 0, 0, 0, 1], [1, 0, 0, 0, 0, 0, 0, 0]], []]: False,
            [[[0, 0, 0, 0, 0, 0, 1, 1], [1, 0, 0, 0, 0, 0, 1, 0]], []]: True}
        """
        return hash((self.parallelogram_polyomino(), self.cut()))

    def __copy__(self):
        r"""
        Copy a directed convex polyomino

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp1 = copy(dcp)
            sage: dcp1 is dcp
            False
            sage: dcp1 == dcp
            True
            sage: dcp1
            [[[0, 0, 1, 0, 1, 0, 1, 1, 1, 1], [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]], [3, 2]]
        """
        return DirectedConvexPolyomino(
            [self.parallelogram_polyomino(), self.cut()]
        )

    def parallelogram_polyomino(self):
        r"""
        Return the minimal parallelogram polyomino that contain the directed
        convex polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: pp == dcp.parallelogram_polyomino()
            True
        """
        return ClonableList.__getitem__(self, 0)

    def pp(self):
        r"""
        Return the minimal parallelogram polyomino that contain the directed
        convex polyomino.

        It is a shortcut for parallelogram_polyomino().

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: pp == dcp.pp()
            True
        """
        return self.parallelogram_polyomino()

    def cut(self):
        """
        Return the partition to remove form the parallelogram polyomino of
        self.parallelogram_polyomino() to obtain the directed convex polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: cut == dcp.cut()
            True
        """
        return ClonableList.__getitem__(self, 1)

    def __init__(self, parent, value, check=True):
        r"""
        Construct a directed convex polyomino.

        The input is a pair of a polyomino parallelogram and a partition.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp
            [[[0, 0, 1, 0, 1, 0, 1, 1, 1, 1], [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]], [3, 2]]

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp
            [[[0, 0, 1, 0, 1, 0, 1, 1, 1, 1], [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]], []]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp
            [[[1], [1]], []]
        """
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise ValueError(
                "Invalid input for directed convex polyomino: %s." % (
                    str(value)
                )
            )
        [pp, cut] = value
        if isinstance(pp, (list, tuple)):
            pp = ParallelogramPolyomino(pp)
        if isinstance(cut, (list, tuple)):
            cut = Partition(cut)
        ClonableList.__init__(self, parent, [pp, cut])
        if check:
            if not isinstance(value, (list, tuple)):
                raise ValueError(
                    "Value %s must be a list or a tuple." % (value)
                )
            self.check()
        self._options = None

    def degree_convexity(self):
        r"""
        Return the degree convexity of a directed convex polyomino.

        A convex polyomino is said to be k-convex if every pair of its cells
        can be connected by a monotone path (path with south and east steps)
        with at most k changes of direction.
        The degree of convexity of a convex polyomino P is the smallest integer
        k such that P is k-convex.

        If the directed convex polyomino is empty, the function return -1.

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1] ,
            ....:         [1, 1, 0, 1, 1, 0, 1, 0, 0]  ]
            ....:
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp.degree_convexity()
            3

            sage: cut = Partition([1])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp.degree_convexity()
            2

        """
        pp = self.parallelogram_polyomino()
        k_pp = pp.degree_convexity()
        if not pp.is_flat():
            return k_pp
        max_degree_cut = self.maximal_degree_cut()
        cut = self.cut()
        if _include(max_degree_cut, cut):
            return k_pp - 1
        else:
            return k_pp

    def is_k_directed(self, k):
        r"""
        Return true if the directed convex polyomino is k-directed.

        A convex polyomino is said to be k-convex if every pair of its cells
        can be connected by a monotone path (path with south and east steps)
        with at most k changes of direction.
        The degree of convexity of a convex polyomino P is the smallest integer
        k such that P is k-convex.

        """
        return self.degree_convexity() <= k

    @cached_method
    def get_array(self):
        r"""
        Return an array of 0s and 1s such that the 1s represent the boxes of
        the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: matrix(dcp.get_array())
            [1 1 0 0 0 0]
            [1 1 1 1 1 1]
            [0 1 1 1 0 0]
            [0 0 1 0 0 0]

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: matrix(dcp.get_array())
            [1 1 0 0 0 0]
            [1 1 1 1 1 1]
            [0 1 1 1 1 1]
            [0 0 1 1 1 1]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: matrix(dcp.get_array())
            []

        """
        array = deepcopy(self.parallelogram_polyomino().get_array())
        cut = self.cut()
        for h in range(len(cut)):
            for w in range(cut[h]):
                array[-1-h][-1-w] = 0
        return array

    def nb_corners(self):
        array = self.get_array()
        def entry( array, i, j ):
            if (
                i<0 or j<0 or
                i >= len(array) or
                len(array) == 0 or j >= len(array[0])
            ):
                return 0
            return array[i][j]
        res = 0
        for h in range( len(array) ):
            for w in range( len(array[0]) ):
                if array[h][w] == 1 :
                    if (entry(array,h,w-1) == 0 and entry(array,h-1,w) == 0):
                        res += 1
                    if (entry(array,h,w-1) == 0 and entry(array,h+1,w) == 0):
                        res += 1
                    if (entry(array,h,w+1) == 0 and entry(array,h-1,w) == 0):
                        res += 1
                    if (entry(array,h,w+1) == 0 and entry(array,h+1,w) == 0):
                        res += 1
        return res

    def maximal_cut(self):
        r"""
        """
        pp = self.parallelogram_polyomino()
        return Partition(_maximal_partition_cut(pp))

    def is_flat(self):
        """
        Return true if the Parallelogram polyomino associated with the directed
        convex polyomino is flat.
        """
        return self.parallelogram_polyomino().is_flat()

    def maximal_degree_cut(self):
        r"""
        If the parallelogram polyomino is flat, return the size of the maximal
        rectangle included in the parallelogram polyomino where cells have a
        maximal degree convexity.
        It the parallelogram polyomino is not flat, return the empty partition.

        A cell of convex polyomino has k as degree convexity if it can be
        connected with the leftmost cell of the top row by a monotone path
        (path with south and east steps) with at most k changes of direction.

        RETURNS: The return value is a list [h, w] where h is the height of the
        parallelogram polyomino and w is the width.

        EXAMPLES::

            sage: dcp = DirectedConvexPolyomino(
            ....:     [
            ....:         [
            ....:             [0, 0, 0, 1, 0, 1, 0, 1, 1],
            ....:             [1, 1, 0, 1, 0, 0, 1, 0, 0]
            ....:         ], []
            ....:     ]
            ....: )
            sage: dcp.maximal_degree_cut()
            [1]

            sage: dcp = DirectedConvexPolyomino(
            ....:     [[[0, 0, 0, 0, 1, 1, 1], [1, 1, 1, 0, 0, 0, 0]], []]
            ....: )
            sage: dcp.maximal_degree_cut()
            [2, 2, 2]

            sage: dcp = DirectedConvexPolyomino(
            ....:     [
            ....:         [
            ....:             [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1],
            ....:             [1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0]
            ....:         ], []
            ....:     ]
            ....: )
            sage: dcp.maximal_degree_cut()
            [3, 3]

            sage: dcp = DirectedConvexPolyomino(
            ....:     [
            ....:         [
            ....:              [0, 0, 1, 0, 1, 0, 0, 1, 1],
            ....:              [1, 1, 0, 1, 0, 1, 0, 0, 0]
            ....:         ], []
            ....:     ]
            ....: )
            sage: dcp.maximal_degree_cut()
            [1, 1]

            sage: dcp = DirectedConvexPolyomino(
            ....:     [[[1], [1]], []]
            ....: )
            sage: dcp.maximal_degree_cut()
            []

        TESTS:

            sage: dcp = DirectedConvexPolyomino(
            ....:     [
            ....:         [
            ....:             [0, 0, 1, 0, 1, 0, 1, 1],
            ....:             [1, 1, 0, 1, 0, 1, 0, 0]
            ....:         ], []
            ....:     ]
            ....: )
            sage: dcp.maximal_degree_cut()
            [1]

            sage: dcp = DirectedConvexPolyomino(
            ....:     [[[0, 1], [1, 0]], []]
            ....: )
            sage: dcp.maximal_degree_cut()
            []

        """
        pp = self.parallelogram_polyomino()
        if pp.size() == 0 or pp.size() == 1:
            return Partition([])
        if not pp.is_flat():
            return Partition([])
        k = pp.degree_convexity()
        if k % 2 == 0:
            direction = 0
        else:
            direction = 1
        h = pp.bounce_path(direction=direction)[-1]
        w = pp.bounce_path(direction=1-direction)[-1]
        return Partition([w for i in range(h)])

    def bounce_path(self, direction=1):
        r"""
        Return the bounce path of directed convex polyomino.

        The bounce path is a path with two steps (1, 0) and (0, 1).

        If 'direction' is 1 (resp. 0), the bounce path is the path
        starting at position position (h=1, w=0) (resp. (h=0, w=1)) with
        initial direction, the vector (0, 1) (resp. (1, 0)), and turning
        each time the path crosses the perimeter of the parallelogram
        polyomino.

        The path is coded by a list of integers. Each integer represents
        the size of the path between two turnings.

        You can visualize the two bounce paths by using the following
        commands:

        """
        raise NotImplementedError()
        return None

    def bounce(self, direction=1):
        r"""
        Return the bounce of the parallelogram polyomino.

        Les p be the bounce path of the parallelogram
        polyomino. (p=self.bounce_path())
        The bounce is defined by:
        sum([(1+floor(i/2))*p[i] for i in range(len(p))])

        """
        raise NotImplementedError()
        return None

    def area(self):
        r"""
        Returns the area of the directed convex polyomino.

        EXAMPLES::

        """
        area = 0
        for line in self.get_array():
            for e in line:
                area += e
        return area

    def _repr_(self):
        r"""
        Return a string representation of the directed convex polyomino.

        EXAMPLES::

        """
        return self.get_options().dispatch(self, '_repr_', 'display')

    def _repr_list(self):
        r"""
        Return a string representation with list style.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp._repr_list()
            '[[[0, 0, 1, 0, 1, 0, 1, 1, 1, 1], [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]], [3, 2]]'
        """
        return ClonableList._repr_(self)

    def geometry(self):
        r"""
        Return a pair [h, w] containing the height and the width of the
        directed convex polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp.geometry()
            [4, 6]

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp.geometry()
            [4, 6]

            sage: pp = ParallelogramPolyomino(
            ....:     [[1], [1]]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp.geometry()
            [0, 1]

        """
        pp = self.parallelogram_polyomino()
        return [pp.height(), pp.width()]

    def size(self):
        r"""
        Return the size of the parallelogram polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp.size()
            9
        """
        return self.parallelogram_polyomino().size()

    def get_options(self):
        r"""
        Return all the options of the object.

        EXAMPLES::

            sage: pp = DirectedConvexPolyomino([[[0, 1], [1, 0]], []])
            sage: pp.get_options()
            options for directed convex polyominoes
        """
        if self._options is None:
            return self.parent().get_options()
        return self._options

    def set_options(self, *get_value, **set_value):
        r"""
        Set new options to the object.

        EXAMPLES::
            TODO
        """
        if self._options is None:
            self._options = deepcopy(self.get_options())
        self._options(*get_value, **set_value)

    def _repr_drawing(self):
        r"""
        Return a string representing a drawing of the directed convex
        polyomino.

        EXAMPLES::

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp._repr_drawing()
            '[1 1 0 0 0 0]\n[1 1 1 1 1 1]\n[0 1 1 1 0 0]\n[0 0 1 0 0 0]'
        """
        return str(matrix(self.get_array()))

    def get_tikz_options(self):
        return self.get_options()['tikz_options']

    def width(self):
        return self.parallelogram_polyomino().width()

    @cached_method
    def widths(self):
        res = []
        for line in self.get_array():
            res.append(sum(line))
        return res

    def height(self):
        return self.parallelogram_polyomino().height()

    @cached_method
    def heights(self):
        width = self.width()
        res = [0 for i in range(width)]
        array = self.get_array()
        for h in range(len(array)):
            for w in range(width):
                res[w] += array[h][w]
        return res

    def _to_tikz_diagram(self):
        tikz_options = self.get_tikz_options()
        grid_width = self.width() + 1
        grid_height = self.height() + 1
        drawing_tool = _drawing_tool(
            tikz_options,
            XY=lambda v: [v[0], grid_height-1-v[1]]
        )
        res = ""
        if self.size() == 0:
            res += drawing_tool.draw_line([0, 0], [1, 0])
            return res

        def carre(x, y):
            res = drawing_tool.draw_polyline(
                [[x, y], [x+1, y], [x+1, y+1], [x, y+1], [x, y]]
            )
            return res
        array = self.get_array()
        for h in range(len(array)):
            for w in range(len(array[h])):
                if array[h][w] == 1:
                    res += carre(w, h)
        return res

    def _to_tikz_bounce(self, directions=[0, 1]):
        pp = self.parallelogram_polyomino()
        pp.get_options()['tikz_options'] = deepcopy(
            self.get_tikz_options()
        )
        t = pp.get_options()['tikz_options']['translation']
        t[0] += -.5
        t[1] += .5
        return pp._to_tikz_bounce(directions=directions)

    def to_tikz(self):
        r"""
        Return the tikz code of the directed convex polyomino.

        This code is the code present inside a tikz latex environment.
        """
        res = ""
        drawing_components = self.get_options()['drawing_components']
        if 'diagram' in drawing_components:
            res += self._to_tikz_diagram()
        directions = []
        if 'bounce_0' in drawing_components:
            directions.append(0)
        if 'bounce_1' in drawing_components:
            directions.append(1)
        if len(directions) != 0:
            res += self._to_tikz_bounce(directions)
        return res

    def _latex_(self):
        r"""
        Return a LaTeX version of ``self``.

        For more on the latex options, see
        :meth:`DirectedConvexPolyominoes.global_options`.
        """
        return self.get_options().dispatch(self, '_latex_', 'latex')

    def _latex_drawing(self):
        r"""
        Return a LaTeX version of ``self`` in a drawing style.
        """
        latex.add_package_to_preamble_if_available("tikz")
        tikz_options = self.get_tikz_options()
        res = "\n\\begin{tikzpicture}[scale=%s]" % (tikz_options['scale'])
        res += self.to_tikz()
        res += "\n\\end{tikzpicture}"
        return res

    def _latex_list(self):
        r"""
        Return a LaTeX version of ``self`` in a list style.
        """
        return "\\[%s\\]" % self._repr_list()


class DirectedConvexPolyominoesFactory(SetFactory):
    r"""
    The directed convex polyominoes factory.
    """
    def __call__(self, size=None, policy=None):
        r"""
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return DirectedConvexPolyominoes_size(size, policy)
        if size is None:
            return DirectedConvexPolyominoes_all(policy)
        raise ValueError("Invalid argument for directed convex polyominoers"
                         "Factory.")

    def add_constraints(self, cons, args_opts):
        r"""
        This function permit to add some enumeration constraint to the 
        factory. The factory make a family using the given constraints.

        :meth:`SetFactory.add_constraints<.set_factories.SetFactory.add_constraints>`.
        """
        args, opts = args_opts
        return cons + args

    @lazy_attribute
    def _default_policy(self):
        return TopMostParentPolicy(self, (), DirectedConvexPolyomino)

    def _repr_(self):
        """
        """
        return "Factory for directed convex polyominoes"

DirectedConvexPolyominoes = DirectedConvexPolyominoesFactory()
DirectedConvexPolyominoes.__doc__ = \
    DirectedConvexPolyominoesFactory.__call__.__doc__


class DirectedConvexPolyominoes_size(
    ParentWithSetFactory, UniqueRepresentation
):
    r"""
    The directed convex polyominoes of size `n`.
    """
    def __init__(self, size, policy):
        r"""
        Construct a set of directed convex polyominoes of a given size.
        """
        self._size = size
        ParentWithSetFactory.__init__(
            self, (size, ), policy, category=FiniteEnumeratedSets()
        )

    def _repr_(self):
        r"""
        Return the string representation of the set of directed convex
        polyominoes

        EXAMPLES::

            sage: DirectedConvexPolyominoes(3)
            directed convex polyominoes of size 3
        """
        return "directed convex polyominoes of size %s" % (self._size)

    def an_element(self):
        r"""
        """
        return next(self.__iter__())

    def check_element(self, el, check):
        r"""
        """
        if el.size() != self.size():
            raise ValueError(
                "The direct convex polyomino has a Wrong size: %s" % (
                    el.size()
                )
            )

    def cardinality(self):
        r"""
        Return the number of directed convex polyominoes.

        EXAMPLES::
            sage: all([
            ....:     DirectedConvexPolyominoes(i).cardinality() ==
            ....:     len(list(DirectedConvexPolyominoes(i)))
            ....:     for i in range(1, 5)
            ....: ])
            True

        """
        n = self.size() - 1
        return binomial(2*n, n)

    def __iter__(self):
        r"""
        Return a directed convex polyomino generator.
        """
        from sage.combinat.parallelogram_polyomino \
            import ParallelogramPolyominoes
        from sage.combinat.partition import Partitions
        for pp in ParallelogramPolyominoes(self.size()):
            [h, w] = _maximal_cut(pp)
            max_cut = [w for i in range(h)]
            for size in range(h*w + 1):
                for partition in Partitions(size, outer=list(max_cut)):
                    yield DirectedConvexPolyomino([pp, partition])

    def get_options(self):
        return self.global_options

    def set_options(self, *get_value, **set_value):
        self.global_options(*get_value, **set_value)

    global_options = DirectedConvexPolyominoesOptions

    def size(self):
        r"""
        Return the size of the convex directed polyominoes generated by this
        parent.

        EXAMPLES::

            sage: DirectedConvexPolyominoes(0).size()
            0
            sage: DirectedConvexPolyominoes(1).size()
            1
            sage: DirectedConvexPolyominoes(5).size()
            5
        """
        return self._size


class DirectedConvexPolyominoes_all(
    ParentWithSetFactory, DisjointUnionEnumeratedSets
):
    r"""
    This class enumerates all the directed convex polyominoes.
    """
    def __init__(self, policy):
        r"""
        Construct the set of all directed convex polyominoes.

        EXAMPLES::

            sage: DCPS = DirectedConvexPolyominoes()
            sage: DCPS
            directed convex polyominoes

            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp in DCPS
            True

            sage: next(DCPS.__iter__()) in DCPS
            True
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category=FiniteEnumeratedSets()
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                NonNegativeIntegers(), self._directed_convex_polyominoes_size
            ),
            facade=True, keepkey=False,
            category=self.category()
        )

    def _directed_convex_polyominoes_size(self, n):
        return DirectedConvexPolyominoes_size(n, policy=self.facade_policy())

    def _repr_(self):
        r"""
        Return a string representation of the set of directed convex
        polyominoes.

        EXAMPLES::

            sage: DCPS = DirectedConvexPolyominoes()
            sage: DCPS
            directed convex polyominoes
        """
        return "directed convex polyominoes"

    def check_element(self, el, check):
        r"""
        Check is a given element `el` is in the set of directed convex
        polyominoes.

        EXAMPLES::

            sage: DCPS = DirectedConvexPolyominoes()
            sage: pp = ParallelogramPolyomino(
            ....:     [
            ....:         [0, 0, 1, 0, 1, 0, 1, 1, 1, 1],
            ....:         [1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
            ....:     ]
            ....: )
            sage: cut = Partition([3, 2])
            sage: dcp = DirectedConvexPolyomino([pp, cut])
            sage: dcp in DCPS
            True
        """
        pass

    def get_options(self):
        r"""
        Return all the options associated with the set of
        directed convex polyomino.

        EXAMPLES::

            sage: DPS = DirectedConvexPolyominoes()
            sage: options = DPS.get_options()
            sage: options
            options for directed convex polyominoes
            sage: options()
            Current options for directed convex polyominoes
              - display:            list
            ...
        """
        return self.global_options

    def set_options(self, *get_value, **set_value):
        self.global_options(*get_value, **set_value)

    global_options = DirectedConvexPolyominoesOptions
