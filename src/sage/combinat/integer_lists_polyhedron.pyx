r"""
Generate lists of integers using polyhedra

AUTHOR:

- Jeroen Demeyer (2015-03-14) implement integer lists using polyhedra,
  see :trac:`17920`
"""

#*****************************************************************************
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include 'sage/ext/interrupt.pxi'
from sage.misc.cachefunc import cached_method
from sage.rings.integer cimport Integer, smallInteger
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent cimport Parent
from sage.structure.list_clone cimport ClonableArray
from sage.rings.infinity import infinity
from sage.misc.lazy_import import LazyImport
Polyhedron = LazyImport('sage.geometry.polyhedron.constructor', 'Polyhedron')


def Polyhedron_inf(ieqs, **kwds):
    """
    Given a list of inequalities, return a :class:`Polyhedron`
    determined by the inequalities.

    Unlike the usual :class:`Polyhedron` constructor, we allow
    infinities as constant terms of the inequalities.

    INPUT:

    - ``ieqs`` -- a list of inequalities for the :class:`Polyhedron`
      constructor.

    - ``**kwds`` -- other arguments passed to :class:`Polyhedron`

    EXAMPLES::

        sage: from sage.combinat.integer_lists_polyhedron import Polyhedron_inf

    We construct the 2-dimensional triangle defined by
    `x >= 0`, `y >= 0` and `x + y <= 1`::

        sage: ieqs = [(0,1,0), (0,0,1), (1,-1,-1)]
        sage: P = Polyhedron_inf(ieqs=ieqs)
        sage: P
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
        sage: P.inequalities()
        (An inequality (-1, -1) x + 1 >= 0,
         An inequality (1, 0) x + 0 >= 0,
         An inequality (0, 1) x + 0 >= 0)

    We add first a trivial inequality, then an impossible inequality::

        sage: ieqs.append((infinity,0,0))
        sage: Polyhedron_inf(ieqs=ieqs)
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices
        sage: ieqs.append((-infinity,0,0))
        sage: Polyhedron_inf(ieqs=ieqs)
        The empty polyhedron in ZZ^2
        sage: Polyhedron_inf(ieqs=ieqs, base_ring=RDF)
        The empty polyhedron in RDF^2

    Finally, we use only the impossible inequality::

        sage: Polyhedron_inf(ieqs=[(-infinity,0,0)], base_ring=RDF)
        The empty polyhedron in RDF^2
    """
    if 'ambient_dim' not in kwds:
        try:
            kwds['ambient_dim'] = len(ieqs[0])-1
        except Exception:
            pass

    # Loop over given inequalities and construct new, checking for
    # infinities.
    cdef list newieqs = []
    for ieq in ieqs:
        c = ieq[0]
        if c == -infinity:
            # Insatisfiable => empty polyhedron
            return Polyhedron(**kwds)
        elif c == infinity:
            # Trivial inequality => do not add it
            pass
        else:
            newieqs.append(ieq)

    return Polyhedron(ieqs=newieqs, **kwds)


def _function_iter(f):
    """
    A generator returning `f(0)`, `f(1)`, ...

    EXAMPLES::

        sage: from sage.combinat.integer_lists_polyhedron import _function_iter
        sage: it = _function_iter(lambda i: i)
        sage: for _ in range(5):
        ....:     print next(it)
        0
        1
        2
        3
        4
    """
    cdef Integer i = smallInteger(0)
    cdef Integer one = smallInteger(1)
    while True:
        yield f(i)
        i = i._add_(one)


def _minimal_rectangle_sum(part, min_length, max_length):
    """
    Given a number of parts, each equal to ``part`` with length
    bounded by ``[min_length, max_length]``, what is the minimal
    possible sum?

    EXAMPLES::

        sage: from sage.combinat.integer_lists_polyhedron import _minimal_rectangle_sum
        sage: _minimal_rectangle_sum(2, 50, 60)
        100
        sage: _minimal_rectangle_sum(-2, 50, 60)
        -120
        sage: _minimal_rectangle_sum(0, 0, infinity)
        0
        sage: _minimal_rectangle_sum(1, 0, infinity)
        0
        sage: _minimal_rectangle_sum(-1, 0, infinity)
        -Infinity
        sage: _minimal_rectangle_sum(infinity, 0, infinity)
        0
        sage: _minimal_rectangle_sum(infinity, 0, 0)
        0
    """
    assert min_length <= max_length
    if not part or max_length <= 0:
        return smallInteger(0)
    elif part < 0:
        return part * max_length
    elif min_length <= 0:
        return smallInteger(0)
    else:
        return part * min_length


def _maximal_rectangle_sum(part, min_length, max_length):
    """
    Given a number of parts, each equal to ``part`` with length
    bounded by ``[min_length, max_length]``, what is the maximal
    possible sum?

    EXAMPLES::

        sage: from sage.combinat.integer_lists_polyhedron import _maximal_rectangle_sum
        sage: _maximal_rectangle_sum(2, 50, 60)
        120
        sage: _maximal_rectangle_sum(-2, 50, 60)
        -100
        sage: _maximal_rectangle_sum(0, 0, infinity)
        0
        sage: _maximal_rectangle_sum(1, 0, infinity)
        +Infinity
        sage: _maximal_rectangle_sum(-1, 0, infinity)
        0
        sage: _maximal_rectangle_sum(infinity, 0, infinity)
        +Infinity
        sage: _maximal_rectangle_sum(infinity, 0, 0)
        0
    """
    assert min_length <= max_length
    if not part or max_length <= 0:
        return smallInteger(0)
    if part > 0:
        return part * max_length
    elif min_length <= 0:
        return smallInteger(0)
    else:
        return part * min_length


cdef inline Integer len_or_0(x):
    try:
        return smallInteger(len(x))
    except Exception:
        return smallInteger(0)


cdef inline signed_infinity(s):
    if not s:
        return smallInteger(0)
    if s > 0:
        return infinity
    else:
        return -infinity


cdef class IntegerList_polyhedron(ClonableArray):
    """
    Element class for :class:`IntegerLists_polyhedron`.
    """
    cpdef check(self):
        """
        Check whether ``self`` is really an element of its parent.

        EXAMPLES::

            sage: C = IntegerLists(max_part=3, min_length=4)
            sage: p = C([1,2,1,1])
            sage: p.parent()
            Integer lists of sum at least 0 satisfying certain constraints
            sage: p.check()
            True
            sage: C([1,2,1,4]).check()
            False
            sage: C([1,2]).check()
            False
        """
        return self in self._parent


class IntegerLists_polyhedron(Parent):
    r"""
    A combinatorial class `C` for integer lists satisfying certain
    sum, length, upper/lower bound and regularity constraints. The
    purpose of this tool is mostly to provide a Constant Amortized
    Time iterator through those lists, in lexicographic order.

    INPUT:

    In the list of input arguments, all arguments starting with
    ``min_`` or ``max_`` take either an integer or a signed infinity.
    Unless otherwise specified, the defaults are always "no condition".

    - ``n`` -- (default: ``None``) the requested sum of parts.
      If ``None``, use ``min_sum`` and ``max_sum`` instead.

    - ``min_sum`` -- minimal sum of parts, only used if ``n is None``

    - ``max_sum`` -- maximal sum of parts, only used if ``n is None``

    - ``min_length`` -- (default: 0) a lower bound for the length,
      must be non-negative

    - ``max_length`` -- an upper bound for the length

    - ``length`` -- overrides both ``min_length`` and ``max_length`` if
      specified

    - ``min_part`` -- (default: 0) a lower bound for all parts

    - ``floor`` -- a function or list to be used for additional lower
      bounds on the parts

    - ``ceiling`` -- a function or list to be used for additional upper
      bounds on the parts

    - ``min_slope`` -- the minimal value of the difference between two
      parts

    - ``max_slope`` -- the maximal value of the difference between two
      parts

    - ``min_part_last`` -- (default: 1) am additional lower bound for
      the last part in the list, only used for lists of length *longer*
      than ``min_length``.

    - ``max_part_last`` -- an additional upper bound for the last part
      in the list, only used for lists of length *longer* than
      ``min_length``.

    An *integer list* is a list `l` of integers, its *parts*. The
    *length* of `l` is the number of its parts; the *sum* of `l` is the
    sum of its parts.

    The constraints on the lists are as follow:

    - Sum: `sum(l) == n`

    - Length: ``min_length <= len(l) <= max_length``

    - Lower and upper bounds: ``floor(i) <= l[i] <= ceiling(i)``, for
      ``i`` from 0 to ``len(l)-1``. If ``floor`` is a list or iterator,
      it should give the values ``floor(0)``, ``floor(1)``, ... where
      the end of the iterator is interpreted as "no more conditions".

    - Slope condition: ``min_slope <= l[i+1] - l[i] <= max_slope``,
      for ``i`` from 0 to ``len(l)-2``

    .. NOTE::

        The arguments ``min_part``, ``floor`` and ``min_part_last`` all
        give lower bounds. All these lower bounds are considered
        together to compute the actual lower bounds. Analogously for
        upper bounds.

    .. NOTE::

        You should not make any assumptions on the *order* of the
        returned list. If you want a lexigraphic ordering, use
        :class:`IntegerListsLex` instead or manually sort the returned
        list.

    .. NOTE::

        The purpose of ``min_part_last`` is to consider two valid
        integer lists equivalent if they only differ by trailing zeros
        (or negative numbers). With the default value
        ``min_part_last=1``, only the list with the least number of
        trailing non-positive numbers will be produced. If you do not
        want this, set ``min_part_last=-infinity``.

    EXAMPLES:

    We create the combinatorial class of lists of length 3 and sum 2::

        sage: C = IntegerLists(2, length=3)
        sage: C
        Integer lists of sum 2 satisfying certain constraints
        sage: C.cardinality()
        6
        sage: [p for p in C]
        [[0, 0, 2], [0, 1, 1], [0, 2, 0], [1, 0, 1], [1, 1, 0], [2, 0, 0]]
        sage: [2, 0, 0] in C
        True
        sage: [2, 0, 1] in C
        False
        sage: "a" in C
        False
        sage: ["a"] in C
        False
        sage: C.first()
        [0, 0, 2]

    One can specify lower and upper bound on each part::

        sage: list(IntegerLists(5, length=3, floor=[1,2,0], ceiling=[3,2,3]))
        [[1, 2, 2], [2, 2, 1], [3, 2, 0]]

    Using the slope condition, one can generate integer partitions
    (but see :mod:`sage.combinat.partition.Partitions`)::

        sage: list(IntegerLists(4, max_slope=0))
        [[4], [2, 2], [3, 1], [2, 1, 1], [1, 1, 1, 1]]

    This is the list of all partitions of `7` with parts at least `2`::

        sage: list(IntegerLists(7, max_slope=0, min_part=2))
        [[7], [4, 3], [5, 2], [3, 2, 2]]

    This is the list of all partitions of `5` and length at most 3
    which are bounded below by [2,1,1]::

        sage: list(IntegerLists(5, max_slope=0, max_length=3, floor=[2,1,1]))
        [[5], [3, 2], [4, 1], [2, 2, 1], [3, 1, 1]]

    Note that ``[5]`` is considered valid, because the lower bound
    constraint only apply to existing positions in the list. To
    obtain instead the partitions containing ``[2,1,1]``, one need to
    use ``min_length``::

        sage: list(IntegerLists(5, max_slope=0, min_length=3, max_length=3, floor=[2,1,1]))
        [[2, 2, 1], [3, 1, 1]]

    This is the list of all partitions of `5` which are contained in
    ``[3,2,2]``::

        sage: list(IntegerLists(5, max_slope=0, max_length=3, ceiling=[3,2,2]))
        [[3, 2], [3, 1, 1], [2, 2, 1]]

    This is the list of all compositions of `4` (but see Compositions)::

        sage: list(IntegerLists(4, min_part=1))
        [[4], [1, 3], [2, 2], [3, 1], [1, 1, 2], [1, 2, 1], [2, 1, 1], [1, 1, 1, 1]]

    This is the list of all integer vectors of sum `4` and length `3`::

        sage: list(IntegerLists(4, length=3))
        [[0, 0, 4],
         [0, 1, 3],
         [0, 2, 2],
         [0, 3, 1],
         [0, 4, 0],
         [1, 0, 3],
         [1, 1, 2],
         [1, 2, 1],
         [1, 3, 0],
         [2, 0, 2],
         [2, 1, 1],
         [2, 2, 0],
         [3, 0, 1],
         [3, 1, 0],
         [4, 0, 0]]

    There are all the lists of sum 4 and length 4 such that l[i] <= i::

        sage: list(IntegerLists(4, length=4, ceiling=lambda i: i))
        [[0, 0, 1, 3], [0, 0, 2, 2], [0, 1, 0, 3], [0, 1, 1, 2], [0, 1, 2, 1]]

    This is the list of all monomials of degree `4` which divide the
    monomial `x^3y^1z^2` (a monomial being identified with its
    exponent vector)::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ....:     return x^exponents[0] * y^exponents[1] * z^exponents[2]
        sage: list( IntegerLists(4, length = len(m), ceiling = m, element_constructor = term) )
        [x^2*z^2, x^3*z, x*y*z^2, x^2*y*z, x^3*y]

    Note the use of the element_constructor feature.

    In the following example, the algorithm breaks down because it
    cannot find any lists and it cannot know for sure that the ceiling
    function always returns zero::

        sage: list(IntegerLists(1, ceiling=lambda i: 0))
        Traceback (most recent call last):
        ...
        RuntimeError: no more lists found, but cannot prove that there are none of length > 9

    Adding a length bound works::

        sage: list(IntegerLists(1, ceiling=lambda i: 0, max_length=30))
        []

    The ``max_part`` argument is more powerful, this works::

        sage: list(IntegerLists(1, max_part=0))
        []

    Negative parts also work (note that you probably want to change
    ``min_part_last`` too in this case!)::

        sage: IntegerLists(4, min_part=-1, max_slope=-1).list()
        [[4], [3, 1]]
        sage: IntegerLists(4, min_part=-1, min_part_last=-infinity, max_slope=-1).list()
        [[4],
         [3, 1],
         [4, 0],
         [5, -1],
         [5, 0, -1],
         [4, 1, -1],
         [3, 1, 0],
         [3, 2, -1],
         [4, 1, 0, -1],
         [3, 2, 0, -1]]

    We allow an infinite iterator::

        sage: L = IntegerLists(ceiling=[0], min_slope=1, max_slope=2)
        sage: it = iter(L)
        sage: for _ in range(20):
        ....:     print next(it)
        []
        [0, 1]
        [0, 2]
        [0, 1, 2]
        [0, 1, 3]
        [0, 2, 3]
        [0, 2, 4]
        [0, 1, 2, 3]
        [0, 1, 2, 4]
        [0, 1, 3, 4]
        [0, 1, 3, 5]
        [0, 2, 3, 4]
        [0, 2, 3, 5]
        [0, 2, 4, 5]
        [0, 2, 4, 6]
        [0, 1, 2, 3, 4]
        [0, 1, 2, 3, 5]
        [0, 1, 2, 4, 5]
        [0, 1, 2, 4, 6]
        [0, 1, 3, 4, 5]

    However, we do not allow infinitely many lists of the same length::

        sage: IntegerLists(4, min_part=-infinity, length=3, min_slope=1).list()
        Traceback (most recent call last):
        ...
        RuntimeError: there seem to be infinitely many lists of length 3

    TESTS::

    All these used to be broken (:trac:`17548`)::

        sage: IntegerLists(4, min_part=1, length=3, min_slope=1).list()
        []
        sage: IntegerLists(6, length=2, ceiling=[4,2], floor=[3,3]).list()
        []
        sage: IntegerLists(6, max_part=3, max_slope=-1).list()
        [[3, 2, 1]]
        sage: IntegerLists(10, min_part=2, max_slope=-1).list()
        [[10], [6, 4], [7, 3], [8, 2], [5, 3, 2]]
        sage: IntegerLists(5, min_part=1, max_part=2, min_slope=1, floor=[2]).list()
        []
        sage: IntegerLists(4, min_slope=0, max_slope=0).list()
        [[4], [2, 2], [1, 1, 1, 1]]
        sage: IntegerLists(6, min_slope=-1, max_slope=-1).list()
        [[6], [3, 2, 1]]
        sage: IntegerLists(min_length=3, max_length=2).list()
        []
        sage: IntegerLists(4, min_part=1, max_part=2, min_slope=1).list()
        []
        sage: IntegerLists(7, min_part=1, max_part=4, floor=[1], min_slope=1).list()
        [[3, 4], [1, 2, 4]]
        sage: IntegerLists(4, min_length=1, max_length=2, floor=[2], ceiling=[2,2], min_slope=0).list()
        [[2, 2]]
        sage: IntegerLists(10, min_length=2, min_part=1, floor=[1,1], max_slope=-1).list()
        [[6, 4],
         [7, 3],
         [8, 2],
         [9, 1],
         [7, 2, 1],
         [6, 3, 1],
         [5, 3, 2],
         [5, 4, 1],
         [4, 3, 2, 1]]
        sage: L = IntegerLists(6)
        sage: it = iter(L)
        sage: while True:
        ....:     if next(it) == L([1,2,3]):
        ....:         break
        sage: IntegerLists(5, floor=[2], ceiling=[2], min_slope=-1, max_slope=1, max_length=6).list()
        [[2, 3],
         [2, 1, 2],
         [2, 2, 1],
         [2, 1, 1, 1],
         [2, 1, 0, 1, 1],
         [2, 1, 1, 0, 1],
         [2, 1, 0, 0, 1, 1],
         [2, 1, 0, 1, 0, 1],
         [2, 1, 1, 0, 0, 1]]

    Older tests::

        sage: g = lambda x: lambda i: x
        sage: list(IntegerLists(0, floor=g(1), min_slope=0))
        [[]]
        sage: list(IntegerLists(0, floor=g(1), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerLists(0, max_length=0, floor=g(1), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerLists(0, max_length=0, floor=g(0), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerLists(0, min_part=1, min_slope=0))
        [[]]
        sage: list(IntegerLists(1, min_part=1, min_slope=0))
        [[1]]
        sage: list(IntegerLists(0, min_length=1, min_part=1, min_slope=0))
        []
        sage: list(IntegerLists(0, min_length=1, min_slope=0))
        [[0]]
        sage: list(IntegerLists(3, max_length=2))
        [[3], [0, 3], [1, 2], [2, 1]]
        sage: partitions = {"min_part": 1, "max_slope": 0}
        sage: partitions_min_2 = {"floor": g(2), "max_slope": 0}
        sage: compositions = {"min_part": 1}
        sage: integer_vectors = lambda l: {"length": l}
        sage: lower_monomials = lambda c: {"length": c, "floor": lambda i: c[i]}
        sage: upper_monomials = lambda c: {"length": c, "ceiling": lambda i: c[i]}
        sage: constraints = { "min_part":1, "min_slope": -1, "max_slope": 0}
        sage: list(IntegerLists(6, **partitions))
        [[6],
         [3, 3],
         [4, 2],
         [5, 1],
         [3, 2, 1],
         [4, 1, 1],
         [2, 2, 2],
         [3, 1, 1, 1],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerLists(6, **constraints))
        [[6],
         [3, 3],
         [2, 2, 2],
         [3, 2, 1],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerLists(1, **partitions_min_2))
        []
        sage: list(IntegerLists(2, **partitions_min_2))
        [[2]]
        sage: list(IntegerLists(3, **partitions_min_2))
        [[3]]
        sage: list(IntegerLists(4, **partitions_min_2))
        [[4], [2, 2]]
        sage: list(IntegerLists(5, **partitions_min_2))
        [[5], [3, 2]]
        sage: list(IntegerLists(6, **partitions_min_2))
        [[6], [3, 3], [4, 2], [2, 2, 2]]
        sage: list(IntegerLists(7, **partitions_min_2))
        [[7], [4, 3], [5, 2], [3, 2, 2]]
        sage: list(IntegerLists(9, **partitions_min_2))
        [[9], [5, 4], [6, 3], [7, 2], [4, 3, 2], [5, 2, 2], [3, 3, 3], [3, 2, 2, 2]]
        sage: list(IntegerLists(10, **partitions_min_2))
        [[10],
         [5, 5],
         [6, 4],
         [7, 3],
         [8, 2],
         [6, 2, 2],
         [5, 3, 2],
         [4, 3, 3],
         [4, 4, 2],
         [4, 2, 2, 2],
         [3, 3, 2, 2],
         [2, 2, 2, 2, 2]]
        sage: list(IntegerLists(4, **compositions))
        [[4], [1, 3], [2, 2], [3, 1], [1, 1, 2], [1, 2, 1], [2, 1, 1], [1, 1, 1, 1]]
        sage: list(IntegerLists(6, min_length=1, floor=[7]))
        []

    ALGORITHM: enumerate integral points in polyhedra. This is a simple
    but relatively slow algorithm (mostly because polyhedra in Sage are
    not very fast).
    """
    Element = IntegerList_polyhedron

    def __init__(self,
                 n=None,
                 length=None, min_length=0, max_length=infinity,
                 floor=None, ceiling=None,
                 min_part=0, max_part=infinity,
                 min_slope=-infinity, max_slope=infinity,
                 *, min_sum=0, max_sum=infinity,
                 min_part_last=1, max_part_last=infinity,
                 name=None,
                 element_constructor=None,
                 element_class=None,
                 **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: C = IntegerLists(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C == loads(dumps(C)) # this did fail at some point, really!
            True
            sage: C is loads(dumps(C)) # todo: not implemented
            True
            sage: C.cardinality().parent() is ZZ
            True
            sage: TestSuite(C).run()
        """
        for k in kwds:
            v = kwds[k]
            if v is not None:
                setattr(self, k, v)

        if n is not None:
            min_sum = n
            max_sum = n
            # Set self.n, which is not used by IntegerLists_polyhedron,
            # but by some derived classes
            self.n = n

        if length is not None:
            min_length = length
            max_length = length

        if min_length < 0:
            min_length = 0

        self.min_sum = min_sum
        self.max_sum = max_sum
        self.min_length = min_length
        self.max_length = max_length
        self.min_part = min_part
        self.max_part = max_part
        self.min_part_last = min_part_last
        self.max_part_last = max_part_last
        self.min_slope = min_slope
        self.max_slope = max_slope

        # Floor/ceiling should either be iterable, or a callable function
        if not floor:
            floor_iter = None
        else:
            try:
                floor_iter = iter(floor)
            except TypeError:
                floor_iter = _function_iter(floor)

        if not ceiling:
            ceil_iter = None
        else:
            try:
                ceil_iter = iter(ceiling)
            except TypeError:
                ceil_iter = _function_iter(ceiling)

        if name is not None:
            self.rename(name)

        self.floor_iter = floor_iter
        self.ceil_iter = ceil_iter
        if element_class is not None:
            self.Element = element_class
        if element_constructor is not None:
            self._element_constructor_ = element_constructor
        Parent.__init__(self, category=FiniteEnumeratedSets())

        # floor_list and ceil_list will contain lower and upper bounds
        # for the i-th part, keeping in mind min/max_part and
        # floor/ceil_iter and also min/max_slope.
        self.floor_list = []
        self.ceil_list = []

        # Dynamic version of max_length: if the conditions imply some
        # upper bound on the length, this attribute will be changed to
        # reflect that.
        if self.is_trivially_zero():
            self.effective_max_length = 0
        else:
            self.effective_max_length = max_length

        # try_length is the minimal length of lists to compute before
        # we worry that there are no more possible lists.
        if self.effective_max_length < infinity:
            self.try_length = self.effective_max_length
        else:
            self.try_length = max([8 + self.min_length,
                len_or_0(floor), len_or_0(ceiling)])
            if self.max_sum < infinity:
                self.try_length = max(self.try_length, self.max_sum)
            if -self.min_sum < infinity:
                self.try_length = max(self.try_length, -self.min_sum)
        assert self.try_length < infinity

    def is_trivially_zero(self):
        """
        Do the conditions trivially exclude a list of length > 0?

        Regardless of what this method returns, the empty list may or
        may not satisfy the constraints.

        EXAMPLES::

            sage: L = IntegerLists(min_part=3, max_part=2)
            sage: L.is_trivially_zero()
            True
            sage: list(L)  # the empty list satisfies
            [[]]
            sage: L = IntegerLists(max_sum=2, min_part=3)
            sage: L.is_trivially_zero()
            True
            sage: list(L)
            [[]]

        Note that we don't look at ``floor`` or ``ceiling``::

            sage: L = IntegerLists(floor=[-2], ceiling=[-3])
            sage: L.is_trivially_zero()
            False
            sage: list(L)
            [[]]
        """
        return (
            self.min_sum == infinity or
            self.min_part == infinity or
            self.max_length <= 0 or
            self.max_sum == -infinity or
            self.max_part == -infinity or
            self.min_sum > self.max_sum or
            self.min_length > self.max_length or
            self.min_part > self.max_part or
            _minimal_rectangle_sum(self.min_part, max(1, self.min_length), self.max_length) > self.max_sum or
            _maximal_rectangle_sum(self.max_part, max(1, self.min_length), self.max_length) < self.min_sum)

    def _element_constructor_(self, lst):
        """
        Construct an element with ``self`` as parent.

        EXAMPLES::

            sage: C = IntegerLists(4)
            sage: C([4])
            [4]
        """
        return self.element_class(self, lst)

    def __cmp__(self, x):
        """
        Compares two different :class:`IntegerLists_polyhedron`.

        For now, the comparison is done just on their repr's which is
        not robust!

        EXAMPLES::

            sage: C = IntegerLists(2, length=3)
            sage: D = IntegerLists(4, length=3)
            sage: repr(C) == repr(D)
            False
            sage: C == D
            False
        """
        return cmp(repr(self), repr(x))

    def _repr_(self):
        """
        Return the name of this combinatorial class.

        EXAMPLES::

            sage: C = IntegerLists(2, length=3)
            sage: C
            Integer lists of sum 2 satisfying certain constraints
            sage: C = IntegerLists(min_sum=-1, max_sum=4)
            sage: C
            Integer lists of sum in [-1, 4] satisfying certain constraints
            sage: C = IntegerLists(min_sum=-infinity, max_sum=4)
            sage: C
            Integer lists of sum at most 4 satisfying certain constraints
            sage: C = IntegerLists()
            sage: C
            Integer lists of sum at least 0 satisfying certain constraints
            sage: C = IntegerLists(min_sum=-infinity, max_sum=infinity)
            sage: C
            Integer lists satisfying certain constraints
            sage: C = IntegerLists(min_sum=1, max_sum=4, name="A given name")
            sage: C
            A given name
        """
        if self.min_sum == self.max_sum:
            s= " of sum {}".format(self.min_sum)
        elif self.min_sum == -infinity:
            if self.max_sum == infinity:
                s = ""
            else:
                s = " of sum at most {}".format(self.max_sum)
        elif self.max_sum == infinity:
            s = " of sum at least {}".format(self.min_sum)
        else:
            s = " of sum in [{}, {}]".format(self.min_sum, self.max_sum)
        return "Integer lists" + s + " satisfying certain constraints"

    def floor(self, i):
        """
        Return the minimum part that can appear at the `i^{th}` position of
        any list produced.

        EXAMPLES::

            sage: C = IntegerLists(4, min_part=1)
            sage: C.floor(0)
            1
            sage: C = IntegerLists(4, floor=[1,2])
            sage: C.floor(0)
            1
            sage: C.floor(1)
            2
            sage: C.floor(2)
            0

            sage: C = IntegerLists(min_part=-1, ceiling=[2], max_slope=1)
            sage: [C.floor(i) for i in range(10)]
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
        """
        return self.get_floor_ceil(i)[0]

    def ceiling(self, i):
        """
        Return the maximum part that can appear in the `i^{th}`
        position in any list produced.

        EXAMPLES::

            sage: C = IntegerLists(4, max_part=3)
            sage: C.ceiling(0)
            3
            sage: C = IntegerLists(4, ceiling=[3,2])
            sage: C.ceiling(0)
            3
            sage: C.ceiling(1)
            2
            sage: C.ceiling(2)
            4

            sage: C = IntegerLists(min_part=-1, ceiling=[2], max_slope=1)
            sage: [C.ceiling(i) for i in range(10)]
            [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        Without conditions, we get infinity::

            sage: C = IntegerLists()
            sage: C.ceiling(0)
            +Infinity
        """
        return self.get_floor_ceil(i)[1]

    def get_floor_ceil(self, Py_ssize_t i):
        """
        Compute ``self.floor_list`` and ``self.ceil_list`` until at
        least position ``i`` and return
        ``(self.floor_list[i], self.ceil_list[i])``.

        EXAMPLES::

            sage: IntegerLists(6).get_floor_ceil(10)
            (0, 6)

        TESTS::

            sage: C = IntegerLists()
            sage: C.get_floor_ceil(-1)
            Traceback (most recent call last):
            ...
            IndexError: part -1 does not exist
        """
        if i < 0:
            raise IndexError("part %r does not exist"%i)

        cdef Py_ssize_t n = len(self.floor_list)
        assert n == len(self.ceil_list)

        while n <= i:
            floor_sum = sum(self.floor_list)
            ceil_sum = sum(self.ceil_list)

            # Compute floor_list[n] and ceil_list[n]
            p = self.min_part
            if self.floor_iter is not None:
                try:
                    p = max(p, next(self.floor_iter))
                except StopIteration:
                    self.floor_iter = None
            if n > 0:
                p = max(p, self.floor_list[n-1] + self.min_slope)
            if self.max_part <= 0:
                # Total sum of *other* parts is at most ceil_sum
                p = max(p, self.min_sum - ceil_sum)
            self.floor_list.append(p)

            p = self.max_part
            if self.ceil_iter is not None:
                try:
                    p = min(p, next(self.ceil_iter))
                except StopIteration:
                    self.ceil_iter = None
            if n > 0:
                p = min(p, self.ceil_list[n-1] + self.max_slope)
            if self.min_part >= 0:
                # Total sum of *other* parts is at least floor_sum
                p = min(p, self.max_sum - floor_sum)
            self.ceil_list.append(p)

            n += 1

        return (self.floor_list[i], self.ceil_list[i])

    def _polyhedron_ieqs(self, Py_ssize_t length, Py_ssize_t sumlength, Py_ssize_t dimension):
        """
        Return a list of inequalities for the :class:`Polyhedron`
        constructor.

        INPUT:

        - ``length`` -- length of sequence to generate inequalities
          for: the first ``length`` variables will have inequalities
          involving floor/ceiling and will have bounds on their slopes.

        - ``sumlength`` -- number of variables to use for sum
          inequalities: the first ``sumlength`` variables will be
          summed.

        - ``dimension`` -- ambient dimension of the target polyhedron,
          i.e. total number of variables to use.

        OUTPUT:

        - if the polyhedron is trivially empty, return ``None``
        
        - otherwise, return a list of inequalities (really, lists of
          length ``dimension+1``)

        EXAMPLES::

            sage: C = IntegerLists(6, min_part=1)
            sage: C._polyhedron_ieqs(3, 3, 3)
            [[-1, 1, 0, 0],
             [6, -1, 0, 0],
             [-1, 0, 1, 0],
             [5, 0, -1, 0],
             [-1, 0, 0, 1],
             [4, 0, 0, -1],
             [+Infinity, -1, 1, 0],
             [+Infinity, 1, -1, 0],
             [+Infinity, 0, -1, 1],
             [+Infinity, 0, 1, -1],
             [-6, 1, 1, 1],
             [6, -1, -1, -1]]

        For length 7, the floor is larger than the ceiling, so the
        inequalities are trivially contradictory and ``None`` is
        returned::

            sage: C.get_floor_ceil(6)
            (1, 0)
            sage: print C._polyhedron_ieqs(7, 7, 7)
            None

        We can ask for more dimensions to manually add conditions::

            sage: C._polyhedron_ieqs(3, 5, 7)
            [[-1, 1, 0, 0, 0, 0, 0, 0],
             [6, -1, 0, 0, 0, 0, 0, 0],
             [-1, 0, 1, 0, 0, 0, 0, 0],
             [5, 0, -1, 0, 0, 0, 0, 0],
             [-1, 0, 0, 1, 0, 0, 0, 0],
             [4, 0, 0, -1, 0, 0, 0, 0],
             [+Infinity, -1, 1, 0, 0, 0, 0, 0],
             [+Infinity, 1, -1, 0, 0, 0, 0, 0],
             [+Infinity, 0, -1, 1, 0, 0, 0, 0],
             [+Infinity, 0, 1, -1, 0, 0, 0, 0],
             [-6, 1, 1, 1, 1, 1, 0, 0],
             [6, -1, -1, -1, -1, -1, 0, 0]]
        """
        assert 0 <= length <= dimension
        assert 0 <= sumlength <= dimension

        cdef list ieqs = []
        cdef list e
        cdef Py_ssize_t i

        # Inequalities for floor/ceiling
        for i in range(length):
            f, c = self.get_floor_ceil(i)
            if f >= c+1:
                # Contradiction: all polyhedra of length > i must be empty.
                # Note: the +1 in the formula above makes this work even
                # if floor and ceiling are both +/- infinity
                self.effective_max_length = min(self.effective_max_length, i)
                return None
            e = [0] * dimension
            e[i] = 1
            ieqs.append([-f] + e)
            e[i] = -1
            ieqs.append([c] + e)

        # Inequalties for slopes
        for i in range(length-1):
            e = [0] * dimension
            e[i+1] = 1
            e[i] = -1
            ieqs.append([-self.min_slope] + e)
            e[i+1] = -1
            e[i] = 1
            ieqs.append([self.max_slope] + e)

        # Inequalities for sum
        e = [1] * sumlength + [0] * (dimension - sumlength)
        ieqs.append([-self.min_sum] + e)
        e = [-1] * sumlength + [0] * (dimension - sumlength)
        ieqs.append([self.max_sum] + e)

        return ieqs

    @cached_method
    def polyhedron(self, Py_ssize_t length):
        """
        Return a polyhedron representing a list of integers of length
        ``length``. The polyhedron has one dimension for every part.

        Note: the length conditions are not checked. It returns a
        possibly non-empty polyhedron for all lengths.

        INPUT:

        - ``length`` -- a non-negative integer: the length of the
          ``IntegerList``, also the dimension of the polyhedron.

        EXAMPLES::

            sage: C = IntegerLists(10, min_part=-2, min_part_last=-2, min_slope=1)
            sage: C.polyhedron(0)
            The empty polyhedron in QQ^0
            sage: C.polyhedron(1)
            A 0-dimensional polyhedron in QQ^1 defined as the convex hull of 1 vertex
            sage: P = C.polyhedron(7)
            sage: P
            A 6-dimensional polyhedron in QQ^7 defined as the convex hull of 7 vertices

        In this general, this is not a lattice polytope::

            sage: P.vertices()
            (A vertex at (-2, -1, 3/5, 8/5, 13/5, 18/5, 23/5),
             A vertex at (-2, -1, 0, 1, 3, 4, 5),
             A vertex at (-2, -1, 0, 7/4, 11/4, 15/4, 19/4),
             A vertex at (-2, -1, 0, 1, 2, 3, 7),
             A vertex at (-2, -1/2, 1/2, 3/2, 5/2, 7/2, 9/2),
             A vertex at (-2, -1, 0, 1, 2, 9/2, 11/2),
             A vertex at (-11/7, -4/7, 3/7, 10/7, 17/7, 24/7, 31/7))

        TESTS::

            sage: C.polyhedron(-1)
            Traceback (most recent call last):
            ...
            ValueError: length must be >= 0
        """
        if length < 0:
            raise ValueError("length must be >= 0")

        cdef list ieqs = self._polyhedron_ieqs(length, length, length)
        if ieqs is None:
            # Empty polyhedron
            return Polyhedron(base_ring=QQ, ambient_dim=length)

        # Extra inequalities for the last part if length more than minimum
        cdef list e
        if length > self.min_length:
            e = [0] * length
            e[length-1] = 1
            ieqs.append([-self.min_part_last] + e)
            e[length-1] = -1
            ieqs.append([self.max_part_last] + e)

        return Polyhedron_inf(ieqs=ieqs, base_ring=QQ)

    @cached_method
    def polyhedron_more(self, Py_ssize_t length):
        """
        Return a polyhedron representing a list of integers of length
        more than ``length``. The polyhedron has ``length+1`` dimensions,
        where the first ``length`` dimensions represent parts as usual,
        and the last dimension represents the sum of all additional
        parts.

        INPUT:

        - ``length`` -- an integer >= 1: the polyhedron has dimension
          ``length+1`` and represents sequences of length more than
          ``length``.

        EXAMPLES::

            sage: C = IntegerLists(10, min_part=-2, min_part_last=-2, min_slope=1)
            sage: P = C.polyhedron_more(1)
            sage: P
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices
            sage: P.vertices()
            (A vertex at (-2, 12), A vertex at (8, 2))
            sage: C.polyhedron_more(5).vertices()
            (A vertex at (-2, 3/4, 7/4, 11/4, 15/4, 3),
             A vertex at (-2, -1, 7/3, 10/3, 13/3, 3),
             A vertex at (-2, -1, 0, 9/2, 11/2, 3),
             A vertex at (-3/5, 2/5, 7/5, 12/5, 17/5, 3),
             A vertex at (-2, -1, 0, 1, 2, 10),
             A vertex at (-2, -1, 0, 1, 9, 3))

        This proves that there are no sequences of length 8 or more::

            sage: C.polyhedron_more(7)
            The empty polyhedron in QQ^8

        The implementation requires ``length >= 1``::

            sage: C.polyhedron_more(0)
            Traceback (most recent call last):
            ...
            ValueError: length must be >= 1
        """
        if length < 1:
            raise ValueError("length must be >= 1")

        # Start with empty polyhedron
        P = Polyhedron(base_ring=QQ, ambient_dim=length+1)

        cdef list ieqs0 = self._polyhedron_ieqs(length, length+1, length+1)
        if ieqs0 is None:
            return P

        # Instead of computing one polyhedron, we compute 3 and take
        # the convex hull: we compute one polyhedron for every choice of
        # sign (<0, =0, >0) of part "length-1".
        fa, ca = self.get_floor_ceil(length-1)
        cdef list e, ieqs
        cdef long sign
        for sign from -1 <= sign <= 1:
            # Floor/ceiling of part "length-1" keeping in mind sign
            f = fa if sign == -1 else max(fa, smallInteger(sign))
            c = ca if sign == 1 else min(ca, smallInteger(sign))
            if f >= c+1:
                continue

            # Limit of floor/ceiling for part index going to infinity
            if f == -infinity:
                flim = self.min_part
            else:
                flim = max(self.min_part, f + signed_infinity(self.min_slope))
            if c == infinity:
                clim = self.max_part
            else:
                clim = min(self.max_part, c + signed_infinity(self.max_slope))

            # Bounds on last part (which has index >= length)
            flast = max([self.min_part, self.min_part_last, min(f+self.min_slope, flim)])
            clast = min([self.max_part, self.max_part_last, max(c+self.max_slope, clim)])
            if flast >= clast+1:
                continue

            # Sum of floors/ceilings
            if f < 0:
                if flim >= 0:
                    # The sum is actually a finite negative number,
                    # try again with a larger length!
                    self.try_length = max(self.try_length, length - f)
                fs = -infinity
            else:
                # Take only last part to get minimal sum
                fs = flast

            if c > 0:
                if clim <= 0:
                    self.try_length = max(self.try_length, length + c)
                cs = infinity
            else:
                cs = clast

            ieqs = ieqs0[:]

            # Inequalities for part "length-1"
            e = [0] * (length+1)
            e[length-1] = 1
            ieqs.append([-f] + e)
            e[length-1] = -1
            ieqs.append([c] + e)

            # Inequalities for parts >= length
            e = [0] * (length+1)
            e[length] = 1
            ieqs.append([-fs] + e)
            e[length] = -1
            ieqs.append([cs] + e)

            P = P.convex_hull(Polyhedron_inf(ieqs=ieqs, base_ring=QQ))

        return P

    def __iter__(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: C = IntegerLists(2, length=3)
            sage: list(C)  # indirect doctest
            [[0, 0, 2], [0, 1, 1], [0, 2, 0], [1, 0, 1], [1, 1, 0], [2, 0, 0]]
        """
        cdef list L = []

        # Iterate by length
        cdef Py_ssize_t length = self.min_length

        # Note: effective_max_length can change during this loop
        while length <= self.effective_max_length:
            sig_check()
            P = self.polyhedron(length)
            if length > self.try_length and P.is_empty():
                # Perhaps all following polyhedra are empty?
                Q = self.polyhedron_more(length)
                if Q.is_empty():
                    self.effective_max_length = min(self.effective_max_length, length-1)
                    return
                elif length > self.try_length:  # Do nothing if try_length was increased
                    raise RuntimeError("no more lists found, but cannot prove that there are none of length > {}".format(length))
            elif not P.is_compact():
                raise RuntimeError("there seem to be infinitely many lists of length {}".format(length))

            for x in P.integral_points():
                yield self._element_constructor_(x.list())

            length += 1

    def count(self):
        """
        Default brute force implementation of count by iteration
        through all the objects.

        EXAMPLES::

            sage: C = IntegerLists(2, length=3)
            sage: C.cardinality() == C.count()
            True
        """
        cdef Py_ssize_t n = 0
        for x in self:
            n += 1
        return n

    def __contains__(self, v):
        """
        Return ``True`` if and only if ``v`` is in ``self``.

        EXAMPLES::

            sage: C = IntegerLists(2, length=3)
            sage: [2, 0, 0] in C
            True
            sage: [2, 0] in C
            False
            sage: [3, 0, 0] in C
            False
            sage: all(v in C for v in C)
            True
        """
        cdef Py_ssize_t length = len(v)
        if length < self.min_length or length > self.effective_max_length:
            return False
        return self.polyhedron(length).__contains__(v)


class IntegerListsLex_polyhedron(IntegerLists_polyhedron):
    """
    Frontend for :class:`IntegerLists_polyhedron` which returns the
    elements sorted in reverse lexicographic order.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ....:     return x^exponents[0] * y^exponents[1] * z^exponents[2]
        sage: C = IntegerListsLex(4, length=len(m), ceiling=m, element_constructor=term)
        sage: list(C)
        [x^3*y, x^3*z, x^2*y*z, x^2*z^2, x*y*z^2]
    """
    def __init__(self, *args, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(7, max_length=3, min_part=2)
            sage: C.list()
            [[7], [5, 2], [4, 3], [3, 4], [3, 2, 2], [2, 5], [2, 3, 2], [2, 2, 3]]
        """
        # If element_constructor is given, don't pass it to
        # IntegerLists_polyhedron but store it in a private attribute.
        ec = kwds.pop('element_constructor', None)
        if ec is not None:
            self.__element_constructor = ec
        else:
            self.__element_constructor = self._element_constructor_
        IntegerLists_polyhedron.__init__(self, *args, **kwds)

    def _repr_(self):
        """
        Return the name of this combinatorial class.

        EXAMPLES::

            sage: C = IntegerListsLex(5, max_length=3)
            sage: C
            Integer lists of sum 5 satisfying certain constraints, in revlex order
        """
        return IntegerLists_polyhedron._repr_(self) + ", in revlex order"

    def __iter__(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C)  # indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        L = [x for x in IntegerLists_polyhedron.__iter__(self)]
        L.sort(key=lambda t: [-a for a in t])
        for t in L:
            yield self.__element_constructor(t)
