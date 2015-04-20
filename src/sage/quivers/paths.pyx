"""
Quiver Paths
"""

#*****************************************************************************
#  Copyright (C) 2012    Jim Stark <jstarx@gmail.com>
#                2013/14 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.data_structures.bounded_integer_sequences cimport *
from cpython.slice cimport PySlice_Check
from cython.operator cimport dereference as deref, preincrement as preinc, predecrement as predec, postincrement as postinc

include "sage/ext/interrupt.pxi"
include 'sage/data_structures/bitset.pxi'

from sage.rings.integer_ring import ZZ

cdef class QuiverPath(MonoidElement):
    r"""
    Class for paths in a quiver.

    A path is given by two vertices, ``start`` and ``end``, and a finite
    (possibly empty) list of edges `e_1, e_2, \ldots, e_n` such that the
    initial vertex of `e_1` is ``start``, the final vertex of `e_i` is
    the initial vertex of `e_{i+1}`, and the final vertex of `e_n` is
    ``end``.  In the case where no edges are specified, we must have
    ``start = end`` and the path is called the trivial path at the given
    vertex.

    INPUT:

    - ``parent`` -- the path semigroup associated with a quiver; this is
      where the path will live
    - ``path`` -- tuple or iterable. If ``path`` is a tuple then it is
      assumed to be of the form ``(vertex, vertex)`` or
      ``(start, end, label)``.  In the first case the trivial path at the
      given vertex is created.  In the second case a path consisting of
      just the given edge is created.  If ``path`` is not a tuple then it
      is assumed to be an iterable variable giving the edges of a path,
      where each edge is in one of the two forms above.
    - ``check`` -- boolean (default: ``True``); if it is ``False``, no
      sanity check will be performed on the given iterable

    OUTPUT:

    - :class:`QuiverPath`

    .. NOTE::

        Do *not* use this constructor directly! Instead, pass the input to the
        path semigroup that shall be the parent of this path.

    EXAMPLES:

    Specify a path by giving a list of edges::

        sage: Q = DiGraph({1:{2:['a','d'], 3:['e']}, 2:{3:['b']}, 3:{1:['f'], 4:['c']}})
        sage: F = Q.path_semigroup()
        sage: p = F([(1, 2, 'a'), (2, 3, 'b')])
        sage: p
        a*b

    Paths are not *unique*, but different representations of "the same" path
    yield *equal* paths::

        sage: q = F([(1, 1)]) * F([(1, 2, 'a'), (2, 3, 'b')]) * F([(3, 3)])
        sage: p is q
        False
        sage: p == q
        True

    The ``*`` operator is concatenation of paths. If the two paths do not
    compose, its result is ``None``::

        sage: print(p*q)
        None
        sage: p*F([(3, 4, 'c')])
        a*b*c
        sage: F([(2,3,'b'), (3,1,'f')])*p
        b*f*a*b

    The length of a path is the number of edges in that path.  Trivial paths
    are therefore length-`0`::

        sage: len(p)
        2
        sage: triv = F([(1, 1)])
        sage: len(triv)
        0

    List index and slice notation can be used to access the edges in a path.
    QuiverPaths can also be iterated over.  Trivial paths have no elements::

        sage: for x in p: print x
        (1, 2, 'a')
        (2, 3, 'b')
        sage: list(triv)
        []

    There are methods giving the initial and terminal vertex of a path::

        sage: p.initial_vertex()
        1
        sage: p.terminal_vertex()
        3
    """
    def __cinit__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1)])  # indirect doctest

        """
        self._path.length = 0

    def __dealloc__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1)]) * Q([(1, 1)])
            sage: del p    # indirect doctest

        """
        biseq_dealloc(self._path)

    cdef QuiverPath _new_(self, int start, int end):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q(['a']) * Q(['b'])    # indirect doctest

        """
        cdef QuiverPath out = PY_NEW(self._parent.element_class)
        out._parent = self._parent
        out._start = start
        out._end = end
        return out

    def __init__(self, parent, start, end, path, check=True):
        """
        Creates a path object.  Type ``QuiverPath?`` for more information.

        INPUT:

        - ``parent``, a path semigroup.
        - ``start``, integer, the label of the initial vertex.
        - ``end``, integer, the label of the terminal vertex.
        - ``path``, list of integers, providing the list of arrows
          occuring in the path, labelled according to the position in
          the list of all arrows (resp. the list of outgoing arrows at
          each vertex).
        - ``check``, whether or not to check the input.

        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 1)]) * Q([(1, 1)])
            sage: Q([(1,3,'x')])
            Traceback (most recent call last):
            ...
            ValueError: (1, 3, 'x') is not in list

        Note that QuiverPath should not be called directly, because
        the elements of the path semigroup associated with a quiver
        may use a sub-class of QuiverPath. Nonetheless, just for test, we
        show that it *is* possible to create a path in a deprecated way::

            sage: p == QuiverPath(Q, 1, 1, [], max(len(Q.quiver().outgoing_edges(v)) for v in Q.quiver().vertices()))
            True
            sage: list(Q([(1, 1)])*Q([(1, 2, 'a')])*Q([(2, 2)])*Q([(2, 3, 'b')])*Q([(3, 3)]))
            [(1, 2, 'a'), (2, 3, 'b')]

        """
        MonoidElement.__init__(self, parent=parent)
        cdef unsigned int l = len(path)
        self._start = start
        self._end   = end
        biseq_init_list(self._path, path, parent._nb_arrows)
        if not check:
            return
        cdef unsigned int n
        cdef tuple arrow
        Q = parent.quiver()
        E = Q.edges()
        if not l:
            if start!=end:
                raise ValueError("Start and endpoint of a path of length 0 must coincide")
            if start not in Q:
                raise ValueError("Vertex {} not in the quiver".format(start))
            return
        V = Q.vertices()
        if start not in V:
            raise ValueError("Startpoint {} should belong to {}".format(start, V))
        if end not in V:
            raise ValueError("Startpoint {} should belong to {}".format(end, V))
        if E[path[0]][0]!=start:
            raise ValueError("First edge should start at vertex {}".format(start))
        if E[path[-1]][1]!=end:
            raise ValueError("Last edge should end at vertex {}".format(end))
        for n from 0<n<l:
            if E[path[n-1]][1]!=E[path[n]][0]:
                raise ValueError("Edge {} ends at {}, but edge {} starts at {}".format(E[path[n-1]][2], E[path[n-1]][1], E[path[n]][2], E[path[n]][0]))

    def __reduce__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q(['a']) * Q(['b'])
            sage: loads(dumps(p)) == p   # indirect doctest
            True
            sage: loads(dumps(p)) is p
            False

        """
        return NewQuiverPath, (self._parent, self._start, self._end,
                               biseq_pickle(self._path))

    def __hash__(self):
        """
        TESTS::

            sage: from sage.quivers.paths import QuiverPath
            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q(['a']) * Q(['b'])
            sage: q = Q([(1, 1)])
            sage: {p:1, q:2}[Q(['a','b'])]    # indirect doctest
            1

        """
        if self._path.length==0:
            return hash(self._start)
        cdef Py_hash_t h = self._start*(<Py_hash_t>1073807360) + biseq_hash(self._path)
        if h==-1:
            return -2
        return h
        ## bitset_hash is not a good hash either
        ## We should consider using FNV-1a hash, see http://www.isthe.com/chongo/tech/comp/fnv/,
        ## Or the hash defined in http://burtleburtle.net/bob/hash/doobs.html
        ## Or http://www.azillionmonkeys.com/qed/hash.html

    def _repr_(self):
        r"""
        Default representation of a path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]) # indirect doctest
            a*b
            sage: Q([(1, 1)]) # indirect doctest
            e_1
        """
        cdef mp_size_t i
        if not self._path.length:
            return 'e_{0}'.format(self._start)
        E = self._parent.quiver().edges()
        return '*'.join([E[biseq_getitem(self._path, i)][2] for i in range(self._path.length)])

    def __len__(self):
        """
        Return the length of the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: len(Q([(1, 2, 'a'), (2, 3, 'b')]))
            2
            sage: len(Q([(1, 1)]))
            0
            sage: len(Q([(1, 2, 'a')]))
            1
        """
        return self._path.length

    def deg(self):
        """
        Return the degree of the path, which is the same as its length.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: Q([(1, 2, 'a'), (2, 3, 'b')]).deg()
            2
            sage: Q([(1, 1)]).deg()
            0
            sage: Q([(1, 2, 'a')]).deg()
            1
        """
        return len(self)

    def __nonzero__(self):
        """
        Implement boolean values for the object.

        .. NOTE::

            The boolean value is always ``True``, since the partial semigroup
            formed by the paths of a quiver does not contain zero.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: a = Q([(1, 2, 'a')])
            sage: b = Q([(2, 3, 'b')])
            sage: bool(a*b)
            True
            sage: bool(Q.idempotents()[0])
            True
        """
        return True

    def __cmp__(left, right):
        return (<Element>left)._cmp(right)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        Comparison for :class:`QuiverPaths`.

        The following data (listed in order of preferance) is used for
        comparison:

        - **Negative** length of the paths
        - initial and terminal vertices of the paths
        - Edge sequence of the paths, by reverse lexicographical ordering.

        .. NOTE::

            This code is used by :class:`CombinatorialFreeModule` to order
            the monomials when printing elements of path algebras.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c'], 4:['d']}}).path_semigroup()
            sage: a = Q([(1, 2, 'a')])
            sage: b = Q([(1, 2, 'b')])
            sage: c = Q([(2, 3, 'c')])
            sage: d = Q([(2, 4, 'd')])
            sage: e = Q.idempotents()[3]
            sage: e < a  # e is shorter than a
            False
            sage: a < e
            True
            sage: d < a*c
            False
            sage: a*c < d
            True
            sage: a < b
            True
            sage: b < a
            False
            sage: a*c < a*d
            True
            sage: a*d < a*c
            False
            sage: a < a
            False

        """
        # Since QuiverPath inherits from Element, it is guaranteed that
        # both arguments are elements of the same path semigroup
        cdef QuiverPath cself, other 
        cself = left
        other = right
        # we want *negative* degree reverse lexicographical order
        cdef int c = cmp(other._path.length, cself._path.length)
        if c!=0:
            return c
        c = cmp((cself._start,cself._end), (other._start,other._end))
        if c or (cself._path.length==0):
            return c
        return biseq_cmp(cself._path, other._path)

    def __getitem__(self, index):
        """
        Implement index notation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: p
            a*b*c
            sage: p[0]
            a
            sage: p[-1]
            c
            sage: p[1:]
            b*c
        """
        cdef list E = self._parent._quiver.edges()
        cdef int start, stop, step
        cdef int init, end
        cdef size_t i,ind
        cdef QuiverPath OUT
        if PySlice_Check(index):
            start,stop,step = index.indices(self._path.length)
            if step!=1:
                raise ValueError("Slicing only possible for step 1")
            if start==0 and stop==self._path.length:
                return self
            init = E[biseq_getitem(self._path, start)][0]
            end   = E[biseq_getitem(self._path, stop)][0]
            OUT = self._new_(init, end)
            biseq_init_slice(OUT._path, self._path, start, stop, step)
            return OUT
        if index<0:
            index = self._path.length+index
        if index<0 or index>=self._path.length:
            raise IndexError("list index out of range")
        init = E[biseq_getitem(self._path, index)][0]
        end = E[biseq_getitem(self._path, index)][1]
        OUT = self._new_(init, end)
        biseq_init_slice(OUT._path, self._path, index, index+1, 1)
        return OUT

    def __iter__(self):
        """
        Iteration over the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b'), (3, 4, 'c')])
            sage: for e in p: print e
            (1, 2, 'a')
            (2, 3, 'b')
            (3, 4, 'c')
        """
        # Return an iterator over an empty tuple for trivial paths, otherwise
        # return an iterator for _path as a list
        cdef mp_size_t i
        E = self._parent.quiver().edges()
        for i from 0<=i<self._path.length:
            yield E[biseq_getitem(self._path, i)]

    cpdef MonoidElement _mul_(self, MonoidElement other):
        """
        Compose two paths.

        .. NOTE::

            ``None`` is returned if the terminal vertex of the first path
            does not coincide with the initial vertex of the second path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}, 3:{4:['c']}, 4:{5:['d']}}).path_semigroup()
            sage: x = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y = Q([(3, 4, 'c'), (4, 5, 'd')])
            sage: print y*x
            None
            sage: x*y
            a*b*c*d
            sage: x*Q([(3, 4, 'c')])
            a*b*c
            sage: x*Q([(3, 4, 'c'), (4, 5, 'd')])
            a*b*c*d
            sage: x*6
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*':
             'Partial semigroup formed by the directed paths of Multi-digraph on 5 vertices'
             and 'Integer Ring'
        """
        # By Sage's coercion model, both paths belong to the same quiver
        # In particular, both are QuiverPath
        cdef QuiverPath right = other
        if self._end != right._start:
            return None
        cdef QuiverPath OUT = self._new_(self._start, right._end)
        biseq_init_concat(OUT._path, self._path,right._path)
        return OUT

    def __mod__(self, other):
        """
        Return ``self`` with ``other`` deleted from the beginning.

        If ``other`` is not the beginning of ``self`` then ``None`` is
        returned.  Deleting the trivial path at vertex `v` from a path that
        begins at `v` does not change the path.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: a = Q([(1, 2, 'a')])
            sage: b = Q([(2, 3, 'b')])
            sage: e1 = Q([(1, 1)])
            sage: e2 = Q([(2, 2)])
            sage: p % a
            b
            sage: print p % b
            None
            sage: p % e1
            a*b
            sage: print p % e2
            None
        """
        cdef QuiverPath right = other
        cdef QuiverPath cself = self
        # Handle trivial case
        if right is None or cself._start!=right._start:
            return None
        if right._path.length==0:
            return self

        # If other is the beginning, return the rest
        cdef QuiverPath OUT
        if (cself._start == right._start) and biseq_startswith(cself._path, right._path):
            OUT = cself._new_(right._end, cself._end)
            biseq_init_slice(OUT._path, cself._path, right._path.length, cself._path.length, 1)
            return OUT
        else:
            return None

    def gcd(self, QuiverPath P):
        """
        Greatest common divisor of two quiver paths, with co-factors.

        INPUT:

        A :class:`QuiverPath` ``P``

        OUTPUT:

        - :class:`QuiverPath`s ``(C1,G,C2)`` such that ``self==C1*G`` and ``P=G*C2``, or
        - ``(None, None, None)``, if the paths do not overlap (or belong to different quivers).

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{1:['b'], 3:['c']}, 3:{1:['d']}}).path_semigroup()
            sage: p1 = Q(['c','d','a','b','a','c','d'])
            sage: p1
            c*d*a*b*a*c*d
            sage: p2 = Q(['a','b','a','c','d','a','c','d','a','b'])
            sage: p2
            a*b*a*c*d*a*c*d*a*b
            sage: S1, G, S2 = p1.gcd(p2)
            sage: S1, G, S2
            (c*d, a*b*a*c*d, a*c*d*a*b)
            sage: S1*G == p1
            True
            sage: G*S2 == p2
            True
            sage: p2.gcd(p1)
            (a*b*a*c*d*a, c*d*a*b, a*c*d)

        We test that a full overlap is detected::

            sage: p2.gcd(p2)
            (e_1, a*b*a*c*d*a*c*d*a*b, e_1)

        The absence of an overlap is detected::

            sage: p2[2:-1]
            a*c*d*a*c*d*a
            sage: p2[1:]
            b*a*c*d*a*c*d*a*b
            sage: print p2[2:-1].gcd(p2[1:])
            (None, None, None)

        """
        if self._parent is not P._parent:
            return (None, None, None)
        cdef size_t i, start
        sig_on()
        i = biseq_startswith_tail(P._path, self._path, 0)
        sig_off()
        if i==-1:
            return (None, None, None)
        return (self[:i], self[i:], P[self._path.length-i:])

    cpdef bint has_subpath(self, QuiverPath subpath) except -1:
        """
        Tells whether ``self`` contains a given sub-path.

        INPUT:

        ``subpath``, a path of positive length in the same path semigroup as
        ``self``.

        EXAMPLES::

            sage: S = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup()
            sage: S.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: (c*b*e*a).has_subpath(b*e)
            1
            sage: (c*b*e*a).has_subpath(b*f)
            0
            sage: (c*b*e*a).has_subpath(e_1)
            Traceback (most recent call last):
            ...
            ValueError: We only consider sub-paths of positive length
            sage: (c*b*e*a).has_subpath(None)
            Traceback (most recent call last):
            ...
            ValueError: The given sub-path is empty

        """
        if subpath is None:
            raise ValueError("The given sub-path is empty")
        if subpath._parent is not self._parent:
            raise ValueError("The two paths belong to different quivers")
        if subpath._path.length == 0:
            raise ValueError("We only consider sub-paths of positive length")
        cdef int v
        cdef object C
        cdef size_t i
        cdef size_t max_i, bitsize
        if self._path.length < subpath._path.length:
            return 0
        if biseq_contains(self._path, subpath._path, 0)==-1:
            return 0
        return 1

    cpdef bint has_initial_segment(self, QuiverPath subpath) except -1:
        """
        Tells whether ``self`` starts with a given sub-path.

        INPUT:

        ``subpath``, a path in the same path semigroup as ``self``.

        OUTPUT:

        ``0`` or ``1``, which stands for ``False`` resp. ``True``.

        EXAMPLES::

            sage: S = DiGraph({0:{1:['a'], 2:['b']}, 1:{0:['c'], 1:['d']}, 2:{0:['e'],2:['f']}}).path_semigroup()
            sage: S.inject_variables()
            Defining e_0, e_1, e_2, a, b, c, d, e, f
            sage: (c*b*e*a).has_initial_segment(b*e)
            0
            sage: (c*b*e*a).has_initial_segment(c*b)
            1
            sage: (c*b*e*a).has_initial_segment(e_1)
            1
            sage: (c*b*e*a).has_initial_segment(e_2)
            0
            
        """
        if subpath._parent is not self._parent:
            raise ValueError("The two paths belong to different quivers")
        if self._start != subpath._start:
            return 0
        if subpath._path.length==0:
            return 1
        if biseq_startswith(self._path, subpath._path):
            return 1
        return 0

    def initial_vertex(self):
        """
        Return the initial vertex of the path.

        OUTPUT:

        - integer, the label of the initial vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.initial_vertex()
            1
        """
        return self._start

    def terminal_vertex(self):
        """
        Return the terminal vertex of the path.

        OUTPUT:

        - integer, the label of the terminal vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: y = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: y.terminal_vertex()
            3
        """
        return self._end

    def reverse(self):
        """
        Return the path along the same edges in reverse order in the
        opposite quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: p = Q([(1, 2, 'a'), (2, 3, 'b')])
            sage: p.reverse()
            b*a
            sage: p.reverse().parent() is Q.reverse()
            True
            sage: e = Q.idempotents()[0]
            sage: e
            e_1
            sage: e.reverse()
            e_1
        """
        Q = self._parent.reverse()
        # Handle trivial paths
        if self._path.length==0:
            return Q.element_class(Q, self._end, self._start, [], check=False)

        # Reverse all the edges in the path, then reverse the path
        cdef mp_size_t i
        return Q.element_class(Q, self._end, self._start, [biseq_getitem(self._path,i) for i in xrange(self._path.length-1,-1,-1)], check=False)

cpdef QuiverPath NewQuiverPath(Q, start, end, biseq_data):
    """
    Return a new quiver path for given defining data.

    INPUT:

    - ``Q``, the path semigroup of a quiver
    - ``start``, an integer, the label of the startpoint
    - ``end``, an integer, the label of the endpoint
    - ``data``, a string: The bitmap encoding the path represented as integer
       at base `32`.
    - ``bitsize``, the number of bits used to store the path
    - ``itembitsize``, the number of bits used to store a single item
    - ``length``, the number of items in the path

    TESTS::

        sage: from sage.quivers.paths import QuiverPath
        sage: Q = DiGraph({1:{2:['a']}, 2:{3:['b']}}).path_semigroup()
        sage: p = Q(['a']) * Q(['b'])
        sage: loads(dumps(p)) == p   # indirect doctest
        True
        sage: p.__reduce__()
        (<...NewQuiverPath>,
         (Partial semigroup formed by the directed paths of Multi-digraph on 3 vertices,
          1,
          3,
          ((0, 4L, 1, ..., (4L,)), 2L, 2)))

    """
    cdef QuiverPath out = PY_NEW(Q.element_class)
    out._parent = Q
    out._start = start
    out._end   = end
    biseq_unpickle(out._path, biseq_data[0], biseq_data[1], biseq_data[2])
    return out
