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

from sage.structure.element cimport MonoidElement, Element
from sage.misc.bounded_integer_sequences cimport biseq_t, allocate_biseq, getitem_biseq, concat_biseq, startswith_biseq, contains_biseq, max_overlap_biseq, slice_biseq, list_to_biseq, biseq_to_list

from sage.rings.integer_ring import ZZ
from cython.operator import dereference as deref

include "sage/ext/stdsage.pxi"
include "sage/libs/ntl/decl.pxi"
include "sage/ext/interrupt.pxi"

cdef extern from "Python.h":
    bint PySlice_Check(PyObject* ob)

cdef extern from "mpz_pylong.h":
    cdef long mpz_pythonhash(mpz_t src)

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
    cdef biseq_t _path
    cdef int _start, _end

    def __dealloc__(self):
        # I tested that this will not crash, even if self._path isn't initialised
        mpz_clear(self._path.data)

    cdef QuiverPath _new_(self, int start, int end, biseq_t data):
        cdef QuiverPath out = PY_NEW(self._parent.element_class)
        out._parent = self._parent
        out._start = start
        out._end = end
        out._path = data
        return out

    def __init__(self, parent, start, end, path, check=True):
        """
        Creates a path object.  Type ``QuiverPath?`` for more information.

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
        use a sub-class of QuiverPath. Nonetheless, just for test, we
        show that it *is* possible to create a path in a deprecated way::

            sage: p == QuiverPath(Q, 1, 1, [])
            True
            sage: list(Q([(1, 1)])*Q([(1, 2, 'a')])*Q([(2, 2)])*Q([(2, 3, 'b')])*Q([(3, 3)]))
            [(1, 2, 'a'), (2, 3, 'b')]
        """
        MonoidElement.__init__(self, parent=parent)
        cdef unsigned int l = len(path)
        cdef mpz_t tmp
        self._start = start
        self._end   = end
        mpz_init_set_ui(tmp, len(parent.quiver().edges())-1)
        self._path = deref(list_to_biseq(deref(allocate_biseq(l, mpz_sizeinbase(tmp, 2))), path))
        mpz_clear(tmp)
        if not check:
            return
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
        cdef unsigned int n
        for n from 0<n<l:
            if E[path[n-1]][1]!=E[path[n]][0]:
                raise ValueError("Edge {} ends at {}, but edge {} starts at {}".format(E[path[n-1]][2], E[path[n-1]][1], E[path[n]][2], E[path[n]][0]))

    def __reduce__(self):
        cdef size_t n
        cdef char *s
        n = mpz_sizeinbase(self._path.data, 32) + 2
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string defining a bounded integer sequence."
        sig_on()
        mpz_get_str(s, 32, self._path.data)
        sig_off()
        data_str = <object> PyString_FromString(s)
        PyMem_Free(s)
        return NewQuiverPath, (self._parent, self._start, self._end, data_str, self._path.bitsize, self._path.itembitsize, self._path.length)

    def __hash__(self):
        if self._path.length==0:
            return hash(self._start)
        return mpz_pythonhash(self._path.data)

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
        if self._path.length:
            E = self._parent.quiver().edges()
            return '*'.join([E[e][2] for e in biseq_to_list(self._path)])
        return 'e_{0}'.format(self._start)

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

        As usual in Sage, the ``__cmp__`` method of a Python sub-class of
        :class:`sage.structure.element.Element` can assume that both arguments
        belong to the same parent.

        If the QuiverPaths are unequal then one of the following data (listed
        in order of preferance) is unequal and used for comparison:

        - **Negative** length of the paths
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
        return mpz_cmp(cself._path.data, other._path.data)

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
        cdef unsigned int init, end
        if PySlice_Check(<PyObject *>index):
            start,stop,step = index.indices(self._path.length)
            if step!=1:
                raise ValueError("Slicing only possible for step 1")
            if start==0 and stop==self._path.length:
                return self
            init = E[getitem_biseq(self._path, start)][0]
            end   = E[getitem_biseq(self._path, stop)][0]
            return self._new_(init, end, deref(slice_biseq(self._path, start, stop, step)))
        if index<0:
            index = self._path.length+index
        if index<0 or index>=self._path.length:
            raise IndexError("list index out of range")
        init = E[getitem_biseq(self._path, index)][0]
        end = E[getitem_biseq(self._path, index)][1]
        return self._new_(init, end, deref(slice_biseq(self._path, index, index+1, 1)))

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
        E = self._parent.quiver().edges()
        for n in biseq_to_list(self._path):
            yield E[n]

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
        return self._new_(self._start, right._end, deref(concat_biseq(self._path,right._path)))

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
        if startswith_biseq(cself._path, right._path):
            return cself._new_(right._end, cself._end, deref(slice_biseq(cself._path, right._path.length, cself._path.length, 1)))
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
        cdef int i
        if startswith_biseq(P._path, self._path):
            i = 0
        else:
            i = max_overlap_biseq(self._path, P._path)
        if i==-1:
            return (None, None, None)
        return (self[:i], self[i:], P[self._path.length-i:])

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
        return Q.element_class(Q, self._end, self._start, biseq_to_list(self._path)[::-1], check=False)

cpdef QuiverPath NewQuiverPath(Q, start, end, data, bitsize, itembitsize, length):
    cdef QuiverPath out = PY_NEW(Q.element_class)
    out._parent = Q
    out._start = start
    out._end   = end
    mpz_init2(out._path.data, bitsize)
    out._path.bitsize = bitsize
    out._path.itembitsize = itembitsize
    out._path.mask_item = ((<unsigned int>1)<<itembitsize)-1
    out._path.length = length
    mpz_set_str(out._path.data, data, 32)
    return out
