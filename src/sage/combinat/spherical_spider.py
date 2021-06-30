# -*- coding: utf-8 -*-
r"""
### Introduction
This is an implementation of finitely generated free spiders.
The account below is intended as an extended introduction.
There are technical details which I skate over and which need to
be dealt with in the implementation.

This is intended to be the first phase of a project to implement
the Bendix-Knuth completion algorithm discussed in [5]_. This algorithm starts
with a finite presentation and iteratively constructs new relations.
If it terminates the final presentation is confluent.

### Spiders

The notion of spiders was introduced in [3]_.
A spider consists of a set (or vector space) with the operations of rotate, join and stitch.
There are axioms for these operations and a spider is equivalent to a strict pivotal
category or, with an extra axiom, to a strict spherical category. The motivation was
to study the category of finite dimensional representations of a quantum group.
These categories give examples of spiders and the main result of the paper was to
give finite confluent presentations for the rank two examples, A2, B2=C2, G2, (the rank one example
was known previously as the skein relation approach to the Jones polynomial).
A rank three example is given in [6]_.

It has been an open problem since then to construct finite confluent presentations
for higher rank examples. It is known (unpublished) that these examples are finitely
generated but not that these are finitely presented. A finite confuent presentation
would give algorithms for computing link polynomials and for doing the caculations
in [1]_.

### Webs.

The first phase of working with finite presentations is to understand the free objects.
The elements of a free spider are called webs.

The data for a free spider is first a (finite) set X with an involution x. Then the objects of
the free strict spherical category with be the free monoid on X with antiinvolution given on the
generators by x.

The usual approach to webs is to define a web as a planar graph with oriented edges labelled by X.
Here the graph is embedded in a disc and the embedding is taken up to isotopy. The benefit of
this approach is that these pictures can be drawn on a blackboard (or similar). However this
approach gives the misleading idea that this is a branch of topology and, more importantly,
is not suitable for a computer implementation.
Here we represent a web by a combinatorial structure (essentially a ribbon graph, constellation)
but we need to allow the web to have boundary points. Our basic operations are rotate and glue.
It is clear that this is equivalent to rotate, join and stitch.

### Operations

Since a spider is a set with operations and we have defined these operations for webs,
it is clear that we have constructed a spider (assuming we have checked that the operations
satisfy the axioms). However we have not explained the sense in which this is a free spider.

One general theory of algebraic theories is given by monads. In this setting we have
an underlying functor which has a left adjoint. In the theory of groups the underlying functor
gives a set and applying the left adjoint to a set gives a free group and we call the original
set the generators of the group.

This is one-dimensional is the sense that free objects are constructed in terms of words.
These are lists of generators and are written on a line. Spiders and webs are two-dimensional
since webs have to be drawn as two dimensional.

In our setting, the set of objects of a spherical category has more structure; namely
the rotation map. Here the set of objects is the free monoid on a set of generators
which is given by a list of generators. Then rotation is rotation of this list.
The underlying functor takes the set of webs, the set of words and the boundary map.
We fix the set of words with rotation. Then this underlying functor gives a set with
rotation together with an equivariant map to words. In order to construct the left adjoint
we start with an object in the underlying category. These are called generators. Then
we take all webs such that every vertex is a generator. Then to justify calling this an
algebraic theory, we need to check that this is the functor of a monad.

To see that this is a monad we change perspective on webs. We take the diagram of a web.
Then we can cut out a small disc around each vertex. This picture has no vertices but
is now a disc with a hole punched out for each vertex. Each of these pictures defines
an operation. Given a spider, we assign an object to each punched out hole so that the
boundaries match. Then we fill in each hole by attaching the object. For example,
glueing has two punched out holes corresponding to the two inputs of the glueing operation.
However we could also fill in a punched out hole by the picture of an operation.
This gives the picture of some operation. This means that these operations have the structure
of an operad. In fact this is a cyclic operad as defined in [2]_ and [4]_. Furthermore a spider
is a cyclic algebra for this cyclic operad. Taking the glueing operation as basic and
building up operations is the componential approach to cyclic operads given in [4]_.

### Implemented

This file has two main classes: the Element class :class:`SphericalWeb` and the Parent class
:class:`SphericalSpider`. Then :meth:`vertex` is the basic construction of a web. Then
:meth:`glue` has input two webs and an integer and output a web; and :meth:`polygon`
has input a list of webs and output a web. These operations build more complicated webs
starting with vertices. Once you have constructed a web you can see a picture using
:meth:`plot`. This takes the combinatorial data for a web and outputs a graphics object.

Related implementations are DiagramAlgebras and #25901

REFERENCES:

.. [1] Predrag Cvitanović
Group Theory: Birdtracks, Lie's, and Exceptional Groups
Princeton University Press, 2008
ISBN    0691118361, 9780691118369

.. [2] E Getzler, MM Kapranov
Cyclic operads and cyclic homology
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.146.2678&rep=rep1&type=pdf

.. [3] G.Kuperberg
Spiders for rank 2 Lie algebras
Commun.Math. Phys. 180, 109–151 (1996)
:arxiv:`q-alg/9712003`

.. [4] Jovana Obradovic
Monoid-like definitions of cyclic operad
Theory and Applications of Categories, Vol. 32, (2017), No. 12, pp 396-436.
http://www.tac.mta.ca/tac/volumes/32/12/32-12.pdf

.. [5] Adam S Sikora and Bruce W Westbury
Confluence theory for graphs
Algebraic & Geometric Topology,
Vol. 7, (2007), pp 439–478
:arxiv:`math/0609832`


.. [6] Bruce W. Westbury
Invariant tensors for the spin representation of so(7)
Mathematical Proceedings of the Cambridge Philosophical Society,
Vol. 144, (2008), pp 217-240
:arxiv:`math/0601209`

.. [7] Colin Scott Hagemeyer
Spiders and Generalized Confluence
Dissertation, UC Davis (2018)
:arxiv:`math/1809.10338`

AUTHORS:

- Bruce Westbury (2021): initial version

"""

#*****************************************************************************
#       Copyright (C) 2021 Bruce Westbury bruce.westbury@gmail.com
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.graphs.graph import Graph
from sage.combinat.baxter_permutations import BaxterPermutations
from sage.structure.richcmp import richcmp, op_EQ, op_NE
from copy import copy
from collections import namedtuple
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras import Algebras
from sage.rings.semirings.all import NN
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.groups.braid import BraidGroup
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.modules.free_module_element import vector
from sage.rings.ring import Algebra

Strand  = namedtuple('Strand', ['oriented','colour'], defaults=[0,'black'])

class halfedge():
    """
    The class of half edges in a surface graph.

    This should probably be an attribute either of SphericalWeb or SphericalSpider
    """
    def __init__(self, st=Strand(), crossing=False):
        """
        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: halfedge()
            <sage.combinat.spherical_spider.halfedge object at ...>
        """
        self.strand = st
        self.crossing = crossing

    def __hash__(self):
        r"""
        A hash function for :class:`halfedge`

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: halfedge().__hash__() # random
        """
        return hash((self.strand, self.crossing))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

        sage: from sage.combinat.spherical_spider import halfedge
        sage: halfedge()._repr_()
        '(0,black,False)'
        sage: halfedge()
        <sage.combinat.spherical_spider.halfedge object at ...>
        """
        return f"({self.strand.oriented},{self.strand.colour},{self.crossing})"

    def dual(self):
        """
        Construct the dual halfedge.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: halfedge().dual()
            <sage.combinat.spherical_spider.halfedge object at ...>
        """
        st = self.strand
        return halfedge(Strand(-st.oriented, st.colour),self.crossing)

class SphericalWeb(Element):
    r"""The class of webs.

    This consists of
    * a set of half-edges
    * a bijection `c` on the set of half-edges with no fixed points
    * an involution `e` on a subset of the half-edges

    The half-edges for which `e` is undefined are the boundary half-edges.
    This set has a total order.

    The only orbits of `c` of order two have both half-edges in the boundary
    or are a loop.
    """

    def __init__(self, c:  dict, e: dict, b: list, check=True):
        r"""
        Initialise an instance of :class:`SphericalWeb`.

        INPUT:

            * `c` a bijection of the set of half-edges
            * `e` a partial involution of the set of half-edges
            * `b` the ordered list of boundary edges

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: b = [halfedge(),halfedge()]
            sage: c = {b[0]:b[1], b[1]:b[0]}
            sage: SphericalWeb(c, dict([]), b)
            The spherical web with c = (1, 0), e = ()
             and edges ().
        """

        bd = [a.strand for a in b]
        parent = SphericalSpider(bd)
        Element.__init__(self, parent)

        self.cp = c
        self.e = e
        self.b = b
        self.normalize()
        if check:
            self.check()

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black')]*3).vertex().__copy__()
            The spherical web with c = (1, 2, 0), e = ()
             and edges ().
        """
        D = {a:halfedge(a.strand,a.crossing) for a in self.cp}
        c = {D[a]:D[self.cp[a]] for a in self.cp}
        e = {D[a]:D[self.e[a]] for a in self.e}
        b = [D[a] for a in self.b]
        return SphericalWeb(c, e, b)

    def check(self):
        r"""
        Implement the abstract method :meth:`check` of :class:`ClonableElement`.

        Check ``self`` is a valid web.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: a = halfedge()
            sage: SphericalWeb({a:a}, dict(), [a], check=False)
            The spherical web with c = (0,), e = ()
             and edges ().
            sage: SphericalWeb({a:a}, dict(), [a])
            Traceback (most recent call last):
            ...
            ValueError: the mapping c has at least one fixed point
        """
        c = self.cp
        e = self.e
        b = self.b
        h = set(c)
        if not all(isinstance(a,halfedge) for a in h):
            raise ValueError("every element must be a half-edge")
        if set(c.values()) != h:
            raise ValueError("the map c is required to be a bijection")
        if any(not e[e[a]] is a for a in e):
            raise ValueError("the map e must be an involution")
        if any(c[a] is a for a in c):
            raise ValueError("the mapping c has at least one fixed point")
        if any(e[a] is a for a in e):
            raise ValueError("the mapping e has at least one fixed point")
        if not set(e.keys()).issubset(h):
            raise ValueError("the domain of e must be a subset of the domain of c")
        if not set(b).issubset(h):
            raise ValueError("the boundary must be a subset of the domain of c")
        if not set(e.keys()).isdisjoint(set(b)):
            raise ValueError("the domain of e must not intersect the boundary")
        #for i,a in enumerate(b):
        #    u = a
        #    while c[u] in e:
        #        u = e[c[u]]
        #    j = b.index(c[u])
        #    if 0 < j < i:
        #       raise ValueError("boundary is inconsistent")

    def normalize(self):
        r"""
        This removes nearly all vertices of degree two.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: h = [halfedge() for i in range(4)]
            sage: c = {h[0]:h[1], h[1]:h[0], h[2]:h[3], h[3]:h[2]}
            sage: e = {h[1]:h[2],h[2]:h[1]}
            sage: b = [h[0],h[3]]
            sage: SphericalWeb(c, e, b) # indirect doctest
            The spherical web with c = (1, 0), e = ()
             and edges ().

        TESTS::

        This should not ever happen.

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black')]*2).vertex() # indirect doctest
            The spherical web with c = (1, 0), e = ()
             and edges ().

        Check loops are not removed.

            sage: SphericalSpider([]).loop(Strand()) # indirect doctest
            A closed spherical web with 1 edges.
        """
        flag = True
        while flag :
            flag = False
            c = self.cp
            e = self.e
            rm = [a for a in e if c[c[a]] == a]
            if len(rm) != 0:
                x = rm[0]
                z = e[x]
                y = c[x]
                if y in e and z != y:
                    flag = True
                    e[z] = e[y]
                    e[e[y]] = z
                    c.pop(x)
                    c.pop(y)
                    e.pop(x)
                    e.pop(y)
                elif y in self.b:
                    flag = True
                    c[y] = c[z]
                    w = [a for a in c if c[a] == z][0]
                    c[w] = y
                    c.pop(x)
                    c.pop(z)
                    e.pop(x)
                    e.pop(z)

    def _repr_(self):
        r"""
        Overload default implementation.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: SphericalSpider([]).polygon([u,u,u])._repr_()
            "The spherical web with c = ..., e = (6, 7, 8, 3, 4, 5)\n and edges ('(0,black,False)',..., '(0,black,False)')."
        """
        if len(self.b) > 0:
            cn, en, hn = self.canonical()
            return f"The spherical web with c = {cn}, e = {en}\n and edges {hn}."
        else:
            return f"A closed spherical web with {int(len(self.e)/2)} edges."

    @cached_method
    def canonical(self):
        r"""
        A canonical labelling of the elements of ``self``.

        This returns two lists of integers and a list of strands such that ``self``
        can be recovered, up to isomorphism, from these two sequences.

        Let ``self`` have `n` elements and `k` boundary elements.
        Then the first list is a bijection on [0,1,...,n-1] and
        the second list is an involution on [k,k+1,...,n-1].

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: cn, en, hn = SphericalSpider([]).polygon([u,u,u,u]).canonical()
            sage: cn
            (4, 6, 8, 10, 5, 0, 7, 1, 9, 2, 11, 3)
            sage: en
            (7, 10, 9, 4, 11, 6, 5, 8)
            sage: hn
            ('(0,black,False)',
             ...
             '(0,black,False)')
        """
        b = self.b
        c = self.cp
        e = self.e
        k = len(b)
        Dp = {a:i for i,a in enumerate(b)}
        Dp.update({a:(k+i) for i,a in enumerate(self._traversal(b))})
        cn = tuple([Dp[c[a]] for a in b]+[Dp[c[a]] for a in self._traversal(b)])
        en = tuple([Dp[e[a]] for a in self._traversal(b)])
        hn = tuple([a._repr_() for a in self._traversal(b)])
        return cn, en, hn

    def __hash__(self):
        r"""
        Overload the :meth:`__hash__`.

        This is needed to put a :class:`SphericalWeb` into a :class:`set`
        or to use it as a key in a :class:`dict'.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: v = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: v.__hash__()  # random

            sage: hash(v) # random

            sage: u = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: set([u,v]) # indirect doctest
            {The spherical web with c = (1, 2, 0), e = ()
             and edges ().}
        """
        return hash((self.parent(),*self.canonical()))

    def _richcmp_(self, other, op):
        """
        Overload :meth:`__eq__` and :meth:`__ne__`.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black')]*4).vertex()
            sage: v = SphericalSpider([Strand(0,'black')]*4).vertex()
            sage: u is v, u == v, u != v # indirect doctest
            (False, True, False)
            sage: u < v # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: '<' not supported between ... and 'SphericalWeb'

        There is a problem with this as :meth:`canonical` does not
        see anything not conected to the boundary.

            sage: S = SphericalSpider([])
            sage: S.loop(Strand()) == S.empty()
            True
        """
        if op == op_EQ or op == op_NE:
            return richcmp(self.canonical(), other.canonical(), op)
        else:
            return NotImplemented

#### End of underscore methods ####

#### Start of methods for working with webs ####

    def _traversal(self, initial):
        """
        A generator for the elements of ``self`` connected to the
        elements in ``initial``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: w = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: w._traversal(w.b)
            <generator object SphericalWeb._traversal at ...>
            sage: w._traversal(w.b[0])
            <generator object SphericalWeb._traversal at ...>
        """
        if isinstance(initial,halfedge):
            initial = tuple([initial])
        else:
            initial = tuple(initial)
        if not set(initial).issubset(self.cp.keys()):
            raise ValueError("initial must be a subset of the set of elements")

        c = self.cp
        e = self.e
        #b = self.b

        visited = list(initial)
        new = list()
        flag = True
        while flag:
            flag = False
            for u in visited:
                v = c[u]
                while v not in visited:
                    new.append(v)
                    yield v
                    flag = True
                    v = c[v]
            visited += new
            for u in visited:
                if u in e:
                    v = e[u]
                    if v not in visited:
                        new.append(v)
                        yield v
                        flag = True
            visited += new
            new = list()
        return

    @staticmethod
    def _stitch(c, e, x, y):
        """
        Connect `x` and `y`.

        This is a low level method that is not intended to be called directly.
        It was written to comply with the DRY principle.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: c,e = u._stitch(u.cp,u.e,u.b[0],u.b[-1])
            sage: len(c),len(e)
            (5, 4)
        """
        st = x.strand
        if y.strand != Strand(-st.oriented, st.colour):
            raise ValueError(f"{x.strand} and {y.strand} must be dual")

        u = halfedge(y.strand)
        v = halfedge(x.strand)
        c[u] = v
        c[v] = u
        e[x] = u
        e[u] = x
        e[y] = v
        e[v] = y

        return c, e

    def rotate(self, k):
        r"""Rotate the boundary anticlockwise `k` steps.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*5).vertex()
            sage: u.rotate(3)
            The spherical web with c = (1, 2, 3, 4, 0), e = ()
             and edges ().
            sage: u.rotate(-1)
            The spherical web with c = (1, 2, 3, 4, 0), e = ()
             and edges ().
        """
        result = self.__copy__()
        b = result.b
        nb = b[k:]+b[:k]
        return SphericalWeb(result.cp, result.e, nb)

    def glue(self, other, n):
        r"""Glue two ribbon graphs together.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: v = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: u.glue(v,1)
            The spherical web with c = (1, 4, 3, 5, 0, 2), e = (5, 4)
             and edges ('(0,black,False)', '(0,black,False)').
            sage: u.glue(v,0)
            The spherical web with c = (1, 2, 0, 4, 5, 3), e = ()
             and edges ().
            sage: u.glue(v,-1)
            Traceback (most recent call last):
            ...
            ValueError: n=-1 cannot be negative
            sage: u.glue(v,4)
            Traceback (most recent call last):
            ...
            ValueError: n=4 is too large
        """
        if n < 0:
            raise ValueError(f"n={n} cannot be negative")

        if n > len(self.b) or n > len(other.b):
            raise ValueError(f"n={n} is too large")

        ns = self.__copy__()
        no = other.__copy__()

        bs = ns.b
        bo = no.b
        if n == 0:
            b = bs+bo
        else:
            b = bs[:-n]+bo[n:]

        c = {**ns.cp, **no.cp}
        e = {**ns.e, **no.e}

        for x,y in zip(reversed(bs[-n:]),bo[:n]):
            c, e =  self._stitch(c,e,x,y)

        return SphericalWeb(c, e, b)

    def mirror_image(self):
        r"""
        Construct the mirror image of ``self``.

        This should be called dual and is not correct.
        It is not currently called.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: u.mirror_image()
            The spherical web with c = (1, 2, 0), e = ()
             and edges ().
            sage: v = SphericalSpider([Strand(0,'black')]*4).vertex()
            sage: u.glue(v,1).mirror_image()
            The spherical web with c = (1, 2, 5, 4, 6, 0, 3), e = (6, 5)
             and edges ('(0,black,False)', '(0,black,False)').
            sage: w = u.glue(u,1).glue(v,1)
            sage: w == w.mirror_image()
            False
        """
        D =  {a:halfedge() for a in self.cp}
        cn = {D[self.cp[a]]:D[a] for a in D}
        en = {D[a]:D[self.e[a]] for a in self.e}
        bn = tuple(reversed([D[a] for a in self.b]))
        return SphericalWeb(cn, en, bn)

    def vertices(self):
        """
        Find the vertices of ``self``.

        These are the orbits of `c`.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: t = SphericalSpider([Strand()]*3).vertex()
            sage: u = SphericalSpider([]).polygon([t]*4)
            sage: [len(a) for a in u.vertices()]
            [3, 3, 3, 3]
        """
        c = self.cp
        he = set(c.keys())
        result = set()
        while len(he) != 0:
            a = he.pop()
            vertex = [a]
            b = c[a]
            while b != a:
                vertex.append(b)
                he.discard(b)
                b = c[b]
            result.add(tuple(vertex))
        return result

    def faces(self):
        """
        Find the faces of ``self``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: len(u.faces())
            3
            sage: w = SphericalSpider([Strand(0,'black')]*4).vertex()
            sage: len(w.faces())
            4
            sage: len(u.glue(u,0).faces())
            6
        """
        c = self.cp
        e = self.e
        he = set(c.keys())
        result = set()

        # First find the external faces.
        for a in self.b:
            u = a
            face = [a]
            he.discard(a)
            while c[u] in e:
                u = e[c[u]]
                face.append(u)
                he.discard(u)
            result.add(tuple(face))

        # Now find the internal faces.
        while len(he) != 0:
            a = he.pop()
            face = [a]
            u = e[c[a]]
            while u != a:
                face.append(u)
                he.discard(u)
                u = e[c[u]]
            result.add(tuple(face))

        return result

    def is_closed(self):
        """
        Return ``True`` if ``self`` is closed.

        Note that this means that the boundary is empty.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand()]*3).vertex().is_closed()
            False
            sage: SphericalSpider([]).loop(Strand()).is_closed()
            True
        """
        return len(self.b) == 0

    def is_connected(self):
        """
        Return ``True`` if ``self`` is connected.

        Note that this means that the diagram including the boundary
        is connected.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand()]*3).vertex().is_connected()
            True
            sage: SphericalSpider([]).loop(Strand()).is_connected()
            False
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: v = SphericalSpider([Strand()]*3).vertex()
            sage: u.glue(v,0).is_connected()
            True
        """
        return len(self.cp) == len(self.canonical()[0])

    def components(self):
        """
        Return the closed components of ``self``.

        This is the complement of the connected component of the boundary.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: u.components()
            (The spherical web with c = (1, 2, 0), e = ()
             and edges ().,
            A closed spherical web with 0 edges.)
            sage: u.glue(u,0)
            The spherical web with c = (1, 2, 0, 4, 5, 3), e = ()
             and edges ().
            sage: u.glue(u,0).components()
            (The spherical web with c = (1, 2, 0, 4, 5, 3), e = ()
             and edges ().,
            A closed spherical web with 0 edges.)
        """
        Dn = {a:halfedge(a.strand) for a in self.b}
        for a in self._traversal(self.b[0]):
            Dn[a] = halfedge(a.strand)

        cn = {Dn[a]:Dn[self.cp[a]] for a in self.cp}
        en = {Dn[a]:Dn[self.e[a]] for a in self.e}
        bn = [Dn[a] for a in self.b]
        wb = SphericalWeb(cn, en, bn)

        Dc = {a:halfedge(a.strand) for a in self.cp if not a in Dn}
        cc = {Dc[a]:Dc[self.cp[a]] for a in Dc}
        ec = {Dc[a]:Dc[self.e[a]] for a in Dc}
        wc = SphericalWeb(cc, ec, [])

        return wb, wc

    def is_decomposable(self):
        """
        Return True if ``self`` is decomposable.

        A web `w` is decomposable if it can be written as `w = u.glue(v,0)`
        where `u` and `v` are non-empty.

        Note that this means that the diagram excluding the boundary
        is connected.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: u.is_decomposable()
            False
            sage: u.glue(u,0).is_decomposable()
            True
        """
        if len(self.b) == 0:
            raise ValueError("not implemented for a closed web")
        return len(self.cp) != len(list(self._traversal(self.b[0])))+1

    def is_separable(self):
        r"""
        Return ``True`` if ``self`` is separable.

        This means each face has distinct vertices.
        Including the boundary faces.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*2).vertex()
            sage: v = SphericalSpider([Strand()]*4).vertex()
            sage: v.glue(u,2).is_separable()
            True
        """
        from itertools import product
        for v,f in product(self.vertices(),self.faces()):
            if len(set(v).intersection(set(f))) > 1:
                return True
        return False

    def is_simple(self):
        """
        Return ``True`` if ``self`` is a simple graph.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*2).vertex()
            sage: v = SphericalSpider([Strand()]*4).vertex()
            sage: v.glue(u,2).is_simple()
            False
            sage: v.glue(v,2).is_simple()
            False
        """
        return all(len(x)>2 for x in self.faces())

    def to_graph(self):
        r"""
        Construct the graph of ``self``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: u.to_graph()
            Graph on 3 vertices
            sage: SphericalSpider([]).polygon([u]*3).to_graph()
            Graph on 9 vertices
        """
        c = self.cp
        e = self.e
        G = Graph({a:[c[a]] for a in c})
        for a in e:
            if not G.has_edge(e[a], a, "e"):
                G.add_edge(a, e[a], "e")
        return G

    def show(self):
        r"""Show the web ``self``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand()]*3).vertex().show()
            Graphics object consisting of 7 graphics primitives
        """
        return self.to_graph().plot(vertex_labels=False, edge_labels=True)

    def _layout(self):
        r"""
        Layout the planar map.

        This uses the barycentric layout. The boundary points are the vertices of
        a regular polygon. Then the vertices are positioned so that each vertex
        is at the centre of mass of its neighbours.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: len(SphericalSpider([]).polygon([u]*3)._layout())
            6

            sage: len(u.glue(u,1)._layout())
            5

        If the graph is not simple the diagram will degenerate.

            sage: v = SphericalSpider([Strand()]*4).vertex()
            sage: w = SphericalSpider([Strand()]*2).vertex()
            sage: v.glue(w,2)._layout()
            {((-1.00000000000000, 0.000000000000000), (0.0, 0.0)),
            ((1.00000000000000, 0.000000000000000), (0.0, 0.0))}


        If there are no boundary points only the boundary circle is drawn.

            sage: SphericalSpider([]).loop(Strand())._layout()
            set()
        """
        from sage.matrix.all import matrix
        from sage.functions.trig import sin, cos
        from sage.all import pi

        vt = list(self.vertices())
        nv = len(vt)
        d = len(self.b)
        M = matrix(nv,nv)
        for i,u in enumerate(vt):
            M[i,i] = len(u)
            x = set(self.e[a] for a in u if a in self.e)
            for j,v in enumerate(vt):
                if i != j:
                    M[i,j] = -len(x.intersection(set(v)))

        U = matrix(nv,1,0.0)
        for i,b in enumerate(self.b):
            x = cos(2*pi*i/d).n()
            for j,v in enumerate(vt):
                if b in v:
                    U[j,0] = U[j,0] + x

        V = matrix(nv,1,0.0)
        for i,b in enumerate(self.b):
            y = sin(2*pi*i/d).n()
            for j,v in enumerate(vt):
                if b in v:
                    V[j,0] = V[j,0] + y

        Mi = M.inverse()
        pos = [(r[0],s[0]) for r,s in zip(Mi*U, Mi*V)]

        result = set()
        for i,u in enumerate(vt):
            for j,b in enumerate(self.b):
                if self.cp[b] in u:
                    x = cos(2*pi*j/d).n()
                    y = sin(2*pi*j/d).n()
                    result.add(((x,y),pos[i],))
            x = set(self.e[a] for a in u if a in self.e)
            for j,v in enumerate(vt):
                if i < j:
                    if any(r in v for r in x):
                        result.add((pos[i],pos[j],))

        return result

    def plot(self):
        r"""
        Plot the planar map.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: SphericalSpider([]).polygon([u]*3).plot()
            Graphics object consisting of 7 graphics primitives

            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: u.glue(u,1).plot()
            Graphics object consisting of 6 graphics primitives

        If the graph is not simple the diagram will degenerate.

            sage: v = SphericalSpider([Strand()]*4).vertex()
            sage: w = SphericalSpider([Strand()]*2).vertex()
            sage: v.glue(w,2).plot()
            Graphics object consisting of 3 graphics primitives

        If there are no boundary points only the boundary circle is drawn.

            sage: SphericalSpider([]).loop(Strand()).plot()
            Graphics object consisting of 1 graphics primitive

        TODO::

        Add colour, direction, under crossing.
        """
        from sage.plot.circle import circle
        from sage.plot.line import line

        lines = self._layout()
        G = circle((0,0),1)
        for a in lines:
                G += line(a,thickness=2,rgbcolor=(1,0,0))
        G.set_aspect_ratio(1)
        G.axes(False)
        return G

    def _latex_(self):
        """
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: SphericalSpider([]).polygon([u]*3)._latex_()
            ...
            sage: u.glue(u,1)._latex_()
            ...

        If the graph is not simple the diagram will degenerate.

            sage: v = SphericalSpider([Strand()]*4).vertex()
            sage: w = SphericalSpider([Strand()]*2).vertex()
            sage: v.glue(w,2)._latex_()
            '\\draw (0,0) circle (1cm);\n\\draw (-1.00000000000000,0.0) -- (0.0,0.0);\n\\draw (1.00000000000000,0.0) -- (0.0,0.0);\n'

        If there are no boundary points only the boundary circle is drawn.

            sage: SphericalSpider([]).loop(Strand())._latex_()
            '\\draw (0,0) circle (1cm);\n'

        TODO::

            Add colour, direction, under crossing.
        """
        lines = self._layout()

        #result = "\\begin{tikzpicture}\n"
        result = "\\draw (0,0) circle (1cm);\n"
        for a in lines:
            result += "\\draw ({},{}) -- ({},{});\n".format(a[0][0],a[1][0],a[1][0],a[1][1])
        #result += "\\end{tikzpicture}\n"

        return result

    def to_snappy():
        r"""
        If SnapPy is installed (see https://snappy.math.uic.edu/installing.html)
        then construct the link from the web. For information on planar diagrams
        in SnapPy see https://snappy.math.uic.edu/spherogram.html
        """
        # The documentation in  mod:`sage.misc.package` says not to use this but to use
        # the framework in :mod:`sage.features` instead. I need some help with this.
        if not is_package_installed('snappy'):
            raise ValueError("This requires the optional package SnapPy.")

#### End of methods for working with webs ####

#### Start of methods for rewriting ####

    def search(self, h):
        r"""
        Find copies of h in ``self``

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: len(list(SphericalSpider([]).polygon([u]*4).search(u)))
            12

        TODO::

        This should be rewritten to use :meth:``_traversal``.
        """

        def test(x):
            flag = True
            while flag:
                flag = False
                newD = dict()
                for a in Dm:
                    while not (c[a] in Dm):
                        if a in Dm:
                            newD[c[a]] = self.cp[Dm[a]]
                        elif a in newD:
                            newD[c[a]] = self.cp[newD[a]]
                        a = c[a]
                for a in Dm:
                    if a in e:
                        if not (e[a] in Dm):
                            if not Dm[a] in self.e:
                                return False
                            newD[e[a]] = self.e[Dm[a]]
                if len(newD) != 0:
                    Dm.update(newD)
                    flag = True
                for a in Dm:
                    if c[a] in Dm:
                        if Dm[c[a]] != self.cp[Dm[a]]:
                            return False
                    if a in e:
                        if e[a] in Dm:
                            if not Dm[a] in self.e:
                                return False
                            if Dm[e[a]] != self.e[Dm[a]]:
                                return False
            return True

        c = h.cp
        e = h.e
        for x in self.cp:
            Dm = {h.b[0]:x}
            if test(x):
                assert [ a for a in c if not (a in Dm) ] == [], "Mapping is not fully defined."
                if len(set(Dm.values())) == len(Dm):
                    yield Dm

    def replace(self, k, D, h):
        r"""
        Replace image of map D:h -> ``self`` by k

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: h = SphericalSpider([Strand()]*3).vertex()
            sage: g = SphericalSpider([]).polygon([h]*3)
            sage: D = next(g.search(h))
            sage: g.replace(g,D,h)
            The spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2, 11, 13, 14, 10, 12, 9), e = (9, 10, 8, 11, 12, 5, 3, 4, 6, 7, 14, 13)
             and edges ('(0,black,False)', ..., '(0,black,False)').
        """
        parent = self.parent()
        if parent != k.parent():
            raise ValueError(f"the two parents {self.parent()} and {k.parent()} are different")
        if parent != h.parent():
            raise ValueError(f"the two parents {self.parent()} and {h.parent()} are different")

        if len(h.b) != len(k.b):
            raise ValueError(f"boundaries of {k} and {h} must have the same length")
        if any( u.dual().strand != v.strand for u,v in zip(h.b, k.b)):
            raise ValueError(f"boundaries of {k} and {h} must match")

        Ds = {a:halfedge(a.strand,a.crossing) for a in self.cp if not a in D.values()}
        Dk = {a:halfedge(a.strand,a.crossing) for a in k.cp}
        c = {Ds[a]:Ds[self.cp[a]] for a in Ds}
        c.update({Dk[a]:Dk[k.cp[a]] for a in Dk})

        e = {Ds[a]:Ds[self.e[a]] for a in Ds if a in self.e and self.e[a] in Ds}
        e.update({Dk[a]:Dk[k.e[a]] for a in k.e})

        Db = {x:y for x,y in zip(h.b,k.b)}
        b = [None]*len(self.b)
        for i,a in enumerate(self.b):
            if a in Ds:
                b[i] = Ds[a]

        for a in h.b:
            if D[a] in self.e:
                x = Ds[self.e[D[a]]]
                y = Dk[Db[a]]
                c, e = self._stitch(c,e,x,y)
            else:
                i = self.b.index(D[a])
                b[i] = Dk[Db[a]]

        return SphericalWeb(c, e, b)

    def replace_linear(self, k, D, h):
        r"""
        Replace image of map D:h -> ``self`` by the linear combination k

        EXAMPLES::

        """
        if h.parent() != k.parent().basis().keys():
            raise ValueError("boundaries are different")

        R = k.parent().base()
        L = FreeSphericalSpider(R, self.parent())

        mc = k.monomial_coefficients()

        return L.sum_of_terms([(self.replace(a, D, h),mc[a]) for a in mc])

    def apply_rule(self, term, replacement):
        r"""
        Find ``term`` in ``self`` and replace with the linear combination ``replacement``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: S = SphericalSpider([])
            sage: term = S.loop(Strand())
            sage: delta = PolynomialRing(ZZ, 'delta').gen()
            sage: L = FreeSphericalSpider(delta.parent(),[])
            sage: replacement = delta * L(S.empty())
            sage: S.loop(Strand()).apply_rule(S.loop(Strand()), replacement)

        """
        try:
            D = next(self.search(term))
            return self.replace_linear(replacement, D, term)
        except StopIteration:
            return self

#### End  of methods for rewriting ####

#### Start of Parent ####

class SphericalSpider(UniqueRepresentation, Parent):
    r"""
    The Parent class for SphericalWeb.

    EXAMPLES::

        sage: from sage.combinat.spherical_spider import Strand
        sage: SphericalSpider([Strand(0,'black'),Strand(0,'black')])
        The spherical spider with boundary [(0, 'black'), (0, 'black')]
        sage: SphericalSpider([])
        The spherical spider with boundary []
    """
    @staticmethod
    def __classcall__(cls, boundary):
        if boundary in NN:
            boundary = [Strand()]*boundary

        return super(SphericalSpider, cls).__classcall__(cls, tuple(boundary))

    def __init__(self, boundary):
        r"""
        Initialise an instance of this class.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black'),Strand(0,'black')])
            The spherical spider with boundary [(0, 'black'), (0, 'black')]
            sage: SphericalSpider(2) == SphericalSpider([Strand(0,'black'),Strand(0,'black')])
            True
        """

        self.boundary = boundary

        Parent.__init__(self)

    def _repr_(self):
        r"""
        Overload the default method.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: P = SphericalSpider([Strand(0,'black')]*3)
            sage: P._repr_()
            "The spherical spider with boundary [(0, 'black'), (0, 'black'), (0, 'black')]"
        """
        return f"The spherical spider with boundary {[(a.oriented,a.colour) for a in self.boundary]}"

    Element = SphericalWeb

    def vertex(self):
        r"""
        Construct a single vertex.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black')]*4).vertex()
            The spherical web with c = (1, 2, 3, 0), e = ()
             and edges ().
        """
        bd = self.boundary
        n = len(bd)
        if n<2:
            raise ValueError(f"n={n} must be at least 2")

        b = [None]*n
        for i, a in enumerate(bd):
            b[i] = halfedge(a)
        c = {b[i-1]:b[i] for i in range(n)}
        e = dict([])
        return SphericalWeb(c,e,b)

    @staticmethod
    def crossing():
        r"""
        Construct an unoriented crossing.

        EXAMPLES::

            sage: s = SphericalSpider([]).crossing(); s
            The spherical web with c = (1, 2, 3, 0), e = ()
             and edges ().
            sage: s.show()
            Graphics object consisting of 9 graphics primitives
        """
        b = [None]*4
        for i in range(2):
            b[2*i] = halfedge(Strand(), False)
            b[2*i+1] = halfedge(Strand(), True)
        c = {b[i-1]:b[i] for i in range(4)}

        return SphericalWeb(c, dict([]), b)

    @staticmethod
    def loop(st: Strand):
        r"""
        Construct a loop.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([]).loop(Strand(0,'black'))
            A closed spherical web with 1 edges.
        """

        h = [halfedge(st),halfedge(st)]
        c = {h[0]:h[1], h[1]:h[0]}
        e = {h[0]:h[1], h[1]:h[0]}
        return SphericalWeb(c,e,[])

    @staticmethod
    def empty():
        """
        Construct the empty diagram.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([]).empty()
            A closed spherical web with 0 edges.
        """

        return SphericalWeb({},{},[])

    @staticmethod
    def polygon(corners):
        """
        Construct a polygon from a list of webs.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black')]*3).vertex()
            sage: SphericalSpider([]).polygon([u,u,u])
            The spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2), e = (6, 7, 8, 3, 4, 5)
             and edges ('(0,black,False)', ..., '(0,black,False)').
            sage: SphericalSpider([]).polygon([])
            A closed spherical web with 0 edges.
        """
        from functools import reduce

        if len(corners) == 0:
            return SphericalSpider(tuple([])).empty()

        # Avoid duplicates.
        cn = copy(corners)
        for i,u in enumerate(corners):
            for v in corners[:i]:
                if u == v:
                    cn[i] = copy(u)
        corners = cn

        c = reduce(lambda r, s: {**r, **s}, [a.cp for a in corners])
        e = reduce(lambda r, s: {**r, **s}, [a.e for a in corners])

        for u,v in zip(corners, corners[1:]):
            x = u.b[-1]
            y = v.b[0]
            c,e = SphericalWeb._stitch(c,e,x,y)

        x = corners[-1].b[-1]
        y = corners[0].b[0]
        c,e = SphericalWeb._stitch(c,e,x,y)

        b = sum([list(a.b[1:-1]) for a in corners],[])
        return SphericalWeb(c, e, b)

    @staticmethod
    def trefoil():
        r"""
        The trefoil as an unoriented long knot.

        EXAMPLES::

            sage: SphericalSpider([]).trefoil()
            The spherical web with c = (2, 5, 3, 4, 0, 6, 7, 1, 11, 8, 9, 10), e = (7, 8, 9, 10, 11, 2, 3, 4, 5, 6)
             and edges ('(0,black,True)', '(0,black,False)',  ..., '(0,black,True)', '(0,black,False)').
            sage: SphericalSpider([]).trefoil().show()
            Graphics object consisting of 35 graphics primitives
        """
        s = SphericalSpider([]).crossing()
        #s.b[0].crossing = True
        #s.b[2].crossing = True

        result = SphericalSpider([]).polygon([s,s,s]).rotate(1)
        cap = SphericalSpider([Strand()]*2).vertex()
        return result.glue(cap.glue(cap,0),4)

    @staticmethod
    def from_permutation(pi, baxter=True):
        r"""
        Construct a planar map from a two stack sorted permutation.

        This implements the algorithm in :arxiv:`math/0805.4180`.
        This algorithm is designed to apply to Baxter permutations.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: S = SphericalSpider([])
            sage: pi = Permutation([5,3,4,9,7,8,10,6,1,2])
            sage: S.from_permutation(pi)
            The spherical web with c = (1, 3, 6, 4, 5, 0, 7, 2, 10, 8, 9), e = (7, 8, 9, 10, 3, 4, 5, 6)
             and edges ('(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)').
            sage: pi = Permutation([2,4,1,3])
            sage: S.from_permutation(pi)
            Traceback (most recent call last):
            ...
            ValueError: [2, 4, 1, 3] is not a Baxter permutation
            sage: S.from_permutation(pi,baxter=False)
            The spherical web with c = (1, 0), e = ()
             and edges ().
        """
        if baxter:
            if not pi in BaxterPermutations():
                raise ValueError(f"{pi} is not a Baxter permutation")
        black = set((i+1,a) for i, a in enumerate(pi))
        n = len(pi)
        white = set([(1/2,1/2),(n+1/2,n+1/2)])
        ascents = [i+1 for i,a in enumerate(zip(pi,pi[1:])) if a[0] < a[1]]
        for a in ascents:
            l = max(pi[i] for i in range(a) if pi[i] < pi[a])
            white.add((a+1/2,l+1/2))

        Dup = {}
        Ddn = {}
        e = {}
        for a in black:
            w = [c for c in white if c[0]>a[0] and c[1]>a[1]]
            w.sort(key=lambda a: a[0]+a[1])
            r = halfedge()
            Dup[(w[0],a)] = r
            w = [c for c in white if c[0]<a[0] and c[1]<a[1]]
            w.sort(key=lambda a: a[0]+a[1])
            s = halfedge()
            Ddn[(w[-1],a)] = s
            e[r] = s
            e[s] = r

        c = {}
        for g in white:
            inco = [a for (w,a) in Ddn if w == g]
            inco.sort(key=lambda t: t[1]-t[0])
            outg = [a for (w,a) in Dup if w == g]
            outg.sort(key=lambda t: t[0]-t[1])
            #inco.reverse()
            for x,y in zip(inco,inco[1:]):
                c[Ddn[(g,x)]] = Ddn[(g,y)]
            for x,y in zip(outg,outg[1:]):
                c[Dup[(g,x)]] = Dup[(g,y)]
            if len(inco) > 0 and len(outg) > 0:
                c[Ddn[(g,inco[-1])]] = Dup[(g,outg[0])]
                c[Dup[(g,outg[-1])]] = Ddn[(g,inco[0])]
            elif len(inco) == 0 and len(outg) > 0:
                c[Dup[(g,outg[-1])]] = Dup[(g,outg[0])]
            elif len(inco) > 0 and len(outg) == 0:
                c[Ddn[(g,inco[-1])]] = Ddn[(g,inco[0])]
            else:
                raise RuntimeError("this can't happen")

        g = (1/2,1/2)
        inco = [a for (w,a) in Ddn if w == g]
        inco.sort(key=lambda t: t[0]-t[1])
        b =  [e[Ddn[(g,x)]] for x in inco]
        for a in inco:
            x = Ddn[(g,a)]
            e.pop(e[x])
            #c.pop(e[x])
            e.pop(x)
            c.pop(x)

        return SphericalWeb(c,e,b)

    @staticmethod
    def from_Gauss_code(G):
        r"""
        Construct the shadow diagram from a Gauss code.

        EXAMPLES::

            sage: SphericalSpider([]).from_Gauss_code([1,2,3,1,2,3])
            A closed spherical web with 6 edges.

            sage: G = [1, 2, 4, 5, 8, 1, 3, 4, 6, 7, 2, 3, 5, 6, 7, 8]
            sage: SphericalSpider([]).from_Gauss_code(G)
            A closed spherical web with 16 edges.
        """
        from sage.combinat.spherical_spider import halfedge
        from sage.knots.gauss_code import recover_orientations
        changed, positive, negative, ori = recover_orientations(G)

        n = len(G)
        h = [None]*(2*n)
        for i in range(2*n):
            h[i] = halfedge()

        e = {h[2*i-1]: h[2*i] for i in range(n)}
        e.update({h[2*i]: h[2*i-1] for i in range(n)})

        c = dict([])
        for r, s in positive:
            c.update({h[2*r]:h[2*s+1], h[2*s+1]:h[2*s], h[2*s]:h[2*r+1], h[2*r+1]:h[2*r]})

        for r, s in negative:
            c.update({h[2*r]:h[2*r+1], h[2*r+1]:h[2*s], h[2*s]:h[2*s+1], h[2*s+1]:h[2*r]})

        return SphericalWeb(c,e,[])

    @staticmethod
    def from_snappy():
        r"""
        If SnapPy is installed (see https://snappy.math.uic.edu/installing.html)
        then construct the web from the link. For information on planar diagrams
        in SnapPy see https://snappy.math.uic.edu/spherogram.html

        EXAMPLES::

            sage: SphericalSpider([]).from_snappy()
            Traceback (most recent call last):
            ...
            ValueError: This requires the optional package SnapPy.
        """
        # The documentation in  mod:`sage.misc.package` says not to use this but to use
        # the framework in :mod:`sage.features` instead. I need some help with this.
        from sage.misc.package import is_package_installed

        if not is_package_installed('snappy'):
            raise ValueError("This requires the optional package SnapPy.")

#### End of Parent ####

class FreeSphericalSpider(CombinatorialFreeModule):
    r"""
    Linear combinations of spherical webs.
    """
    @staticmethod
    def __classcall__(cls, base, boundary):
        if boundary in NN:
            boundary = [Strand()]*boundary

        return super(FreeSphericalSpider, cls).__classcall__(cls, base, tuple(boundary))

    def __init__(self, base, boundary):
        r"""
        Initialise ``self``.`

        EXAMPLES::FreeSphericalSpider

            sage: FreeSphericalSpider(QQ,[])
            Free module generated by The spherical spider with boundary [] over Rational Field
            sage: FreeSphericalSpider(QQ,[]).an_element()
            0
            sage: F = FreeSphericalSpider(QQ,[])
            sage: isinstance(0, F.element_class)
            False
            sage: isinstance(F(0), F.element_class)
            True
        """
        CombinatorialFreeModule.__init__(self, base, SphericalSpider(boundary))

    def boundary(self):
        return self.basis().keys().boundary

    def vertex(self):
        r"""
        Return the vertex on the boundary as a monomial.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: FreeSphericalSpider(QQ,[Strand()]*3).vertex()
            B[The spherical web with c = (1, 2, 0), e = ()
             and edges ().]
            """
        v = self.basis().keys().vertex()
        return self.monomial(v)

    def rotate(self, k):
        r"""
        Extend :meth:'rotate' by linearity

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: L = FreeSphericalSpider(QQ,[Strand()]*3)
            sage: rot = L.rotate(1); rot
            Generic endomorphism of Free module generated by The spherical spider ...
            sage: rot(L.vertex())
            B[The spherical web with c = (1, 2, 0), e = ()
             and edges ().]
        """
        b = self.boundary()
        codomain = FreeSphericalSpider(self.base(), b[k:]+b[:k])
        on_basis = lambda x : codomain(x.rotate(k))
        return self.module_morphism(codomain=codomain, on_basis=on_basis)

    class Element(CombinatorialFreeModule.Element):
        r"""
        A class for elements of FreeSphericalSpider.

            TESTS::

        """

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            EXAMPLES::
                sage: from sage.combinat.spherical_spider import SL2_relations
                sage: (SL2_relations()[0][1])._latex_()
                '\\delta\\draw (0,0) circle (1cm);\n'
            """
            mc = self.monomial_coefficients()
            return ''.join([(mc[a])._latex_()+a._latex_() for a in mc])

        def rotate(self, k):
            r"""
            Extend :meth:'rotate' by linearity

            EXAMPLES::

                sage: from sage.combinat.spherical_spider import Strand
                sage: L = FreeSphericalSpider(QQ,[Strand()]*3)
                sage: L.vertex().rotate(1)
                B[The spherical web with c = (1, 2, 0), e = ()
                 and edges ().]
            """
            b = self.parent().boundary()
            L = FreeSphericalSpider(self.parent().base(), b[k:]+b[:k])
            mc = self.monomial_coefficients()
            return L.sum_of_terms([(a.rotate(k), mc[a]) for a in mc])

        def _glue_right(self, diagram, k):
            r"""
            Extend :method`glue` by linearity.

            """
            bs = self.boundary
            bo = diagram.boundary
            if k == 0:
                bd = bs+bo
            else:
                bd = bs[:-k]+bo[k:]
            codomain = FreeSphericalSpider(self.base(), bd)
            on_basis = lambda x : codomain(x.glue(diagram, k))
            return self.module_morphism(codomain=codomain, on_basis=on_basis)

        def glue(self, other, k):
            r"""
            Extend :meth:`glue` by bilinearity.

            """
            from itertools import product
            bs = self.parent().boundary
            bo = other.parent().boundary
            if k == 0:
                bd = bs+bo
            else:
                bd = bs[:-k]+bo[k:]
            L = FreeSphericalSpider(self.parent().base(), bd)
            ms = self.monomial_coefficients()
            mo = other.monomial_coefficients()
            return L.sum((x.glue(y, k), ms[x]*mo[y]) for x,y in product(ms, mo))

        def apply_rule(self, term, replacement):
            r"""
            Extend :method'apply_rule' by linearity

            EXAMPLES::

            """
            codomain = self
            on_basis = lambda x : codomain(x.apply_rule(term, replacement))
            return self.module_morphism(codomain=codomain, on_basis=on_basis)

        def simplify(self, term, replacement):
            r"""
            Simplify by repeatedly applying :meth:'apply_rule'

            EXAMPLES::

            """
            new = self
            finished = False
            while not finished:
                old = new
                new = old.apply_rule(term, replacement)
                finished = new == old

            return old

class LinearSphericalSpider(FreeSphericalSpider):
    r"""
    Linear combinations of spherical webs.
    """
    @staticmethod
    def __classcall__(cls, base, boundary, rewrite_rules=None):
        if boundary in NN:
            boundary = [Strand()]*boundary

        return super(FreeSphericalSpider, cls).__classcall__(cls, base, tuple(boundary), rewrite_rules)

    def __init__(self, base, boundary, rewrite_rules):
        r"""
        Initialise ``self``.`

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import SL2_relations
            sage: LinearSphericalSpider(QQ, [], SL2_relations())
            Free module generated by The spherical spider with boundary [] over Rational Field
            sage: LinearSphericalSpider(QQ,[]).an_element()
            0
            sage: F = LinearSphericalSpider(QQ,[])
            sage: isinstance(F(0), F.element_class)
            True
        """
        FreeSphericalSpider.__init__(self, base, boundary)
        self.rewrite_rules = rewrite_rules

    def normal_form(self):
        r"""
        Reduce to normal form.

        """

    def glue(self, other, k):
        return lambda self, other : super(self).glue(other, k).normal_form()
        r"""
        The version of glue which  applies the rewrite rules.

        """

def SL2_relations(delta=None):
    r"""

    The quantum `SL(2)` skein relations.

    The only relation replaces a loop by `\delta`.

    EXAMPLES::

        sage: sage.combinat.spherical_spider.SL2_relations()
        ((A closed spherical web with 1 edges.,
         delta*B[A closed spherical web with 0 edges.]),)
    """

    if delta == None:
        delta = PolynomialRing(ZZ, 'delta').gen()

    S = SphericalSpider([])
    L = FreeSphericalSpider(delta.parent(), [])
    return tuple([tuple([S.loop(Strand()), delta*L(S.empty())])])

def A2_relations(delta=None):
    """
    The quantum `A_2` skein relations.

    EXAMPLES::

        sage: sage.combinat.spherical_spider.A2_relations()
        ((A closed spherical web with 1 edges.,
        (delta^2-1)*B[A closed spherical web with 0 edges.]),
        (A closed spherical web with 1 edges.,
        (delta^2-1)*B[A closed spherical web with 0 edges.]),
        (The spherical web with c = (2, 4, 3, 0, 5, 1), e = (5, 4, 3, 2)
        and edges ('(1,black,False)', '(1,black,False)', '(-1,black,False)', '(-1,black,False)').,
        -delta*B[The spherical web with c = (1, 0), e = ()
        and edges ().]),
        (The spherical web with c = (4, 6, 8, 10, 5, 0, 7, 1, 9, 2, 11, 3), e = (7, 10, 9, 4, 11, 6, 5, 8)
        and edges ('(-1,black,False)', ..., '(1,black,False)').,
        B[The spherical web with c = (1, 0, 3, 2), e = ()
        and edges ().] + B[The spherical web with c = (3, 2, 1, 0), e = ()
        and edges ().]))
    """
    if delta == None:
        delta = PolynomialRing(ZZ, 'delta').gen()

    sp = Strand(oriented=1, colour='black')
    sm = Strand(oriented=-1, colour='black')
    vm = SphericalSpider([sm]*3).vertex()
    vp = SphericalSpider([sp]*3).vertex()

    relations = []

    # Loop relations
    S = SphericalSpider([])
    L = FreeSphericalSpider(delta.parent(), [])
    relations.append(tuple([S.loop(sp), (delta**2-1)*L(S.empty())]))
    relations.append(tuple([S.loop(sm), (delta**2-1)*L(S.empty())]))

    # Digon relation
    S = SphericalSpider([sp, sm])
    L = FreeSphericalSpider(delta.parent(), [sp, sm])
    relations.append(tuple([vp.glue(vm,2), -delta*L(SphericalSpider([sp,sm]).vertex())]))

    # Square relation
    S = SphericalSpider([sm, sp, sm, sp])
    square = S.polygon([vm, vp, vm, vp])
    edgep = SphericalSpider([sm,sp]).vertex()
    U = edgep.glue(edgep, 0)
    edgem = SphericalSpider([sp,sm]).vertex()
    I = edgem.glue(edgem, 0).rotate(1)
    L = FreeSphericalSpider(delta.parent(), [sm, sp, sm, sp])
    relations.append(tuple([square, L(U)+L(I)]))

    return tuple(relations)

def B2_relations(delta=None):
    """
    The quantum `B_2` skein relations.

    EXAMPLES::

        sage: sage.combinat.spherical_spider.B2_relations()
        ((A closed spherical web with 1 edges.,
        (delta^2+delta-2)*B[A closed spherical web with 0 edges.]),
        (A closed spherical web with 1 edges.,
        (delta^3-2*delta+1)*B[A closed spherical web with 0 edges.]),
        (The spherical web with c = (1, 2, 0), e = (2, 1)
        and edges ('(0,black,False)', '(0,black,False)').,
        0),
        (The spherical web with c = (2, 4, 3, 0, 5, 1), e = (5, 4, 3, 2)
        and edges ('(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)').,
        (-delta-2)*B[The spherical web with c = (1, 0), e = ()
        and edges ().]),
        (The spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2), e = (6, 7, 8, 3, 4, 5)
        and edges ('(0,green,False)', ..., '(0,green,False)').,
        0))
    """
    if delta == None:
        delta = PolynomialRing(ZZ, 'delta').gen()

    sy = Strand(oriented=0, colour='black')
    so = Strand(oriented=0, colour='green')
    u = SphericalSpider([so,sy,sy]).vertex()
    v = SphericalSpider([sy,sy,so]).vertex()
    K = v.glue(u, 1)
    H = K.rotate(1)

    relations = []

    # Loop relations
    S = SphericalSpider([])
    L = FreeSphericalSpider(delta.parent(), [])
    relations.append(tuple([S.loop(sy), (delta**2+delta-2)*L(S.empty())]))
    relations.append(tuple([S.loop(so), (delta**3-2*delta+1)*L(S.empty())]))

    # Tadpole relation
    S = SphericalSpider([so])
    L = FreeSphericalSpider(delta.parent(), [so])
    tadpole = u.glue(SphericalSpider([sy,sy]).vertex(),2)
    relations.append(tuple([tadpole, 0*L(tadpole)]))

    # Digon relation
    S = SphericalSpider([so,so])
    L = FreeSphericalSpider(delta.parent(), [so,so])
    relations.append(tuple([u.glue(v,2),-(delta+2)*L(S.vertex())]))

    # Triangle relation
    S = SphericalSpider([so,so,so])
    L = FreeSphericalSpider(delta.parent(), [so,so,so])
    v = SphericalSpider([so,so,so]).vertex()
    triangle = S.polygon([v]*3)
    relations.append(tuple([triangle, 0*L(triangle)]))

    return tuple(relations)

def G2_relations(delta=None):
    """
    The quantum `G_2` skein relations.

    EXAMPLES::

        sage: sage.combinat.spherical_spider.G2_relations()
        ((A closed spherical web with 1 edges.,
        (delta^5+delta^4-5*delta^3-4*delta^2+6*delta+3)*B[A closed spherical web with 0 edges.]),
        (The spherical web with c = (1, 2, 0), e = (2, 1)
        and edges ('(0,black,False)', '(0,black,False)').,
        0),
        (The spherical web with c = (2, 4, 3, 0, 5, 1), e = (5, 4, 3, 2)
        and edges ('(0,black,False)', '(0,black,False)', '(0,black,False)', '(0,black,False)').,
        (-delta^3-delta^2+2*delta+2)*B[The spherical web with c = (1, 0), e = ()
        and edges ().]),
        (The spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2), e = (6, 7, 8, 3, 4, 5)
        and edges ('(0,black,False)', '(0,black,False)', ..., '(0,black,False)').,
        (delta^2-1)*B[The spherical web with c = (1, 2, 0), e = ()
        and edges ().]),
        (The spherical web with c = (4, 6, 8, 10, 5, 0, 7, 1, 9, 2, 11, 3), e = (7, 10, 9, 4, 11, 6, 5, 8)
        and edges ('(0,black,False)', ..., '(0,black,False)', '(0,black,False)').,
        -delta*B[The spherical web with c = (1, 4, 3, 5, 0, 2), e = (5, 4)
        and edges ('(0,black,False)', '(0,black,False)').] - delta*B[The spherical web with c = (4, 2, 5, 0, 3, 1), e = (5, 4)
        and edges ('(0,black,False)', '(0,black,False)').] + (delta^2-1)*B[The spherical web with c = (1, 0, 3, 2), e = ()
        and edges ().] + (delta^2-1)*B[The spherical web with c = (3, 2, 1, 0), e = ()
        and edges ().]),
        (The spherical web with c = (5, 7, 9, 11, 13, 6, 0, 8, 1, 10, 2, 12, 3, 14, 4), e = (8, 13, 10, 5, 12, 7, 14, 9, 6, 11)
        and edges ('(0,black,False)', ..., '(0,black,False)').,
        B[The spherical web with c = (1, 2, 0, 4, 3), e = ()
        and edges ().] + B[The spherical web with c = (1, 4, 3, 2, 0), e = ()
        ...
        and edges ('(0,black,False)', ..., '(0,black,False)').]))
    """
    if delta == None:
        delta = PolynomialRing(ZZ, 'delta').gen()

    s = Strand(oriented=0, colour='black')
    edge = SphericalSpider([s]*2).vertex()
    v = SphericalSpider([s]*3).vertex()

    relations = []

    # Loop relations
    S = SphericalSpider([])
    L = FreeSphericalSpider(delta.parent(), [])
    relations.append(tuple([S.loop(s), (delta**5+delta**4-5*delta**3-4*delta**2+6*delta+3)*L(S.empty())]))

    # Tadpole relation
    S = SphericalSpider([s])
    L = FreeSphericalSpider(delta.parent(), [s])
    tadpole = v.glue(edge,2)
    relations.append(tuple([tadpole, 0*L(tadpole)]))

    # Digon relation
    S = SphericalSpider([s]*2)
    L = FreeSphericalSpider(delta.parent(), [s]*2)
    relations.append(tuple([v.glue(v,2),-(delta+1)*(delta**2-2)*L(edge)]))

    # Triangle relation
    S = SphericalSpider([s]*3)
    L = FreeSphericalSpider(delta.parent(), [s]*3)
    triangle = S.polygon([v]*3)
    relations.append(tuple([triangle, (delta**2-1)*L(v)]))

    # Square relation
    S = SphericalSpider([s]*4)
    L = FreeSphericalSpider(delta.parent(), [s]*4)
    square = S.polygon([v]*4)
    H = v.glue(v, 1)
    K = H.rotate(1)
    U = edge.glue(edge, 0)
    I = U.rotate(1)
    relations.append(tuple([square, -delta*(L(H)+L(K))+(delta**2-1)*(L(U)+L(I))]))

    # Pentagon relation
    S = SphericalSpider([s]*5)
    L = FreeSphericalSpider(delta.parent(), [s]*5)
    pentagon = S.polygon([v]*5)
    a = v.glue(edge, 0)
    ar = sum([L(a.rotate(i)) for i in range(5)], L(0))
    b = v.glue(v.glue(v, 1), 1)
    br = sum([L(b.rotate(i)) for i in range(5)], L(0))
    relations.append(tuple([pentagon, ar-br]))

    return tuple(relations)

"""
def F4():

    The `F_4` skein relations.

"""

class WebAlgebra(CombinatorialFreeModule, Algebra):
    r"""
    Linear combinations of spherical webs.
    """
    @staticmethod
    def __classcall__(cls, base, boundary):
        if boundary in NN:
            boundary = [Strand()]*boundary

        bp = list(boundary)
        bn = [Strand(-a.oriented, a.colour) for a in bp]
        bd = tuple(bp + list(reversed(bn)))
        return super(WebAlgebra, cls).__classcall__(cls, base, bd)

    def __init__(self, base, boundary):
        r"""
        Initialise ``self``.`

        EXAMPLES::

            sage: WebAlgebra(QQ,[])
            Free module generated by The spherical spider with boundary [] over Rational Field
            sage: WebAlgebra(QQ,[]).an_element()
            0
        """
        CombinatorialFreeModule.__init__(self, base, SphericalSpider(boundary),
            category=Algebras(base).WithBasis().Unital())

    def boundary(self):
        r"""
        Return the boundary.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: A = WebAlgebra(QQ,[Strand()]*3); A
            Free module generated by The spherical spider with boundary ... over Rational Field
            sage: len(A.basis().keys().boundary)
            6
            sage: len(A.boundary())
            3
        """
        b = list(self.basis().keys().boundary)
        n = len(b)//2
        return tuple(b[:n])

    @cached_method
    def one_basis(self):
        r"""
        Return the identity element, as per ``AlgebrasWithBasis.ParentMethods.one_basis``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: WebAlgebra(QQ, []).one_basis()
            A closed spherical web with 0 edges.
            sage: WebAlgebra(QQ, 3).one_basis()
            The spherical web with c = (5, 4, 3, 2, 1, 0), e = ()
             and edges ().
            sage: WebAlgebra(QQ, 3).one()
            B[The spherical web with c = (5, 4, 3, 2, 1, 0), e = ()
             and edges ().]
        """
        from sage.combinat.spherical_spider import halfedge
        b = list(self.basis().keys().boundary)

        h = [ halfedge(a) for a in b ]
        c = {x: y for x, y in zip(h,reversed(h))}

        return SphericalWeb(c, {}, h)

    def product_on_basis(self, X, Y):
        r"""
        Return the product of two basis elements, as per
        ``AlgebrasWithBasis.ParentMethods.product_on_basis``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: A = WebAlgebra(QQ,[Strand()]*3)
            sage: A.one() * A.one() == A.one()
            True
        """
        n = len(self.boundary())
        return self(X.glue(Y, n))

    def markov_trace(self):
        r"""
        Construct the Markov trace on ``self``.

        EXAMPLES::

        """
        from sage.combinat.spherical_spider import halfedge
        b = list(self.basis().keys().boundary)
        n = len(b)
        h = [ halfedge(a) for a in b ]
        c = {h[i]: h[n-i-1] for i in range(n)}
        caps = SphericalWeb(c, {}, b)

        codomain = FreeSphericalSpider(self.base(),[])
        on_basis = lambda x : codomain(x.glue(caps, n))
        return self.module_morphism(codomain=codomain, on_basis=on_basis)

    def U(self, k):
        r"""
        Construct the element `U_i`

        EXAMPLES::

            sage: WebAlgebra(QQ,2).U(1)
            B[The spherical web with c = (1, 0, 3, 2), e = ()
             and edges ().]
            sage: A = WebAlgebra(QQ,3)
            sage: A.U(1) * A.U(2) * A.U(1) == A.U(1)
            True
            sage: A.U(2) * A.U(1) * A.U(2) == A.U(2)
            True
            sage: A = WebAlgebra(QQ,4)
            sage: A.U(1) * A.U(3) == A.U(3) * A.U(1)
            True
        """
        from sage.combinat.spherical_spider import halfedge
        b = list(self.basis().keys().boundary)
        n = len(b)
        if not 0 < k < n//2:
            return NotImplemented
        if b[k-1] != Strand(-b[k].oriented, b[k].colour):
            return NotImplemented


        h = [ halfedge(a) for a in b ]
        c = {h[i]: h[n-i-1] for i in range(n)}
        c[h[k]] = h[k-1]
        c[h[k-1]] = h[k]
        c[h[n-k-1]] = h[n-k]
        c[h[n-k]] = h[n-k-1]
        return self(SphericalWeb(c, {}, h))

def temperley_lieb(n, R=ZZ, q=None):
    r"""
    Construct the algebra homomorphism from the group algebra of the `n`-string braid group
    to the `n`-string Temperley-Lieb algebra.

    EXAMPLES::

        sage: from sage.combinat.spherical_spider import temperley_lieb
        sage: temperley_lieb(2)
            Generic morphism:
                From: Algebra of Braid group on 2 strands over Integer Ring
                To:   Free module generated by The spherical spider with boundary [...] over ... in q over Integer Ring
        sage: braid = BraidGroup(2).algebra(ZZ)(BraidGroup(2)([1]))
        sage: s = temperley_lieb(2)(braid); s
        q*B[The spherical web with c = (3, 2, 1, 0), e = ()
        and edges ().] - B[The spherical web with c = (1, 0, 3, 2), e = ()
        and edges ().]
        sage: braid = BraidGroup(2).algebra(ZZ)(BraidGroup(2)([-1]))
        sage: t = temperley_lieb(2)(braid); t
        (q^-1)*B[The spherical web ...] - B[The spherical web ...]
        sage: s * t
        B[The spherical web ...] - (q^-1-1+q)*B[The spherical web ...]
        sage: s1 = temperley_lieb(3)(BraidGroup(3).algebra(ZZ)(BraidGroup(2)([1])))
        sage: s2 = temperley_lieb(3)(BraidGroup(3).algebra(ZZ)(BraidGroup(2)([1])))
        sage: s1 * s2 * s1 == s2 *s1 * s2
        True
    """
    if not n in NN:
        raise ValueError(f"{n} must be nonnegative")

    domain = BraidGroup(n).algebra(R)
    if q == None:
        q = LaurentPolynomialRing(ZZ, 'q').gen()

    codomain = WebAlgebra(q.parent(), n)

    def on_basis(x):
        def sign(k):
            if k > 0:
                return 1
            return -1

        result = codomain.one()
        for k in x.Tietze():
            result *= q**sign(k) - codomain.U(abs(k))

        return result

    return domain.module_morphism(codomain=codomain, on_basis=on_basis)

def from_diagram(n, delta):
    r"""
    Construct the algebra morhism from the Temperley-Lieb algebra as a diagram algebra
    to the Temperley-Lieb algebra as a web algebra.

    EXAMPLES::

        sage: from sage.combinat.spherical_spider import from_diagram
        sage: from sage.combinat.diagram_algebras import TemperleyLiebAlgebra
        sage: t = TemperleyLiebAlgebra(2, 2).an_element()
        sage: from_diagram(2,2)
        Generic morphism:
        From: Temperley-Lieb Algebra of rank 2 with parameter 2 over Integer Ring
        To:   Free module generated by The spherical spider with boundary ... over Integer Ring
        sage: from_diagram(2,2)(t)
        2*B[The spherical web with c = (1, 0, 3, 2), e = ()
        and edges ().] + 2*B[The spherical web with c = (3, 2, 1, 0), e = ()
        and edges ().]
        """
    from sage.combinat.diagram_algebras import TemperleyLiebAlgebra
    from sage.combinat.spherical_spider import halfedge

    domain = TemperleyLiebAlgebra(n, delta, delta.parent())
    codomain = WebAlgebra(delta.parent(), n)
    def on_basis(x):
        h = [halfedge() for _ in range(2*n)]
        m = {i+1: i for i in range(n)}
        m.update({-i-1: 2*n-i-1 for i in range(n)})
        c = dict([])
        for r, s in x.diagram():
            c[h[m[r]]] = h[m[s]]
            c[h[m[s]]] = h[m[r]]

        return codomain(SphericalWeb(c, {}, h))

    return domain.module_morphism(codomain=codomain, on_basis=on_basis)

class Path():
    r"""
    An instance of this class is a generalisation of a highest weight word in a tensor product of (finite) crystals.
    """
    def __init__(self, steps, weight_path, Dyck_path):
        r"""
        A path is a sequence of steps and each step is labelled as up or down.
        The sequence of up/down labels is a Dyck path.
        """
        spaces = {a.parent() for a in steps}
        if len(spaces) > 1:
            raise ValueError("All steps must be in the same vector space.")
        if len(spaces) == 0:
            raise ValueError("The set of steps must be non-empty")
        space = spaces.pop()
        self.space = space
        self.rank = space.dimension()

        n = len(weight_path)
        if len(Dyck_path) != n:
            raise ValueError("The weight_path and the Dyck_path must have the same length.")
        self.length = n

        if any(abs(x) != 1 for x in Dyck_paths):
            raise ValueError("Dyck_path must be a sequence of +1 or -1.")
        self.Dyck_path = Dyck_path

        if any(a not in steps for a in weight_path):
            raise ValueError("Every element of the weight_path must be an element of steps.")
        self.steps = steps
        self.weight_path = weight_path

        heights = [None]*(n+1)
        heights[0] = space(0)
        for i in range(n):
            heights[i+1] = heights[i] + Dyck_path[i]*weight_path[i]

        if any( any(i<0 for i in a) for a in heights):
            raise ValueError("The path is not dominant.")

        self.heights = heights

    def switch(self, i):
        r"""
        Change a down-up to an up-down,

        """
        if self.Dyck_path[i] != -1:
            raise RuntimeError(f"Position {i} of {self.Dyck_path} must be -1")
        if self.Dyck_path[i+1] != 1:
            raise RuntimeError(f"Position {i+1} of {self.Dyck_path} must be +1")

        left  = [max(x-y,0) for x, y in zip(weight_path[i], weight_path[i+1])]
        right = [max(0,y-x) for x, y in zip(weight_path[i], weight_path[i+1])]

        wp = copy(self.weight_path)
        wp[i] = vector(QQ, left)
        wp[i+1] = vector(QQ, right)
        Dp = copy(self.Dyck_path)
        Dp[i] = 1
        Dp[i+1] = -1

        return Path(self.steps, wp, Dp)
