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
#from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.graphs.graph import Graph
#from sage.combinat.permutation import Permutation
from sage.combinat.baxter_permutations import BaxterPermutations
from sage.structure.richcmp import richcmp, op_EQ, op_NE
from copy import copy
from collections import namedtuple
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.homset import Hom

Strand  = namedtuple('Strand', ['oriented','colour','crossing'], defaults=[0,'black',False])

class halfedge():
    """
    The class of half edges in a surface graph.

    This should probably be an attribute either of SphericalWeb or SphericalSpider
    """
    def __init__(self, st=Strand()):
        """
        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: halfedge()
            <sage.combinat.spherical_spider.halfedge object at ...>
        """
        self.strand = st

    def __hash__(self):
        return hash(self.strand)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

        sage: from sage.combinat.spherical_spider import halfedge
        sage: halfedge()._repr_()
        '(0,black,False)'
        sage: halfedge() # indirect test
        <sage.combinat.spherical_spider.halfedge object at ...>
        """
        return f"({self.strand.oriented},{self.strand.colour},{self.strand.crossing})"

    def dual(self):
        """
        Construct the dual halfedge.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: halfedge().dual()
            <sage.combinat.spherical_spider.halfedge object at ...>
        """
        st = self.strand
        return (halfedge(Strand(-st.oriented, st.colour, st.crossing)))

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
        Initialise an instance of :class`SphericalWeb`.

        INPUT:

            * `c` a bijection of the set of half-edges
            * `e` a partial involution of the set of half-edges
            * `b` the ordered list of boundary edges

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import halfedge
            sage: b = [halfedge(),halfedge()]
            sage: c = {b[0]:b[1], b[1]:b[0]}
            sage: SphericalWeb(c,{},b)
            The spherical web with c = (1, 0) and e = ().
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
        Implement the abstract method :meth:`__copy__` of :class:`ClonableElement`.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black',False)]*3).vertex().__copy__()
            The spherical web with c = (1, 2, 0) and e = ().
        """
        D = {a:halfedge(a.strand) for a in self.cp}
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
            The spherical web with c = (0,) and e = ().
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
            The spherical web with c = (1, 0) and e = ().

        TESTS::

        This should not ever happen.

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black',False)]*2).vertex() # indirect doctest
            The spherical web with c = (1, 0) and e = ().

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
            'The spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2) and e = (6, 7, 8, 3, 4, 5).'
        """
        if len(self.b) > 0:
            cn, en = self.canonical()
            return f"The spherical web with c = {cn} and e = {en}."
        else:
            return f"A closed spherical web with {int(len(self.e)/2)} edges."

    @cached_method
    def canonical(self):
        r"""
        A canonical labelling of the elements of ``self``.

        This returns two lists of integers such that ``self``
        can be recovered, up to isomorphism, from these two sequences.

        Let ``self`` have `n` elements and `k` boundary elements.
        Then the first list is a bijection on [0,1,...,n-1] and
        the second list is an involution on [k,k+1,...,n-1].

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand()]*3).vertex()
            sage: cn, en = SphericalSpider([]).polygon([u,u,u,u]).canonical()
            sage: cn
            (4, 6, 8, 10, 5, 0, 7, 1, 9, 2, 11, 3)
            sage: en
            (7, 10, 9, 4, 11, 6, 5, 8)
        """
        b = self.b
        c = self.cp
        e = self.e
        k = len(b)
        gen = self._traversal(b)
        gl = list(gen) # This defeats the purpose of using a generator.
        Dp = {a:i for i,a in enumerate(b)}
        Dp.update({a:(k+i) for i,a in enumerate(gl)})
        cn = tuple([Dp[c[a]] for a in b]+[Dp[c[a]] for a in gl])
        en = tuple([Dp[e[a]] for a in gl])
        return cn, en

    def __hash__(self):
        r"""
        Overload the :method:`__hash__`.

        This is needed to put a :class:`SphericalWeb` into a :class:`set`
        or to use it as a key in a :class:`dict'.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: v = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: v.__hash__()  # random

            sage: hash(v) # random

            sage: u = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: set([u,v]) # indirect doctest
            {The spherical web with c = (1, 2, 0) and e = ().}
        """
        return hash((self.parent(),*self.canonical()))

    def _richcmp_(self, other, op):
        """
        Overload :meth:`__eq__` and :meth:`__ne__`.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black',False)]*4).vertex()
            sage: v = SphericalSpider([Strand(0,'black',False)]*4).vertex()
            sage: u is v, u == v, u != v # indirect doctest
            (False, True, False)
            sage: u < v # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: '<' not supported between ... and 'SphericalWeb'

        TODO::

            This should take the parent and/or type into account.
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
            sage: w = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: w._traversal(w.b)
            <generator object SphericalWeb._traversal at ...>
            sage: w._traversal(w.b[0])
            <generator object SphericalWeb._traversal at ...>
        """
        if isinstance(initial,halfedge):
            initial = tuple([initial])
        else:
            initial = tuple(initial)
        if not set(initial).issubset(self.cp):
            raise ValueError("initial must be a subset of the set of elements")

        c = self.cp
        e = self.e
        b = self.b

        visited = list(initial)
        new = list()
        flag = True
        while flag:
            flag = False
            for a in visited:
                b = c[a]
                while b not in visited:
                    new.append(b)
                    yield b
                    flag = True
                    b = c[b]
            visited += new
            for a in visited:
                if a in e:
                    b = e[a]
                    if b not in visited:
                        new.append(b)
                        yield b
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
        xd = Strand(-st.oriented, st.colour, st.crossing)
        if y.strand != xd:
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
            The spherical web with c = (1, 2, 3, 4, 0) and e = ().
            sage: u.rotate(-1)
            The spherical web with c = (1, 2, 3, 4, 0) and e = ().
        """
        result = self.__copy__()
        b = result.b
        nb = b[k:]+b[:k]
        return SphericalWeb(result.cp, result.e, nb)

    def glue(self, other, n):
        r"""Glue two ribbon graphs together.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: v = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: u.glue(v,1)
            The spherical web with c = (1, 4, 3, 5, 0, 2) and e = (5, 4).
            sage: u.glue(v,0)
            The spherical web with c = (1, 2, 0, 4, 5, 3) and e = ().
            sage: u.glue(v,4)
            Traceback (most recent call last):
            ...
            ValueError: n=4 is too large
        """
        #if n < 0:
        #    raise ValueError(f"n={n} cannot be negative")
        #parent = self.parent()
        #if parent != other.parent():
        #    raise ValueError(f"the two parents {self.parent()} and {other.parent()} are different")
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

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: u = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: u.mirror_image()
            The spherical web with c = (1, 2, 0) and e = ().
            sage: v = SphericalSpider([Strand(0,'black',False)]*4).vertex()
            sage: u.glue(v,1).mirror_image()
            The spherical web with c = (1, 2, 5, 4, 6, 0, 3) and e = (6, 5).
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
            sage: u = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: len(u.faces())
            3
            sage: w = SphericalSpider([Strand(0,'black',False)]*4).vertex()
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
            (The spherical web with c = (1, 2, 0) and e = ().,
            A closed spherical web with 0 edges.)
            sage: u.glue(u,0)
            The spherical web with c = (1, 2, 0, 4, 5, 3) and e = ().
            sage: u.glue(u,0).components()
            (The spherical web with c = (1, 2, 0, 4, 5, 3) and e = ().,
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
        # This adds each edge twice.
        for a in e:
            G.add_edge(a,e[a],'e')
        return G

    def show(self):
        r"""Show the web ``self``.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand()]*3).vertex().show()
            Graphics object consisting of 4 graphics primitives
        """
        return self.to_graph().plot(vertex_labels=False)

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
            '\\begin{tikzpicture}\n\\draw (0,0) circle (1cm);\n\\draw (-1.00000000000000,0.0) -- (0.0,0.0);\n\\draw (1.00000000000000,0.0) -- (0.0,0.0);\n\\end{tikzpicture}\n'

            If there are no boundary points only the boundary circle is drawn.

            sage: SphericalSpider([]).loop(Strand())._latex_()
            '\\begin{tikzpicture}\n\\draw (0,0) circle (1cm);\n\\end{tikzpicture}\n'

        TODO::

            Add colour, direction, under crossing.
        """
        lines = self._layout()

        result = "\\begin{tikzpicture}\n"
        result += "\\draw (0,0) circle (1cm);\n"
        for a in lines:
            result += "\\draw ({},{}) -- ({},{});\n".format(a[0][0],a[1][0],a[1][0],a[1][1])
        result += "\\end{tikzpicture}\n"

        return result

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
            The spherical web with c = (3, 5, ... 12, 9) and e = (9, 10, ... 14, 13).
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

        Ds = {a:halfedge(a.strand) for a in self.cp if not a in D.values()}
        Dk = {a:halfedge(a.strand) for a in k.cp}
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

#### End  of methods for rewriting ####

#### Start of Parent ####

class SphericalSpider(UniqueRepresentation, Parent):
    r"""
    The Parent class for SphericalWeb.

    EXAMPLES::

        sage: from sage.combinat.spherical_spider import Strand
        sage: SphericalSpider([Strand(0,'black',False),Strand(0,'black',False)])
        The spherical spider with boundary (Strand(oriented=0, colour='black', crossing=False), ...)
        sage: SphericalSpider([])
        The spherical spider with boundary ()
    """
    @staticmethod
    def __classcall__(cls, boundary):
        return super(SphericalSpider, cls).__classcall__(cls, tuple(boundary))


    def __init__(self, boundary):
        r"""
        Initialise an instance of this class.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black',False),Strand(0,'black',False)])
            The spherical spider with boundary (Strand(oriented=0, colour='black', crossing=False), ...)
        """

        self.boundary = boundary

        Parent.__init__(self)

    def _repr_(self):
        r"""
        Overload the default method.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: P = SphericalSpider([Strand(0,'black',False)]*3)
            sage: P._repr_()
            "The spherical spider with boundary (Strand(oriented=0, colour='black', crossing=False),  ...)"
        """
        return f"The spherical spider with boundary {self.boundary}"

    Element = SphericalWeb

    def vertex(self):
        r"""
        Construct a single vertex.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([Strand(0,'black',False)]*4).vertex()
            The spherical web with c = (1, 2, 3, 0) and e = ().
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
    def loop(st: Strand):
        r"""
        Construct a loop.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: SphericalSpider([]).loop(Strand(0,'black',False))
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
            sage: u = SphericalSpider([Strand(0,'black',False)]*3).vertex()
            sage: SphericalSpider(tuple([])).polygon([u,u,u])
            The spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2) and e = (6, 7, 8, 3, 4, 5).
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

        for u,v in zip(corners,corners[1:]):
            x = u.b[-1]
            y = v.b[0]
            c,e = SphericalWeb._stitch(c,e,x,y)

        x = corners[-1].b[-1]
        y = corners[0].b[0]
        c,e = SphericalWeb._stitch(c,e,x,y)

        b = sum([list(a.b[1:-1]) for a in corners],[])
        return SphericalWeb(c, e, b)

    @staticmethod
    def from_permutation(pi, baxter=True):
        r"""
        Construct a planar map from a two stack sorted permutation.

        This implements the algorithm in :arxiv:`math/0805.4180`.
        This algorithm is designed to apply to Baxter permutations.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: S = SphericalSpider(tuple({}))
            sage: pi = Permutation([5,3,4,9,7,8,10,6,1,2])
            sage: S.from_permutation(pi)
            The spherical web with c = (...) and e = (...).
            sage: pi = Permutation([2,4,1,3])
            sage: S.from_permutation(pi)
            Traceback (most recent call last):
            ...
            ValueError: [2, 4, 1, 3] is not a Baxter permutation
            sage: S.from_permutation(pi,baxter=False)
            The spherical web with c = (1, 0) and e = ().
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

#### End of Parent ####

class LinearSphericalSpider(CombinatorialFreeModule):
    r"""
    Linear combinations of spherical webs.
    """
    @staticmethod
    def __classcall__(cls, base_ring, boundary):
        return super(LinearSphericalSpider, cls).__classcall__(cls, base_ring, tuple(boundary))

    def __init__(self, base_ring, boundary):
        r"""
        Initialise ``self``.`

        EXAMPLES::LinearSphericalSpider

            sage: LinearSphericalSpider(QQ,[])
            Free module generated by The spherical spider with boundary () over Rational Field
            sage: LinearSphericalSpider(QQ,[]).an_element()
            0
            sage: F = LinearSphericalSpider(QQ,[])
            sage: isinstance(0, F.element_class)
            False
            sage: isinstance(F(0), F.element_class)
            True
        """
        CombinatorialFreeModule.__init__(self, base_ring, SphericalSpider(boundary))

    def boundary(self):
        return self.basis().keys().boundary

    def vertex(self):
        r"""
        Return the vertex on the boundary as a monomial.

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: LinearSphericalSpider(QQ,[Strand()]*3).vertex()
            B[The spherical web with c = (1, 2, 0) and e = ().]
            """
        v = self.basis().keys().vertex()
        return self.monomial(v)

    def rotate(self, k):
        r"""
        Extend :method'rotate' by linearity

        EXAMPLES::

            sage: from sage.combinat.spherical_spider import Strand
            sage: LinearSphericalSpider(QQ,[Strand()]*3).rotate(1)

        """
        b = self.boundary()
        codomain = LinearSphericalSpider(self.base_ring, b[k:]+b[:k])
        on_basis = lambda x : x.leading_support().rotate(k)
        return self._module_morphism(codomain=codomain, on_basis=on_basis)

