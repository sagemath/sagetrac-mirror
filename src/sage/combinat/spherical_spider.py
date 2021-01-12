r"""

### Introduction
This is an implementation of finitely generated free spiders.
The account below is intended as an extended introduction.
There are technical details which I skate over and which need to
be dealt with in the implementation.

This is intended to be the first phase of a project to implement
the Bendix-Knuth completion algorithm discussed in [5]. This algorithm starts
with a finite presentation and iteratively constructs new relations.
If it terminates the final presentation is confluent.

### Spiders

The notion of spiders was introduced in [3].
A spider consists of a set (or vector space) with the operations of rotate, join and stitch.
There are axioms for these operations and a spider is equivalent to a strict pivotal
category or, with an extra axiom, to a strict spherical category. The motivation was
to study the category of finite dimensional representations of a quantum group.
These categories give examples of spiders and the main result of the paper was to
give finite confluent presentations for the rank two examples, A2, B2=C2, G2, (the rank one example
was known previously as the skein relation approach to the Jones polynomial).

It has been an open problem since then to construct finite confluent presentations
for higher rank examples. It is known (unpublished) that these examples are finitely
generated but not that these are finitely presented. A finite confuent presentation
would give algorithms for computing link polynomials and for doing the caculations
in [1].

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
of an operad. In fact this is a cyclic operad as defined in [2]. Furthermore a spider
is a cyclic algebra for this cyclic operad. Taking the glueing operation as basic and
building up operations is the componential approach to cyclic operads given in

REFERENCES:

[1] Predrag Cvitanović
Group Theory: Birdtracks, Lie's, and Exceptional Groups
Princeton University Press, 2008
ISBN    0691118361, 9780691118369

[2] E Gretzler, MM Kapranov
Cyclic operads and cyclic homology
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.146.2678&rep=rep1&type=pdf

[3] G.Kuperberg
Spiders for rank 2 Lie algebras
Commun.Math. Phys. 180, 109–151 (1996)
https://arxiv.org/abs/1103.3519.

[4] Jovana Obradovic
Monoid-like definitions of cyclic operad
Theory and Applications of Categories, Vol. 32, 2017, No. 12, pp 396-436.
http://www.tac.mta.ca/tac/volumes/32/12/32-12.pdf

[5] Adam S Sikora and Bruce W Westbury
Confluence theory for graphs
Algebraic & Geometric Topology 7 (2007) 439–478
https://arxiv.org/abs/math/0609832

AUTHORS:

- Bruce Westbury (2021): initial version

"""

#*****************************************************************************
#       Copyright (C) 2021 Bruce Westbury bruce.westbury@gmail.com
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from dataclasses import dataclass
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableElement
from sage.graphs.graph import Graph

class halfedge():
    """
    The class of half edges in a surface graph.

    This should probably be an attribute either of SphericalWeb or SphericalSpider
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: halfedge()
            <sage.combinat.spherical_spider.halfedge object at ...>
        """

class SphericalWeb(ClonableElement):
    r"""The class of mutable ribbon graphs.

    This consists of
    * a set of half-edges
    * a bijection `c` on the set of half-edges with no fixed points
    * an involution `e` on a subset of the half-edges

    The half-edges for which `e` is undefined are the boundary half-edges.
    This set has a total order.

    The only orbits of `c` of order two have both half-edges in the boundary.
    """

    def __init__(self, c: dict, e: dict, b: list, parent: Parent):
        r"""
        Initialise a ``SphericalWeb``.

        INPUT::

            * `c` a bijection of the set of half-edges
            * `e` a partial involution of the set of half-edges
            * `b` the ordered list of boundary edges

        EXAMPLES::

            sage: b = [halfedge(),halfedge()]
            sage: c = {b[0]:b[1], b[1]:b[0]}
            sage: SphericalSpider('plain')(c,{},b)
            The plain spherical web with c = (1, 0) and e = ().
        """
        ClonableElement.__init__(self,parent=parent)
        self.cp = c
        self.e = e
        self.boundary = tuple(b)
        self.normalize()
        self.set_immutable()
        self.check()

    def __copy__(self):
        r"""
        Implement the abstract method of :class:`ClonableElement`.

        EXAMPLES::

            sage: SphericalSpider('plain').vertex(3).__copy__()
            The plain spherical web with c = (1, 2, 0) and e = ().
        """
        D = {a:halfedge() for a in self.cp}
        c = {D[a]:D[self.cp[a]] for a in self.cp}
        e = {D[a]:D[self.e[a]] for a in self.e}
        b = [D[a] for a in self.boundary]
        result = self.parent()(c,e,b)
        return result

    def check(self):
        r"""
        Implement the abstract method of :class:`ClonableElement`.

        Check ``self`` is a valid web.
        """
        c = self.cp
        e = self.e
        b = self.boundary
        h = set(c.keys())
        if not all(isinstance(a,halfedge) for a in h):
            raise ValueError("every element must be a half-edge")
        if set(c.values()) != h:
            raise ValueError("the map c is required to be a bijection")
        if any(e[e[a]] != a for a in e):
            raise ValueError("the map e must be an involution")
        if any(c[a] == a for a in c):
            raise ValueError("the mapping c has at least one fixed point")
        if any(e[a] == a for a in e):
            raise ValueError("the mapping e has at least one fixed point")
        if not set(e.keys()).issubset(h):
            raise ValueError("the domain of e must be a subset of the domain of c")
        if not set(b).issubset(h):
            raise ValueError("the boundary must be a subset of the domain of c")
        if not set(e.keys()).isdisjoint(set(b)):
            raise ValueError("the domain of e must not intersect the boundary")
        for i,a in enumerate(b):
            u = a
            while c[u] in e:
                u = e[c[u]]
            j = b.index(c[u])
            if 0 < j < i:
                raise ValueError("boundary is inconsistent")

    def normalize(self):
        r"""
        Overload the :method:`normalize` of :class:`ClonableElement`.

        This removes nearly all vertices of degree two.

        EXAMPLES::

            sage: h = [halfedge() for i in range(4)]
            sage: c = {h[0]:h[1], h[1]:h[0], h[2]:h[3], h[3]:h[2]}
            sage: e = {h[1]:h[2],h[2]:h[1]}
            sage: b = [h[0],h[3]]
            sage: SphericalSpider('plain')(c,e,b) # indirect doctest
            The plain spherical web with c = (1, 0) and e = ().

            This should not ever happen.

            sage: SphericalSpider('plain').vertex(2) # indirect doctest
            The plain spherical web with c = (1, 0) and e = ().

            Check loops are not removed.

            sage: SphericalSpider('plain').loop() # indirect doctest
            A closed plain spherical web with 1 edges.
        """
        self._set_mutable()
        flag = True
        while(flag):
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
                elif y in self.boundary:
                    flag = True
                    c[y] = c[z]
                    w = [a for a in c if c[a] == z][0]
                    c[w] = y
                    c.pop(x)
                    c.pop(z)
                    e.pop(x)
                    e.pop(z)
        self.set_immutable()

    def _repr_(self):
        r"""
        Overload default implementation.

        EXAMPLES::

            sage: polygon_web(3)._repr_()
            'The plain spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2) and e = (6, 7, 8, 3, 4, 5).'
        """
        if len(self.boundary) > 0:
            cn, en = self.canonical()
            return f"The {self.parent()._name} spherical web with c = {cn} and e = {en}."
        else:
            return f"A closed {self.parent()._name} spherical web with {int(len(self.e)/2)} edges."

    def __str__(self):
        r"""
        Overload default implementation.

        EXAMPLES::

            sage: str(polygon_web(4))
            'A plain spherical web with 4 internal edges and 4 boundary edges.'
        """
        return f"A {self.parent()._name} spherical web with {int(len(self.e)/2)} internal edges and {len(self.boundary)} boundary edges."

    @cached_method
    def canonical(self):
        r"""
        A canonical labelling of the elements of `self``.

        This returns two lists of integers such that ``self``
        can be recovered, up to isomorphism, from these two sequences.

        Let ``self`` have `n` elements and `k` boundary elements.
        Then the first list is a bijection on [0,1,...,n-1] and
        the second list is an involution on [k,k+1,...,n-1].

        EXAMPLES::

            sage: cn, en = polygon_web(4).canonical()
            sage: cn
            (4, 6, 8, 10, 5, 0, 7, 1, 9, 2, 11, 3)
            sage: en
            (7, 10, 9, 4, 11, 6, 5, 8)
        """
        b = self.boundary
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

    def _hash_(self):
        r"""
        Implement the :method:`_hash_` of :class:`ClonableElement`.

        EXAMPLES::

            sage: polygon_web(4)._hash_()  # random

            sage: hash(polygon_web(4)) # random
        """
        return hash((self.parent(),*self.canonical()))

    __hash__ = _hash_
    """
    Overload the :method:`__hash__`.

    This is needed to put a :class:`SphericalWeb` into a :class:`set`
    or to use it as a key in a :class:`dict'.

    EXAMPLES::

        sage: u = vertex(3)
        sage: v = vertex(3)
        sage: set([u,v]) # indirect doctest
        {The plain spherical web with c = (1, 2, 0) and e = ().}
    """

    def __eq__(self,other):
        """
        Overload :method:`__eq__`.

        EXAMPLES::

            sage: u = polygon_web(4)
            sage: v = polygon_web(4)
            sage: u is v, u == v # indirect doctest
            (False, True)
        """
        if type(self) != type(other):
            return False
        if self.parent() != other.parent():
            return False
        return self.canonical() == other.canonical()

    def __ne__(self,other):
        """
        Overload :method:`__ne__`.

        EXAMPLES::

            sage: u = polygon_web(4)
            sage: v = polygon_web(4)
            sage: u != v # indirect doctest
            False
        """
        return not self == other

#### End of underscore methods ####

#### Start of methods for working with webs ####

    def _traversal(self,initial):
        """
        A generator for the elements of ``self`` connected to the
        elements in ``initial``.

        EXAMPLES::

            sage: w = SphericalSpider('plain').vertex(3)
            sage: w._traversal(w.boundary)
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

    def rotate(self,k: int):
        r"""Rotate the boundary anticlockwise `k` steps.

        EXAMPLES::

            sage: polygon_web(5).rotate(3)
            The plain spherical web with c = (5, 7, 9, 11, 13, 6, 0, 8, 1, 10, 2, 12, 3, 14, 4) and e = (8, 13, 10, 5, 12, 7, 14, 9, 6, 11).
            sage: polygon_web(5).rotate(-1)
            The plain spherical web with c = (5, 7, 9, 11, 13, 6, 0, 8, 1, 10, 2, 12, 3, 14, 4) and e = (8, 13, 10, 5, 12, 7, 14, 9, 6, 11).
        """
        result = self.__copy__()
        b = result.boundary
        result._set_mutable()
        result.boundary = b[k:]+b[:k]
        result.set_immutable()
        return result

    def glue(self, other, n: int):
        r"""Glue two ribbon graphs together.

        EXAMPLES::

            sage: u = SphericalSpider('plain').vertex(3)
            sage: v = SphericalSpider('plain').vertex(3)
            sage: u.glue(v,1)
            The plain spherical web with c = (1, 4, 3, 5, 0, 2) and e = (5, 4).
        """
        if n < 0:
            raise ValueError(f"n={n} cannot be negative")
        parent = self.parent()
        if parent != other.parent():
            raise ValueError(f"the two parents {self.parent()} and {other.parent()} are different")
        if n > len(self.boundary) or n > len(other.boundary):
            raise ValueError(f"n={n} is too large")

        Ds = {a:halfedge() for a in self.cp}
        c = {Ds[a]:Ds[self.cp[a]] for a in self.cp}
        e = {Ds[a]:Ds[self.e[a]] for a in self.e}
        bs = [Ds[a] for a in self.boundary]

        Do = {a:halfedge() for a in other.cp}
        c.update({Do[a]:Do[other.cp[a]] for a in other.cp})
        e.update({Ds[a]:Ds[other.e[a]] for a in other.e})
        bo = [Do[a] for a in other.boundary]

        b = bs[:-n]+bo[n:]

        topL = [halfedge() for i in range(n)]
        botL = [halfedge() for i in range(n)]
        topR = [halfedge() for i in range(n)]
        botR = [halfedge() for i in range(n)]

        e.update({topR[i]:botL[i] for i in range(n)})
        e.update({botL[i]:topR[i] for i in range(n)})
        e.update({x:y for x,y in zip(bs[-n:],topL)})
        e.update({y:x for x,y in zip(bs[-n:],topL)})
        e.update({x:y for x,y in zip(bo[:n],reversed(botR))})
        e.update({y:x for x,y in zip(bo[:n],reversed(botR))})

        c.update({topR[i]:topL[i] for i in range(n)})
        c.update({topL[i]:topR[i] for i in range(n)})
        c.update({botR[i]:botL[i] for i in range(n)})
        c.update({botL[i]:botR[i] for i in range(n)})

        return self.parent()(c,e,b)

    def mirror_image(self):
        r"""
        Construct the mirror image of ``self``.
        """
        D =  {a:halfedge() for a in self.cp}
        cn = {D[self.cp[a]]:D[a] for a in D}
        en = {D[a]:D[self.e[a]] for a in self.e}
        bn = reversed([D[a] for a in self.boundary])
        return self.parent()(cn,en,bn)

    def vertices(self):
        """
        Find the vertices of ``self``.

        These are the orbits of `c`.

        EXAMPLES::

            sage: [len(a) for a in polygon_web(4).vertices()]
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

            sage: len(SphericalSpider('plain').vertex(4).faces())
            4
        """
        c = self.cp
        e = self.e
        he = set(c.keys())
        result = set()

        # First find the external faces.
        for i,a in enumerate(self.boundary):
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
        Return True if ``self`` is closed.

        Note that this means that the boundary is empty.

        EXAMPLES::

            sage: SphericalSpider('plain').vertex(3).is_closed()
            False
            sage: SphericalSpider('plain').loop().is_closed()
            True
        """
        return len(self.boundary) == 0

    def is_connected(self):
        """
        Return True if ``self`` is connected.

        Note that this means that the diagram including the boundary
        is connected.

        EXAMPLES::

            sage: SphericalSpider('plain').vertex(3).is_connected()
            True
            sage: SphericalSpider('plain').loop().is_connected()
            False
            sage: u = SphericalSpider('plain').vertex(3)
            sage: v = SphericalSpider('plain').vertex(3)
            sage: u.glue(v,0).is_connected()
            False
        """
        return len(self.cp) == len(self.canonical()[0])

    def components(self):
        """
        Return the closed components of ``self``.

        This is the complement of the connected component of the boundary.
        """
        Dn = {a:halfedge() for a in self.boundary}
        for a in self._traversal(self.boundary):
            Dn[a] = halfedge()

        cn = {Dn[a]:Dn[self.cp[a]] for a in self.cp}
        en = {Dn[a]:Dn[self.e[a]] for a in self.e}
        bn = [Dn[a] for a in self.boundary]
        wb = SphericalWeb(cn,en,bn)

        Dc = {a:halfedge() for a in self.cp if not a in Dn}
        cc = {Dc[a]:Dc[self.cp[a]] for a in self.cp}
        ec = {Dc[a]:Dc[self.e[a]] for a in self.e}
        wc = SphericalWeb(cc,ec,[])

        return wb, wc

    def is_decomposable(self):
        """
        Return True if ``self`` is decomposable.

        A web `w` is decomposable if it can be written as `w = u.glue(v,0)`
        where `u` and `v` are non-empty.

        Note that this means that the diagram excluding the boundary
        is connected.
        """
        if len(self.boundary) == 0:
            raise ValueError("not implemented for a closed web")
        return len(self.cp) == len(list(self._traversal(self.boundary[0])))

    def is_separable(self):
        r"""
        Return ``True`` if ``self`` is separable.

        This means each face has distinct vertices.
        Including the boundary faces.

        EXAMPLES::

            sage: S = SphericalSpider('plain')
            sage: S.vertex(4).glue(S.vertex(2),2).is_separable()
            True
        """
        from itertools import product
        for v,f in product(self.vertices(),self.faces()):
            if len(set(v).intersection(set(f))) > 1:
                return True
        return False

    def to__permutation(self):
        """
        If ``self`` is non separable, encode ``self`` as a
        two stack sortable permutation.

        There should also be a from_permutation
        """
        # Construct rooted non-separable planar map.
        ## Add new vertex and connect to boundary, root connects to first boundary point.
        ## For inverse, delete vertex containing root and reconstruct boundary.
        # Construct bipolar oriented planar map.
        ## This means edges round a vertex are partitioned into two ordered subsets.
        # Construct the two ordered trees.
        # Traverse the trees and construct the permutation.

    def to_graph(self):
        r"""
        Construct the graph of ``self``.

        EXAMPLES::

            sage: SphericalSpider('plain').vertex(3).to_graph()
            Graph on 3 vertices
            sage: polygon_web(3).to_graph()
            Graph on 9 vertices
        """
        c = self.cp
        e = self.e
        G = Graph({a:[c[a]] for a in c})
        for a in e:
            G.add_edge(a,e[a],'e')
        return G.copy(immutable=True)

    def show(self):
        r"""Show the web ``self``.

        EXAMPLES::

            sage: SphericalSpider('plain').vertex(3).show()
            Graphics object consisting of 4 graphics primitives
        """
        return self.to_graph().plot(vertex_labels=False)

#### End of methods for working with webs ####

#### Start of methods for rewriting ####

    def search(self, h):
        r"""
        Find copies of h in ``self``

        EXAMPLES::

            sage: len(list(polygon_web(4).search(SphericalSpider('plain').vertex(3))))
            12

        TODO::

            This should be rewritten to use :meth:``_traversal``.
        """
        if self.parent() != h.parent():
            raise ValueError(f"the two parents {self.parent()} and {other.parent()} are different")

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
            Dm = {h.boundary[0]:x}
            if test(x):
                assert [ a for a in c if not (a in Dm) ] == [], "Mapping is not fully defined."
                if len(set(Dm.values())) == len(Dm):
                    yield Dm

    def replace(self,k,D: dict,h):
        r"""
        Replace image of map D:h -> ``self`` by k

        EXAMPLES::

            sage: g = polygon_web(4)
            sage: h = SphericalSpider('plain').vertex(3)
            sage: D = next(g.search(h))
            sage: g.replace(polygon_web(3),D,h)
            The plain spherical web with c = (4, 6, 8, ... 13, 15, 12) and e = (12, 13, 9, ... 10, 17, 16).
        """
        parent = self.parent()
        if parent != k.parent():
            raise ValueError(f"the two parents {self.parent()} and {k.parent()} are different")
        if parent != h.parent():
            raise ValueError(f"the two parents {self.parent()} and {h.parent()} are different")
        n = len(k.boundary)
        if len(h.boundary) != n:
            raise ValueError("boundaries of k and h must have the same length")

        Ds = {a:halfedge() for a in self.cp if not a in D.values()}
        Dk = {a:halfedge() for a in k.cp}
        c = {Ds[a]:Ds[self.cp[a]] for a in Ds}
        c.update({Dk[a]:Dk[k.cp[a]] for a in Dk})

        e = {Ds[a]:Ds[self.e[a]] for a in Ds if a in self.e and self.e[a] in Ds}
        e.update({Dk[a]:Dk[k.e[a]] for a in k.e})

        Db = {x:y for x,y in zip(h.boundary,k.boundary)}
        b = [None]*len(self.boundary)
        for i,a in enumerate(self.boundary):
            if a in Ds:
                b[i] = Ds[a]

        for a in h.boundary:
            if D[a] in self.e:
                s = halfedge()
                t = halfedge()
                c[s] = t
                c[t] = s
                x = Ds[self.e[D[a]]]
                y = Dk[Db[a]]
                e[x] = s
                e[s] = x
                e[y] = t
                e[t] = y
            else:
                i = self.boundary.index(D[a])
                b[i] = Dk[Db[a]]

        return parent(c,e,b)

#### End  of methods for rewriting ####

#### Start of Parent ####

class SphericalSpider(Parent,UniqueRepresentation):
    r"""
    The Parent class for SphericalWeb.

    EXAMPLES::

        sage: SphericalSpider('plain')
        The plain spherical spider.
    """
    def __init__(self,name: str):
        Parent.__init__(self)

        self._name = name

        def name(self):
            return self._name

    def _repr_(self):
        r"""
        Overload the default method.

        EXAMPLES::

            sage: P = SphericalSpider('plain')
            sage: P._repr_()
            'The plain spherical spider.'
        """
        return f"The {self._name} spherical spider."

    def _element_constructor_(self,c,e,b):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

        """
        return self.element_class(c,e,b,self)

    def _an_element_(self):
        """

        EXAMPLES::

            sage: SphericalSpider('plain')._an_element_()
            The plain spherical web with c = (1, 0) and e = ().
        """
        b = [halfedge(),halfedge()]
        c = {b[0]:b[1], b[1]:b[0]}
        return self.element_class(c,{},b,self)

    Element = SphericalWeb

    def vertex(self,n: int):
        r"""Construct a single vertex of valency `n`.

        EXAMPLES::

            sage: SphericalSpider('plain').vertex(4)
            The plain spherical web with c = (1, 2, 3, 0) and e = ().
        """
        if n<2:
            raise ValueError(f"n={n} must be at least 2")
        b = [ halfedge() for i in range(n) ]
        c = {b[i-1]:b[i] for i in range(n)}
        e = dict([])
        return self.element_class(c,e,b,self)

    def loop(self):
        r"""
        Construct a loop.

        EXAMPLES::

            sage: SphericalSpider('plain').loop()
            A closed plain spherical web with 1 edges.
        """
        h = [halfedge(),halfedge()]
        c = {h[0]:h[1], h[1]:h[0]}
        e = {h[0]:h[1], h[1]:h[0]}
        return self.element_class(c,e,[],self)

    def empty(self):
        """
        Construct the empty diagram.

        EXAMPLES::

            sage: SphericalSpider('plain').empty()
            A closed plain spherical web with 0 edges.
        """
        return self.element_class({},{},[],self)

#### End of Parent ####

#### Start of generators (work in progress) ####

############################################################
# I am not clear on how this should be implemented.
#
# The intuition comes from the spiders which arise from the
# category of finite dimensional representations of a quantised
# enveloping algebra. Here a strand is modelled on an irreducible
# representation has a 'name', and probably
# a 'grading' which is either 'odd' or 'even'. Each instance of
# the parent class SphericalSpider should have an attribute
# which defines a (finite) set of strands and this set has an involution.
# This involution is modelled on taking the dual of a representaion.
# A strand is self-dual if it is fixed by this involution and a
# self-dual representation is either 'orthogonal' or 'symplectic'.
#
# Each halfedge() of a web will have an attribute which
# is an instance of the Decoration class. One of the entries of the
# Decoration class is a Strand. We require that this is an element
# of the strand set of the parent of the web. The mechanism for
# enforcing this is that we only construct vertices that satisfy
# this condition and we only construct webs using rotate and glue.
# For example, :func:``polygon_web`` should be rewritten to become a
# method that uses rotate and glue to build the polygon.
#
# In python we have @dataclass, nametuple, NamedTuple
# which all seem similar.
#
# In Sage we have set, Set, FiniteEnumeratedSet, FiniteFamily, ...
# and I am not clear what benefits there are to each of these.
#
# I should probably steal ideas from IndexedFreeMonoid
###############################################################

@dataclass
class Strand:
    name: str
    grading: bool
    parity: bool

class Strands():
    """
    This defines a set with involution.
    """
    def __init__(self,objs: set, duals: dict=None):
        if any(notinstance(a,Strand) for a in objs):
            raise ValueError(f"entries of {objs} must all be a Strand")
        duals = {str(x):str(y) for x,y in duals.items()}
        objs = set(str(a) for a in objs)
        if not duals.keys().issubset(objs):
            raise ValueError("domain of dual map is not a subset of halfedges")
        if not duals.values().issubset(objs):
            raise ValueError("range of dual map is not a subset of halfedges")
        if any(duals[duals[a]] != a):
            raise ValueError("the dual map must be an involution")

        self.objects = objs
        duals.update({a:a for a in objs if not a in duals})
        self.duals = duals

@dataclass
class Decoration:
    r"""
    This class that should probably follow IndexedGenerators.

    This should record the following information.

    - a direction, either 'in' or 'out'.
    - a label, an instance of :class:``Strand``.
    - possibly some marking to deal with symmetry
    """
    label: Strand
    direction: bool
    marking: int

#### End of generators (work in progress) ####

#### Start of functions ####

############################################################
# These should all be moved to become become parent methods.
# I have not done this because glue is not working properly;
# these functions use glue and are not used elsewhere.

def polygon_web(n: int):
    r"""Construct a polygon with n sides.

    EXAMPLES::

        sage: polygon_web(3)
        The plain spherical web with c = (3, 5, 7, 4, 0, 6, 1, 8, 2) and e = (6, 7, 8, 3, 4, 5).
        sage: polygon_web(1)
        The plain spherical web with c = (1, 2, 0) and e = (2, 1).
    """
    if n<1: raise ValueError
    a =  [halfedge() for i in range(n)]
    b1 = [halfedge() for i in range(n)]
    b2 = [halfedge() for i in range(n)]
    c = {a[i]:b1[i] for i in range(n)}
    c.update({b1[i]:b2[i] for i in range(n)})
    c.update({b2[i]:a[i] for i in range(n)})
    e = {b1[i-1]:b2[i] for i in range(n)}
    e.update({b2[i]:b1[i-1] for i in range(n)})
    return SphericalSpider('plain')(c,e,a)



def enumerate_diagrams(generators,max: int,forbidden=[],connected=False):
    """
    Enumerate diagrams.
    """
    old = set(generators)
    new = set()
    diagrams = old
    while len(old) != 0:
        for a in old:
            for b in generators:
                la = len(a.boundary)
                lb = len(b.boundary)
                r = 1 if connected else 0
                s = max(r,(la+lb-max)/2)
                p = min(la,lb)
                for k in range(r,s):
                    c = a.glue(b,k)
                    # Check c
                    new.add(c)
        diagrams.update(new)
        old = new
        new = set()
    return diagrams
