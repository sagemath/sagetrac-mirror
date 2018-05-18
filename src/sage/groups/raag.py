r"""
Right-Angled Artin Groups

A *right-angled Artin group* (often abbreviated as RAAG) is a group which
has a presentation whose only relations are commutators between generators.
These are also known as graph groups, since they are (uniquely) encoded by
(simple) graphs, or partially commutative groups.

AUTHORS:

- Travis Scrimshaw (2013-09-01): Initial version
- Travis Scrimshaw (2018-02-05): Made compatible with
  :class:`~sage.groups.artin.ArtinGroup`
"""

#****************************************************************************
#       Copyright (C) 2013,2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import division, absolute_import, print_function
import six

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp
from sage.groups.finitely_presented import FinitelyPresentedGroup, FinitelyPresentedGroupElement
from sage.groups.free_group import FreeGroup
from sage.groups.artin import ArtinGroup, ArtinGroupElement
from sage.graphs.graph import Graph
from sage.combinat.root_system.coxeter_matrix import CoxeterMatrix
from sage.combinat.root_system.coxeter_group import CoxeterGroup

class RightAngledArtinGroup(ArtinGroup):
    r"""
    The right-angled Artin group defined by a graph `G`.

    Let `\Gamma = \{V(\Gamma), E(\Gamma)\}` be a simple graph.
    A *right-angled Artin group* (commonly abbreviated as RAAG) is the group

    .. MATH::

        A_{\Gamma} = \langle g_v : v \in V(\Gamma)
        \mid [g_u, g_v] \text{ if } \{u, v\} \notin E(\Gamma) \rangle.

    These are sometimes known as graph groups or partially commutative groups.
    This RAAG's contains both free groups, given by the complete graphs,
    and free abelian groups, given by disjoint vertices.

    .. WARNING::

        This is the opposite convention of some papers.

    Right-angled Artin groups contain many remarkable properties and have a
    very rich structure despite their simple presentation. Here are some
    known facts:

    - The word problem is solvable.
    - They are known to be rigid; that is for any finite simple graphs
      `\Delta` and `\Gamma`, we have `A_{\Delta} \cong A_{\Gamma}` if and
      only if `\Delta \cong \Gamma` [Dro1987]_.
    - They embed as a finite index subgroup of a right-angled Coxeter group
      (which is the same definition as above except with the additional
      relations `g_v^2 = 1` for all `v \in V(\Gamma)`).
    - In [BB1997]_, it was shown they contain subgroups that satisfy the
      property `FP_2` but are not finitely presented by considering the
      kernel of `\phi : A_{\Gamma} \to \ZZ` by `g_v \mapsto 1` (i.e. words of
      exponent sum 0).
    - `A_{\Gamma}` has a finite `K(\pi, 1)` space.
    - `A_{\Gamma}` acts freely and cocompactly on a finite dimensional
      `CAT(0)` space, and so it is biautomatic.
    - Given an Artin group `B` with generators `s_i`, then any subgroup
      generated by a collection of `v_i = s_i^{k_i}` where `k_i \geq 2` is a
      RAAG where `[v_i, v_j] = 1` if and only if `[s_i, s_j] = 1` [CP2001]_.

    The normal forms for RAAG's in Sage are those described in [VW1994]_ and
    gathers commuting groups together.

    INPUT:

    - ``G`` -- a graph
    - ``names`` -- a string or a list of generator names

    EXAMPLES::

        sage: Gamma = Graph(4)
        sage: G = RightAngledArtinGroup(Gamma)
        sage: a,b,c,d = G.gens()
        sage: a*c*d^4*a^-3*b
        v0^-2*v1*v2*v3^4

        sage: Gamma = graphs.CompleteGraph(4)
        sage: G = RightAngledArtinGroup(Gamma)
        sage: a,b,c,d = G.gens()
        sage: a*c*d^4*a^-3*b
        v0*v2*v3^4*v0^-3*v1

        sage: Gamma = graphs.CycleGraph(5)
        sage: G = RightAngledArtinGroup(Gamma)
        sage: G
        Right-angled Artin group of Cycle graph
        sage: a,b,c,d,e = G.gens()
        sage: d*b*a*d
        v1*v3^2*v0
        sage: e^-1*c*b*e*b^-1*c^-4
        v2^-3

    We create the previous example but with different variable names::

        sage: G.<a,b,c,d,e> = RightAngledArtinGroup(Gamma)
        sage: G
        Right-angled Artin group of Cycle graph
        sage: d*b*a*d
        b*d^2*a
        sage: e^-1*c*b*e*b^-1*c^-4
        c^-3

    REFERENCES:

    - [Cha2006]_
    - [BB1997]_
    - [Dro1987]_
    - [CP2001]_
    - [VW1994]_

    - :wikipedia:`Artin_group#Right-angled_Artin_groups`
    """
    @staticmethod
    def __classcall_private__(cls, G, names=None):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: G1 = RightAngledArtinGroup(graphs.CycleGraph(5))
            sage: Gamma = Graph([(0,1),(1,2),(2,3),(3,4),(4,0)])
            sage: G2 = RightAngledArtinGroup(Gamma)
            sage: G3 = RightAngledArtinGroup([(0,1),(1,2),(2,3),(3,4),(4,0)])
            sage: G4 = RightAngledArtinGroup(Gamma, 'v')
            sage: G1 is G2 and G2 is G3 and G3 is G4
            True

        Handle the empty graph::

            sage: RightAngledArtinGroup(Graph())
            Traceback (most recent call last):
            ...
            ValueError: the graph must not be empty
        """
        if not isinstance(G, Graph):
            G = Graph(G, immutable=True)
        else:
            G = G.copy(immutable=True)
        if G.num_verts() == 0:
            raise ValueError("the graph must not be empty")
        if names is None:
            names = 'v'
        if isinstance(names, six.string_types):
            if ',' in names:
                names = [x.strip() for x in names.split(',')]
            else:
                names = [names + str(v) for v in G.vertices()]
        names = tuple(names)
        if len(names) != G.num_verts():
            raise ValueError("the number of generators must match the"
                             " number of vertices of the defining graph")
        return super(RightAngledArtinGroup, cls).__classcall__(cls, G, names)

    def __init__(self, G, names):
        """
        Initialize ``self``.

        TESTS::

            sage: G = RightAngledArtinGroup(graphs.CycleGraph(5))
            sage: TestSuite(G).run()
        """
        self._graph = G
        F = FreeGroup(names=names)
        CG = Graph(G).complement()  # Make sure it's mutable
        CG.relabel()  # Standardize the labels
        cm = [[-1]*CG.num_verts() for _ in range(CG.num_verts())]
        for i in range(CG.num_verts()):
            cm[i][i] = 1
        for u,v in CG.edge_iterator(labels=False):
            cm[u][v] = 2
            cm[v][u] = 2
        self._coxeter_group = CoxeterGroup(CoxeterMatrix(cm, index_set=G.vertices()))
        rels = tuple(F([i + 1, j + 1, -i - 1, -j - 1])
                     for i, j in CG.edge_iterator(labels=False))  # +/- 1 for indexing
        FinitelyPresentedGroup.__init__(self, F, rels)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: RightAngledArtinGroup(graphs.CycleGraph(5))
            Right-angled Artin group of Cycle graph
        """
        return "Right-angled Artin group of {}".format(self._graph)

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.gen(2)
            v2
        """
        return self.element_class(self, ([i, 1],))

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.gens()
            (v0, v1, v2, v3, v4)
            sage: Gamma = Graph([('x', 'y'), ('y', 'zeta')])
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.gens()
            (vx, vy, vzeta)
        """
        return tuple(self.gen(i) for i in range(self._graph.num_verts()))

    def ngens(self):
        """
        Return the number of generators of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.ngens()
            5
        """
        return self._graph.num_verts()

    def graph(self):
        """
        Return the defining graph of ``self``.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.graph()
            Cycle graph: Graph on 5 vertices
        """
        return self._graph

    @cached_method
    def one(self):
        """
        Return the identity element `1`.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G.one()
            1
        """
        return self.element_class(self, ())

    one_element = one

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.

        TESTS::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: elt = G([[0,3], [3,1], [2,1], [1,1], [3,1]]); elt
            v0^3*v3*v2*v1*v3
            sage: G(elt)
            v0^3*v3*v2*v1*v3
            sage: G(1)
            1
        """
        if isinstance(x, RightAngledArtinGroup.Element):
            if x.parent() is self:
                return x
            raise ValueError("there is no coercion from {} into {}".format(x.parent(), self))
        if x == 1:
            return self.one()
        verts = self._graph.vertices()
        x = [[verts.index(s[0]), s[1]] for s in x]
        return self.element_class(self, self._normal_form(x))

    def _normal_form(self, word):
        """
        Return the normal form of the word ``word``. Helper function for
        creating elements.

        EXAMPLES::

            sage: Gamma = graphs.CycleGraph(5)
            sage: G = RightAngledArtinGroup(Gamma)
            sage: G._normal_form([[0,2], [3,1], [2,1], [0,1], [1,1], [3,1]])
            ([0, 3], [3, 1], [2, 1], [1, 1], [3, 1])
            sage: a,b,c,d,e = G.gens()
            sage: a^2 * d * c * a * b * d
            v0^3*v3*v2*v1*v3
            sage: a*b*d == d*a*b and a*b*d == a*d*b
            True
            sage: a*c*a^-1*c^-1
            1
            sage: (a*b*c*d*e)^2 * (a*b*c*d*e)^-2
            1
        """
        pos = 0
        G = self._graph
        v = G.vertices()
        w = [list(x) for x in word]  # Make a (2 level) deep copy
        while pos < len(w):
            comm_set = [w[pos][0]]
            # The current set of totally commuting elements
            i = pos + 1

            while i < len(w):
                letter = w[i][0]  # The current letter
                # Check if this could fit in the commuting set
                if letter in comm_set:
                    # Try to move it in
                    if any(G.has_edge(v[w[j][0]], v[letter])
                           for j in range(pos + len(comm_set), i)):
                        # We can't, so go onto the next letter
                        i += 1
                        continue
                    j = comm_set.index(letter)
                    w[pos + j][1] += w[i][1]
                    w.pop(i)
                    i -= 1  # Since we removed a syllable
                    # Check cancellations
                    if w[pos + j][1] == 0:
                        w.pop(pos + j)
                        comm_set.pop(j)
                        i -= 1
                        if not comm_set:
                            pos = 0
                            # Start again since cancellation can be pronounced effects
                            break
                elif all(not G.has_edge(v[w[j][0]], v[letter])
                         for j in range(pos, i)):
                    j = 0
                    for x in comm_set:
                        if x > letter:
                            break
                        j += 1
                    w.insert(pos + j, w.pop(i))
                    comm_set.insert(j, letter)

                i += 1
            pos += len(comm_set)
        return tuple(w)

    class Element(ArtinGroupElement):
        """
        An element of a right-angled Artin group (RAAG).

        Elements of RAAGs are modeled as lists of pairs ``[i, p]`` where
        ``i`` is the index of a vertex in the defining graph (with some
        fixed order of the vertices) and ``p`` is the power.
        """
        def __init__(self, parent, lst):
            """
            Initialize ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: elt = G.prod(G.gens())
                sage: TestSuite(elt).run()
            """
            self._data = lst
            elt = []
            for i, p in lst:
                if p > 0:
                    elt.extend([i + 1] * p)
                elif p < 0:
                    elt.extend([-i - 1] * -p)
            FinitelyPresentedGroupElement.__init__(self, parent, elt)

        def __reduce__(self):
            """
            Used in pickling.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: elt = G.prod(G.gens())
                sage: loads(dumps(elt)) == elt
                True
            """
            P = self.parent()
            V = P._graph.vertices()
            return (P, ([[V[i], p] for i, p in self._data],))

        def _repr_(self):
            """
            Return a string representation of ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: a * b^2 * e^-3
                v0*v1^2*v4^-3
                sage: Gamma = Graph([('x', 'y'), ('y', 'zeta')])
                sage: G = RightAngledArtinGroup(Gamma)
                sage: x,y,z = G.gens()
                sage: z * y^-2 * x^3
                vzeta*vy^-2*vx^3
                sage: G.<a,b,c> = RightAngledArtinGroup(Gamma)
                sage: c * b^-2 * a^3
                c*b^-2*a^3
            """
            if not self._data:
                return '1'
            v = self.parent().variable_names()

            def to_str(name, p):
                if p == 1:
                    return "{}".format(name)
                else:
                    return "{}^{}".format(name, p)

            return '*'.join(to_str(v[i], p) for i, p in self._data)

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: latex(a*b*e^-4*d^3)
                \sigma_{0}\sigma_{1}\sigma_{4}^{-4}\sigma_{3}^{3}
                sage: latex(G.one())
                1
                sage: Gamma = Graph([('x', 'y'), ('y', 'zeta')])
                sage: G = RightAngledArtinGroup(Gamma)
                sage: x,y,z = G.gens()
                sage: latex(x^-5*y*z^3)
                \sigma_{\text{\texttt{x}}}^{-5}\sigma_{\text{\texttt{y}}}\sigma_{\text{\texttt{zeta}}}^{3}
            """
            if not self._data:
                return '1'

            from sage.misc.latex import latex
            latexrepr = ''
            v = self.parent()._graph.vertices()
            for i, p in self._data:
                latexrepr += "\\sigma_{{{}}}".format(latex(v[i]))
                if p != 1:
                    latexrepr += "^{{{}}}".format(p)
            return latexrepr

        def _mul_(self, y):
            """
            Return ``self`` multiplied by ``y``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: a * b
                v0*v1
                sage: b * a
                v1*v0
                sage: a*b*c*d*e
                v0*v1*v2*v3*v4
                sage: a^2*d*c*a*b*d
                v0^3*v3*v2*v1*v3
                sage: e^-1*a*b*d*c*a^-2*e*d*b^2*e*b^-3
                v4^-1*v0*v3*v1*v0^-2*v2*v1^-1*v4*v3*v4
            """
            P = self.parent()
            lst = self._data + y._data
            return self.__class__(P, P._normal_form(lst))

        def __pow__(self, n):
            """
            Implement exponentiation.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: elt = G.prod(G.gens())
                sage: elt**3
                v0*v1*v2*v3*v4*v0*v1*v2*v3*v4*v0*v1*v2*v3*v4
                sage: elt^-2
                v4^-1*v3^-1*v2^-1*v1^-1*v0^-1*v4^-1*v3^-1*v2^-1*v1^-1*v0^-1
                sage: elt^0
                1
            """
            P = self.parent()
            if not n:
                return P.one()

            if n < 0:
                lst = sum((self._data for i in range(-n)), ())  # Positive product
                lst = [[x[0], -x[1]] for x in reversed(lst)]  # Now invert
                return self.__class__(P, P._normal_form(lst))

            lst = sum((self._data for i in range(n)), ())
            return self.__class__(self.parent(), P._normal_form(lst))

        def __invert__(self):
            """
            Return the inverse of ``self``.

            TESTS::

                sage: Gamma = graphs.CycleGraph(5)
                sage: G = RightAngledArtinGroup(Gamma)
                sage: a,b,c,d,e = G.gens()
                sage: (a * b)^-2
                v1^-1*v0^-1*v1^-1*v0^-1
            """
            P = self.parent()
            lst = [[x[0], -x[1]] for x in reversed(self._data)]
            return self.__class__(P, P._normal_form(lst))

        def _richcmp_(self, other, op):
            """
            Compare ``self`` and ``other``.

            TESTS::

                sage: A = ArtinGroup(['B',3])
                sage: x = A([1, 2, 1])
                sage: y = A([2, 1, 2])
                sage: x == y
                True
                sage: x < y^(-1)
                True
                sage: A([]) == A.one()
                True
                sage: x = A([2, 3, 2, 3])
                sage: y = A([3, 2, 3, 2])
                sage: x == y
                True
                sage: x < y^(-1)
                True
            """
            return richcmp(self._data, other._data, op)

