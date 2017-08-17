r"""
Gammoids

Let `D` be a directed graph and let `E` and `T` be not necessarily disjoint
sets of vertices of `D`. Say a subset `X` of `E` is in a collection `I` if
`X` can be linked into a subset of `T`. This defines a gammoid `(E,I)`,
where `E` is the ground set of a matroid and `I` is its independent sets.

Some authors use a reverse convention, where instead of a set `T` of roots,
they have a starting set `S` that is linked into subsets of `E`. Here we use
the convention that the vertices `T` are at the end of the directed paths,
not the beginning.
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2017 Zachary Gershkoff <zgersh2@lsu.edu>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from .matroid import Matroid
from .utilities import sanitize_contractions_deletions, setprint_s

from sage.graphs.digraph import DiGraph
from .minor_matroid import MinorMatroid

from copy import copy, deepcopy

class Gammoid(Matroid):
    r"""
    The gammoid class.

    INPUT:

    - ``D`` -- A DiGraph representing the gammoid.
    - ``roots`` -- A subset of the vertices.
    - ``groundset`` -- (optional) A subset of the vertices. If not specified,
      the entire vertex set is used (and the gammoid will be strict).

    OUTPUT:

    An instance of ``Gammoid``.

    EXAMPLES::

        sage: from sage.matroids.gammoid import Gammoid
        sage: D = digraphs.TransitiveTournament(5)
        sage: M = Gammoid(D, roots=[3,4]); M
        Gammoid of rank 2 on 5 elements
        sage: M.is_isomorphic(matroids.Uniform(2,5))
        True
        sage: D.add_vertex(6)
        sage: N = Gammoid(D, roots=[3,4])
        sage: N.loops()
        frozenset({6})
        sage: O = Gammoid(D, roots=[3,4,6])
        sage: O.coloops()
        frozenset({6})
        sage: O.full_rank()
        3
        sage: P = Gammoid(D, roots=[3,4], groundset=[0,2,3]); P
        Gammoid of rank 2 on 3 elements
    """

    def __init__(self, D, roots, groundset=None):
        """
        See class definition for full documentation.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: edgedict = {1:[2], 2:[3], 4:[1,5], 5:[2,3,8], 6:[4,7], 7:[5,8]}
            sage: D = DiGraph(edgedict)
            sage: M = Gammoid(D, roots=[1,2,3])
            sage: N1 = Gammoid(D, groundset=range(1,8), roots=[1,2,3])
            sage: N2 = M.delete(8)
            sage: N1 == N2
            True
        """
        self._roots = frozenset(roots)
        vertices = frozenset(D.vertices())
        if not self._roots.issubset(vertices):
            raise ValueError("roots must be a subset of the vertices")

        if groundset is None:
            self._groundset = vertices
        else:
            self._groundset = frozenset(groundset)
            if not self._groundset.issubset(vertices):
                raise ValueError("ground set must be a subset of the vertices")

        # discard edge labels
        self._D = DiGraph(D.edges(labels=False))
        self._prune_vertices()
        self._G = self._D.copy(immutable=True)
        self._rootv = self._D.add_vertex()
        for v in roots:
            self._D.add_edge(v, self._rootv)

    def _prune_vertices(self):
        """
        Remove irrelevant vertices from the internal graph.

        This will remove vertices that are not part of the ground set and
        cannot be used in a valid path between an element and a root.
        However, this will not remove a cycle of such vertices.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: D = digraphs.TransitiveTournament(5)
            sage: D.add_vertex(6)
            sage: M = Gammoid(D, roots=[3,4], groundset=[0,1,4])
            sage: M.digraph().vertices()
            [0, 1, 2, 3, 4]
        """
        vertices_c = self._roots.union(self._groundset)
        extra_sources = set(self._D.sources()).difference(vertices_c)
        extra_sinks = set(self._D.sinks()).difference(vertices_c)
        while extra_sources or extra_sinks:
            self._D.delete_vertices(set(extra_sources).union(extra_sinks))
            extra_sources = set(self._D.sources()).difference(vertices_c)
            extra_sinks = set(self._D.sinks()).difference(vertices_c)

    def __hash__(self):
        r"""
        Return an invariant of the matroid.

        This function is called when matroids are added to a set. It is very
        desirable to override it so it can distinguish matroids on the same
        groundset, which is a very typical use case!

        .. WARNING::

            This method is linked to __richcmp__ (in Cython) and __cmp__ or
            __eq__/__ne__ (in Python). If you override one, you should (and in
            Cython: MUST) override the other!

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: D = digraphs.TransitiveTournament(5)
            sage: M1 = Gammoid(D, roots=[3,4])
            sage: M2 = Gammoid(D, roots=[2,3])
            sage: hash(M1) == hash(M2)
            False
            sage: N1 = M1.delete(2)
            sage: N2 = Gammoid(D, roots=[3,4], groundset=[0,1,3,4])
            sage: hash(N1) == hash(N2)
            True
            sage: D.delete_edge(0,4)
            sage: M3 = Gammoid(D, roots=[3,4])
            sage: hash(M1) == hash(M3)
            False
        """
        return hash((self._G, self._groundset, self._roots))

    def __eq__(self, other):
        """
        Compare two matroids.

        For two graphic matroids to be equal, all attributes of the underlying
        graphs must be equal.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        ``True`` if ``self`` and ``other`` have the same digraph, roots, and
        ground set; ``False`` otherwise.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: D = digraphs.TransitiveTournament(5)
            sage: M1 = Gammoid(D, roots=[3,4])
            sage: M2 = Gammoid(D, roots=[2,3])
            sage: M1 == M2
            False
            sage: N1 = M1.delete(2)
            sage: N2 = Gammoid(D, roots=[3,4], groundset=[0,1,3,4])
            sage: N1 == N2
            True
        """
        if not isinstance(other, Gammoid):
            return False
        # The roots are implied by self._D
        return (self._D == other._D and self._groundset == other._groundset)

    def __ne__(self, other):
        """
        Compare two matroids.

        INPUT:

        - ``other`` -- A matroid.

        OUTPUT:

        ``False`` if ``self`` and ``other`` have the same digraphs, roots, and
        ground set; ``True`` otherwise.

        EXAMPLES::

        sage: from sage.matroids.gammoid import Gammoid
        sage: edgedict = {1:[2], 2:[3], 4:[1,5], 5:[2,3,8], 6:[4,7], 7:[5,8]}
        sage: D = DiGraph(edgedict)
        sage: M = Gammoid(D, roots=[1,2,3])
        sage: D.relabel([1,2,3,4,5,6,7,'a'])
        sage: N = Gammoid(D, roots=[1,2,3])
        sage: M == N
        False
        sage: M.delete(8) == N.delete('a')
        True
        """
        return (not self == other)

    def digraph(self):
        """
        Return the DiGraph associated with the gammoid.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: edgelist = [(0, 4), (0, 5), (1, 0), (1, 4), (2, 0),
            ....: (2, 1), (2, 3), (2, 6), (3, 4), (3, 5), (4, 0), (5, 2), (6, 5)]
            sage: D = DiGraph(edgelist)
            sage: M = Gammoid(D, roots=[4,5,6])
            sage: M.digraph() == D
            True
        """
        return self._G.copy(data_structure='sparse')

    def digraph_plot(self):
        """
        Plot the graph with color-coded vertices.

        Vertices that are elements but not roots will be shown as blue. Vertices that
        are roots but not elements are red. Vertices that are both are pink, and vertices
        that are neither are grey.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: D = digraphs.TransitiveTournament(4)
            sage: M = Gammoid(D, roots=[2,3])
            sage: M.digraph_plot()
            Graphics object consisting of 11 graphics primitives
        """
        self._inter = frozenset(self._roots.intersection(self._groundset))
        self._buckets = frozenset(self._roots.difference(self._inter))
        self._ending = frozenset(self._groundset.difference(self._inter))
        self._therest = set(self._G.vertices()).difference(self._roots)
        self._therest = frozenset(self._therest.difference(self._groundset))
        # Vertices just in buckets set are red "#D55E00"
        # Vertices just in the ending set are blue "#0072B2"
        # Vertices in both are pink "#CC79A7"
        # Vertices in neither are grey "#999999"
        d = {"#D55E00": list(self._buckets), "#CC79A7": list(self._inter),
            "#0072B2": list(self._ending), "#999999": list(self._therest)}
        return self._G.plot(vertex_colors = d)

    def _rank(self, X):
        """
        Return the rank of a set.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: D = digraphs.TransitiveTournament(4)
            sage: M = Gammoid(D, roots=[2,3])
            sage: M.rank([2])
            1
            sage: M.full_rank()
            2
            sage: M.rank(M.groundset())
            2
        """
        source = self._D.add_vertex()
        for x in X:
            self._D.add_edge(source, x)

        rank = len(self._D.vertex_disjoint_paths(source, self._rootv))
        self._D.delete_vertex(source)
        return rank

    def groundset(self):
        """
        Return the ground set of the matroid.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: edgelist = [(0, 4), (0, 5), (1, 0), (1, 4), (2, 0),
            ....: (2, 1), (2, 3), (2, 6), (3, 4), (3, 5), (4, 0), (5, 2), (6, 5)]
            sage: D = DiGraph(edgelist)
            sage: M = Gammoid(D, roots=[4,5,6])
            sage: sorted(M.groundset())
            [0, 1, 2, 3, 4, 5, 6]
        """
        return self._groundset

    def _repr_(self):
        """
        Returns a string representation of the matroid.

        EXAMPLES:

            sage: from sage.matroids.gammoid import Gammoid
            sage: edgelist = [(0, 4), (0, 5), (1, 0), (1, 4), (2, 0),
            ....: (2, 1), (2, 3), (2, 6), (3, 4), (3, 5), (4, 0), (5, 2), (6, 5)]
            sage: D = DiGraph(edgelist)
            sage: M = Gammoid(D, roots=[4,5,6]); M
            Gammoid of rank 3 on 7 elements
        """
        self._mrank = str(self._rank(self._groundset))
        self._elts = str(len(self._groundset))

        return "Gammoid of rank " + self._mrank + " on " + self._elts + " elements"

    def _minor(self, contractions=frozenset([]), deletions=frozenset([])):
        """
        Return a minor.

        INPUT:

        - ``contractions`` -- frozenset, subset of ``self.groundset()`` to be contracted
        -  ``deletions`` -- frozenset, subset of ``self.groundset()`` to be deleted

        Assumptions: contractions are independent, deletions are coindependent,
        contractions and deletions are disjoint.

        OUTPUT:

        If there are contractions, a MinorMatroid. If not, a Gammoid.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: edgedict = {1:[2], 2:[3], 4:[1,5], 5:[2,3,8], 6:[4,7], 7:[5,8]}
            sage: D = DiGraph(edgedict)
            sage: M = Gammoid(D, roots=[1,2,3]); M
            Gammoid of rank 3 on 8 elements
            sage: N1 = M.delete([3,4])
            sage: N1.groundset()
            frozenset({1, 2, 5, 6, 7, 8})
        """
        # Unlike TransversalMatroid, we don't check for coloops to delete
        # because that would mean altering the roots
        if deletions:
            new_groundset = self.groundset().difference(deletions)
            N = Gammoid(self.digraph(), self._roots, new_groundset)
        else:
            N = self

        if contractions:
            return MinorMatroid(N, contractions=contractions, deletions=frozenset([]))
        else:
            return N

    def gammoid_extension(self, vertex, neighbors=None):
        """
        Return a gammoid extended by an element.

        The new element can be a vertex of the digraph that is not in the starting set,
        or it can be a new source vertex.

        INPUT::

        - ``vertex`` -- A vertex of the gammoid's digraph that is not already in the
          ground set, or a new vertex.
        - ``neighbors`` -- (optional) If ``vertex`` is not already in the graph,
          its neighbors will be specified. The new vertex will have in degree `0`
          regardless of this option.
        """
        pass
