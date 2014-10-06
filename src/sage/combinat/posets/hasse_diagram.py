r"""
Hasse diagrams of posets
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
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

from copy import copy
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.misc.misc import uniq
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

class HasseDiagram(DiGraph):
    """
    The Hasse diagram of a poset. This is just a transitively-reduced,
    directed, acyclic graph without loops or multiple edges.

    .. note::

       We assume that ``range(n)`` is a linear extension of the poset.
       That is, ``range(n)`` is the vertex set and a topological sort of
       the digraph.

    This should not be called directly, use Poset instead; all type
    checking happens there.

    EXAMPLES::

        sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
        sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]}); H
        Hasse diagram of a poset containing 4 elements
        sage: TestSuite(H).run()
    """
    def _repr_(self):
        r"""
        TESTS::
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H._repr_()
            'Hasse diagram of a poset containing 4 elements'
        """
        return "Hasse diagram of a poset containing %s elements"%self.order()

    def linear_extension(self):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extension()
            [0, 1, 2, 3]
        """
        # Recall: we assume range(n) is a linear extension.
        return range(len(self))

    def linear_extensions(self):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.linear_extensions()
            [[0, 1, 2, 3], [0, 2, 1, 3]]
        """
        return self.topological_sort_generator()

    def is_linear_extension(self,lin_ext=None):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_linear_extension(range(4))
            True
            sage: H.is_linear_extension([3,2,1,0])
            False
        """
        if lin_ext is None or lin_ext == range(len(self)):
            for x,y in self.cover_relations_iterator():
                if not x < y:
                    return False
            return True
        else:
            for x,y in self.cover_relations_iterator():
                if not lin_ext.index(x) < lin_ext.index(y):
                    return False
            return True

    # Could this be achieved by adding some options to
    # GenericGraph.plot, and just overriding graphics_array_defaults?

    def plot(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black', layout = "acyclic", **kwds):
        """
        Returns a Graphics object corresponding to the Hasse diagram.

        EXAMPLES::

            sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
            sage: elm_lbls = Permutations(3).list()
            sage: P = Poset(uc,elm_lbls)
            sage: H = P._hasse_diagram
            sage: levels = H.level_sets()
            sage: heights = dict([[i, levels[i]] for i in range(len(levels))])
            sage: type(H.plot(label_elements=True))
            <class 'sage.plot.graphics.Graphics'>

        ::

            sage: P = Posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [3,4,1,2])
            sage: P._hasse_diagram.plot()
            Graphics object consisting of 42 graphics primitives
        """
        # Set element_labels to default to the vertex set.
        if element_labels is None:
            element_labels = range(self.num_verts())

        # Create the underlying graph.
        graph = DiGraph(self)
        graph.relabel(element_labels)

        return graph.plot(layout = layout, **kwds)

    def show(self, label_elements=True, element_labels=None,
            label_font_size=12,label_font_color='black',
            vertex_size=300, vertex_colors=None,**kwds):
        """
        Shows the Graphics object corresponding to the Hasse diagram.
        Optionally, it is labelled.

        INPUT:


        -  ``label_elements`` - whether to display element
           labels

        -  ``element_labels`` - a dictionary of element
           labels


        EXAMPLES::

            sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
            sage: elm_lbls = Permutations(3).list()
            sage: P = Poset(uc,elm_lbls)
            sage: H = P._hasse_diagram
            sage: levels = H.level_sets()
            sage: heights = dict([[i, levels[i]] for i in range(len(levels))])
            sage: H.show(label_elements=True)
        """
        self.plot(label_elements=label_elements, element_labels=element_labels,
            label_font_size=label_font_size,label_font_color=label_font_color,
            vertex_size=vertex_size, vertex_colors=vertex_colors).show(**kwds)

    def cover_relations_iterator(self):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: list(H.cover_relations_iterator())
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        for u,v,l in self.edge_iterator():
            yield (u,v)

    def cover_relations(self,element=None):
        r"""
        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2,3], 1:[3,4], 2:[5], 3:[5], 4:[5]})
            sage: H.cover_relations()
            [(0, 2), (0, 3), (1, 3), (1, 4), (2, 5), (3, 5), (4, 5)]
        """
        return [c for c in self.cover_relations_iterator()]

    def is_lequal(self, i, j):
        """
        Returns True if i is less than or equal to j in the poset, and
        False otherwise.

        .. note::

            If the :meth:`lequal_matrix` has been computed, then this method is
            redefined to use the cached matrix (see :meth:`_alternate_is_lequal`).

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0, 1, 4
            sage: H.is_lequal(x,y)
            False
            sage: H.is_lequal(y,x)
            False
            sage: H.is_lequal(x,z)
            True
            sage: H.is_lequal(y,z)
            True
            sage: H.is_lequal(z,z)
            True
        """
        return i == j or \
                (i < j and j in self.breadth_first_search(i))

    def is_less_than(self, x, y):
        r"""
        Returns True if ``x`` is less than or equal to ``y`` in the
        poset, and False otherwise.

        TESTS::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0, 1, 4
            sage: H.is_less_than(x,y)
            False
            sage: H.is_less_than(y,x)
            False
            sage: H.is_less_than(x,z)
            True
            sage: H.is_less_than(y,z)
            True
            sage: H.is_less_than(z,z)
            False
        """
        if x == y:
            return False
        else:
            return self.is_lequal(x,y)

    def is_gequal(self, x, y):
        r"""
        Returns ``True`` if ``x`` is greater than or equal to ``y``, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: Q = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0,1,4
            sage: Q.is_gequal(x,y)
            False
            sage: Q.is_gequal(y,x)
            False
            sage: Q.is_gequal(x,z)
            False
            sage: Q.is_gequal(z,x)
            True
            sage: Q.is_gequal(z,y)
            True
            sage: Q.is_gequal(z,z)
            True
        """
        return self.is_lequal(y,x)

    def is_greater_than(self, x, y):
        """
        Returns ``True`` if ``x`` is greater than but not equal to
        ``y``, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: Q = HasseDiagram({0:[2], 1:[2], 2:[3], 3:[4], 4:[]})
            sage: x,y,z = 0,1,4
            sage: Q.is_greater_than(x,y)
            False
            sage: Q.is_greater_than(y,x)
            False
            sage: Q.is_greater_than(x,z)
            False
            sage: Q.is_greater_than(z,x)
            True
            sage: Q.is_greater_than(z,y)
            True
            sage: Q.is_greater_than(z,z)
            False
        """
        return self.is_less_than(y,x)

    def minimal_elements(self):
        """
        Returns a list of the minimal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P(0) in P.minimal_elements()
            True
            sage: P(1) in P.minimal_elements()
            True
            sage: P(2) in P.minimal_elements()
            True
        """
        indegs = self.in_degree(labels=True)
        return [x for x in indegs if indegs[x]==0]

    def maximal_elements(self):
        """
        Returns a list of the maximal elements of the poset.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.maximal_elements()
            [4]
        """
        outdegs = self.out_degree(labels=True)
        return [x for x,d in outdegs.iteritems() if d==0]

    def bottom(self):
        """
        Returns the bottom element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.bottom() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.bottom()
            0
        """
        min_elms = self.minimal_elements()
        if len(min_elms) == 1: return min_elms[0]
        return None

    def has_bottom(self):
        """
        Returns True if the poset has a unique minimal element.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4],4:[]})
            sage: P.has_bottom()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_bottom()
            True
        """
        if self.bottom() is not None: return True
        return False

    def top(self):
        """
        Returns the top element of the poset, if it exists.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.top() is None
            True
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.top()
            1
        """
        max_elms = self.maximal_elements()
        if len(max_elms) == 1: return max_elms[0]
        return None

    def has_top(self):
        """
        Returns ``True`` if the poset contains a unique maximal element, and
        ``False`` otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.has_top()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.has_top()
            True
        """
        if not self.top() is None: return True
        return False

    def is_bounded(self):
        """
        Returns True if the poset contains a unique maximal element and a
        unique minimal element, and False otherwise.

        EXAMPLES::

            sage: P = Poset({0:[3],1:[3],2:[3],3:[4,5],4:[],5:[]})
            sage: P.is_bounded()
            False
            sage: Q = Poset({0:[1],1:[]})
            sage: Q.is_bounded()
            True
        """
        return self.has_top() and self.has_bottom()

    def is_chain(self):
        """
        Returns True if the poset is totally ordered, and False otherwise.

        EXAMPLES::

            sage: L = Poset({0:[1],1:[2],2:[3],3:[4]})
            sage: L.is_chain()
            True
            sage: V = Poset({0:[1,2]})
            sage: V.is_chain()
            False

        TESTS:

        Check :trac:`15330`::

            sage: p = Poset(DiGraph({0:[1],2:[1]}))
            sage: p.is_chain()
            False
        """
        # Quick check: for |V| verts there must be |V|-1 edges.
        if self.num_edges()+1 != self.num_verts():
            return False
        # There is one minimum and all other vertices have out-degree 1
        seen_0 = False
        for d in self.out_degree():
            if d == 1:
                pass
            elif d == 0:
                if seen_0:
                    return False
                seen_0 = True
            else:
                return False

        # Maximum in-degree is 1
        return all(d<=1 for d in self.in_degree())

    def dual(self):
        """
        Returns a poset that is dual to the given poset.

        EXAMPLES::

            sage: P = Posets.IntegerPartitions(4)
            sage: H = P._hasse_diagram; H
            Hasse diagram of a poset containing 5 elements
            sage: H.dual()
            Hasse diagram of a poset containing 5 elements

        TESTS::

            sage: H = Posets.IntegerPartitions(4)._hasse_diagram
            sage: H.is_isomorphic( H.dual().dual() )
            True
            sage: H.is_isomorphic( H.dual() )
            False
        """
        H = self.reverse()
        H.relabel(perm=range(H.num_verts()-1,-1,-1), inplace=True)
        return HasseDiagram(H)

    def interval(self, x, y):
        """
        Return a list of the elements `z` of ``self`` such that
        `x \leq z \leq y`. The order is that induced by the
        ordering in ``self.linear_extension``.

        INPUT:

        -  ``x`` -- any element of the poset

        -  ``y`` -- any element of the poset

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: I = set([2,5,6,4,7])
            sage: I == set(H.interval(2,7))
            True
        """
        return [z for z in range(self.order())[x:y+1] if
                self.is_lequal(x,z) and self.is_lequal(z,y)]

    closed_interval = interval

    def open_interval(self, x, y):
        """
        Return a list of the elements `z` of ``self`` such that
        `x < z < y`. The order is that induced by the ordering in
        ``self.linear_extension``.

        EXAMPLES::

            sage: uc = [[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]]
            sage: dag = DiGraph(dict(zip(range(len(uc)),uc)))
            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram(dag)
            sage: set([5,6,4]) == set(H.open_interval(2,7))
            True
            sage: H.open_interval(7,2)
            []
        """
        ci = self.interval(x,y)
        if len(ci) == 0:
            return []
        else:
            return ci[1:-1]

    def rank_function(self):
        r"""
        Return the (normalized) rank function of the poset,
        if it exists.

        A *rank function* of a poset `P` is a function `r`
        that maps elements of `P` to integers and satisfies:
        `r(x) = r(y) + 1` if `x` covers `y`. The function `r`
        is normalized such that its minimum value on every
        connected component of the Hasse diagram of `P` is
        `0`. This determines the function `r` uniquely (when
        it exists).

        OUTPUT:

        - a lambda function, if the poset admits a rank function
        - ``None``, if the poset does not admit a rank function

        EXAMPLES::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: P.rank_function() is not None
            True
            sage: P = Poset(([1,2,3,4],[[1,4],[2,3],[3,4]]), facade = True)
            sage: P.rank_function() is not None
            True
            sage: P = Poset(([1,2,3,4,5],[[1,2],[2,3],[3,4],[1,5],[5,4]]), facade = True)
            sage: P.rank_function() is not None
            False
            sage: P = Poset(([1,2,3,4,5,6,7,8],[[1,4],[2,3],[3,4],[5,7],[6,7]]), facade = True)
            sage: f = P.rank_function(); f is not None
            True
            sage: f(5)
            0
            sage: f(2)
            0

        TESTS::

            sage: P = Poset([[1,3,2],[4],[4,5,6],[6],[7],[7],[7],[]])
            sage: r = P.rank_function()
            sage: for u,v in P.cover_relations_iterator():
            ...    if r(v) != r(u) + 1:
            ...        print "Bug in rank_function!"

        ::

            sage: Q = Poset([[1,2],[4],[3],[4],[]])
            sage: Q.rank_function() is None
            True

        test for ticket :trac:`14006`::

            sage: H = Poset()._hasse_diagram
            sage: s = dumps(H)
            sage: f = H.rank_function()
            sage: s = dumps(H)
        """
        if(self._rank_dict is None):
            return None
        return self._rank_dict.__getitem__ #the rank function is just the getitem of the dict

    def is_meet_semilattice(self):
        r"""
        Returns ``True`` if ``self`` has a meet operation, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_meet_semilattice()
            True

            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[]})
            sage: H.is_meet_semilattice()
            True

            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.is_meet_semilattice()
            False
        """
        try:
            self.meet_matrix()
        except ValueError:
            return False
        else:
            return True

    def is_join_semilattice(self):
        r"""
        Returns ``True`` if ``self`` has a join operation, and
        ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_join_semilattice()
            True
            sage: H = HasseDiagram({0:[2,3],1:[2,3]})
            sage: H.is_join_semilattice()
            False
            sage: H = HasseDiagram({0:[2,3],1:[2,3],2:[4],3:[4]})
            sage: H.is_join_semilattice()
            False
        """
        try:
            self.join_matrix()
        except ValueError:
            return False
        else:
            return True

    def is_distributive_lattice(self): # still a dumb algorithm...
        r"""
        Returns ``True`` if ``self`` is the Hasse diagram of a
        distributive lattice, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,3,2],1:[4],2:[4,5,6],3:[6],4:[7],5:[7],6:[7],7:[]})
            sage: H.is_distributive_lattice()
            False
            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3]})
            sage: H.is_distributive_lattice()
            True
            sage: H = HasseDiagram({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: H.is_distributive_lattice()
            False
        """
        try:
            jn = self.join_matrix()
            mt = self.meet_matrix()
        except ValueError:
            return False
        n = jn.ncols()
        for x in range(n):
            for y in range(n):
                for z in range(n):
                    if mt[x][jn[y][z]]!=jn[mt[x][y]][mt[x][z]]: return False
        return True

    def is_complemented_lattice(self):
        r"""
        Returns ``True`` if ``self`` is the Hasse diagram of a
        complemented lattice, and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: H.is_complemented_lattice()
            True

            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[4]})
            sage: H.is_complemented_lattice()
            False
        """
        try:
            jn = self.join_matrix()
            mt = self.meet_matrix()
        except ValueError:
            return False
        n = self.cardinality()
        c = [-1 for x in range(n)]
        for x in range(n):
            for y in range(x,n):
                if jn[x][y]==n-1 and mt[x][y]==0:
                    c[x]=y
                    c[y]=x
        return all([c[x]!=-1 for x in range(n)])

    def complements(self):
        r"""
        Return a list ``l`` such that ``l[i]`` is a complement of
        ``i`` in ``self``, or ``None`` if no such complement exists.

        A complement of ``x`` is an element ``y`` such that the meet
        of ``x`` and ``y`` is the bottom element of ``self`` and the
        join of ``x`` and ``y`` is the top element of ``self``.

        EXAMPLES::

            sage: from sage.combinat.posets.hasse_diagram import HasseDiagram
            sage: H = HasseDiagram({0:[1,2,3],1:[4],2:[4],3:[4]})
            sage: H.complements()
            [4, 3, 3, 2, 0]

            sage: H = HasseDiagram({0:[1,2],1:[3],2:[3],3:[4]})
            sage: H.complements()
            [4, None, None, None, 0]
        """
        jn = self.join_matrix()
        mt = self.meet_matrix()
        n = self.cardinality()
        c = [None for x in range(n)]
        for x in range(n):
            for y in range(x,n):
                if jn[x][y]==n-1 and mt[x][y]==0:
                    c[x]=y
                    c[y]=x
        return c

