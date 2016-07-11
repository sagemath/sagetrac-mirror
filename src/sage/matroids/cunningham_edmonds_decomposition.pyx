r"""
Cunningham Edmons decomposition

These guys are characterized by an aborescence, ``T``, with nodes which are
graphs, and edge names that correspond to the parent marker edges. 

Construction
============



EXAMPLES:

REFERENCES
==========

..  [BW88] \R. E. Bixby, D. K. Wagner, An Almost Linear-Time Algorithm for Graph Realization. In Mathematics of Operations Research (????), 99-123

AUTHORS:

- 

TESTS::



Methods
=======
"""
#*****************************************************************************
#       Copyright (C) 2016 Tara Fife <fi.tara@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from matroid cimport Matroid            # We'll need this for later.
from sage.matrix.matrix import Matrix
from sage.matrix.constructor import matrix
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.structure.sage_object import SageObject



cdef class Node(sage.structure.sage_object.SageObject):
    """
    

    INPUT:

    - ``G`` -- a graph with edge lables and multiedges allowed
    - ``pm`` -- an edge lable of an edge in ``G``
    - ``f`` -- and integer representing the way that ``self`` is glued to it's parrent marker.


    EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.get_graph()
            Multi-graph on 10 vertices
            sage: N.get_graph().edges()
            [(0, 1, 1),
             (0, 4, 2),
             (0, 5, 3),
             (1, 2, 4),
             (1, 6, 5),
             (2, 3, 6),
             (2, 7, 7),
             (3, 4, 8),
             (3, 8, 9),
             (4, 9, 10),
             (5, 7, 11),
             (5, 8, 12),
             (6, 8, 13),
             (6, 9, 14),
             (7, 9, 15)]
            sage: N.get_parent_marker()
            3
            sage: N.get_parent_marker_edge()
            (0, 5, 3)
            sage: N.get_f()
            -1
            sage: N.set_f(1)
            sage: N.get_f()
            1
            sage: N.is_polygon()
            False
            sage: N.is_path({1, 2, 3, 4})
            False
            sage: N.is_path({1, 3, 4, 7, 11})
            False
            sage: N.is_path({1, 2, 9, 12})
            False
            sage: N.is_path({1, 2, 4, 6})
            True
            sage: N.is_cycle({1, 2, 4, 6})
            False
            sage: N.is_cycle({1, 2, 4, 6, 8})
            True
            sage: N.typing({1, 4, 7, 11})
            1
            sage: N.typing({1, 2, 4, 7, 11})
            2
            sage: N.typing({1, 4, 6})
            3
            sage: N.typing({1, 4, 6, 12})
            4
            sage: N.typing({1, 2, 4, 8})
            5


            sage: G= Graph([(0, 1, 1), (1, 2, 2), (2, 3, 3), (3, 4, 4), (4, 5, 5), (5, 6, 6), (6, 7, 7), (7, 8, 8), (8, 9, 9), (9, 10, 10), (10, 11, 11), (11, 0, 0)], multiedges=True)
            sage: N = Node(G, 0)
            sage: N.is_polygon()
            True
            sage: N.__relink1({5, 7}, {2, 3, 5, 7, 11}).get_graph().edges()
            [(0, 1, 3),
             (0, 11, 0),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 5),
             (4, 5, 1),
             (5, 6, 4),
             (6, 7, 6),
             (7, 8, 8),
             (8, 9, 9),
             (9, 10, 10),
             (10, 11, 7)]
            sage: N.__relink1({5}, {2, 3, 5, 7, 11}).get_graph().edges()
            [(0, 1, 3),
             (0, 11, 0),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 7),
             (4, 5, 5),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
            sage: N.__relink1({}, {2, 3, 5, 7, 11}).get_graph().edges()
            [(0, 1, 3),
             (0, 11, 0),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 5),
             (4, 5, 7),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]


            sage: N.__relink2({5, 7}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 5),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 7),
             (4, 5, 0),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
            sage: N.__relink2({5}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 5),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 7),
             (4, 5, 0),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
            sage: N.__relink2({}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 10),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 5),
             (4, 5, 7),
             (5, 6, 0),
             (6, 7, 1),
             (7, 8, 4),
             (8, 9, 6),
             (9, 10, 8),
             (10, 11, 9)]
            sage: N.__relink2({5, 2}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 2),
             (1, 2, 11),
             (2, 3, 7),
             (3, 4, 5),
             (4, 5, 0),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]

    """

    def __init__(self, G, pm, int f = -1):
        r"""
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.get_graph()
            Multi-graph on 10 vertices
            sage: N.get_graph().edges()
            [(0, 1, 1),
             (0, 4, 2),
             (0, 5, 3),
             (1, 2, 4),
             (1, 6, 5),
             (2, 3, 6),
             (2, 7, 7),
             (3, 4, 8),
             (3, 8, 9),
             (4, 9, 10),
             (5, 7, 11),
             (5, 8, 12),
             (6, 8, 13),
             (6, 9, 14),
             (7, 9, 15)]
            sage: N.get_parent_marker()
            3
            sage: N.get_parent_marker_edge()
            (0, 5, 3)
            sage: N.get_f()
            -1
        """
        self.graph = G
        self.parent_marker = pm
        self.f = f


    cpdef get_graph(self):
        r"""
        EXAMPLES:
            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.get_graph()
            Multi-graph on 10 vertices
            sage: N.get_graph().edges()
            [(0, 1, 1),
             (0, 4, 2),
             (0, 5, 3),
             (1, 2, 4),
             (1, 6, 5),
             (2, 3, 6),
             (2, 7, 7),
             (3, 4, 8),
             (3, 8, 9),
             (4, 9, 10),
             (5, 7, 11),
             (5, 8, 12),
             (6, 8, 13),
             (6, 9, 14),
             (7, 9, 15)]
        """
        return self.graph


    cpdef get_parent_marker(self):
        r"""
        EXAMPLES:
    
            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.get_parent_marker()
            3
        """
        return self.parent_marker


    cpdef get_parent_marker_edge(self):
        r"""
        EXAMPLES:
            
            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.get_parent_marker_edge()
            (0, 5, 3)
        """
        for e in self.graph.edges():
            if e[2] == self.parent_marker:
                return e


    cpdef get_f(self):
        r"""
        EXAMPLES:
        
            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.get_f()
            -1
        """
        return self.f


    cpdef set_f(self, int n):
        r"""
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            
            sage: N.get_f()
            -1
            sage: N.set_f(1)
            sage: N.get_f()
            1
        """
        self.f = n
        return None


    cpdef is_polygon(self):
        r"""
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.is_polygon()
            False
        """
        if not self.graph.size() == self.graph.order():
            return False
        if not len(self.graph.connected_components()) == 1:
            return False
        if self.graph.size() == 2:
            return False
        return True

    cpdef is_path(self, P):
        r"""
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.is_path({1, 2, 3, 4})
            False
            sage: N.is_path({1, 3, 4, 7, 11})
            False
            sage: N.is_path({1, 2, 9, 12})
            False
            sage: N.is_path({1, 2, 4, 6})
            True
        """
        H = self.get_graph().copy()
        # S = set(H.edge_labels()).difference(P)
        for e in H.edges():
            if not e[2] in P:
                H.delete_edge(e)

        num_leaves=0
        for i in H.degree():
            if i > 2:
                return False
            if i == 1:
                num_leaves = num_leaves + 1
        if num_leaves > 2 or num_leaves == 0:
            return False
        return True


    cpdef is_cycle(self, P):
        r"""
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.is_cycle({1, 2, 4, 6})
            False
            sage: N.is_cycle({1, 2, 4, 6, 8})
            True 
        """
        if len(set(P)) == 0:
            return True

        H = self.get_graph().copy()
        for e in H.edges():
            if not e[2] in P:
                H.delete_edge(e)

        num_leaves=0
        for i in H.degree():
            if i > 2:
                return False
            if i == 1:
                return False
        return True


    cpdef typing(self, P):
        r"""
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (0, 4, 2), (0, 5, 3), (1, 2, 4), (1, 6, 5), (2, 3, 6), (2, 7, 7), (3, 4, 8), (3, 8, 9), (4, 9, 10), (5, 7, 11), (5, 8, 12), (6, 8, 13), (6, 9, 14), (7, 9, 15)], multiedges=True)
            sage: N = Node(G, 3)
            sage: N.typing({1, 4, 7, 11})
            1
            sage: N.typing({1, 2, 4, 7, 11})
            2
            sage: N.typing({1, 4, 6})
            3
            sage: N.typing({1, 4, 6, 12})
            4
            sage: N.typing({1, 2, 4, 8})
            5     
        """
        C = set(P) | {self.get_parent_marker()}
        if self.is_cycle(C):
            return 1
        if self.is_path(P):
            if self.is_path(C):
                return 3

            H = self.get_graph().copy()
            for e in H.edges():
                if not e[2] in C:
                    H.delete_edge(e)

            num_leaves=0
            for i in H.degree():
                if i == 1:
                    num_leaves = num_leaves + 1
            if num_leaves == 1:
                return 2
        if self.is_path(C):
            return 4
        return 5

    cpdef __relink1(self, Z=set(()), WQ=set(())):
        """
        reorders ``Q`` for use in __hyperpath.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path and, one end incident to the parent marker, the other incident to ``m_i`` (if it exists), and ``m_2`` (if it exitst) incident to ``m``.
    
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (1, 2, 2), (2, 3, 3), (3, 4, 4), (4, 5, 5), (5, 6, 6), (6, 7, 7), (7, 8, 8), (8, 9, 9), (9, 10, 10), (10, 11, 11), (11, 0, 0)], multiedges=True)
            sage: N = Node(G, 0)

            sage: N.__relink1({5, 7}, {2, 3, 5, 7, 11}).get_graph().edges()
            [(0, 1, 3),
             (0, 11, 0),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 5),
             (4, 5, 1),
             (5, 6, 4),
             (6, 7, 6),
             (7, 8, 8),
             (8, 9, 9),
             (9, 10, 10),
             (10, 11, 7)]
            sage: N.__relink1({5}, {2, 3, 5, 7, 11}).get_graph().edges()
            [(0, 1, 3),
             (0, 11, 0),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 7),
             (4, 5, 5),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
            sage: N.__relink1({}, {2, 3, 5, 7, 11}).get_graph().edges()
            [(0, 1, 3),
             (0, 11, 0),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 5),
             (4, 5, 7),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
        """
        m_1 = {}
        m_2 = {}
        S = set(WQ).difference(set(Z))
        if len(Z) > 0:
            m_1 = {Z.pop()}
            if len(Z) > 0:
                m_2 = {Z.pop()}
        E = set(self.graph.edge_labels()).difference(set(WQ)).difference({self.parent_marker})
        n = len(self.graph.edges())
        G = Graph()
        V = self.graph.vertices()
        i = 0
        for e in S:
            G.add_edge(V[i], V[i+1], e)
            i = i + 1
        for e in m_1:
            G.add_edge(V[i], V[i+1], e)
            i = i + 1
        for e in E:
            G.add_edge(V[i], V[i+1], e)
            i = i + 1
        for e in m_2:
            G.add_edge(V[i], V[i+1], e)
            i = i + 1
        G.add_edge(V[0], V[i], self.parent_marker)

        return Node(G, self.parent_marker, self.f)


    cpdef __relink2(self, Z=set(()), WQ=set(())):
        """
        reorders ``Q`` for use in __typing.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path an, if ``m_i`` exists (`i=1,2`), so that ``m_1`` is incident to an end of the path ``WQ`` - ``Z``
        
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: G= Graph([(0, 1, 1), (1, 2, 2), (2, 3, 3), (3, 4, 4), (4, 5, 5), (5, 6, 6), (6, 7, 7), (7, 8, 8), (8, 9, 9), (9, 10, 10), (10, 11, 11), (11, 0, 0)], multiedges=True)
            sage: N = Node(G, 0)
            
            sage: N.__relink2({5, 7}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 5),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 7),
             (4, 5, 0),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
            sage: N.__relink2({5}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 5),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 7),
             (4, 5, 0),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]
            sage: N.__relink2({}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 10),
             (1, 2, 2),
             (2, 3, 11),
             (3, 4, 5),
             (4, 5, 7),
             (5, 6, 0),
             (6, 7, 1),
             (7, 8, 4),
             (8, 9, 6),
             (9, 10, 8),
             (10, 11, 9)]
            sage: N.__relink2({5, 2}, {2, 3, 5, 7, 11}).edges()
            [(0, 1, 3),
             (0, 11, 2),
             (1, 2, 11),
             (2, 3, 7),
             (3, 4, 5),
             (4, 5, 0),
             (5, 6, 1),
             (6, 7, 4),
             (7, 8, 6),
             (8, 9, 8),
             (9, 10, 9),
             (10, 11, 10)]

        """
        i = 0
        G = Graph()
        V = self.graph.vertices()
        n = len(self.graph.edges())
        E = set(self.graph.edge_labels()).difference(set(WQ))
        S = set(WQ).difference(set(Z))
        if len(Z) > 0:
            G.add_edge(V[0], V[n - 1], Z.pop())
        for e in S:
            G.add_edge(V[i], V[(i + 1) % n], e)
            i = i + 1
        for e in Z:
            G.add_edge(V[i], V[i+1], e)
            i = i + 1
        for e in E:
            G.add_edge(V[i], V[(i + 1) % n], e)
            i = i + 1

        return G







cdef class CunninghamEdmondsDecomposition(sage.structure.sage_object.SageObject):
    """
    

    INPUT:

    - ``M`` -- (default: ``None``) a ``r`` by ``c`` matrix.

    OUTPUT:

    - If the input is a matrix ``M``, return a ``CunninghamEdmondsDecomposition``
      instance representing ``M``.

    EXAMPLES:

        sage: from sage.matroids.advanced import *

    """

    # NECESSARY
    # We are assuming that if ``T`` is specified, then so are ``nodes``, ``next_edge``, and ``next_vertex``.
    def __init__(self, T=None, nodes = None, int next_edge=0, next_vertex=None, K_1=-1, K_2=-1, u_1=-1, u_2=-1):
        """
        Initialization of the decomposition.

        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: D = CunninghamEdmondsDecomposition(next_edge=17)
            sage: D.get_arborescence()
            Digraph on 1 vertex
            sage: D.get_nodes()[0].get_graph()
            Multi-graph on 2 vertices
            sage: D.get_nodes()[0].get_graph().edges()
            [(0, 1, 0), (0, 1, 17)]
            sage: D.get_nodes()[0].get_graph().vertices()
            [0, 1]


            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: T = DiGraph([(0, 1)])
            sage: nodes = [Node(Graph([(0, 1, 1), (0, 1, 10), (0, 1, 18)], multiedges = True), 10, 1), Node(Graph([(2, 3, 18), (3, 4, 7), (4, 5, 8), (5, 6, 9), (6, 2, 11)], multiedges = True), 10, 1)]
            sage: D = CunninghamEdmondsDecomposition(T, nodes, Integer(19), Integer(7))


        """
        # This first option is mostly set up for testing right now. We are assuming that ``T`` is a digraph wich is a tree, whose
        # vertices are labled `0, \ldots, n'. And that ``nodes`` is an ordered list of the nodes which belong to the arboresence,
        # so that ``nodes[i]`` coresponds to the i'th vertex in ``T``, and that the root o ``T`` is 0.
        if not T == None:
            self.arborescence = T
            self.nodes = nodes
            self.root = self.arborescence.sources()[0]          # This should work, because our ``arborescence`` should only have one source. We might want to check for that.
            
            self.next_arborescence_vertex = T.order()
            self.next_vertex = next_vertex
            self.next_edge = next_edge
            self.K_1 = K_1
            self.K_2 = K_2
            self.u_1 = u_1
            self.u_2 = u_2

        else:
            self.arborescence = DiGraph()
            self.arborescence.add_vertex(0)
            N = Node(Graph([(0, 1, 0), (0, 1, next_edge)], multiedges=True), next_edge, 0)
            self.nodes = [N]
            self.root = 0                       # When we do define the aborescence to be something other than the empty DiGraph, we will use ``self.root = self.arborescence.sources()[0]``.
                                                # This should work, because our ``arborescence`` should only have one source. We might want to check for that.
            self.next_arborescence_vertex = 1
            self.next_vertex = 2
            self.next_edge = next_edge + 1
            self.K_1 = -1
            self.K_2 = -1
            self.u_1 = -1
            self.u_2 = -1


    cpdef relink1(self, Q, Z=set(()), WQ=set(())):
        """
        reorders ``Q`` for use in __hyperpath.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path an, one end incident to ``m``, the other incident to ``m_i`` (if it exists), and ``m_2`` (if it exitst) incident to ``m``.
        
        EXAMPLES:

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: T = DiGraph([(0, 1)])
            sage: nodes = [Node(Graph([(0, 1, 1), (0, 1, 10), (0, 1, 18)], multiedges = True), 10, 1), Node(Graph([(2, 3, 18), (3, 4, 7), (4, 5, 8), (5, 6, 9), (6, 2, 11)], multiedges = True), 18, 1)]
            sage: D = CunninghamEdmondsDecomposition(T, nodes, Integer(19), Integer(7))
            sage: D.relink1(1, {11}, {7, 8, 11})
            sage: D.get_nodes()[1].get_graph().edges()
            [(2, 3, 8), (2, 6, 18), (3, 4, 7), (4, 5, 11), (5, 6, 9)]
        """
        self.nodes[Q] = self.nodes[Q].__relink1(Z, WQ)
        return None


    cpdef __typing(self, P, pi):
        """
        INPUT:

        - ``D_hat`` -- a reduced `t`-decomposition, corresponding to ``D`` and ``P``
        - ``P`` -- a nonempty subset of the edges of ``m(D)``
        - ``pi`` -- a depth partition ``(pi_0, ... , pi_s)`` of ``D_hat``
          where ``G`` is in ``pi_k`` if the path from ``G`` to the root of
          ``D_hat`` has ``k`` arcs.

        OUTPUT:

        A reordering of ``D_hat``, and hence ``D`` such that each complete child
        ``H_i`` of ``root(D_hat)`` is good and `T(H_i)` is known.

        EXAMPLES:

        sage: from sage.matroids.cunningham_edmonds_decomposition import *
        sage: D = CunninghamEdmondsDecomposition(next_edge=17)
        sage: Dp = D.__hypopath(set([17]))
        sage: Dp.get_arborescence()
        Digraph on 1 vertex
        sage: Dp.get_nodes()[0].get_graph()
        Multi-graph on 2 vertices
        sage: Dp.get_nodes()[0].get_graph().edges()
        [(0, 1, 0), (0, 1, 17)]
        """
        D_hat = self
        for v in self.arborescence.vertices():
            if len(self.nodes(v).get_graph().edge_labels() & P) == 0:
                D_hat.arborescence.delete_vertex(v)
        s = len(pi) - 1
        if s == 0:
                                                                                    #throw error
            return False                                                            # Fix this.
        #T1
        # For each polygon `H\in \pi_s`, apply ``Relink1(H)``. If `T(H) = 5`
        # for some `H\in \pi_s`, stop--Lemma 4.4 implies ``P`` is not a hypopath
        # of `m(D)`. (Note that each `H\in \pi_s` is a polygon, bond, or prime,
        # and so the computation of `T(H)` is easy. The data structures required
        # to do this computation in time linear in `P\cap E(H)` are discussed in
        # section 6) set i = s-1.
        for H in pi(s):
            if D_hat.nodes[H].is_polygon():
                D_hat.nodes[H].__relink1(P)
                if D_hat.nodes[H].typing == 5:
                    return False
            i = s - 1
        #T2
        # If `i = 0`, stop--the desired ``D_hat`` has been found.
        while i > 0:
            #T3
            # Let `Q\in \pi_i`, let `H_i, \ldots, H_t` be the complete children of ``Q``,
            # and let `H\in Q_f[H_1, ldots, H_t]` where `f` is the current orientation of
            # ``D_hat``. If more than two of the `H_i` are not type 1, stop--Lemma 4.3 implies
            # ``P`` is not a hypopath of `m(D)`. If ``Q`` is a polygon, apply ``Relink1(Q)``.
            # Find, if possible, an orientation ``f'`` of and a corresponding new ``D_hat``
            # such that `Q_{f'}[H_1,\ldots , H_t]` is good--see (4.7). If `T(H) = 5`, stop.

            # Insert note about how things are slightly different with `Q = rood(D_hat)`.
            for Q in pi[i]:
                children = set(D_hat.arborescence.neighbors_out(Q))                                                       #Needs to be fixed
                num_twos = {}
                num_threes = {}
                num_fours = {}
                num_fives = {}
                for H in children:
                    n = 1                                                           #Change to n=T(H), when I figure out how to do this.
                    if 2:
                        num_twos.add(H)
                    elif n == 3:
                        num_threes.add(H)
                    elif n == 4:
                        num_fours.add(H)
                    elif n == 5:
                        num_fives.add(H)
                    if len(num_twos) + len(num_threes) + len(num_fours) + len(num_fives) > 2:
                        return False
                    if D_hat.nodes[Q].is_polygon():
                        D_hat.nodes[Q].__relink1()
                    # R1
                    # If `T(H_i) = 1` for all `i` and `T(H) = 2` or 3, either set `K_1 = Q`,
                    # if ``Q`` hasn't already been assigned, or set `K_2 = Q`.
                    if len(num_fives) + len(num_fours) + len(num_threes) + len(num_twos) == 0:
                        if D_hat.nodes[H].typing() == 2 or D_hat.nodes[H].typing() == 3:                                  # I need to figure out how to compute T(H)
                            if D_hat.K_1 == -1:
                                D_hat.K_1 = Q
                            else:
                                D_hat.K_2 = Q
                        # R2
                        # If `T(H_i) = 1` for all `i` and `T(H) = 4`, set ``K_1`` and ``K_2``
                        # equal to ``Q``.
                        if D_hat.nodes[H].typing() == 4:
                            D_hat.K_1 = Q
                            D_hat.K_2 = Q
                    # R3
                    # Suppose `T(H_k)` equals 2 or 3 for some `k`, `T(H_i) = 1` for all `i \not = k`, and `T(H) = 4`.
                    # Let ``K_2`` be the node on the unique path in ``D`` between ``Q`` and ``K_1, that is nearest
                    # ``K_1``  and contains the same end of ``P`` as ``Q``
                    if len(num_threes) + len(num_twos) == 1 and len(num_fours) == 0 and D_hat.nodes[H].typing() == 4:
                        D_hat.K_2 = Q
                        # for i in range[D_hat.arborescence.order()]:
                        # I just realized that I missed the last line "and contains the same end of ''P'' as ''Q''."
                        # This means that I have to think more deeply about what I'm doing.
                        # D_hat.arborescence.neighbors_in(D_hat.K_1)[0]   # This is the parent of ``K_1``

                    # Here, we need to find (if possible an orientation `f'` and a corresponding new ``D_hat`` such that F_f'[H1, ..., Ht] is good--see (4.7)
                    if D_hat.nodes[H].typing() == 5:                                                   #Again, figure out how to do this!
                        return False
            i = i -1


        return D_hat


    cpdef __relink2(Q, Z=set(()), WQ=set(())):
        """
        reorders ``Q`` for use in __typing.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path an, if ``m_i`` exists (`i=1,2`), so that ``m_1`` is incident to an end of the path ``WQ`` - ``Z``
        """
        Q.__relink1(Z=None, WQ=None)
        return None


    cpdef __hypopath(self, P):
        """
        Finds a hyperpath in a graph decomposition, or the conclusion that none exists
        INPUT:

        - ``D`` -- a `t`-decomposition
        - ``P`` -- a nonempty subset of the edges of ``m(D)``

        OUTPUT:

        A new `t`-decomposition ``D`` such that ``P`` is a path of `Qf[H_1, ... , H_t],
        or the conclusion that ``P`` is not a hypopath of ``D``.

        """
        # H1
        # Compute the reduced t-decomposition ``D_hat`` for ``P``. If ``D_hat`` has
        # just one member, let this member be ``Q`` and go to step H3.
        D_hat = self
        for v in D_hat.arborescence.vertices():
            if len(set(D_hat.nodes[v].get_graph().edge_labels()) & P) == 0:
                D_hat.arborescence.delete_vertex(v)
        if D_hat.arborescence.order() == 1:
            Q = D_hat.root
        else:
            # H2
            # Apply ``Typing``. If this results in the conclusion that ``P`` is not
            # a hypopath, stoop: otherwise, set ``Q`` equal to the root of ``D_hat``.
            pi = D_hat.__get_pi()
            worked = D_hat.__typing(P, pi)
            if worked == False:
                return False
            Q = D_hat.root

        #H3
        # Let `H_1, ldots, H_t` be the complete children of ``Q`` in ``D_hat``. If
        # more than two `H_i` are not of type 1, stop--``P`` is not a hypopath.
        # Apply ``Relink2(Q) if ``Q`` is a polygon. Find, if possible, and orientation
        # `f'` and corresponding ``D_hat`` such that ``P`` is a path of
        # `Q_{f'}[H_1, ldots, H_t]`.--see(4.7)
        num_twos = 0
        num_threes = 0
        num_fours = 0
        num_fives = 0
        children = D_hat.arborescence.neighbors_out(D_hat.root)
        for H in children:
            n = D_hat.__typing(P | D_hat.nodes[H].edges, D_hat.__get_pi())                                                                  # Change to n = T(Hi)
            if n == 2:
                num_twos = num_twos + 1
                # # This is meant to be used to figure out how to merge in a correct way.
                # G = graph()
                # for e in D_hat.nodes[H].get_graph().edges():
                #     if D_hat.nodes[H].get_graph().edge_label(e) in P:
                #         G.add_edge(e)
                # pm = D_hat.nodes[H].get_parent_marker_edge()
                
                # # We want f to be the vertex of the parent marker edge that has degree one in ``P``.
                # if G.degree(e[0]) == 1:
                #     f = e[0]
                # else:
                #     f = e[1]
                # D_hat.nodes[H].set_f(f)
            elif n == 3:
                num_threes = num_threes + 1
            elif n == 4:
                num_fours = num_fours + 1
            elif n == 5:
                num_fives = num_fives + 1
            if num_twos + num_threes + num_fours + num_fives > 2:
                return False
        if D_hat.nodes[Q].get_graph().size() == D_hat.nodes[Q].get_graph().order() and D_hat.nodes[Q].get_graph().size() > 2:
            D_hat.nodes[Q].__relink2(P)


        # # Find, if possible, an oriantation ``f'`` and corresponding ``D_hat`` such that ``P`` is a path of `Qf[H1, ... , Ht]`.--see (4.7).
        # for V in D_hat.arborescence.vertices():
        #     if D_hat.nodes[V]:
        # R4
        # If ``K_1`` has not been assigned, set ``K_1`` and ``K_2`` equal to ``Q``. Otherwise,
        # let ``K_1`` be the node on the unique path in ``D`` between ``Q`` and ``K_1`` and
        # contains the same end of ``P`` as ``Q``.
        if D_hat.K_1 == -1:
            D_hat.K_1 = Q
            D_hat.K_2 = Q
        else:
            D_hat.K_2 = D_hat.arborescence.neighbors_in(D_hat.K_1)[0] # This is also wrong!
        #R 5
        # Suppose ``K_1`` and ``K_2`` have been selected and `K_1 = K_2`. If the ends of ``P``
        # are the ends of `pm(K_1)` and ``K_1`` is not a bond, set ``K_1`` and ``K_2`` equal
        # to the parent of ``K_1``. If the ends of ``P`` are the ends of another marker ``m_i``
        # of ``K_1``, ``K_1`` is a polygon, and `m_i` is an edge of a bond `G`, set ``K_1`` and
        # ``K_2`` equal to ``G``.
        if not D_hat.K_1 == -1 and D_hat.K_1 == D_hat.K_2 and D_hat.nodes[D_hat.K_1].get_graph().order() > 2:
            if {D_hat.u_1, D_hat.u_2} == {D_hat.nodes[D_hat.K_1].get_parent_marker()[0], D_hat.nodes[D_hat.K_1].get_parent_marker()[1]}:
                D_hat.K_1 = D_hat.arborescence.neighbors_in(self.K_1)[0]
                D_hat.K_2 = D_hat.K_1
            if D_hat.nodes[D_hat.K_1].is_polygon():
                for v in D_hat.arborescence.neighbors_out(D_hat.K_1):
                    u1 = D_hat.nodes[v][0]
                    u2 = D_hat.nodes[v][1]
                    if {u1, u2} == {D_hat.nodes[D_hat.K_1].get_parent_marker()[0], D_hat.nodes[D_hat.K_1].get_parent_marker()[1]} and D_hat.nodes[v].get_graph().order() > 2:
                        D_hat.K_1 = v
                        D_hat.K_2 = v
        return D_hat


    # # This is used in ``__update``.
    # cpdef __squeeze(self, N, L):
    #     """
    #     Finds a hyperpath in a graph decomposition, or the conclusion that none exists
    #     INPUT:

    #     - ``L`` -- a 

    #     OUTPUT:

    #     A `t`-decomposition ``D*`` of the graph obtained from ``m(D)`` by adding the edges of `C-P` so that ``C`` is a cycle, and `C-P` is incident to ``m(D)`` at exactly two nodes.

    #     """

    cpdef __update(self, P, C):
        """
        Finds a hyperpath in a graph decomposition, or the conclusion that none exists
        INPUT:

        - ``D`` -- a `t`-decomposition
        - ``P`` -- a path of ``m(D)``
        - ``C`` -- a set such that `C\cap E(m(D)) = P`, and `C-P\not = \emptyset`
        - ``K_1`` -- node of ``D``
        - ``K_2`` -- node of ``D``
        - ``u_1`` -- node of ``K_1``
        - ``u_2`` -- node of ``K_2`` We are assuming that the conditions of Lemma 5.1 hold and that (R5) has been applied #Then, why don't we just apply (R5) here?

        OUTPUT:

        A `t`-decomposition ``D*`` of the graph obtained from ``m(D)`` by adding the edges of `C-P` so that ``C`` is a cycle, and `C-P` is incident to ``m(D)`` at exactly two nodes.

        """
        # U0
        # Set ``D_star`` equal to ``D``. If `|C-P| = 1`, set `\{f\} = (C - P)`; otherwise, let ``f`` be
        # a new element, form a polygon with edge-set `\{f\} \cup (C - P)`(the order of the edges is
        # irrelevent) and add this polygon to ``D_star``
        D_star = self
        

        # # This is all ugly and gross. The goal is just to find which verticies in ``P`` are the end vertices.
        # for N in D_hat.arborescence.vertices():
        #     if R.degree(D_hat.nodes[N].get_parent_marker()[0]) == 1:
        #         D_hat.nodes[N].set_f(D_hat.nodes[N].get_parent_marker()[0])
        #     else:
        #         D_hat.nodes[N].set_f(D_hat.nodes[N].get_parent_marker()[1])
        # G = D_hat.merge()
        # for v in D_hat.arborescence().vertices():
        #     for e in D_hat.nodes[v].get_graph().edges():
        #         if D_hat.nodes[v].get_graph().edge_label(e) in P:
        #             G.add_edge(e)
        # G.add_edges(P.edges())
        # V = G.vertices()
        # for v in V:
        #     if G.degree(v) == 2:
        #         V.remove(v)
        # D_hat.u_1 = V[0]
        # D_hat.u_2 = V[1]


        num_edges = 0
        for N in D_star.arborescence.vertices():
            num_edges = num_edges + D_star.nodes[N].get_graph().order()

        if len(C)-len(P) == 1:
            f = (set(C).difference(set(P)))[0]
        else:
            f = num_edges
            num_edges = num_edges + 1
            D_star.__add_cycle(({f} | set(C)).difference(set(P)))

        K_1 = D_star.K_1
        K_2 = D_star.K_2
        u_1 = D_star.u_1
        u_2 = D_star.u_2



        if K_1 == K_2:
            # U1
            # (Case `K_1 = K_2`) Aply the appropriate one of (U1.1)-(U1.3) and stop
            # U1.1
            # If ``K_1`` is not a polygon, join ``u_1`` and ``u_2`` by ``f`` in ``K_1``.
            if D_star.nodes[K_1].get_graph().size() > D_star.nodes[K_1].get_graph().order():
                D_star.nodes[K_1].get_graph().add_edge(u_2, u_2, num_edges)
                num_edges = num_edges + 1
            # U1.2
            # Suppose ``K_1`` is a polygon, and ``u_1`` and ``u_2`` are joined in ``K_1``
            # by ``f'`` in ``K_1``. Let ``f''`` be a new element. replace ``f'`` by f''``
            # in ``K_1`` and add a bond with edge-set `\{f, f', f''\}' to ``D_star``.
            elif D_star.nodes[K_1].is_polygon() and D_star.nodes[K_1].get_graph().has_edge(u_1, u_2):
                f_p = D_star.nodes[K_1].get_graph().edge_label(u_1, u_2)
                f_pp = num_edges
                num_edges = num_edges + 2
                D_star.arborescence.add_edge(K_1, len(D_star.nodes))
                G = Graph([(0,1,f), (0,1,f_p), (0,1,f_pp)], multiedges=True)
                N = Node (G, f, 1)
                D_star.nodes.append(N)
            # U1.3
            # Suppose ``K_1`` is a polygon and ``u_1`` and ``u_2`` are not adjacent. Then there are two
            # distinct paths in ``K_1``, say ``L_1`` and ``L_2`` joining ``u_1`` and ``u_2``. Let ``f_1``
            # and ``f_2`` be new elements. Delete ``K_1`` from ``D_star``, add polgyongs formed by
            # joining the ends of ``L_i`` with ``f_i`` (`i = 1, 2), and add a bond with edge set
            # `\{f, f_1, f_2\}`. 
            else:
                H = Graph(D_star.nodes[K_1].get_graph(), multiedges = False)
                L_1 = H.edge_disjoint_paths(u_2, u_1)[0]
                L_2 = H.edge_disjoint_paths(u_2, u_1)[1]
                f_1 = num_edges
                f_2 = num_edges + 1
                num_edges = num_edges + 2

                G_1 = Graph()
                G_1.add_edges(L_1.edges())
                G_1.add_edge((u_2, u_1, f_1))
                G_2 = Graph()
                G_2.add_edges(L_2.edges())
                G_2.add_edge((u_2, u_1, f_1))


                D_star.arborescence.add_edge(K_1, len(D_star.nodes))
                D_star.arborescence.add_edge(len(D_star.nodes), len(D_star.nodes) + 1)
                D_star.arborescence.add_edge(len(D_star.nodes), len(D_star.nodes) + 2)
                N = Node(Graph([(0,1,f), (0,1,f_1), ], multiedges=True), f, 1)
                D_star.nodes.append(N)
                D_star.nodes.append(Node(G_1, f_1, 1))
                D_star.nodes.append(Node(G_2, f_2, 1))



        else:
            # U2
            # (Case `K_1 \not = K_2) Let ``R`` be the unique path `K_1 = J_0, \ldots, J_s = K_2`
            # in ``D_star``. Let `\{m_j} = E(J_j) \cap E(J_{j+1})` (`0 \leq j \leq s - 1).
            # Apply (U2.1)-(U2.4)
            temp = Graph(D_star.arborescence.edges())
            R = D_star.arborescence.all_paths(K_1, K_2)[0]
            s = R.length()
            m = []
            J = []
            J[0] = K_1
            for j in range(1, s + 1):
                R.delete_edge(J[j-1], J[j])
                J[j] = D_star.neighbors(J[j - 1])[0]
                m[j] = (set(D_star.nodes[J[j - 1]].edge_labels()) & set(D_star.nodes[J[j]].edge_labels())).pop()          #There should only be one element in this intersection.


                # U2.1
                # Suppose `0 \leq j \leq s - 1`, `J_j` is prime and `\{m_{j - 1}, m_j\}` is a
                # cycle of ``J_j``, or ``J_j`` is a bond on at least 4 edges and `p(J_j)` is not
                # in ``R``. Let ``f'`` be a new element. Let ``J_j'`` be ``J_j`` with
                # `\{m_{j - 1}, m_j\}` deleted and ``f'`` added (with the same ends), and let ``B``
                # be a bond with edge-set `\{m_{j - 1}, f, m_j\}`. Replace ``J_j`` in ``R`` by
                # ``B``; delete ``J_j`` from ``D_star`` and replace it by ``J_j`` and ``B``.
                # (Note: There can be at most one ``J_j`` as above.)
            for j in range(1, s):
                if not self.node(J[j]).is_polygon() and not self.node(J[j]).order() == 2 and set(m[j - 1][0], m[j - 1][1]) == set(m[j][0], m[j][1]) or self.node(J[j]).order() == 2 and self.node(J[j]).size() > 3 and not R.has_vertex(self.get_parent(J[j])):
                    f_p = num_edges
                    J_p = D_star.nodes[J[j]]
                    for e in J_p.edges():
                        if J_p.edge_label(e) == m[j-1] or J_p.edge_label(e) == m[j]:
                            J_p.add_edge(e[0], e[1], f)
                            J_p.delete_edge(e)
                            #Finish this later

                    G = Graph([(0,1,f), (0,1,m[j - 1]), (0,1,m[j])], multiedges=True)
                    D_star.nodes[J[j]] = Node(G, f, 1)
                    D_star.arborescence.add_edge(J[j], D_star.next_arborescence_vertex)
                    D_star.nodes[D_star.next_arborescence_vertex] = Node(J_p, f, 1)
                    D_star.next_arborescence_vertex = D_star.next_arborescence_vertex + 1
            # U2.2
            # If ``K_1`` is a polygon, let ``L_1`` and ``L_2`` be the two paths joining
            # ``m_1`` and `u_1`. Apply ``Squeeze(L_i)`` (`i = 1, 2`). Do the same for
            # ``K_2``, if ``K_2`` is a polygon.
            if D_star.nodes[K_1].is_polygon():
                G = Graph(D_star.nodes[K_1])
                G.delete_edge(D_star.nodes[K_1].get_parent_marker_edge())
                L_1 = G.edge_disjoint_paths(D_star.nodes[K_1].get_parent_marker_edge()[0], u_1)[0]
                L_1 = set(L_1.edges())
                L_2 = set(G.edges()).difference(L_1)
                D_star.__squeeze(K_1, L_1)
                D_star.__squeeze(K_1, L_2)
            if D_star.nodes[K_2].is_polygon():
                G = Graph(D_star.nodes[K_2])
                G.delete_edge(D_star.nodes[K_2].get_parent_marker_edge())
                L_1 = G.edge_disjoint_paths(D_star.nodes[K_2].get_parent_marker_edge()[0], u_2)[0]
                L_1 = set(L_1.edges())
                L_2 = set(G.edges()).difference(L_1)
                D_star.__squeeze(K_2, L_1)
                D_star.__squeeze(K_2, L_2)
            # U2.3
            # For each internal ``J_j'' in ``R`` that is a polygon, let ``L_1`` and ``L_2`` be the
            # two components of `J_j - \{m_1, m_2,\}` and apply ``Squeeze(L_i)`` ('i = 1, 2').
            for j in range(1, s + 1):
                if D_star.nodes[J[j]].is_polygon():
                    L_1 = set(Graph(D_star.nodes[J[j]].get_graph().delete_edges(m[j], m[j + 1])).connected_components[0].edges())
                    L_2 = set(D_star.nodes[J[j]]).difference(L_1 | {m[j], m[j + 1]})
                    D_star.__squeeze(K_2, L_1)
                    D_star.__squeeze(K_2, L_2)
            # U2.4
            # Let ``G`` be `m(R)` with ``f`` joining ``u_1`` and ``u_2``. Delete R from ``D_star``
            # and add ``G``.
        return None

    cpdef __is_graphic(self):
        """
        Test if a binary matroid is graphic via Bixby and Wagner's paper

        INPUT:

        - ``M`` -- a tottally nonseparable `(0,1)`-matrix, with first column a singleton.

        OUTPUT:

        The conclusion that ``M`` is not realizable, or a realizing graph ``G`` for ``M``.        


        """
        # G1``M_1`` is realized by a bond ``G_1` with two edges. Set `D_1 = \{G_i\}` where
        # for `pm(G_1)` we choose the nontree edge of `G_i`. (Note that ``D_1`` is a
        # decomposition, but not a t-decomposition because ``G_i`` has only two edges. All
        # further ``D_j``, `j \geq 2`, will be t-decompositions.) If `c = 1` stop (``M`` is
        # graphic). Otherwise, set `j = 2`.
        
        # G2
        # Set `P = P_j` and `D = D_{j - 1}`. Apply Hypopath. If ``P`` is not a hypopath of
        # `m(D)`, stop.
       

        # G3
        # Set `C = C_j` and apply Update. (Note taht the quantities ``K_1``, ``K_2``,
        # ``u_1``, ``u_2`` required for input to UPDATE are calculated in HYPOPATH, and the
        # associated calls to TYPING, through application of rules (R1)-(R5).) Set
        # ``D_j = D_star``


        # G4
        # If `j = c`, set `G = m(D_j)` and stop; otherwise set `j = j + 1` and go to step G2.
        return None



    cpdef merge_with_parent(self, N, N_vertex=-1, P_vertex=-1):
        G_N = self.nodes[N].graph
        m = self.nodes[N].get_parent_marker_edge()
        P = self.get_parent(N)
        G_P = self.nodes[P].graph

        for e in G_P.edges():
            if m[2] == e[2]:            # If the lables are the same.
                f = e

        G = G_N.union(G_P)
        if N_vertex == -1:
            u = f[0]                    # Which vertices we are merging with which makes no difference.
            v = f[1]
            x = m[0]
            y = m[1]
        else:
            u = N_vertex
            v = {f[0], f[1]}.difference({N_vertex})
            x = P_vertex
            y = {f[0], f[1]}.difference({P_vertex})
        G.merge_vertices(u,x)           # The order of the vertices determines the name of the resulting vertex.
        G.merge_vertices(v,y)

        self.arborescence.merge_vertices(P, N)
        return Node(G, P.parent_marker, P.get(f))

        return None


    cpdef merge_branch(self, N):
        root = self.root
        n = self.arborescence.longest_path(self.root).size()

        pi = self.__get_pi()
        for i in range (1, n + 1):
            for N in pi[i]:
                self.__merge(N)

        return None

    cpdef __add_cycle(self, cycle):
        return None

    cpdef get_arborescence(self):
        return self.arborescence

    cpdef get_nodes(self):
        return self.nodes

    cpdef get_root(self):
        return self.root

    cpdef __get_pi(self):
        pi = []
        n = self.arborescence.longest_path(self.root).size()
        for i in range(n + 1):
            pi.append([])

        for v in self.arborescence:
            l = self.arborescence.shortest_path_length(self.root, v)
            pi[l].append(v)
        return pi

    cdef get_parent(self, N):
        return N.neighbors_in(N)[0]

    cpdef branch(self, N):
        T = self.arborescence.copy()
        for u in self.arborescence.vertices():
            if T.all_paths(N, u).size() == 0:
                T.delete_vertex(u)
        return T