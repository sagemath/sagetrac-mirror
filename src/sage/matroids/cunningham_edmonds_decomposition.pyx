r"""
Cunningham Edmons decomposition

These guys are characterized by an abhorence, ``T``, with nodes which are
graphs, and edge names that correspond to the parent marker edges. 

Construction
============



EXAMPLES::

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


    def __init__(self, G, pm, int f = -1):
        self.graph = G
        self.parent_marker = pm
        self.f = f


    cpdef get_graph(self):
        return self.graph


    cdef int get_parent_marker(self):
        return self.parent_marker


    cpdef get_parent_marker_edge(self):
        for e in self.graph.edges():
            if self.graph.edge_label(e) == self.parent_marker:
                return e


    cdef int get_f(self):
        return self.f


    cpdef set_f(self, int n):
        self.f = n
        return None


    cpdef is_polygon(self):
        if not self.graph.size() == self.graph.order():
            return False
        if not len(self.graph.connected_components()) == 1:
            return False
        if self.graph.size() == 2:
            return False
        return True

    cpdef is_path(self, P):
        H = self.get_graph().copy()
        S = set(H.edges()).difference(P)
        H.delete_edges(S)

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
        if set(P).size() == 0:
            return True

        H = self.get_graph().copy()
        S = set(H.edges()).difference(P)
        H.delete_edges(S)

        num_leaves=0
        for i in H.degree():
            if i > 2:
                return False
            if i == 1:
                return False
        return True


    cpdef typing(self, P):
        C = set(P) | set (self.get_parent_marker())
        if self.is_cycle(C):
            return 1
        if self.is_path(P):
            if self.is_path(C):
                return 3

            H = self.get_graph().copy()
            S = set(H.edges()).difference(P)
            H.delete_edges(S)

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
        """
        m_1 = set(())
        m_2 = set(())
        if len(Z) > 0:
            m_1 = set(Z[0])
            if len(Z) > 1:
                m_2 = set(Z[1])
        E = set(((self.edges().difference(set(WQ))).difference(Z)).difference(set(self.parent_marker)))
        n = self.graph.edges().size()
        G = Graph()
        V = self.vertices()
        i = 0
        for e in WQ:
            G.add_edge(V(i), V(i+1), (self.graph.edge_label(e)))
            i = i + 1
        for e in m_1:
            G.add_edge(V(i), V(i+1), (self.graph.edge_label(e)))
            i = i + 1
        for e in E:
            G.add_edge(V(i), V(i+1), (self.graph.edge_label(e)))
            i = i + 1
        for e in m_2:
            G.add_edge(V(i), V(i+1), (self.graph.edge_label(e)))
            i = i + 1
        G.add_edge(V(0), V(i), G.edge_label(self.graph.edge_label(self.parent_marker)))

        return G


    cpdef __relink2(self, Z=set(()), WQ=set(())):
        """
        reorders ``Q`` for use in __typing.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path an, if ``m_i`` exists (`i=1,2`), so that ``m_1`` is incident to an end of the path ``WQ`` - ``Z``
        """
        m_2 = set(())
        i = 0
        G = Graph()
        V = self.vertices()
        n = self.graph.edges().size()
        E = set(((self.edges()).difference(set(WQ))).difference(Z))
        if len(Z) > 0:
            i = 1
            G.add_edge(V(0), V(n - 1), self.graph.edge_label(Z[0]))
            if len(Z) > 1:
                m_2 = set(Z[1])
        for e in WQ:
            G.add_edge(V(i), V((i + 1) % n), (self.graph.edge_label(e)))
            i = i + 1
        for e in m_2:
            G.add_edge(V(i), V(i+1), (self.graph.edge_label(e)))
            i = i + 1
        for e in E:
            G.add_edge(V(i), V((i + 1) % n), (self.graph.edge_label(e)))
            i = i + 1

        return G







cdef class CunninghamEdmondsDecomposition(sage.structure.sage_object.SageObject):
    """
    

    INPUT:

    - ``M`` -- (default: ``None``) a ``r`` by ``c`` matrix.

    OUTPUT:

    - If the input is a matrix ``M``, return a ``CunninghamEdmondsDecomposition``
      instance representing ``M``.

    EXAMPLES::

        sage: from sage.matroids.advanced import *

    """

    # NECESSARY
    # Idealy, we would like to be able to construct this from either a matroid or a matrix, but we'll start with the matrix.
    def __init__(self, int next_edge=0, M=None):
        """
        Initialization of the decomposition.

        EXAMPLES::

            sage: from sage.matroids.cunningham_edmonds_decomposition import *
            sage: D = CunninghamEdmondsDecomposition(17)
            sage: D.get_arborescence()
            Digraph on 1 vertex
            sage: D.get_nodes()[0].get_graph()
            Multi-graph on 2 vertices
            sage: D.get_nodes()[0].get_graph().edges()
            [(0, 1, 0), (0, 1, 17)]
            sage: D.get_nodes()[0].get_graph().vertices()
            [0, 1]
        """
        if not M == None:
            2 + 2
            #Do something.

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


    cpdef __relink1(Q, Z=set(()), WQ=set(())):
        """
        reorders ``Q`` for use in __hyperpath.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path an, one end incident to ``m``, the other incident to ``m_i`` (if it exists), and ``m_2`` (if it exitst) incident to ``m``.
        """
        Q.__relink1(Z=set(()), WQ=set(()))
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
        sage: D = CunninghamEdmondsDecomposition(17)
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
        for H in pi(s):
            if D_hat.nodes[H].is_polygon():
                D_hat.nodes[H].__relink1(P)
                if D_hat.nodes[H].typing == 5:                                                     #This is not code, but  a rather poor stub.
                    return False
            i = s - 1
        #T2
        while i > 0:
            #T3
            for Q in pi[i]:
                children = {}                                                       #Needs to be fixed
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
                    #R1
                    if len(num_fives) + len(num_fours) + len(num_threes) + len(num_twos) == 0:
                        if D_hat.nodes[H].typing() == 2 or D_hat.nodes[H].typing() == 3:                                  # I need to figure out how to compute T(H)
                            if D_hat.K_1 == -1:
                                D_hat.K_1 = Q
                            else:
                                D_hat.K_2 = Q
                        #R2
                        if D_hat.nodes[H].typing() == 4:
                            D_hat.K_1 = Q
                            D_hat.K_2 = Q
                    #R3
                    if len(num_threes) + len(num_twos) == 1 and len(num_fours) == 0 and D_hat.nodes[H].typing() == 4:
                        D_hat.K_2 = D_hat.arborescence.neighbors_in(D_hat.K_1)[0]   # This is the parent of ``K_1``

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
        #H1
        D_hat = self
        for v in D_hat.arborescence.vertices():
            if len(set(D_hat.nodes[v].get_graph().edge_labels()) & P) == 0:
                D_hat.arborescence.delete_vertex(v)
        if D_hat.arborescence.order() == 1:
            Q = D_hat.root
        else:
            #H2
            # for G in D_hat:                                                       # Somehow, I need to build pi
            pi = D_hat.__get_pi()
            worked = D_hat.__typing(P, pi)
            if worked == False:
                return False
            Q = D_hat.root

        #H3
        num_twos = 0
        num_threes = 0
        num_fours = 0
        num_fives = 0
        children = D_hat.arborescence.neighbors_out(D_hat.root)
        for H in children:
            n = D_hat.__typing(P | D_hat.nodes[H].edges, D_hat.__get_pi())                                                                  # Change to n = T(Hi)
            if n == 2:
                num_twos = num_twos + 1
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


        # Find, if possible, an oriantation ``f'`` and corresponding ``D_hat`` such that ``P`` is a path of `Qf[H1, ... , Ht]`.--see (4.7).
        # for V in D_hat.arborescence.vertices():
        #     if D_hat.nodes[V]:
        # R4
        if D_hat.K_1 == -1:
            D_hat.K_1 = Q
            D_hat.K_2 = Q
        else:
            D_hat.K_2 = D_hat.arborescence.neighbors_in(D_hat.K_1)[0]
        #R5
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
        D_star = self
        

        # # This is all ugly and gross. The goal is just to find which verticies in ``P`` are the end vertices.
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
            # U1.1
            if D_star.nodes[K_1].get_graph().size() > D_star.nodes[K_1].get_graph().order():
                D_star.nodes[K_1].get_graph().add_edge(u_2, u_2, num_edges)
                num_edges = num_edges + 1
            # U1.2
            elif D_star.nodes[K_1].is_polygon() and D_star.nodes[K_1].get_graph().has_edge(u_1, u_2):
                f_p = D_star.nodes[K_1].get_graph().edge_label(u_1, u_2)
                f_pp = num_edges
                num_edges = num_edges + 2
                D_star.arborescence.add_edge(K_1, len(D_star.nodes))
                G = Graph([(0,1,f), (0,1,f_p), (0,1,f_pp)], multiedges=True)
                N = Node (G, f, 1)
                D_star.nodes.append(N)
            # U1.3
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


            for j in range(1, s + 1):
                # U2.1
                if 1 < j and j < s:
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
                    # U2.3
                    # U2.4
        return None

    cpdef __is_graphic(self):
        """
        Test if a binary matroid is graphic via Bixby and Wagner's paper

        INPUT:

        - ``M`` -- a tottally nonseparable `(0,1)`-matrix, with first column a singleton.

        OUTPUT:

        The conclusion that ``M`` is not realizable, or a realizing graph ``G`` for ``M``.        


        """
        # G1
        
        # G2
        # G3
        # G4
        return None



    cpdef __merge(self, N):
        G_N = self.nodes[N].graph
        m = self.nodes[N].parent_marker
        P = self.get_parent(self.nodes[N])
        G_P = self.nodes[P].graph

        for e in G_P.edges():
            if G_N.edge_label(m) == G_P.edge_label(e):
                f = e
        u = f[0]
        v = f[1]
        x = m[0]
        y = m[1]

        G = G_N.union(G_P)
        if self.nodes[N].f == -1:
            G.merge_vertices(u,x)           # The order of the vertices determines the name of the resulting vertex.
            G.merge_vertices(v,y)
        else:                               # This needs to change.
            G.merge_vertices(u,y)
            G.merge_vertices(v,x)


        return Node(G, P.parent_marker, P.get(f))

        return None


    cpdef merge(self):
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