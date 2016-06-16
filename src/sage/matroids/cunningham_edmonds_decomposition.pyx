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


    def __init__(self, G, pm, int f):
        self.G = G
        self.pm = pm
        self.f = f


    cdef get_graph(self):
        return self.G


    cdef get_parent_marker(self):
        return self.parent_marker


    cdef get_f(self):
        return self.f


    cdef set_f(self, int n):
        self.f = n
        return None


    cdef is_polygon(self):
        if not self.G.size() == self.G.order():
            return False
        if not self.G.num_components() == 1:
            return False
        return True


    cdef __relink1(self, Z=None, WQ=None, m=None):
        if len(Z) > 0:
            m_1 = Z[0]
            if len(Z) > 1:
                m_2 = Z[1]
                if len(Z) > 2:
                    #Throw error.
                    return False
        return None


    cdef __relink2(self, Z=None, WQ=None):
        if len(Z) > 0:
            m_1 = Z[0]
            if len(Z) > 1:
                m_2 = Z[1]
                if len(Z) > 2:
                    #Throw error.
                    return False
        return None



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
    def __init__(self, M=None):
        """
        Initialization of the decomposition.

        EXAMPLES::

            sage: from sage.matroids.advanced import *
        """
        if not M == None:
            2 + 2
            #Do something.

        else:
            self._arborescence = DiGraph()
            self._nodes = {}
            self._root = None                   # When we do define the aborescence to be something other than the empty DiGraph, we will use ``self._root = self._arborescence.sources()[0]``.
                                                # This should work, because our ``_arborescence`` should only have one source. We might want to check for that.



    cdef __relink1(Q, Z=None, WQ=None, m=None):
        """
        reorders ``Q`` for use in __hyperpath.

        INPUT:

        - ``Q`` -- a cycle
        - ``WQ`` -- a subset of the edges ``Q``
        - ``Z`` -- a subset of ``WQ``, containing at most two elements

        OUTPUT:

        A reordering of ``Q`` so that ``WQ`` - ``Z`` is a path an, one end incident to ``m``, the other incident to ``m_i`` (if it exists), and ``m_2`` (if it exitst) incident to ``m``.
        """
        Q.__relink1(Z=None, WQ=None, m=None)
        return None


    cdef __typing(D_hat, P, pi):
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

        """
        s = len(pi) - 1
        if s == 0:
                                                                                    #throw error
            return False                                                            # Fix this.
        #T1
        for H in pi(s):
            if D_hat._nodes[H].size() == D_hat._nodes[H].order():   # we are checking that H is a cycle,
                                        #H should never have any isolated verticies.
                D_hat._nodes[H].__relink1(P)
                # if T(H) == 5:                                                     #This is not code, but  a rather poor stub.
                #     return False
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
                    elif 3:
                        num_threes.add(H)
                    elif 4:
                        num_fours.add(H)
                    elif 5:
                        num_fives.add(H)
                    # if num_fives > 0:
                    #     return False
                    # if num_fours > 1:                                             # I think that this should be here, but I'm not sure, and it doesn't say so in the paper right here.
                    #     return False
                    if len(num_twos) + len(num_threes) + len(num_fours) + len(num_fives) > 2:
                        return False
                    if D_hat._nodes[Q].is_polygon():
                        D_hat._nodes[Q].__relink1()
                    # #R1
                    # if len(num_fives) + len(num_fours) + len(num_threes) + len(num_twos) == 0:
                    #     if T(H) == 2 or T(H) == 3:                                  # I need to figure out how to compute T(H)
                    #         if D_hat.K_1 == None:
                    #             D_hat.K_1 = Q
                    #         else:
                    #             D_hat.K_2 = Q
                    #     #R2
                    #     if T(H) == 4:
                    #         D_hat.K_1 = Q
                    #         D_hat.K_2 = Q
                    # #R3
                    # if len(num_threes) + len(num_twos) == 1 and len(num_fours) == 0 and T(H) == 4:
                    #     D_hat.K_2 = D_hat._arborescence.neighbors_in(D.K_1)[0]                    # I don't understand this, because I don't know why we are assuming that ``K_1`` is already assigned.
                    #                                                                 # This is the parent of ``K_1``

                    # # Here, we need to find (if possible an orientation `f'` and a corresponding new ``D_hat`` such that F_f'[H1, ..., Ht] is good--see (4.7)
                    # if T(H) == 5:                                                   #Again, figure out how to do this!
                    #     return False
            i = i -1


        return None


    cdef __relink2(Q, Z=None, WQ=None):
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


    cdef __hypopath(D, P):
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
        D_hat = D                                                                   # This needs to change so that ``D_hat`` is the reduced `t`-decomposition for ``P``.
        if D.size() == 1:
            Q = D.random_vertex()                                                   # This might not work, but we would need ``Q`` to be the only member of ``D``.
        if D.size() > 1:
            #H2
            # for G in D_hat:                                                       # Somehow, I need to build pi
            pi = {}
            worked = D_hat.__typing(P, pi)
            if worked == False:
                return False
            # Q = root D_hat                                                        # I need to decide how to find the root. This will also be needed to build pi.
            Q = []                                                                  # This is garbage, but I want a placeholder.
        #H3
        num_twos = 0
        num_threes = 0
        num_fours = 0
        num_fives = 0
        # for H in root(D_hat).children():                                            # update this.
        children = {}
        for H in children:
            n = 1                                                                   # Change to n = T(Hi)
            if 2:
                num_twos = num_twos + 1
            elif 3:
                num_threes = num_threes + 1
            elif 4:
                num_fours = num_fours + 1
            elif 5:
                num_fives = num_fives + 1
            if num_twos + num_threes + num_fours + num_fives > 2:
                return False
        if Q.size() == Q.order():
            Q.__relink2(P)
        # 
        # Find, if possible, an oriantation ``f'`` and corresponding ``D_hat`` such that ``P`` is a path of `Qf[H1, ... , Ht]`.--see (4.7).
        #R4
        if D.K_1 == None:
            D.K_1 = Q
            D.K_2 = Q
        else:
            D.K_2 = D._arborescence.neighbors_in(D.K_1)[0]
        # #R5
        # if not D.K_1 == None and D.K_1 == D.K_2:
        #     if P.ends() == D._nodes[K_1].parent_marker().ends():              This is what I want it to be, except I'm not sure yet how to get the vertices on the ends of the path and edge.
        #         D.K_1 = D._arborescence.neighbors_in(D.K_1)[0]
        #         D.K_2 = D.K_1
        #     There is another if statement that needs to be added as part of (R5), but I'm not sure which ``if`` if any it goes under, and its complicated, so I'm leaving it off for now.
        #     This question can probably be answered by reading Lemma 5.4, since the paper says that it is needed for the proof of Lemma 5.4.
        return None


    cdef __update(self, P, C):
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
        if C.size()-P.size() == 1:
            f = (set(C).difference(set(P)))[0]
        else:
            f = self._arborescence.size() + 1
            D_star.__add_cycle(({f} | set(C)).difference(set(P)))

        if D_star.K_1 == D_star.K_2:
            # U1
            # U1.1
            if D_star._node(D_star.K_1).size() > D_star._node(D_star.K_1).order():
                #add f so that it is joined to ``u_1`` and ``u_2`` in K_1
                2 + 2
            # U1.2
            elif D_star._node(D_star.K_1).size() == D_star._node(D_star.K_1).order():
                # Somehow, I need to figure out how to grab an edge between ``u_1`` and ``u_2``.
                # There's more to do here, but I'm going to leave it for now.
                2 + 2
            # U1.3
            else:
                f_1 = f + 1
                f_2 = f + 2
                # Do stuff
        else:
            # U2
            # Define R
            # for j in [s]:
                # Define m_j 
            j = 0
            s = 0                   # s and j should be defined when we define R.
            # U2.1
            # if 1 < j and j < s:
            #     if not self._node(J_j).is_polygon() and not self._node(J_j).order() == 2:   # In other words, if ``self._node(J_1)`` is prime 
            #         #Finish this later
            # U2.2
            # U2.3
            # U2.4
        return None

    cdef __is_graphic(self):
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



    cdef __merge(self, G):
        return None


    cdef merge(self):
        return None



    cdef __add_cycle(self, cycle):
        return None