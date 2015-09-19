r"""
Interface with TdLib (algorithms for tree decompositions)

This module defines functions based on TdLib, a
library that implements algorithms for tree
decompositions written by Lukas Larisch. 

**Definition** :

A `tree decomposition` of a graph `G` is a pair `(T, \beta)` consisting of a 
tree T and a function `\beta: V(T) \rightarrow 2^{V(G)}` associating with each node `t \in V(T)` a set
of vertices `\beta (t) \subseteq V(G)` such that 

(T1) for every edge `e \in E(G)` there is a node `t \in V(T)` with `e \subseteq \beta (t)`, and

(T2) for all `v \in V(G)` the set `\beta^{-1} := \{t \in V(T): v \in \beta (t)\}` is non-empty and connected in T.

The width of `(T, \beta)` is defined as `max\{|\beta (t)|-1: t \in V(T) \}`.
The treewidth of G is defined as the minimum width over all tree decompositions of G.

**Some known results** :

    - Trees have treewidth 1
    
    - Cycles have treewidth 2
    
    - Series-parallel graphs have treewidth at most 2
    
    - Cliques must be contained in some bag of a tree decomposition


Computing the treewidth or a tree decomposition of a given graph is NP-hard in general.

**This module containes the following functions** :

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`preprocessing_MD` | Applies save reduction rules followed by the minDegree-heuristic to a given graph
    :meth:`preprocessing_FI_TM` | Applies save reduction rules followed by the minFill-heuristic followed by triangulation minimization to a given graph

    :meth:`lower_bound` | Computes a lower bound with respect to treewidth of a given graph

    :meth:`CR_greedy_decomp` | Computes a tree decomposition of exact width of a given graph (faster than the dynamic version in practice)
    :meth:`CR_dynamic_decomp` | Computes a tree decomposition of exact width of a given graph (asymptotically faster than the greedy version)

    :meth:`seperator_algorithm` | Computes a 4-approximate tree decomposition of a given graph

    :meth:`MSVS` | Possibly reduces the width of a given tree decomposition with help of minimal seperating vertex sets
    :meth:`minimalChordal` | Possibly reduces the width of a given tree decomposition by triangulation minimization

    :meth:`is_valid_decomposition` | Checks, if a tree decomposition is valid with respect to a given graph
    :meth:`treedec_to_ordering` | Computes an elimination ordering out of a tree decomposition
    :meth:`ordering_to_treedec` | Computes a tree decomposition out of an elimination ordering
    :meth:`get_width` | Returns the width of a given tree decomposition


AUTHOR: Lukas Larisch (10-19-2015): Initial version

REFERENCE:

.. [1] P. D. Seymour and Robin Thomas. 1993. Graph searching and a min-max theorem for tree-width. J. Comb. Theory Ser. B 58, 1 (May 1993), 22-33.
      [`<http://dl.acm.org/citation.cfm?id=158568>`_]
.. [2] 


Methods
-------
"""

from libcpp.vector cimport vector

from tdlib cimport sage_preprocessing_MD
from tdlib cimport sage_preprocessing_FI_TM

from tdlib cimport sage_deltaC_min_d
from tdlib cimport sage_deltaC_max_d
from tdlib cimport sage_deltaC_least_c

from tdlib cimport sage_LBN_deltaC
from tdlib cimport sage_LBNC_deltaC
from tdlib cimport sage_LBP_deltaC
from tdlib cimport sage_LBPC_deltaC

from tdlib cimport sage_CR_greedy_decomp
from tdlib cimport sage_CR_dynamic_decomp

from tdlib cimport sage_seperator_algorithm

from tdlib cimport sage_is_valid_decomposition
from tdlib cimport sage_ordering_to_treedec
from tdlib cimport sage_treedec_to_ordering
from tdlib cimport sage_get_width

from sage.sets.set import Set
from sage.graphs.graph import Graph

include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'


#!!!!!!   NOTICE   !!!!!!!!
#Sage vertices have to be named by unsigned integers
#Sage bags of decompositions have to be lists of unsigned integers
#!!!!!!!!!!!!!!!!!!!!!!!!!!

##############################################################
############ GRAPH/DECOMPOSITION ENCODING/DECODING ###########
#the following will be used implicitly do the translation
#between Sage graph encoding and TdLib graph encoding,
#which is based on BGL

class TreeDecomposition(Graph):
    #This is just for the repr-message.
 
    def __repr__(self):
        r"""
        Returns a short string representation of self.

        EXAMPLE::

            sage: T = TreeDecomposition()
            sage: T
            Treedecomposition of width -1 on 0 vertices 
        """
        return "Treedecomposition of width " + str(get_width(self)) + " on " + str(self.order()) + " vertices"

cdef cython_make_tdlib_graph(G, vector[unsigned int] &V, vector[unsigned int] &E):
    V_python = G.vertices()
    for i in range(0, len(V_python)):
        V.push_back(V_python[i])

    E_python = G.edges()
    for i in range(0, len(E_python)):
        v,w,l = E_python[i]
        E.push_back(v)
        E.push_back(w)

cdef cython_make_tdlib_decomp(T, vector[vector[int]] &V, vector[unsigned int] &E):
    V_python = T.vertices()
    for i in range(0, len(V_python)):
        V.push_back(V_python[i])

    E_python = T.edges()
    for i in range(0, len(E_python)):
        v,w,l = E_python[i]
        E.push_back(V_python.index(v))
        E.push_back(V_python.index(w))


cdef cython_make_sage_graph(G, vector[unsigned int] &V, vector[unsigned int] &E):
    for i in range(0, len(V)):
        G.add_vertex(V[i])

    for i in range(0, len(E), 2):
        G.add_edge(V[E[i]], V[E[i+1]])


cdef cython_make_sage_decomp(G, vector[vector[int]] &V, vector[unsigned int] &E):
    for i in range(0, len(V)):
        G.add_vertex(Set(V[i]))

    for i in range(0, len(E), 2):
        G.add_edge(Set(V[E[i]]), Set(V[E[i+1]]))


##############################################################
############ PREPROCESSING ###################################

cdef cython_preprocessing_MD(G):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)

    cdef int c_lb = -1

    py_lb = sage_preprocessing_MD(V_G, E_G, V_T, E_T, c_lb)

    T = TreeDecomposition()
    cython_make_sage_decomp(T, V_T, E_T)

    print("Tree decomposition of width " + str(get_width(T)) + " computed")

    return T


def preprocessing_MD(G):
    """
    Returns a tree decomposition of exact width iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced 
    instance G' of G will be processed by the minDegree heuristic, which
    successivly eliminates a vertex of minimal degree. The returned tree
    decomposition then may be of non-optimal width. 

    INPUT:

    - ``G`` -- a generic graph

    OUTPUT:

    - A tree decomposition ``T`` of G

    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.show()
        sage: t = sage.graphs.tdlib.preprocessing_MD(g)
        Tree decomposition of width 4 computed
        sage: t
        Treedecomposition of width 4 on 10 vertices
        sage: t.show()
        sage: sage.graphs.tdlib.is_valid_decomposition(g, t)
        Valid tree decomposition
        sage: t.delete_vertices([sage.sets.set.Set([0,1,4,5])])
        sage: sage.graphs.tdlib.is_valid_decomposition(g, t)
        Invalid tree decomposition: not all vertices/edges covered
    """
    return cython_preprocessing_MD(G)


cdef cython_preprocessing_FI_TM(G):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)

    cdef int c_lb = -1

    sig_on()
    py_lb = sage_preprocessing_FI_TM(V_G, E_G, V_T, E_T, c_lb)
    sig_off()

    T = TreeDecomposition()
    cython_make_sage_decomp(T, V_T, E_T)

    print("Tree decomposition of width " + str(get_width(T)) + " computed")

    return T


def preprocessing_FI_TM(G):
    """
    Returns a tree decomposition of exact width, iff the treewidth of
    the given graph G is not greater than 3, otherwise the reduced 
    instance G' of G will be processed by the fillIn heuristic, which
    successivly eliminates a vertex, that will cause least new edges
    within the elimination process. The resulting tree decomposition 
    will be postprocessed by the minimalChordal algorithm, that may
    reduce the width of the tree decomposition.

    INPUT:

    - ``G`` -- a generic graph

    OUTPUT:

    - A tree decomposition ``T`` of G

    EXAMPLES::

        sage: g = graphs.FruchtGraph()
        sage: g.show()
        sage: t = sage.graphs.tdlib.preprocessing_FI_TM(g)
        Tree decomposition of width 3 computed
        sage: t
        Treedecomposition of width 3 on 12 vertices
        sage: t.show()
        sage: t.vertices()
        [{9, 2, 11, 5},
         {1, 2, 6, 7},
         {8, 9, 2, 11},
         {10, 11},
         {2, 11, 6, 7},
         {2, 11, 5, 6},
         {11},
         {10, 11, 5, 6},
         {0, 1, 6, 7},
         {10, 11, 6},
         {9, 2, 4, 5},
         {9, 2, 3, 4}]
        sage: sage.graphs.tdlib.is_valid_decomposition(g,t)
        Valid tree decomposition
    """
    return cython_preprocessing_FI_TM(G)



##############################################################
############ LOWER BOUNDS ####################################

cdef cython_lower_bound(G, algorithm):
    cdef vector[unsigned int] V, E
    cython_make_tdlib_graph(G, V, E)
    cdef int c_lb = 0

    sig_on()

    if(algorithm == "deltaC_min_d"):
        c_lb = sage_deltaC_min_d(V, E)
    elif(algorithm == "deltaC_max_d"):
        c_lb = sage_deltaC_max_d(V, E)
    elif(algorithm == "deltaC_least_c"):
        c_lb = sage_deltaC_least_c(V, E)
    elif(algorithm == "LBN_deltaC"):
        c_lb = sage_LBN_deltaC(V, E)
    elif(algorithm == "LBNC_deltaC"):
        c_lb = sage_LBNC_deltaC(V, E)
    elif(algorithm == "LBP_deltaC"):
        c_lb = sage_LBP_deltaC(V, E)
    elif(algorithm == "LBPC_deltaC"):
        c_lb = sage_LBPC_deltaC(V, E)
    else:
        print("Invalid lower bound algorithm")
        return -2

    sig_off()

    py_lb = c_lb

    return py_lb

def lower_bound(G, algorithm = "deltaC_least_c"):
    """
    Calls one of the following algorithms to compute a
    lower bound on the treewidth of a given graph:

        - deltaC_min_d
        - deltaC_max_d
        - deltaC_least_c
        - LBN_deltaC
        - LBNC_deltaC
        - LBP_deltaC
        - LBPC_deltaC

    INPUTS:

    - ``G`` -- a generic graph

    - ``algorithm`` -- (default: ``'deltaC_least_c'``) specifies the algorithm to use for computing a lower bound
                       on the treewidth of G. The algorithms have to be choosen among the above mentioned 

    OUTPUT:

    - A lower bound on the treewidth of ``G``

    EXAMPLES::
        sage: g = graphs.DurerGraph()
        sage: tw = 4
        sage: lb = -1

        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "deltaC_min_d"))
        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "deltaC_max_d"))
        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "deltaC_least_c"))
        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "LBN_deltaC"))
        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "LBNC_deltaC"))
        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "LBP_deltaC"))
        sage: lb = max(lb, sage.graphs.tdlib.lower_bound(g, "LBPC_deltaC"))

        sage: if(lb > tw):
        ....:     print("error")
    """

    return cython_lower_bound(G, algorithm)


##############################################################
############ EXACT ALGORITHMS ################################

cdef cython_CR_greedy_decomp(G, lb):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)

    cdef int c_lb = lb 

    sig_on()
    sage_CR_greedy_decomp(V_G, E_G, V_T, E_T, c_lb)
    sig_off()

    T = TreeDecomposition()
    cython_make_sage_decomp(T, V_T, E_T)

    print("tree decomposition of width " + str(get_width(T)) + " computed")

    return T


def CR_greedy_decomp(G, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound 
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed. This algorithm is faster than CR_dynamic_decomp in
    practical, but asymptotical slower.

    INPUTS:

    - ``G`` -- a generic graph

    - ``lb`` -- a lower bound to the treewidth of G, e.g. computed by lower_bound (default: ``'-1'``)

    OUTPUT:

    - A tree decomposition of G of tw(G), if the lower bound was not greater than tw(G), otherwise a tree decomposition of width = lb. 

..  WARNING::

    The computation can take a lot of time for graphs on more than about 15 vertices

EXAMPLES::

        sage: g = graphs.HouseGraph()
        sage: t = sage.graphs.tdlib.CR_greedy_decomp(g)
        tree decomposition of width 2 computed
        sage: t.show(vertex_size=2000)

TEST::

        sage: g = graphs.HouseGraph()
        sage: t = sage.graphs.tdlib.CR_greedy_decomp(g)
        tree decomposition of width 2 computed
    """

    return cython_CR_greedy_decomp(G, lb)


cdef cython_CR_dynamic_decomp(G, lb):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)

    cdef int c_lb = lb 

    sig_on()
    sage_CR_dynamic_decomp(V_G, E_G, V_T, E_T, c_lb)
    sig_off()

    T = TreeDecomposition()
    cython_make_sage_decomp(T, V_T, E_T)

    print("tree decomposition of width " + str(get_width(T)) + " computed")

    return T

def CR_dynamic_decomp(G, lb=-1):
    """
    Computes a tree decomposition of exact width, iff the given lower bound 
    is not greater than the treewidth of the input graph. Otherwise
    a tree decomposition of a width than matches the given lower bound
    will be computed.

    INPUTS:

    - ``G`` -- a generic graph

    - ``lb`` -- a lower bound to the treewidth of G, e.g. computed by lower_bound (default: ``'-1'``)

    OUTPUT:

    - A tree decomposition of G of tw(G), if the lower bound was not greater than tw(G), otherwise a tree decomposition of width = lb. 

..  WARNING::

    The computation can take a lot of time for graphs on more than about 15 vertices

EXAMPLES::

        sage: g = graphs.HouseGraph()
        sage: t = sage.graphs.tdlib.CR_greedy_decomp(g)
        tree decomposition of width 2 computed
        sage: t.show(vertex_size=2000)
    """

    return cython_CR_dynamic_decomp(G, lb)

##############################################################
############ APPROXIMATIVE ALGORITHMS ########################

cdef cython_seperator_algorithm(G):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)

    sig_on()
    sage_seperator_algorithm(V_G, E_G, V_T, E_T);
    sig_off()

    T = TreeDecomposition()
    cython_make_sage_decomp(T, V_T, E_T)

    print("Tree decomposition of width " + str(get_width(T)) + " computed")

    return T


def seperator_algorithm(G):
    """
    Computes a tree decomposition of a given graph using
    nearly balanced seperators. The returned width is at most 
    `4 \cdot tw(G)`.

    INPUT:

    - ``G`` -- a generic graph

    OUTPUT:

    - A tree decomposition of width at most `4 \cdot tw(G)`

EXAMPLES::

        sage: g = graphs.McGeeGraph()
        sage: t1 = sage.graphs.tdlib.seperator_algorithm(g)
        Tree decomposition of width 10 computed
        sage: t1
        Treedecomposition of width 10 on 14 vertices
        sage: t2 = sage.graphs.tdlib.preprocessing_MD(g)
        Tree decomposition of width 8 computed
        sage: t1.show(vertex_size=2000)
    """

    return cython_seperator_algorithm(G)


##############################################################
############ POSTPROCESSING ##################################


cdef cython_MSVS(G, T):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)
    cython_make_tdlib_decomp(T, V_T, E_T)

    old_width = get_width(T)

    sig_on()
    sage_MSVS(V_G, E_G, V_T, E_T)
    sig_off()

    T_ = TreeDecomposition()
    cython_make_sage_decomp(T_, V_T, E_T)

    new_width = get_width(T_)

    print("MSVS reduced the width by " + str(old_width-new_width) + ", new width: " + str(new_width))

    return T_


def MSVS(G, T):
    """
    This may reduce the maximal bag size of a tree decomposition.

    INPUTS:

    - ``G`` -- a generic graph

    - ``T`` -- a tree decomposition of ``G``

    OUTPUT:

    - A tree decomposition of ``G`` with possibly smaller width than ``T``

EXAMPLES::

        sage: g = graphs.RandomGNP(15, 0.1)
        sage: g.show()
        sage: t = sage.graphs.tdlib.trivial_decomposition(g)
        sage: t.show(vertex_size=2000)
        sage: t_ = sage.graphs.tdlib.MSVS(g, t)
        MSVS reduced the width by 12, new width: 2
        sage: t_.show(vertex_size=2000)
    """
    return cython_MSVS(G, T)


cdef cython_minimalChordal(G, O):
    cdef vector[unsigned int] V, E, old_elim_ordering, new_elim_ordering
    cython_make_tdlib_graph(G, V, E)

    for l in range(0, len(O)):
        old_elim_ordering.push_back(O[l])

    sig_on()
    sage_minimalChordal(V, E, old_elim_ordering, new_elim_ordering)
    sig_off()
    
    py_new_elim_ordering = []
    cdef i;
    for i in range(0, len(new_elim_ordering)):
        py_new_elim_ordering.append(new_elim_ordering[i])

    return py_new_elim_ordering

def minimalChordal(G, O):
    """
    Returns an alternativ elimination ordering to the given one, that
    may cause a lower width than the given one, when applied to the
    input graph for computing a tree decomposition.

    INPUTS:

    - ``G`` -- a generic graph

    - ``O`` -- an elimination ordering on ``G``

    OUTPUT:

    - An elimination ordering on ``G`` that may cause a lower width of the tree decomposition, that can be made out of
      it, than the width, that ``O`` will cause.

EXAMPLES::

        sage: g = graphs.HigmanSimsGraph()
        sage: g.show()
        sage: o1 = g.vertices()
        sage: t1 = sage.graphs.tdlib.ordering_to_treedec(g, o1)
        sage: t1.show(vertex_size=100)
        sage: o2 = sage.graphs.tdlib.minimalChordal(g, o1)
        sage: t2 = sage.graphs.tdlib.ordering_to_treedec(g, o2)
        sage: t2.show(vertex_size=100)
    """
    return cython_minimalChordal(G, O)


##############################################################
############ MISC ############################################

cdef cython_ordering_to_treedec(G, O):
    cdef vector[unsigned int] V_G, E_G, E_T, elim_ordering
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)

    for i in range(0, len(O)):
        elim_ordering.push_back(O[i])

    sage_ordering_to_treedec(V_G, E_G, V_T, E_T, elim_ordering)

    T = TreeDecomposition()
    cython_make_sage_decomp(T, V_T, E_T)

    return T

def ordering_to_treedec(G, O):
    """
    Applies an elimination ordering to a graph and returns 
    the resulting tree decomposition.

    INPUTS:

    - ``G`` -- a generic graph

    - ``O`` -- an elimination ordering on ``G``

    OUTPUT:

    - A tree decomposition, that has been made by applying ``O`` on ``G``

EXAMPLES::

        sage: g = graphs.RandomGNP(10, 0.05)
        sage: ordering = sage.graphs.tdlib.fillIn_ordering(g)
        sage: t = sage.graphs.tdlib.ordering_to_treedec(g, ordering)
        sage: t.show(vertex_size=2000)
    """

    return cython_ordering_to_treedec(G, O)

cdef cython_treedec_to_ordering(T):
    cdef vector[unsigned int] E_T, elim_ordering
    cdef vector[vector[int]] V_T

    cython_make_tdlib_decomp(T, V_T, E_T)

    sage_treedec_to_ordering(V_T, E_T, elim_ordering)

    py_elim_ordering = []
    cdef i;
    for i in range(0, len(elim_ordering)):
        py_elim_ordering.append(elim_ordering[i])

    return py_elim_ordering

def treedec_to_ordering(T):
    """
    Converts a tree decomposition to an elimination ordering.

    INPUT:

    - ``T`` -- a tree decomposition

    OUTPUT:

    - An elimination ordering, that fits T

EXAMPLES::

        sage: g = graphs.RandomGNP(10, 0.05)
        sage: t = sage.graphs.tdlib.seperator_algorithm(g)
        Tree decomposition of width 3 computed
        sage: ordering = sage.graphs.tdlib.treedec_to_ordering(t)
    """

    return cython_treedec_to_ordering(T)

cdef cython_is_valid_decomposition(G, T):
    cdef vector[unsigned int] V_G, E_G, E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_graph(G, V_G, E_G)
    cython_make_tdlib_decomp(T, V_T, E_T)

    cdef c_status;
    c_status = sage_is_valid_decomposition(V_G, E_G, V_T, E_T);
    
    py_status = c_status

    if(py_status == 0):
        print("Valid tree decomposition")
    elif(py_status == -1):
        print("Invalid tree decomposition: tree decomposition is not a tree")
    elif(py_status == -2):
        print("Invalid tree decomposition: not all vertices/edges covered")
    elif(py_status == -3):
        print("Invalid tree decomposition: some coded vertices are not connected in the tree decomposition")
    else:
        pass

def is_valid_decomposition(G, T):
    """
    Checks, if the definition of a tree decomposition holds for
    a tree decomposition and a graph.

    INPUTS:

    - ``G`` -- a generic graph

    - ``T`` -- a tree decomposition

EXAMPLES::

        sage: g = graphs.FranklinGraph()
        sage: t = sage.graphs.tdlib.seperator_algorithm(g)
        Tree decomposition of width 6 computed
        sage: sage.graphs.tdlib.is_valid_decomposition(g, t)
        Valid tree decomposition
    """

    cython_is_valid_decomposition(G, T)


cdef cython_get_width(T):
    cdef vector[unsigned int] E_T
    cdef vector[vector[int]] V_T

    cython_make_tdlib_decomp(T, V_T, E_T)

    return sage_get_width(V_T);

def get_width(T):
    """
    Returns the width (maximal size of a bag minus one) of a given tree decomposition.

    INPUT:

    - ``T`` -- a tree decomposition

    OUTPUT:

    - The width of ``T``

EXAMPLES::

        sage: g = graphs.RandomGNP(10, 0.05)
        sage: t = sage.graphs.tdlib.trivial_decomposition(g)
        sage: sage.graphs.tdlib.get_width(t)
        9
    """

    return cython_get_width(T)

