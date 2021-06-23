r"""
Elser Numbers and U-Nucleus Complexes

This module computes the Elser numbers and U-Nucleus complexes of a graph, as described in the paper [(DB)HLMNVW2021]_.

**This module contains the following methods:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`nuclei_by_size` | Computes all nuclei of the graph and groups them by size.
    :meth:`elser_number` | Computes the kth Elser number of the graph
    :meth:`nucleus_complex` | For a graph G and subset of vertices U, computes the U-Nucleus complex of G
    :meth:`all_nucleus_complexes` | For a graph G, computes all U-nucleus complexes associated to that graph
    :meth:`nucleus_complex_bettis` | Computes the Betti numbers of the U-nucleus complex, if no U is specified, computes all such complexes


Definition of Elser Numbers
-----------

Let 'G' be a graph. A nucleus of 'G' is a connected subgraph 'N' of 'G' such that N is a vertex cover; that is, every edge of G has at least one
endpoint in 'V(N)'. Let '\mathcal{N}(G)' denote the set of all nuclei of G. The 'k'th Elser number of G is

.. MATH::

    \text{els}_k(G) = (-1)^{\#V(G) + 1} \sum_{N \in \mathcal{N}(G)} (-1)^{\#E(N)} V(N)^k.

For example, if 'G' is the cycle graph on three elements, then the kth Elser number is

.. MATH::

    \text{els}_k(C_3) = 6(3^{k-1} - 2^{k-1}).

These numbers were originally introduced by Veit Elser [Els84]_, who conjectured that they were nonnegative for all graphs 'G'
and integers 'k\geq 2'. In [(DB)HLMNVW2021]_, Galen Dorpalen-Barry, Cyrus Hettle, David C. Livingston, Jeremy L. Martin, George Nasr, Julianne Vega, and Hays Whitlatch proved Elser's conjuecture using a chain complex, which they call the 'U'-nucleus complex.
This module allows you to compute the Elser numbers as well as the 'U'-nucleus complexes.

REFERENCES:

 - [(DB)HLMNVW2021]_
 - [Els84]_

Examples
-----------
In this section, we guide the reader through some initial examples of the functions in this module.
More details on the individual functions can be found in the Methods section.

We can check that the cycle graph on three vertices has positive Elser numbers for '2\leq k\geq 9' using::

    sage: G = graphs.CycleGraph(3); [G.elser_number(k) for k in range(10)]
    [-1, 0, 6, 30, 114, 390, 1266, 3990, 12354, 37830]


For each `U\in V(G)`, the `U`-nucleus complex `\Delta_U^G` is a simplicial complex whose vertices are the edges of `G`; its faces are complements of nuclei.
We can compute the U-nucleus complex using the nucleus_complex method::

    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: G.nucleus_complex([3,4])
    Chain complex with at most 4 nonzero terms over Rational Field

There is also a slightly faster implementation, for when you want to compute the nucleus complexes for all U for a given G::

    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: G.all_nucleus_complexes()
    {(): Chain complex with at most 5 nonzero terms over Rational Field,
     (1,): Chain complex with at most 5 nonzero terms over Rational Field,
     (1, 2): Chain complex with at most 5 nonzero terms over Rational Field,
     (1, 2, 3): Chain complex with at most 4 nonzero terms over Rational Field,
     (1, 2, 3, 4): Chain complex with at most 3 nonzero terms over Rational Field,
     (1, 2, 4): Chain complex with at most 4 nonzero terms over Rational Field,
     (1, 3): Chain complex with at most 4 nonzero terms over Rational Field,
     (1, 3, 4): Chain complex with at most 4 nonzero terms over Rational Field,
     (1, 4): Chain complex with at most 4 nonzero terms over Rational Field,
     (2,): Chain complex with at most 5 nonzero terms over Rational Field,
     (2, 3): Chain complex with at most 4 nonzero terms over Rational Field,
     (2, 3, 4): Chain complex with at most 4 nonzero terms over Rational Field,
     (2, 4): Chain complex with at most 4 nonzero terms over Rational Field,
     (3,): Chain complex with at most 4 nonzero terms over Rational Field,
     (3, 4): Chain complex with at most 4 nonzero terms over Rational Field,
     (4,): Chain complex with at most 4 nonzero terms over Rational Field}

From this, we can extract the nucleus complex of any U. For example, if `G` is the cycle graph on four vertices, applying all_nucleus_complexes will compute all the nucleus complexes for all subsets of the vertices. Taking the Uth element of that dictionary will return the Uth nucleus complex. We illustrate this below with `U = \{0,1\}`::

    sage: from sage.graphs.elser_numbers import nucleus_complex
    sage: G = graphs.CycleGraph(4)
    sage: G.all_nucleus_complexes()[(0,1)] == nucleus_complex(graphs.CycleGraph(4),[0,1])
    True

We can also compute the Betti numbers of these complexes. This method allows you the option of specifying a vertex subset `U` or not. If a U is given, then only the relevant Betti numbers will be returned. Otherwise, the Betti numbers of all nucleus complexes are given::

    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: G.nucleus_complex_bettis([1])
    {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: G.nucleus_complex_bettis()[(1,2)]
    {1: 0, 2: 0, 3: 1, 4: 0, 5: 0}

Authors
-------

- Galen Dorpalen-Barry, Cyrus Hettle, David C. Livingston, Jeremy L. Martin, George Nasr, Julianne Vega, and Hays Whitlatch (2021-06-21), initial implementation

Methods
---------
"""

# ****************************************************************************
#       Copyright (C) 2021 Galen Dorpalen-Barry <dorpa003@umn.edu>
#       Copyright (C) 2021 Cyrus Hettle <chettle3@gatech.edu>
#       Copyright (C) 2021 David C. Livingston <dliving5@uwyo.edu>
#       Copyright (C) 2021 Jeremy L. Martin <jlmartin@ku.edu>
#       Copyright (C) 2021 George Nasr <george.nasr@huskers.unl.edu>
#       Copyright (C) 2021 Julianne Vega <jvega30@kennesaw.edu>
#       Copyright (C) 2021 Hays Whitlatch <whitlatch@gonzaga.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.rest_index_of_methods import doc_index, gen_thematic_rest_table_index
from sage.misc.misc import powerset
from sage.graphs.all import graphs
from sage.rings.rational_field import QQ
from sage.homology.chain_complex import ChainComplex
from sage.matrix.constructor import matrix

@doc_index("Basic methods")
def is_vertex_cover(G,L):
    '''
    Determines if a subset of vertices is a vertex cover of the graph.
    That is, checks if every edge of the graph has at least one vertex in L.

    INPUT:

    - ``G`` -- a graph
    - ``L`` -- a subset of the vertices of G (as a list)

    OUTPUT:

    True if L is a vertex cover of G. False otherwise.

    EXAMPLES::

        sage: G = graphs.CycleGraph(4); G.is_vertex_cover([0,1])
        False
        sage: G = graphs.CycleGraph(4); G.is_vertex_cover([0,1,2])
        True
        sage: G = graphs.StarGraph(5); G.is_vertex_cover([0])
        True
        sage: G = graphs.StarGraph(5); G.is_vertex_cover([1])
        False
        sage: G = graphs.PetersenGraph(); G.is_vertex_cover([6,7,8,9,0,1,3,5])
        True

    TESTS::

        sage: G = graphs.PathGraph(3); G.is_vertex_cover([0,1,'hat'])
        True
        sage: G = graphs.PathGraph(1); G.is_vertex_cover([])
        True
        sage: G = graphs.PathGraph(2); G. is_vertex_cover([])
        False

    '''
    for e in G.edges():
        if e[0] not in L and e[1] not in L:
            return False
    return True

@doc_index("Basic methods")
def connected_subgraphs_of(G):
    '''
    Computes all connected subgraphs of a given graph G.

    INPUT:

    - ``G`` -- a graph

    OUTPUT:

    The collection of connected subgraphs of G.

    EXAMPLES::

        sage: G = graphs.CycleGraph(3); G.connected_subgraphs_of()
        [Graph on 0 vertices, Graph on 2 vertices, Graph on 2 vertices, Graph on 3 vertices, Graph on 2 vertices, Graph on 3 vertices, Graph on 3 vertices, Graph on 3 vertices]
        sage: G = graphs.PathGraph(3); G.connected_subgraphs_of()
        [Graph on 0 vertices, Graph on 2 vertices, Graph on 2 vertices, Graph on 3 vertices]

    Behavior with disconnected graphs::

        sage: G1 = graphs.CycleGraph(3); G2 = graphs.PathGraph(3)
        sage: G = G1.disjoint_union(G2); G.connected_subgraphs_of()
        [Graph on 0 vertices, Graph on 2 vertices, Graph on 2 vertices, Graph on 3 vertices, Graph on 2 vertices, Graph on 3 vertices, Graph on         3 vertices, Graph on 3 vertices, Graph on 2 vertices, Graph on 2 vertices, Graph on 3 vertices]

    TESTS::

        sage: G = graphs.CompleteGraph(4); len(G.connected_subgraphs_of()) == 61
        True

    '''
    from sage.graphs.graph import Graph
    LIST=[]
    for s in powerset(G.edges()):
        H = Graph(s)
        if H.is_connected():
            LIST.append(H)
    return LIST

@doc_index("Leftovers")
def nuclei_by_size(G,U=[]):
    '''
    A nucleus 'N' of a graph 'G' is a connected subgraph of 'G' whose vertex set 'V(N)' is a vertex cover of 'G.'
    This function computes all nuclei of 'G' and groups them by size.

    INPUT:

    - ``G`` -- a graph
    - ``U`` -- a subset of the vertices of G (if no U is specified, assumes U is the emptyset)

    OUTPUT:

    A dictionary whose (i)th element is a list of nuclei with i edges whose vertex set is contained in U.

    EXAMPLES::

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
        sage: G.nuclei_by_size([3,4])
        {2: [Graph on 3 vertices, Graph on 3 vertices],
         3: [Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices],
         4: [Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices,
          Graph on 4 vertices],
         5: [Graph on 4 vertices]}
         sage: G = graphs.CycleGraph(3); G.nuclei_by_size()
         {1: [Graph on 2 vertices, Graph on 2 vertices, Graph on 2 vertices], 2: [Graph on 3 vertices, Graph on 3 vertices, Graph on 3 vertices], 3: [Graph on 3 vertices]}
         sage: G = graphs.PathGraph(3); G.nuclei_by_size()[2][0].edges()
         [(0, 1, None), (1, 2, None)]

    .. SEEALSO::

        * :meth:`~sage.graphs.elser_number` --
          computes the kth Elser number of a graph

    TESTS::

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]); N = G.nuclei_by_size([3,4]); N[3][1].is_isomorphic(graphs.PathGraph(4))
        True
        sage: G = graphs.CompleteGraph(4); N = G.nuclei_by_size([0,2]); N[4][1] == Graph([[0,1],[0,3],[1,3],[0,2]])
        True
        sage: G = graphs.CompleteGraph(2); G.nuclei_by_size()[1][0] == graphs.CompleteGraph(2)
        True
    '''
    output = {}
    UU = set(U)
    for H in connected_subgraphs_of(G):
        S = H.vertices()
        if set(S).issuperset(UU) and is_vertex_cover(G,S):
            m = H.num_edges()
            if m in output:
                output[m].append(H)
            else:
                output[m] = [H]
    return output

@doc_index("Leftovers")
def elser_number(G,k):
    '''
    Computes the kth elser number of G.

    INPUT:

    - ``G`` -- a graph
    - ``k`` -- an integer

    OUTPUT:

    The kth Elser number of ``G``

    EXAMPLES::

        sage: G = graphs.CycleGraph(3); [G.elser_number(k) for k in range(10)]
        [-1, 0, 6, 30, 114, 390, 1266, 3990, 12354, 37830]

        sage: G = graphs.PathGraph(4); [G.elser_number(k) for k in range(10)]
        [0, 0, 2, 18, 110, 570, 2702, 12138, 52670, 223290]

        sage: G = graphs.PathGraph(5); [G.elser_number(k) for k in range(10)]
        [0, 0, 2, 24, 194, 1320, 8162, 47544, 266114, 1448520]

    .. SEEALSO::

        * :meth:`~sage.graphs.nuclei_by_size` --
          computes the nuclei of a graph and groups them by size

    TESTS::

        sage: G = graphs.CycleGraph(3); [G.elser_number(k) for k in range(10)] == [6*(3**(k-1) - 2**(k-1)) for k in range(10)]
        True

        sage: G = graphs.PathGraph(4); [G.elser_number(k) for k in range(10)] == [4**k - 2*(4-1)**k + (4-2)**k for k in range(10)]
        True

        sage: G = graphs.PathGraph(5); [G.elser_number(k) for k in range(10)] == [5**k - 2*(5-1)**k + (5-2)**k for k in range(10)]
        True

        sage: G = graphs.PathGraph(6); [G.elser_number(k) for k in range(10)] == [6**k - 2*(6-1)**k + (6-2)**k for k in range(10)]
        True
    '''
    nuclei = nuclei_by_size(G)
    return (-1)**(len(G.vertices())+1) * sum([(-1)**(i) * sum([len(N.vertices())**k for N in nuclei[i]]) for i in nuclei])

def coboundary_entry(Nyuck,k,i,j):
    '''
    Returns (i,j) entry of the kth coboundary map for the U-nucleus complex associated to Nyuck.

    INPUT

    - ``Nyuck`` -- a dictionary with keys: integers and entries lists of lists
    - ``k`` -- one of the keys of ``Nyuck``
    - ``i`` -- the index of an element in ``Nyuck[k]``
    - ``j`` -- the index of an element in ``Nyuck[k-1]``

    OUTPUT:

    The (i,j)th entry of the kth coboundary map of a U-nucleus complex.

    EXAMPLES::

        sage: from sage.graphs.elser_numbers import coboundary_entry
        sage: Nyuck = {2: [[0, 1]], 3: [[0, 1, 2], [0, 1, 3], [0, 1, 4]], 4: [[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4]], 5: [[0, 1, 2, 3, 4]]}
        sage: coboundary_entry(Nyuck,3,2,0)
        -1
        sage: coboundary_entry(Nyuck,4,2,0)
        0

    '''
    ColSet = Nyuck[k-1][j]
    RowSet = Nyuck[k][i]
    if set(ColSet).issubset(set(RowSet)):
        x = list(set(RowSet)-set(ColSet))[0]
        return -(-1)**(list(RowSet).index(x))
    else:
        return 0

@doc_index("Leftovers")
def nucleus_complex(G, U, ring=QQ):
    '''
    For a graph G and vertex subset U, computes the U-nucleus comples of G defined in [(DB)HLMNVW 2021].

    INPUT:

    - ``G`` -- a graph
    - ``U`` -- a subset of the vertices of G
    - ``ring`` -- coefficient ring (if no ring is specified, the rational numbers are used)

    OUTPUT:

    the U-nucleus complex of G

    EXAMPLES::

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]); G.nucleus_complex([3,4])
        Chain complex with at most 4 nonzero terms over Rational Field

        sage: G = graphs.StarGraph(5); NC = G.nucleus_complex([1,2]); bett = NC.betti(); [bett[deg] for deg in bett]
        [0, 0, 0, 0]

        sage: G = graphs.CycleGraph(5); NC = G.nucleus_complex([0,1]); bett = NC.betti(); [bett[deg] for deg in bett]
        [0, 1, 0]

    .. SEEALSO::

        * :meth:`~sage.graphs.nuclei_by_size` --
          computes the nuclei of a graph and groups them by size
        * :meth:`~sage.graphs.elser_number` --
          computes the Elser number of a graph

    '''
    EDGES = list(G.edges())
    nuclei = nuclei_by_size(G,U)
    Nyuck = { i : [sorted([EDGES.index(e) for e in list(N.edges())]) for N in nuclei[i]] for i in nuclei}
    Cooboundary = lambda k: matrix(ring,len(Nyuck[k]),len(Nyuck[k-1]),lambda i,j: coboundary_entry(Nyuck,k,i,j))
    HomSupport_Indices = list(Nyuck.keys())
    HomSupport_Indices.remove(min(HomSupport_Indices))
    return ChainComplex({k-1: Cooboundary(k) for k in HomSupport_Indices})

@doc_index("Leftovers")
def all_nucleus_complexes(G, ring=QQ):
    '''
    INPUT:

    - ``G`` -- a graph
    - ``ring`` -- a ring (if no ring is specified, the field of rational numbers is used)

    OUTPUT:

    all summands of the Elser complex, as a dictionary

    Note: This is a slight improvement over the function nucleus_complex, as we compute the complete list of nuclei by size at first, then refer to it to compute each summand, one U at a time. Its faster if we want to compute all the Elser summands.

    EXAMPLES::

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]); G.all_nucleus_complexes()[(3,4)]
        Chain complex with at most 4 nonzero terms over Rational Field
        sage: from sage.graphs.elser_numbers import nucleus_complex
        sage: G = graphs.CycleGraph(4); G.all_nucleus_complexes()[(0,1)] == nucleus_complex(graphs.CycleGraph(4),[0,1])
        True
        sage: G = graphs.PathGraph(5); G.all_nucleus_complexes()[(0,1)] == nucleus_complex(graphs.PathGraph(5),[0,1])
        True
        sage: G = graphs.PathGraph(3); G.all_nucleus_complexes()[(0,1)] == nucleus_complex(G,[0,1])
        True

    TESTS::

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]); NC = G.all_nucleus_complexes(); bet = NC[(2,3,4)].betti(); [bet[deg] for deg in bet]
        [0, 3, 0, 0]

    '''
    EDGES = list(G.edges())
    nuclei = nuclei_by_size(G)
    Summands = {}
    for U in powerset(G.vertices()):
        Nyuck = {}
        for i in nuclei:
            U_compatible = []
            for N in nuclei[i]:
                if set(U).issubset(N.vertices()):
                    U_compatible.append(N)
            Nyuck[i] = [sorted([EDGES.index(e) for e in list(N.edges())]) for N in U_compatible]
        Cooboundary = lambda k: matrix(ring,len(Nyuck[k]),len(Nyuck[k-1]),lambda i,j: coboundary_entry(Nyuck,k,i,j))
        HomSupport_Indices = list(Nyuck.keys())
        HomSupport_Indices.remove(min(HomSupport_Indices))
        Summands[tuple(U)] = ChainComplex({k-1: Cooboundary(k) for k in HomSupport_Indices})
    return Summands

@doc_index("Leftovers")
def nucleus_complex_bettis(G, U=None, ring=QQ):
    '''
    Computes the betti numbers of the nucleus complex of a graph.
    If a U is given, then it will compute the Betti numbers of that U-nucleus complex.
    If no U is given, then it will compute the Betti numbers of the nucleus complex of every subset of the vertices.

    INPUT:

    - ``G`` -- a graph
    - ``U`` -- a subset of the vertices of G (if no U is specified, then the function will iterate over all U in the power set of the vertices)
    - ``ring`` -- a coefficient ring (if no ring is specified, the rational numbers are used)

    OUTPUT:

    The betti numbers of U-nucleus comples of G.
    If no U is specified, a dictionary for which the Uth entry is the Betti table of the Uth Elser summand.

    EXAMPLES::

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]); G.nucleus_complex_bettis([1])
        {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

        sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]); G.nucleus_complex_bettis()[(1,2)]
        {1: 0, 2: 0, 3: 1, 4: 0, 5: 0}

    '''
    if U == None:
        summands = all_nucleus_complexes(G)
        return {V : summands[V].betti() for V in summands}
    else:
        return nucleus_complex(G,U).betti()
