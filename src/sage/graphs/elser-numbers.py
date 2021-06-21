def is_vertex_cover(H,L):
    '''
    INPUT:  a graph H and subset of the vertices L
    OUTPUT: A boolean. True if L is a vertex cover of H. False otherwise.

    This is used for computing the Nuclei of a graph
    '''
    value=True
    for e in H.edges():
        if e[0] not in L and e[1] not in L:
            value=False
            break
    return value

def connected_subgraphs_of(G):
    '''
    INPUT:  a graph G
    OUTPUT: The collection of connected subgraphs of G

    This is used for computing the Nuclei of a graph
    '''
    LIST=[]
    for s in powerset(G.edges()):
        H=Graph(s)
        if H.is_connected():
            LIST.append(H)
    return LIST

def NucleiBySize(G,U=[]):
    '''
    INPUT:  a graph G and subset of the vertices of G, which we call U (if no U is specified, assumes U is the emptyset)
    OUTPUT: dictionary whose ith element is a list of nuclei with i edges whose vertex set is contained in U.

    EXAMPLE
    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: NucleiBySize(G,[3,4])
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
    '''
    output={}
    UU=Set(U)
    for H in connected_subgraphs_of(G):
        S = H.vertices()
        if Set(S).issuperset(UU) and is_vertex_cover(G,S):
            m = H.num_edges()
            if m in output:
                output[m].append(H)
            else:
                output[m] = [H]
    return output


def elser_number(G,k,ring=QQ):
    '''
    INPUT:  A graph, an integer, and a ring (default is the rationals)
    OUTPUT: the kth Elser number of G

    EXAMPLE
    sage: [elser_number(graphs.CycleGraph(3),k) for k in range(10)]
    [-1, 0, 6, 30, 114, 390, 1266, 3990, 12354, 37830]

    sage: [elser_number(graphs.PathGraph(4),k) for k in range(10)]
    [0, 0, 2, 18, 110, 570, 2702, 12138, 52670, 223290]

    sage: [elser_number(graphs.PathGraph(5),k) for k in range(10)]
    [0, 0, 2, 24, 194, 1320, 8162, 47544, 266114, 1448520]

    sage: [elser_number(graphs.PathGraph(6),k) for k in range(10)] == [6**k - 2*(6-1)**k + (6-2)**k for k in range(10)]
    True

    '''
    nuclei_by_size = NucleiBySize(G,U=[])
    return (-1)**(len(G.vertices())+1) * sum([(-1)**(i) * sum([len(N.vertices())**k for N in nuclei_by_size[i]]) for i in nuclei_by_size])



def _coboundary_entry(Nyuck,k,i,j):
    '''
    Returns (i,j) entry of the kth coboundary map for the Elser-complex associated to Nyuck.
    This is used in the Elser Summand computation.
    '''
    ColSet = Nyuck[k-1][j]
    RowSet = Nyuck[k][i]
    if set(ColSet).issubset(set(RowSet)):
        x = list(set(RowSet)-set(ColSet))[0]
        return -(-1)^(list(RowSet).index(x))
    else:
        return 0

def ElserSummand(G,U, ring=QQ):
    '''
    INPUT:  a graph G, a vertex set U, and a coefficient ring (if no ring is specified, the rational numbers are used)
    OUTPUT: the U-summand of the Elser cochain complex of G

    EXAMPLE
    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: ElserSummand(G,[3,4])
    Chain complex with at most 4 nonzero terms over Rational Field
    '''
    EDGES = list(G.edges())
    nuclei_by_size = NucleiBySize(G,U)
    Nyuck = { i : [sorted([EDGES.index(e) for e in list(N.edges())]) for N in nuclei_by_size[i]] for i in nuclei_by_size}
    Cooboundary = lambda k: Matrix(ring,len(Nyuck[k]),len(Nyuck[k-1]),lambda i,j: _coboundary_entry(Nyuck,k,i,j))
    HomSupport_Indices = list(Nyuck.keys())
    HomSupport_Indices.remove(min(HomSupport_Indices))
    return ChainComplex({k-1: Cooboundary(k) for k in HomSupport_Indices})


def all_ElserSummands(G,ring=QQ):
    '''
    INPUT:  a graph G
    OUTPUT: all summands of the Elser complex, as a dictionary

    Note: This is a slight improvement over the "ElserSummand function," as we compute the complete list of nuclei by size at first, then refer to it to compute each summand, one U at a time. Its faster if we want to compute all the Elser summands.

    EXAMPLES
    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: all_ElserSummands(G)
    {(): Chain complex with at most 5 nonzero terms over Rational Field, (1,): Chain complex with at most 5 nonzero terms over Rational Field, (2,): Chain complex with at most 5 nonzero terms over Rational Field, (1, 2): Chain complex with at most 5 nonzero terms over Rational Field, (3,): Chain complex with at most 4 nonzero terms over Rational Field, (1, 3): Chain complex with at most 4 nonzero terms over Rational Field, (2, 3): Chain complex with at most 4 nonzero terms over Rational Field, (1, 2, 3): Chain complex with at most 4 nonzero terms over Rational Field, (4,): Chain complex with at most 4 nonzero terms over Rational Field, (1, 4): Chain complex with at most 4 nonzero terms over Rational Field, (2, 4): Chain complex with at most 4 nonzero terms over Rational Field, (1, 2, 4): Chain complex with at most 4 nonzero terms over Rational Field, (3, 4): Chain complex with at most 4 nonzero terms over Rational Field, (1, 3, 4): Chain complex with at most 4 nonzero terms over Rational Field, (2, 3, 4): Chain complex with at most 4 nonzero terms over Rational Field, (1, 2, 3, 4): Chain complex with at most 3 nonzero terms over Rational Field}

    sage: all_ElserSummands(graphs.CycleGraph(4))[(0,1)] == ElserSummand(graphs.CycleGraph(4),[0,1])
    True

    sage: all_ElserSummands(graphs.PathGraph(5))[(0,1)] == ElserSummand(graphs.PathGraph(5),[0,1])
    True

    '''
    EDGES = list(G.edges())
    nuclei_by_size = NucleiBySize(G)
    Summands = {}
    for U in powerset(G.vertices()):
        Nyuck = {}
        for i in nuclei_by_size:
            U_compatible = []
            for N in nuclei_by_size[i]:
                if set(U).issubset(N.vertices()):
                    U_compatible.append(N)
            Nyuck[i] = [sorted([EDGES.index(e) for e in list(N.edges())]) for N in U_compatible]
        Cooboundary = lambda k: Matrix(ring,len(Nyuck[k]),len(Nyuck[k-1]),lambda i,j: _coboundary_entry(Nyuck,k,i,j))
        HomSupport_Indices = list(Nyuck.keys())
        HomSupport_Indices.remove(min(HomSupport_Indices))
        Summands[tuple(U)] = ChainComplex({k-1: Cooboundary(k) for k in HomSupport_Indices})
    return Summands


def ElserComplex_Bettis(G,U=None,ring=QQ):
    '''
    INPUT:  a graph G, a subset of vertices U (optional), and a coefficient ring (if no ring is specified, the rational numbers are used)
    OUTPUT: A dictionary.
    If no U is specified, then a dictionary for which the Uth entry is the Betti table of the Uth Elser summand.
    If U is specified, then the betti table of U is given.

    EXAMPLES

    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: ElserSummand_Bettis(G,[1])
    {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    sage: G = Graph([[1,2],[1,3],[1,4],[2,3],[2,4]]);
    sage: ElserSummand_Bettis(G)[(1,2)]
    {1: 0, 2: 0, 3: 1, 4: 0, 5: 0}

    '''
    if U == None:
        summands = all_ElserSummands(G)
        return {V : summands[V].betti() for V in summands}
    else:
        return ElserSummand(G,U).betti()
