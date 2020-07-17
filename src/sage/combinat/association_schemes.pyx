# -*- coding: utf-8 -*-

def is_from_association_scheme(list arr):
    r"""
    Return (n,r,i) if graph_from_association_scheme(n,r,i) has the given intersection array
    It returns False if this is not the case
    """
    from sage.arith.misc import is_prime_power

    if len(arr) != 6: return False
    n = arr[0]
    mu = arr[4]

    if n <= 0 or  mu <= 0 or (n-1)%mu != 0: return False
    r = (n-1) // mu

    if r == 1:#complete graph
        return False
    if arr != [n,(r-1)*mu,1,1,mu,n]:
        return False

    # Check we can build the appropriate association scheme
    schemeType = -1  # indicates the type of association scheme
    if is_prime_power(n) and (n % 2) * (((n-1)//r) % 2) == 0:
        # cyclotomic scheme
        schemeType = 1

    if schemeType != -1:
        return (n, r, schemeType)
    else:
        return False

def _find_quasigroup(scheme):
    r"""
    scheme is the relationship matrix of an association scheme
    try to find a commutative binary operation on I = {1, ..., r} that make I a quasigroup
    with p_{op(i,l), op(j,l)}^{op(k,l)} = p_{i,j}^k for all i,j,k,l
    """
    # not sure how to do it
    # hopefully not needed


def graph_from_association_scheme(const int n, const int r, const int schemeType):
    r"""
    build an antipodal r-cover of K_{n+1} using association schemes with the given type
    """
    from sage.libs.gap.libgap import libgap
    from sage.matrix.constructor import matrix
    
    if schemeType == 1:
        # cyclotomic scheme
        libgap.LoadPackage("AssociationSchemes")
        s = libgap.CyclotomicScheme(n,r)
        S = matrix(s.RelationMatrix())
        def inv(c,a):
            b = (c-a) % r
            if b == 0:
                b = r
            return b
        op = None

    return association_scheme_graph(S, op, inv)

def association_scheme_graph(scheme, op=None, inv=None):
    r"""
    scheme is the relation matrix representing the association scheme
    op is a binary function on I = {1,...,r} where r is the number of classes
    inv is the partial inverse of op, i.e.:
    op(a,b) = op(b,a) = c => inv(c,a) = b and inv(c,b) = a

    if inv is None, then inv is deduced by op;
    if inv is not None, then op is useless
    """
    from sage.graphs.graph import Graph

    r = max(map(max, scheme.rows()))
    I = list(range(1,r+1))
    n = scheme.ncols()
    inf = n

    if inv is None:
        _inv = {}
        for i in I:
            for j in I:
                c = op(i, j)
                _inv[(c,i)] = j
                _inv[(c,j)] = i

        inv = lambda c, a: _inv[(c,a)]

    edges = []
    for x in range(n):
        for i in I:
            edges.append(( (inf,i),(x,i) ))

    for x in range(n):
        for y in range(x+1, n):
            ij = scheme[x,y]
            for i in I:
                j = inv(ij,i)
                edges.append(( (x,i),(y,j) ))
                edges.append(( (y,i),(x,j) ))
    
    G = Graph(edges,format="list_of_edges")
    return G


def distance_regular_graph(list array, check=True):

    t = is_from_association_scheme(array)
    if t is not False:
        G = graph_from_association_scheme(*t)

    if check:
        t = G.is_distance_regular(True)
        if t is False or (t[0][:-1] +  t[1][1:]) != array:
            raise RuntimeError("Sage built the wrong graph")

    return G
