"""
Matching Polynomial

This module contains the following methods:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`matching_polynomial` | Computes the matching polynomial of a given graph
    :meth:`complete_poly` | Compute the matching polynomial of the complete graph on `n` vertices.

    :meth:`matching_generating_poly` | Computes the matching generating polynomial

AUTHORS:

- Robert Miller, Tom Boothby - original implementation

REFERENCE:

.. [Godsil93] Chris Godsil (1993) Algebraic Combinatorics.


Methods
-------
"""

#*****************************************************************************
#                       Copyright (C) 2010 Robert Miller
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.integer_ring import ZZ
from sage.rings.integer cimport Integer
from sage.misc.all import prod
include 'sage/ext/interrupt.pxi'
include 'sage/ext/cdefs.pxi'
include 'sage/ext/stdsage.pxi'

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_poly cimport *

x = polygen(ZZ, 'x')


def matching_polynomial(G, complement=True, name=None, algorithm='ButeraPernici'):
    """
    Computes the matching polynomial of the graph `G`.

    If `p(G, k)` denotes the number of `k`-matchings (matchings with `k` edges)
    in `G`, then the matching polynomial is defined as [Godsil93]_:

    .. MATH::

        \mu(x)=\sum_{k \geq 0} (-1)^k p(G,k) x^{n-2k}

    INPUT:

    - ``complement`` - (default: ``True``) whether to use Godsil's duality
      theorem to compute the matching polynomial from that of the graphs
      complement (see ALGORITHM).

    - ``name`` - optional string for the variable name in the polynomial

    .. NOTE::

        The ``complement`` option uses matching polynomials of complete graphs,
        which are cached. So if you are crazy enough to try computing the
        matching polynomial on a graph with millions of vertices, you might not
        want to use this option, since it will end up caching millions of
        polynomials of degree in the millions.

    ALGORITHM:

    The `Godsil` algorithm used is a recursive one, based on the following observation
    [Godsil93]_:

    - If `e` is an edge of `G`, `G'` is the result of deleting the edge `e`, and
      `G''` is the result of deleting each vertex in `e`, then the matching
      polynomial of `G` is equal to that of `G'` minus that of `G''`.

      (the algorithm actually computes the *signless* matching polynomial, for
      which the recursion is the same when one replaces the substraction by an
      addition. It is then converted into the matching polynomial and returned)

    Depending on the value of ``complement``, Godsil's duality theorem
    [Godsil93]_ can also be used to compute `\mu(x)` :

    .. MATH::

        \mu(\overline{G}, x) = \sum_{k \geq 0} p(G,k) \mu( K_{n-2k}, x)


    Where `\overline{G}` is the complement of `G`, and `K_n` the complete graph
    on `n` vertices.

    The `ButeraPernici` algorithm is faster on large graphs; for a
    description of the algorithm see :meth:`matching_generating_poly`


    EXAMPLES::

        sage: g = graphs.PetersenGraph()
        sage: g.matching_polynomial()
        x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
        sage: g.matching_polynomial(complement=False)
        x^10 - 15*x^8 + 75*x^6 - 145*x^4 + 90*x^2 - 6
        sage: g.matching_polynomial(name='tom')
        tom^10 - 15*tom^8 + 75*tom^6 - 145*tom^4 + 90*tom^2 - 6
        sage: g = Graph()
        sage: L = [graphs.RandomGNP(8, .3) for i in range(1, 6)]
        sage: prod([h.matching_polynomial() for h in L]) == sum(L, g).matching_polynomial()  # long time (up to 10s on sage.math, 2011)
        True

    ::

        sage: for i in range(1, 12):  # long time (10s on sage.math, 2011)
        ....:     for t in graphs.trees(i):
        ....:         if t.matching_polynomial() != t.characteristic_polynomial():
        ....:             raise RuntimeError('bug for a tree A of size {0}'.format(i))
        ....:         c = t.complement()
        ....:         if c.matching_polynomial(complement=False) != c.matching_polynomial():
        ....:             raise RuntimeError('bug for a tree B of size {0}'.format(i))

    ::

        sage: from sage.graphs.matchpoly import matching_polynomial
        sage: matching_polynomial(graphs.CompleteGraph(0))
        1
        sage: matching_polynomial(graphs.CompleteGraph(1))
        x
        sage: matching_polynomial(graphs.CompleteGraph(2))
        x^2 - 1
        sage: matching_polynomial(graphs.CompleteGraph(3))
        x^3 - 3*x
        sage: matching_polynomial(graphs.CompleteGraph(4))
        x^4 - 6*x^2 + 3
        sage: matching_polynomial(graphs.CompleteGraph(5))
        x^5 - 10*x^3 + 15*x
        sage: matching_polynomial(graphs.CompleteGraph(6))
        x^6 - 15*x^4 + 45*x^2 - 15
        sage: matching_polynomial(graphs.CompleteGraph(7))
        x^7 - 21*x^5 + 105*x^3 - 105*x
        sage: matching_polynomial(graphs.CompleteGraph(8))
        x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
        sage: matching_polynomial(graphs.CompleteGraph(9))
        x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
        sage: matching_polynomial(graphs.CompleteGraph(10))
        x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
        sage: matching_polynomial(graphs.CompleteGraph(11))
        x^11 - 55*x^9 + 990*x^7 - 6930*x^5 + 17325*x^3 - 10395*x
        sage: matching_polynomial(graphs.CompleteGraph(12))
        x^12 - 66*x^10 + 1485*x^8 - 13860*x^6 + 51975*x^4 - 62370*x^2 + 10395
        sage: matching_polynomial(graphs.CompleteGraph(13))
        x^13 - 78*x^11 + 2145*x^9 - 25740*x^7 + 135135*x^5 - 270270*x^3 + 135135*x

    ::

        sage: G = Graph({0:[1,2], 1:[2]})
        sage: matching_polynomial(G)
        x^3 - 3*x
        sage: G = Graph({0:[1,2]})
        sage: matching_polynomial(G)
        x^3 - 2*x
        sage: G = Graph({0:[1], 2:[]})
        sage: matching_polynomial(G)
        x^3 - x
        sage: G = Graph({0:[], 1:[], 2:[]})
        sage: matching_polynomial(G)
        x^3

    ::

        sage: matching_polynomial(graphs.CompleteGraph(0), complement=False)
        1
        sage: matching_polynomial(graphs.CompleteGraph(1), complement=False)
        x
        sage: matching_polynomial(graphs.CompleteGraph(2), complement=False)
        x^2 - 1
        sage: matching_polynomial(graphs.CompleteGraph(3), complement=False)
        x^3 - 3*x
        sage: matching_polynomial(graphs.CompleteGraph(4), complement=False)
        x^4 - 6*x^2 + 3
        sage: matching_polynomial(graphs.CompleteGraph(5), complement=False)
        x^5 - 10*x^3 + 15*x
        sage: matching_polynomial(graphs.CompleteGraph(6), complement=False)
        x^6 - 15*x^4 + 45*x^2 - 15
        sage: matching_polynomial(graphs.CompleteGraph(7), complement=False)
        x^7 - 21*x^5 + 105*x^3 - 105*x
        sage: matching_polynomial(graphs.CompleteGraph(8), complement=False)
        x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105
        sage: matching_polynomial(graphs.CompleteGraph(9), complement=False)
        x^9 - 36*x^7 + 378*x^5 - 1260*x^3 + 945*x
        sage: matching_polynomial(graphs.CompleteGraph(10), complement=False)
        x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
        sage: matching_polynomial(graphs.CompleteGraph(11), complement=False)
        x^11 - 55*x^9 + 990*x^7 - 6930*x^5 + 17325*x^3 - 10395*x
        sage: matching_polynomial(graphs.CompleteGraph(12), complement=False)
        x^12 - 66*x^10 + 1485*x^8 - 13860*x^6 + 51975*x^4 - 62370*x^2 + 10395
        sage: matching_polynomial(graphs.CompleteGraph(13), complement=False)
        x^13 - 78*x^11 + 2145*x^9 - 25740*x^7 + 135135*x^5 - 270270*x^3 + 135135*x

    TESTS:

    Non-integer labels should work, (:trac:`15545`)::

        sage: G = Graph(10);
        sage: G.add_vertex((0,1))
        sage: G.add_vertex('X')
        sage: G.matching_polynomial()
        x^12
    """

    cdef int nverts, nedges, i, j, cur
    cdef int *edges1, *edges2, *edges_mem, **edges
    cdef fmpz_poly_t pol
    cdef Integer c_ZZ

    if G.has_multiple_edges():
        raise NotImplementedError

    nverts = G.num_verts()

    # Using Godsil's duality theorem when the graph is dense

    if complement and G.density() > 0.5:  # this cutoff could probably be tuned
        f_comp = matching_polynomial(G.complement()).list()
        f = x.parent().zero()
        for i from 0 <= i <= nverts / 2:  # implicit floor
            j = nverts - 2 * i
            f += complete_poly(j) * f_comp[j] * (-1)**i
        return f

    nedges = G.num_edges()

    # Relabelling the vertices of the graph as [0...n-1] so that they are sorted
    # in increasing order of degree

    L = []
    for v, d in G.degree_iterator(labels=True):
        L.append((d, v))
    L.sort()
    d = {}
    for i from 0 <= i < nverts:
        d[L[i][1]] = i
    G = G.relabel(d, inplace=False)
    G.allow_loops(False)

    if algorithm == 'Godsil':
        # Initialization of pol, edges* variables.

        # The edges_mem table is of size (2 * nedges * nedges), and is to be read as
        # nedges blocks of size (2 * nedges). These blocks of size (2 * nedges) are
        # themselves to be read as two blocks of length nedges, addressed as edges1
        # and edges2

        # Only the first block of size (2*nedges) is here filled. The function
        # delete_and_add will need the rest of the memory.

        fmpz_poly_init(pol)  # sets to zero
        edges_mem = <int *> sage_malloc(2 * nedges * nedges * sizeof(int))
        edges = <int **> sage_malloc(2 * nedges * sizeof(int *))
        if edges_mem is NULL or edges is NULL:
            if edges_mem is not NULL:
                sage_free(edges_mem)
            if edges is not NULL:
                sage_free(edges)
            raise MemoryError("Error allocating memory for matchpoly.")

        for i from 0 <= i < 2 * nedges:
            edges[i] = edges_mem + i * nedges

        edges1 = edges[0]
        edges2 = edges[1]

        cur = 0
        for i, j in sorted(map(sorted, G.edges(labels=False))):
            edges1[cur] = i
            edges2[cur] = j
            cur += 1

        # Computing the signless matching polynomial

        sig_on()
        delete_and_add(edges, nverts, nedges, nverts, 0, pol)
        sig_off()
        # Building the actual matching polynomial

        coeffs_ZZ = []
        for i from 0 <= i <= nverts:
            c_ZZ = Integer(0)
            fmpz_poly_get_coeff_mpz(c_ZZ.value, pol, i)
            coeffs_ZZ.append(c_ZZ * (-1)**((nverts - i) / 2))

        f = x.parent()(coeffs_ZZ)
        sage_free(edges_mem)
        sage_free(edges)
        fmpz_poly_clear(pol)
        if name is not None:
            return f.change_variable_name(name)
        return f
    elif algorithm == 'ButeraPernici':
        p = matching_generating_poly(G)
        a = p.coefficients()
        n = G.order()
        b = [0]*(n + 1)
        for i in range(len(a)):
            b[n - 2*i] = a[i]*(-1)**i
        if name is None:
            K = x.parent()
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            K = PolynomialRing(ZZ, name)
        p = K(b)
        return p
    else:
        raise ValueError('algorithm must be one of "Godsil" or "ButeraPernici".')




# The following is a cache of complete graph matching polynomials.

cdef list complete_matching_polys = [x.parent().one(), x]


def complete_poly(n):
    """
    Compute the matching polynomial of the complete graph on `n` vertices.

    INPUT:

    - ``n`` -- order of the complete graph

    .. TODO::

        This code could probably be made more efficient by using FLINT
        polynomials and being written in Cython, using an array of
        fmpz_poly_t pointers or something...  Right now just about the
        whole complement optimization is written in Python, and could
        be easily sped up.

    EXAMPLES::

        sage: from sage.graphs.matchpoly import complete_poly
        sage: f = complete_poly(10)
        sage: f
        x^10 - 45*x^8 + 630*x^6 - 3150*x^4 + 4725*x^2 - 945
        sage: f = complete_poly(20)
        sage: f[8]
        1309458150
        sage: f = complete_poly(1000)
        sage: len(str(f))
        406824

    TESTS:

    Checking the numerical results up to 20::

        sage: from sage.functions.orthogonal_polys import hermite
        sage: p = lambda n: 2^(-n/2)*hermite(n, x/sqrt(2))
        sage: all(p(i) == complete_poly(i) for i in range(2, 20))
        True
    """
    # global complete_matching_polys # if we do eventually make it a C array...
    cdef int lcmp = len(complete_matching_polys)

    if n < lcmp:
        return complete_matching_polys[n]
    lcmp -= 2
    a = complete_matching_polys[lcmp]
    b = complete_matching_polys[lcmp + 1]
    while lcmp + 2 <= n:
        a, b, lcmp = b, x * b - (lcmp + 1) * a, lcmp + 1
        complete_matching_polys.append(b)
    return b

cdef void delete_and_add(int **edges, int nverts, int nedges, int totverts, int depth, fmpz_poly_t pol):
    """
    Add matching polynomial to pol via recursion.

    NOTE : at the end of this function, pol represents the *SIGNLESS*
    matching polynomial.
    """
    cdef int i, j, k, edge1, edge2, new_edge1, new_edge2, new_nedges
    cdef int *edges1, *edges2, *new_edges1, *new_edges2
    cdef fmpz * coeff

    if nverts == 3:
        coeff = fmpz_poly_get_coeff_ptr(pol, 3)
        if coeff is NULL:
            fmpz_poly_set_coeff_ui(pol, 3, 1)
        else:
            fmpz_add_ui(coeff, coeff, 1)
        coeff = fmpz_poly_get_coeff_ptr(pol, 1)
        if coeff is NULL:
            fmpz_poly_set_coeff_ui(pol, 1, nedges)
        else:
            fmpz_add_ui(coeff, coeff, nedges)
        return

    if nedges == 0:
        coeff = fmpz_poly_get_coeff_ptr(pol, nverts)
        if coeff is NULL:
            fmpz_poly_set_coeff_ui(pol, nverts, 1)
        else:
            fmpz_add_ui(coeff, coeff, 1)
        return

    edges1 = edges[2 * depth]
    edges2 = edges[2 * depth + 1]
    new_edges1 = edges[2 * depth + 2]
    new_edges2 = edges[2 * depth + 3]

    nedges -= 1

    # The last edge is (edge1, edge2)
    edge1 = edges1[nedges]
    edge2 = edges2[nedges]
    new_nedges = 0

    # The new edges are all the edges that are not incident with (edge1, edge2)

    for i from 0 <= i < nedges:
        if edge1 == edges1[i]:
            break  # since the rest of the edges are incident to edge1
                   # (the edges
                   # are sorted by increasing order of their first component)

        if edge1 != edges2[i] and edge2 != edges2[i]:
            # since edge2 > edge1 > edges1[i], only need to check edges2[i]
            new_edges1[new_nedges] = edges1[i]
            new_edges2[new_nedges] = edges2[i]
            new_nedges += 1

    delete_and_add(edges, nverts - 2, new_nedges, totverts, depth + 1, pol)
    delete_and_add(edges, nverts, nedges, totverts, depth, pol)

class Hobj(object):
    """
    class used to compute products of nilpotent quantities associated
    to hard objects, like non-adjacent edges or non-adjacent vertices,
    used respectively to compute the matching and the independence
    polynomial.
    """
    def __init__(self):
        self.links = []
        self.dt = {}
        self.freedt = list(range(1000, -1, -1))

    def iadd_object(hb, p, val, obj, free):
        """
        Multiply ``p`` by ``(1 + t*val*eta_i*eta_j)``, for `obj=(i,j)`.

        INPUT:
            - ``p`` - polynomial

            - ``val`` - value

            - ``obj`` - a tuple of object labels

            - ``free`` - list of free indices

        ..NOTE::

            ``p`` is changed only if ``free`` is empty

            ``free`` is the list of indices of ``eta`` elements which
            are integrated (that is, put to ``1`` after performing the product).

            For example let `p = 1 + eta_0*eta_1`; multiply it by
            `1 + t*eta_1*eta_2`, knowing that after performing the product
            `eta_1` can be set to `1`, then
            p*(1 + t*eta_1*eta_2) = 1 + t*eta_0 + t*eta_2
            `val = `; obj = (1, 2); free = [1];

        EXAMPLES::

            sage: from sage.graphs.matchpoly import Hobj
            sage: p = {0: 1}
            sage: hb = Hobj()
            sage: x = polygen(ZZ, 'x')
            sage: p = hb.iadd_object(p, x, (0, 1), []); p
            {0: 1, 3: x}
            sage: p = hb.iadd_object(p, x, (1, 2), [1]); p
            {0: 1, 1: x, 4: x}

        """
        links = hb.links
        dt = hb.dt
        freedt = hb.freedt
        exp2 = 0
        # if an element is already present, get its pointer;
        # else assign it a new pointer from freedt.
        # Encode in binary the object
        for i in obj:
            if i in dt:
                j = dt[i]
            else:
                j = freedt.pop()
                dt[i] = j
            exp2 += 1 << j
        free = [dt[i] for i in free]
        # check that the object is not considered twice; save it in links
        t = tuple(sorted(obj))
        if t in links:
            raise ValueError('%s in %s' %(t, links))
        links.append(t)
        #valx = x*val
        if free:
            p1 = p
            p = {}

            # encode free in binary
            mask_free = 0
            for i in free:
                mask_free += 1 << i
            get = p.get
            for exp1, v1 in p1.iteritems():
                # multiply by 1 + t*val*prod_i eta_i
                for ii in range(2):
                    if not ii:
                        exp = exp1
                        v = v1
                    else:
                        if exp1 & exp2:
                            continue
                        exp = exp1 | exp2
                        v = v1*val
                    # eliminate the free elements from exp
                    if exp & mask_free:
                        for i in free:
                            if exp & (1 << i):
                                exp = exp.__xor__(1 << i)
                    p[exp] = get(exp, 0) + v
            # return the pointers not any more used to freedt
            for exp in free:
                freedt.append(exp)
            return p

        else:
            a = []
            for exp1, v1 in p.iteritems():
                if exp1 & exp2:
                    continue
                exp = exp1 | exp2
                v = v1*val
                try:
                    p[exp] = p[exp] + v
                except KeyError:
                    a.append((exp, v))
            # since `exp1 & exp2 = 0`, for `exp = exp1 | exp2`
            # there is a unique `exp1`; therefore the entries `(exp, v)`
            # in a have all different `exp`
            for exp, v in a:
                p[exp] = v

        return p

def obj_free(objects):
    """
    Return a list of tuples ``(obj, free)``

    INPUT:

        - ``objects`` - is a list of tuples (e.g 2-tuples in the case of dimers)

    .. NOTE::

        `free` is the list of elements in the tuple `obj`  which do not occur
        in the following objects.

        The objects are required to be in `range(n)`

    EXAMPLES::

        sage: from sage.graphs.matchpoly import obj_free
        sage: obj_free([(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)])
        [((0, 1), []), ((1, 2), [1]), ((2, 3), [2]), ((3, 4), [3]), ((4, 0), [0, 4])]
    """
    # nvars number of elements in objects
    s = set()
    for t in objects:
        for i in t:
            s.add(i)
    nvars = len(s)
    b = list(sorted(s))
    b.sort()
    if b != list(range(nvars)):
         raise ValueError('elements should be labelled in 0,...,nvars-1')

    v = []
    done = set()
    n = len(objects)
    for i in range(n):
        free = []
        obj = objects[i]
        for j in range(nvars):
            if j in done:
                continue
            hit = 0
            for k in range(i + 1, n):
                if j in objects[k]:
                    hit = 1
                    break
            if not hit:
                free.append(j)
                done.add(j)
        v.append((obj, free))
    return v

def gen_count_hobj(objects, values=None):
    """
    Counting polynomial for hard object from a list of edges

    INPUT:

        - ``objects`` - list of tuples of element indices

    .. NOTE::

        This function is used in `matching_generating_poly`
        and `independence_poly`.

    EXAMPLES::

        sage: from sage.graphs.matchpoly import gen_count_hobj
        sage: gen_count_hobj([(1, 0), (2, 1), (3, 2), (4, 0), (4, 3)])
        5*x^2 + 5*x + 1
    """
    # relabel objects so that its elements are in 0,...,n-1
    a = []
    for obj in objects:
        for i in obj:
            if i not in a:
                a.append(i)
    dt = dict(zip(a, range(len(a))))
    links = []
    for obj in objects:
        obj1 = tuple(dt[i] for i in obj)
        links.append(obj1)
    objects = links
    # add to each object the list of free indices
    obfr = obj_free(objects)
    hb = Hobj()
    if values:
        K = values.values()[0].parent()
    else:
        K = ZZ
    x = polygen(K, 'x')
    p = {0: K.one()}
    if not values:
        for obj, free in obfr:
            p = hb.iadd_object(p, x, obj, free)
    else:
        di = {}
        for k, v in dt.items():
            di[v] = k
        for obj, free in obfr:
            p = hb.iadd_object(p, x*values[(di[obj[0]], di[obj[1]])], obj, free)

    assert len(p) == 1
    return p[0]

def matching_generating_poly(g, links=None, labels=None):
    r"""
    Return the matching generating polynomial `M`.

    INPUT:

    - ``g`` - graph
    - ``links`` - list of the edges of the graph
    - ``labels`` - flag for computing `M` for the weighted graph.

    The matching generating polynomial of a simple graph is given by

    .. MATH::

        M(x)=\sum_{k \geq 0} p(G,k) x^{2k}

    where `p(G, k)` denotes the number of `k`-matchings in `G`,
    that is the number of configurations of non-intersecting edges (dimers).

    If the graph has edge weights, `p(G, k)` denotes the sum of the product
    of the weights of the edges of the `k`-matchings in `G`.

    ALGORITHM:

        The algorithm is described in [ButPer].
        Consider the algebra `R = K[\eta_1, \eta_2, \ldots, \eta_n]`
        where the `\eta_i` are commuting, nilpotent of order `2`
        (i.e. `\eta_i^2 = 0`), and where `n` is the order of the graph `G`.
        We will mostly consider the ring `R[x]` of polynomials over `R`.

        Introduce an "integration" operation `\langle p \rangle` over `R` and
        `R[x]` consisting in the sum of the coefficients of the non-vanishing
        monomials in `\eta_i`

        A useful simplification property of `\langle . \rangle`
        is the following: let
        `p_1(\eta_1, \ldots, \eta_n)` and `p_2(\eta_j, \ldots, \eta_n)` be
        polynomials in `R[x]` where `j \ge 1`.  Then one has

        .. MATH::

            \langle p_1(\eta_1, \ldots, \eta_n) p_2 \rangle =
            \langle p_1(1, \ldots, 1, \eta_j, \ldots, \eta_n) p_2 \rangle

        where `\eta_1,..,\eta_{j-1}` are replaced by `1` in `p_1`. Informally,
        we can "integrate" these variables *before* performing the product. More
        generally, if a monomial `\eta_i` is missing in one of the terms of a
        product of two terms, then it can be integrated in the other term.

        Associate to a dimer with endpoints `i, j` the quantity
        ``O_{i,j} = 1 + x \eta_i \eta_j``; `1` represents the case in which
        the dimer is absent; the other term represent the case in which the
        dimer is present.
        The matching generating polynomial is given by

        .. MATH::

            M(x) = \langle \prod_{(i,j) \in E(G)} (1 + x \eta_i \eta_j) \rangle

        where `E(G)` is the set of edges of G.
        Using the semplification property mentioned above, one can perform
        `\prod_{(i,j) \in E_k}`, where `E_k` is included in `E(G)`,
        replacing `\eta_i` with `1` as soon as the node `i` is absent from
        `E(G) - E_k`.

        The result of performing these product is independent from the
        ordering, but the efficiency of the algorithm depends exponentially
        from the number of `active nodes`, that is the nodes present in
        `G_k`, the graph defined by `E_k`, and which have in `G`
        degree greater than in `G_k` (if a node has the same degree in
        `G_k` as in `G` it means that one can set `eta_k=1`.

        A simple greedy algorithm tries to find an efficient ordering of
        edges to compute the matching polynomial.

        There is no guarantee that an efficient ordering of the edges is found.
        Alternatively can provide explicitly the list of edges `links`.

    EXAMPLES::

        sage: from sage.graphs.matchpoly import matching_generating_poly
        sage: g = Graph({0:[1,4], 1:[0,2], 2:[1,3], 3:[2,4], 4:[0,3]})
        sage: matching_generating_poly(g)
        5*x^2 + 5*x + 1
        sage: g = Graph()
        sage: g.add_edges([(0,2,1),(0,3,2),(1,2,3),(1,3,4)])
        sage: matching_generating_poly(g, labels=True)
        10*x^2 + 10*x + 1
    """
    from sage.graphs.graph_decompositions.active_nodes import ordered_links, _d_relabel
    d = g.to_dictionary()
    for k in list(d.keys()):
        if len(d[k]) == 0:
            del d[k]
    if not d:
        return x.parent()(1)
    # compute ord_links
    if not links:
        k0 = d.keys()[0]
        links = [k0, d[k0][0]]
        ord_links = ordered_links(d, *links)
    else:
        ord_links = links

    # check that ord_links has the right number of edges
    num_edges = sum([len(v) for v in d.values()]) // 2
    if num_edges != len(ord_links):
        raise ValueError('wrong number of links')

    if not labels:
        p = gen_count_hobj(ord_links)
    else:
        values = {}
        for i, j, v in g.edges():
            if i < j and i in d:
                values[(i, j)] = v
        p = gen_count_hobj(ord_links, values)
    return p

