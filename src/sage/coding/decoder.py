"""
Decoding methods for linear error-correcting codes.
Methods implemented:

 * nearest neighbor
 * syndrome
 * Groebner basis
 * Groebner representation

AUTHOR:
    -- David Joyner (2009-02-01): initial version
    -- Veronica Suaste (2013-09-21): add new decoding methods from [Marquez2013]_.

REFERENCES:

    .. [Marquez2013] Irene Marquez-Corbella, "Combinatorial Commutative Algebra Approach to
    Complete Decoding", PhD Thesis, University of Valladolid, 2013.
    
TODO:
  Add lots more methods!
"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.decorators import rename_keyword
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.constructor import FiniteField as GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.misc.misc import prod

def coset_leader(C, v):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns an
    element of the syndrome of v of lowest weight.

    EXAMPLES:
        sage: C = codes.HammingCode(2,GF(3)); C
        Linear code of length 4, dimension 2 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: from sage.coding.decoder import coset_leader
        sage: coset_leader(C, v)
        ((0, 0, 1, 0), 1)
        sage: coset_leader(C, v)[0]-v in C
        True

    """
    coset = [[c + v, (c + v).hamming_weight()] for c in C]
    wts = [x[1] for x in coset]
    min_wt = min(wts)
    s = C[0]  # initializing
    w = v.hamming_weight()  # initializing
    for x in coset:
        if x[1] == min_wt:
            w = x[1]
            s = x[0]
            break
    return s, w

@rename_keyword(deprecation=6094, method="algorithm")
def decode(C, v, algorithm="syndrome"):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns an
    element in C which is closest to v in the Hamming
    metric.

    Methods implemented include "nearest neighbor" (essentially
    a brute force search), "syndrome", "groebner_representation" and
    "groebner_basis".

    Methods "groebner_representation" and "groebner_basis" were implemented
    following the work in [Marquez2013]_.

    EXAMPLES::

        sage: C = codes.HammingCode(2,GF(3))
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: v in C
        False
        sage: from sage.coding.decoder import decode
        sage: c = decode(C, v);c
        (0, 2, 2, 1)
        sage: c in C
        True
        sage: c = decode(C, v, algorithm="nearest neighbor");c
        (0, 2, 2, 1)
        sage: C = codes.HammingCode(3,GF(3)); C
        Linear code of length 13, dimension 10 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 13)
        sage: v = V([2]+[0]*12)
        sage: decode(C, v)  # long time (9s on sage.math, 2011)
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        sage: C = codes.HammingCode(3,GF(2))
        sage: v = vector(GF(2),(0,0,1,1,0,1,0))
        sage: decode(C,v,"groebner_representation")
        (1, 0, 1, 1, 0, 1, 0)

        sage: C = codes.ReedSolomonCode(5,4,GF(5))
        sage: v = vector(GF(5),(4,0,3,3,4))
        sage: decode(C,v,"groebner_basis")
        (0, 0, 1, 2, 2)
    """
    V = C.ambient_space()
    if not isinstance(v, list):
        v = v.list()
    v = V(v)
    if algorithm == "nearest neighbor":
        diffs = [[c - v, (c - v).hamming_weight()] for c in C]
        diffs.sort(lambda x, y:  x[1] - y[1])
        return diffs[0][0] + v
    if algorithm == "syndrome":
        return -V(syndrome(C, v)[0]) + v
    if algorithm == "groebner_representation":
        return decode_groebner_representation(C,v)
    if algorithm == "groebner_basis":
        return decode_groebner_basis(C, v)

def decode_groebner_basis(C, y):
    r"""
    Gradient descent decoding algorithm: decodes the received word ``y`` to an element
    ``c`` in this code using the test-set of the code ``C``.

    The algorithm is described in Algorithm 20 in page 145 of [Marquez2013]_.

    INPUT:

    - ``C`` -- a :class:`~sage.coding.linear_code.LinearCode` instance.

    - ``y`` --Vector of the same length as a codeword

    OUTPUT:

    - Vector representing a word in this code closest to ``y``

    EXAMPLES::

        sage: C = codes.WalshCode(3)
        sage: v = vector(GF(2),(0, 1, 0, 1, 1, 0, 0, 0))
        sage: from sage.coding.decoder import decode_groebner_basis
        sage: dec_word = decode_groebner_basis(C,v)
        sage: dec_word
        (0, 1, 0, 1, 1, 0, 1, 0)
        sage: dec_word in C
        True

        sage: G = matrix(GF(2),[[1,0,0,0,1,1,1,1,1,1],[0,1,0,0,0,0,1,1,1,1],[0,0,1,0,0,1,0,1,1,1],[0,0,0,1,0,1,1,0,1,1]])
        sage: C = LinearCode(G)
        sage: v = vector(GF(2),(0, 1, 0, 1, 0, 1, 1, 0, 0, 1))
        sage: decode_groebner_basis(C,v)
        (0, 1, 0, 1, 0, 1, 0, 1, 0, 0)

        sage: C = codes.HammingCode(4,GF(2))
        sage: v = vector(GF(2),(0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0))
        sage: decode_groebner_basis(C,v)
        (0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0)

        sage: C = codes.BCHCode(8,3,GF(3))
        sage: v = vector(GF(3),(2, 2, 1, 2, 1, 2, 0, 0))
        sage: decode_groebner_basis(C,v)
        (1, 2, 1, 2, 2, 2, 2, 0)

        sage: F = GF(4,'a')
        sage: a = F.primitive_element()
        sage: C = codes.HammingCode(2,F)
        sage: v = vector(F,(a, 0, a + 1, a + 1, a + 1))
        sage: decode_groebner_basis(C,v)
        (a + 1, 0, a + 1, a + 1, a + 1)

        sage: C = codes.ReedSolomonCode(6,4,GF(7))
        sage: v = vector(GF(7),(5, 5, 4, 3, 6, 2))
        sage: decode_groebner_basis(C,v)
        (0, 6, 3, 5, 5, 3)
    """
    t_s = test_set_groebner(C)
    c = vector(C.base_ring(),C.length())
    y = vector(C.base_ring(),y)
    for t in t_s:
        if (y-t).hamming_weight() < y.hamming_weight():
            c = c+t
            y = y+t
    return c

def decode_groebner_representation(C, y):
    r"""
    Gradient descent decoding algorithm: decodes the received word ``y`` to an element
    ``c`` in the code ``C`` using its groebner representation.

    The algorithm is described in Algorithm 12 in page 72 of [Marquez2013]_.

    INPUT:

    - ``C`` -- a :class:`~sage.coding.linear_code.LinearCode` instance.

    - ``y`` --Vector of the same length as a codeword

    OUTPUT:

    - Vector representing a word in this code closest to ``y``

    EXAMPLES::

        sage: C = codes.HammingCode(4,GF(2))
        sage: v = vector(GF(2),(0,1,1,1,0,1,1,1,0,1,0,1,0,1,0))
        sage: from sage.coding.decoder import decode_groebner_representation
        sage: decode_groebner_representation(C,v)
        (0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0)

        sage: G = matrix(GF(2),[[1,0,0,1,1,1],[0,1,0,0,1,1],[0,0,1,1,0,1]])
        sage: C = LinearCode(G)
        sage: v = vector(GF(2),(1,1,1,1,1,0))
        sage: decode_groebner_representation(C,v)
        (0, 1, 1, 1, 1, 0)
        sage: v = vector(GF(2),(1,0,1,0,1,1))
        sage: decode_groebner_representation(C,v)
        (1, 0, 1, 0, 1, 0)
    """
    if not C.base_ring().order() == 2:
        raise NotImplementedError("The 'decode_groebner_representation' algorithm is only implemented for binary codes")
    GR_dic = groebner_representation(C)[1]
    s = y.support()
    y_t = (0,)*C.length()
    for i in range (C.length()):
        if i in s:
            y_t = GR_dic[(y_t,i)]
    return vector(GF(2),y_t)+y

def groebner_basis_fglm(C):
    """
    This function computes the groebner basis of the ideal associated to
    code ``C``, using an adapted fglm algorithm for this case, and
    a graduated order implicitly.
    In this algorithm we use vectors of length's code dimension, the value
    entry ``x`` in position ``i`` represents the variable ``x_{ij}`` where
    ``j`` is given by ``x = \alpha^j``, ``\alpha`` is a generator of the
    finite field of ``C``.

    The algorithm is described in Algorithm 21 in page 148  of [Marquez2013]_.

    INPUT:

    - ``C`` -- a :class:`~sage.coding.linear_code.LinearCode` instance.

    OUTPUT:

    - Generator iterable object with groebner basis elements as vectors.

    EXAMPLES::

        sage: C = codes.HammingCode(2,GF(3))
        sage: from sage.coding.decoder import groebner_basis_fglm
        sage: list(groebner_basis_fglm(C))
        [((0, 0, 2, 2), (1, 0, 0, 0)),
         ((0, 0, 1, 2), (0, 2, 0, 0)),
         ((0, 0, 2, 1), (0, 1, 0, 0)),
         ((0, 0, 1, 1), (2, 0, 0, 0)),
         ((0, 2, 0, 2), (2, 0, 0, 0)),
         ((0, 1, 0, 2), (0, 0, 2, 0)),
         ((0, 2, 0, 1), (0, 0, 1, 0)),
         ((0, 1, 0, 1), (1, 0, 0, 0)),
         ((0, 2, 2, 0), (0, 0, 0, 2)),
         ((0, 1, 2, 0), (0, 2, 0, 2)),
         ((0, 2, 1, 0), (0, 0, 2, 2)),
         ((0, 1, 1, 0), (0, 0, 0, 1)),
         ((2, 0, 0, 2), (0, 0, 1, 0)),
         ((1, 0, 0, 2), (0, 0, 2, 1)),
         ((2, 0, 0, 1), (0, 2, 0, 0)),
         ((1, 0, 0, 1), (0, 0, 2, 0)),
         ((2, 0, 2, 0), (0, 0, 0, 1)),
         ((1, 0, 2, 0), (2, 0, 0, 1)),
         ((2, 0, 1, 0), (0, 0, 2, 1)),
         ((1, 0, 1, 0), (0, 0, 0, 2)),
         ((1, 2, 0, 0), (0, 0, 0, 1)),
         ((2, 1, 0, 0), (0, 0, 0, 2)),
         ((1, 1, 0, 0), (0, 2, 0, 1))]
    """
    from sage.combinat.combination import Combinations
    from sage.combinat.cartesian_product import CartesianProduct
    n = C.length()
    maxdegree = n - C.dimension()+2
    Fq = C.base_ring()
    v1 = vector(Fq,n)
    genMat = [v1]
    alpha = Fq.primitive_element()
    Fqstar = Fq.list()[1:]
    One = Fq(1)
    Fqstar_wmax = [Fqstar]*(n-C.dimension())
    if Fq.is_prime_field():
        for g in C.gen_mat():
            genMat.append(g)
    else:
        for g in C.gen_mat():
            genMat.extend([a*g for a in Fqstar])
    #stores leader terms of groebner basis
    #in a convenient way to check for multiples
    grob_bb=[]
    w2 = {}
    for g in genMat:
        w2[tuple(g)]= v1
    for wt in xrange(1,maxdegree):
        w1 = []
        for one_pos in xrange(n-wt+1):
            for c in Combinations(xrange(one_pos+1,n),wt-1):
                v = v1.__copy__()
                v[one_pos] = One
                for values in CartesianProduct(*(Fqstar_wmax[:wt-1])):
                    for vi,ci in zip(values,c):
                        v[ci] = vi
                    if not multiple_fq(v,grob_bb):
                        for qe in Fqstar:
                            v[one_pos] = qe
                            w1.append(v.__copy__())
        if not w1:
            break
        while w1:
            v = w1.pop()
            for g in genMat:
                temp = tuple(v + g)
                if temp in w2:
                    grob_bb.append([set(v.support()),v.list_from_positions(v.support())])
                    yield (v,w2[temp])
                    break
                w2[temp] = v

def groebner_basis_singular(C, order="degrevlex"):
    r"""
    Computes the grobner basis of the ideal associated to linear code ``C`` w.r.t
    monomial ordering ``order``. Using algorithm fglm from singular.
    We use this function only for binary linear codes.

    INPUT:

    - ``C`` -- a :class:`~sage.coding.linear_code.LinearCode` instance.

    - ``order`` --string (default:``"degrevlex"``) -- a degree ordering
    See :mod:`~sage.rings.polynomial.term_order` for the orderings.

    OUTPUT:

    - Generator iterable object of polynomials representing a reduced groebner basis.

    EXAMPLES::

        sage: C = codes.WalshCode(2)
        sage: from sage.coding.decoder import groebner_basis_singular
        sage: list(groebner_basis_singular(C))
        [x0^2 + 1, x3^2 + 1, x1 + x3, x2 + x3]

        sage: C = codes.HammingCode(3,GF(2))
        sage: list(groebner_basis_singular(C))
        [x0^2 + 1,
         x0*x1 + x2,
         x1^2 + 1,
         x0*x2 + x1,
         x1*x2 + x0,
         x2^2 + 1,
         x0*x3 + x4,
         x1*x3 + x5,
         x2*x3 + x6,
         x3^2 + 1,
         x0*x4 + x3,
         x1*x4 + x6,
         x2*x4 + x5,
         x3*x4 + x0,
         x4^2 + 1,
         x0*x5 + x6,
         x1*x5 + x3,
         x2*x5 + x4,
         x3*x5 + x1,
         x4*x5 + x2,
         x5^2 + 1,
         x0*x6 + x5,
         x1*x6 + x4,
         x2*x6 + x3,
         x3*x6 + x2,
         x4*x6 + x1,
         x5*x6 + x0,
         x6^2 + 1]

        sage: G = Matrix(GF(2),[[1,0,1,1],[0,1,1,0]])
        sage: C = LinearCode(G)
        sage: list(groebner_basis_singular(C))
        [x0^2 + 1, x0*x2 + x3, x2^2 + 1, x0*x3 + x2, x2*x3 + x0, x3^2 + 1, x1 + x2]
    """
    if not C.base_ring().order() == 2:
        raise NotImplementedError("The groebner_basis_singular function for a code is only implemented for binary codes")
    R = PolynomialRing(GF(2),C.length(),'x',order = order)
    Rgens = R.gens()
    gens = []
    for g in C.gen_mat():
        p = prod(Rgens[i] for i in g.support())
        gens.append(p-1)
    I = R.ideal([_**2 -1 for _ in R.gens()]+ gens)
    return I.groebner_basis('libsingular:stdfglm')

def groebner_representation(C, order="degrevlex"):
    r"""
    Returns a Groebner representation of the code ``C``.
    A Groebner representation of an `[n,k]` binary linear code
    ``C`` is a pair `(N,\Phi)` such that:

    #. `N` is a transversal of the cosets in `F_2^n/C`
    (i.e one element of each coset) verifying that `0` belongs to `N` and for
    each `n` in N`\`0` there exists `e_i`
    with `i` in `\{1,...,n\}` such that `n = n' + e_i` with `n'` in  `N`
    #. `\phi` is a function called ``Matphi function`` that maps each pair `(n,e_i)`
    with `n` in `N` to the element of `N` representing the coset that contains `n+e_i`.

    The algorithm is described in Algorithm 6 in page 52 of [Marquez2013]_.

    INPUT:

    - ``C`` -- a :class:`~sage.coding.linear_code.LinearCode` instance.

    - ``order`` -- String (default:``"degrevlex"``) -- a degree ordering
    See :mod:`~sage.rings.polynomial.term_order` for the orderings.

    OUTPUT:

    - List with one vector of each coset of the code ``C`` `(N)`

    - Dictionary representing the ``Matphi function``  `(\Phi)`

    EXAMPLE::

    sage: H = matrix(GF(2),[[1,0,0,1,1,1],[0,1,0,1,0,1],[0,0,1,0,1,1]])
    sage: C = codes.LinearCodeFromCheckMatrix(H)
    sage: from sage.coding.decoder import groebner_representation
    sage: GR = groebner_representation(C)
    sage: dic = GR[1]
    sage: dic[((0,1,0,0,0,0),1)]
    (0, 0, 0, 0, 0, 0)
    sage: dic[((0,1,0,0,0,0),2)]
    (1, 0, 0, 0, 0, 1)
    """
    if not C.base_ring().order() == 2:
        raise NotImplementedError("The 'groebner_representation' function is only implemented for binary codes")
    from sage.rings.polynomial.polydict import ETuple
    from sage.rings.polynomial.term_order import TermOrder
    t_order = TermOrder(order)
    H = C.check_mat().transpose()
    B = C.ambient_space().basis()
    S = []
    N = []
    List = [ETuple(vector(GF(2),C.length()).list())]
    PHI = {}
    while List:
        w = vector(GF(2),List.pop())
        s = w*H
        j = 0
        if s in S:
            j = S.index(s)
        if j:
            for k in w.support():
                if (w-B[k]) in N:
                    PHI[(tuple(w-B[k]),k)] = tuple(N[j])
        else:
            N.append(w)
            S.append(s)
            List =insert_nextnew(w,List,t_order)
            for k in w.support():
                if (w-B[k]) in N:
                    PHI[(tuple(w-B[k]),k)] = tuple(w)
                    PHI[(tuple(w),k)] = tuple(w-B[k])
    return [N,PHI]

def insert_nextnew(v, List, order):
    r"""
    Inserts vectors ``v + e_i``, with ``i`` not in support of ``v`` and ``e_i`` standard basis vector,
    in List according the order ``order``.
    Function for internal use only.

    INPUT:

    - ``v``--Vector
    - ``List``-- List which stores vectors. Not necessarily empty.
    - ``order``--string (default:``"degrevlex"``) -- a degree ordering
    See :mod:`~sage.rings.polynomial.term_order` for the orderings.

    EXAMPLE::

        sage: C = codes.HammingCode(3,GF(2))
        sage: v = vector(GF(2),C.length())
        sage: List = []
        sage: t = TermOrder("degrevlex")
        sage: from sage.coding.decoder import insert_nextnew
        sage: insert_nextnew(v,List,t)
        [(1, 0, 0, 0, 0, 0, 0),
         (0, 1, 0, 0, 0, 0, 0),
         (0, 0, 1, 0, 0, 0, 0),
         (0, 0, 0, 1, 0, 0, 0),
         (0, 0, 0, 0, 1, 0, 0),
         (0, 0, 0, 0, 0, 1, 0),
         (0, 0, 0, 0, 0, 0, 1)]
    """
    s = [i for i in range(v.degree()) if i not in v.support()]
    from sage.rings.polynomial.polydict import ETuple
    from sage.matrix.constructor import identity_matrix
    B = identity_matrix(GF(2),v.degree())
    comp = getattr(order,'compare_tuples_'+order.name())
    if not List:
        List.append(ETuple((v+B[s[0]]).list()))
    for i in s:
        new = ETuple((v+B[i]).list())
        if comp(List[-1],new)==1:
            List.append(new)
            continue
        j = 0
        while comp(List[j],new)==1:
            j+=1
        if new != List[j]:
            List.insert(j,new)
    return List

def multiple_fq(w, groebner_basis):
    r"""
    This function checks if a polynomial term represented(exponents) by ``w`` is
    multiple of any leader term of the Groebner basis ``groebner_basis``.
    For internal use only.

    INPUT:

    - ``w`` -- Vector representing the exponents of a polynomial term
    - ``groebner_basis`` -- List representing the groebner basis elements so far.
    The variables we are working with are in the form `x_{ij}`.
    So the first entry of each element in ``groebner_basis`` must be the support
    of a vector in which each entry indicate the index `i` of the variable
    that is present in the leader term of the groebner basis element.
    And the second entry are the values `j` for each entry `i` in the support.

    OUTPUT:

    - ``True`` if polynomial term represented by ``w`` is multiple of
    any leader term of the groebner representation. ``False`` otherwise.

    EXAMPLES::

        sage: G=matrix(GF(3),[[1,2,1],[0,1,1]])
        sage: C=LinearCode(G)
        sage: from sage.coding.decoder import groebner_basis_fglm
        sage: gb=list(groebner_basis_fglm(C))
        sage: gb
        [((0, 2, 0), (0, 0, 1)), ((0, 1, 0), (0, 0, 2)), ((1, 0, 0), (0, 2, 0))]
        sage: GB=[[set(g[0].support()),g[0].list_from_positions(g[0].support())] for g in gb]
        sage: w = vector(GF(3),(1,2,0))
        sage: from sage.coding.decoder import multiple_fq
        sage: multiple_fq(w,GB)
        True
    """
    ws = set(w.support())
    for g in groebner_basis:
        if g[0].issubset(ws):
            wsl = [w[i] for i in g[0]]
            if g[1]==wsl:
                return True
    return False

def syndrome(C, v):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns the
    elements in V (including v) which belong to the
    syndrome of v (ie, the coset v+C, sorted by weight).

    EXAMPLES::

        sage: C = codes.HammingCode(2,GF(3)); C
        Linear code of length 4, dimension 2 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: from sage.coding.decoder import syndrome
        sage: syndrome(C, v)
         [(0, 0, 1, 0), (0, 2, 0, 1), (2, 0, 0, 2), (1, 1, 0, 0), (2, 2, 2, 0), (1, 0, 2, 1), (0, 1, 2, 2), (1, 2, 1, 2), (2, 1, 1, 1)]

    """
    V = C.ambient_space()
    if not isinstance(v, list):
        v = v.list()
    v = V(v)
    coset = [[c + v, (c + v).hamming_weight()] for c in C]
    coset.sort(lambda x, y: x[1] - y[1])
    return [x[0] for x in coset]

def test_set_groebner(C):
    r"""
    Let `G = \{g_1,...,g_s\}` be a reduced Grobner basis of the ideal
    associated to the code ``C``. Where `g_i = X^{g_i^+} - X^{g_i^-}`
    The set `T = \{ g_i^+ - g_i^- \vert i = 1,...,s\}` is a test-set for
    code ``C``.
    This function returns the test-set for ``C``.
    For internal use only.

    This algorithm is discribed in Proposition 4.30 in page 142 of [Marquez2013]_.

    .. Note::

        The test set doesn't contain duplicates, so cardinality of test-set is less equal
        than groebner basis cardinality.

    OUTPUT:

    - List of vectors representing the test-set.

    EXAMPLES::

        sage: C = codes.HammingCode(3,GF(2))
        sage: from sage.coding.decoder import groebner_basis_singular
        sage: len(groebner_basis_singular(C))
        28
        sage: from sage.coding.decoder import test_set_groebner
        sage: test_set_groebner(C)
        [(0, 0, 0, 0, 0, 0, 0),
         (1, 1, 1, 0, 0, 0, 0),
         (1, 0, 0, 1, 1, 0, 0),
         (0, 1, 0, 1, 0, 1, 0),
         (0, 0, 1, 1, 0, 0, 1),
         (0, 1, 0, 0, 1, 0, 1),
         (0, 0, 1, 0, 1, 1, 0),
         (1, 0, 0, 0, 0, 1, 1)]

        sage: C = codes.BCHCode(8,3,GF(3))
        sage: test_set_groebner(C)
        [(0, 0, 0, 0, 2, 0, 2, 2),
         (0, 0, 0, 0, 1, 0, 1, 1),
         (0, 0, 0, 2, 0, 2, 2, 0),
         (0, 0, 0, 1, 0, 1, 1, 0),
         (0, 0, 2, 0, 2, 2, 0, 0),
         (0, 0, 1, 0, 1, 1, 0, 0),
         (0, 0, 0, 1, 2, 1, 0, 2),
         (0, 0, 0, 2, 1, 2, 0, 1),
         (0, 2, 0, 2, 2, 0, 0, 0),
         (0, 1, 0, 1, 1, 0, 0, 0),
         (0, 0, 1, 0, 0, 1, 2, 2),
         (0, 0, 2, 0, 0, 2, 1, 1),
         (0, 0, 1, 2, 1, 0, 2, 0),
         (2, 0, 2, 2, 0, 0, 0, 0),
         (1, 0, 1, 1, 0, 0, 0, 0),
         (0, 1, 0, 1, 0, 0, 2, 2),
         (0, 1, 0, 0, 1, 2, 2, 0),
         (0, 2, 0, 0, 2, 1, 1, 0),
         (0, 1, 2, 1, 0, 2, 0, 0),
         (0, 2, 1, 2, 0, 1, 0, 0),
         (1, 0, 1, 0, 0, 2, 2, 0),
         (1, 0, 0, 1, 2, 2, 0, 0),
         (2, 0, 0, 2, 1, 1, 0, 0),
         (1, 2, 1, 0, 2, 0, 0, 0),
         (2, 1, 2, 0, 1, 0, 0, 0)]
    """
    if C._get_groebner_test_set() is None:
        test_set =[]
        if C.base_ring().order() == 2:
            GB = groebner_basis_singular(C)
            t_s = [vector(GF(2),gb.exponents()[0].eadd(gb.exponents()[1])) for gb in GB]
        else:
            GB = groebner_basis_fglm(C)
            t_s = [g[0]-g[1] for g in GB]
        for t in t_s:
            if t not in test_set:
                test_set.append(t)
        C._set_groebner_test_set(test_set)
    return C._get_groebner_test_set()
