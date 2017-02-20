from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.matrix.special import companion_matrix, block_matrix


def min_orb(X, v):
    """
    Return the minimal polynomial of a matrix `X` acting on a vector `v`.

    This computes the minimal polynomial Min_X(v) and a basis for the
    space Orb_X(v).

    INPUT:

    - `X` -- matrix acting on a vector space `V`

    - `v` -- vector in the vector space `V`

    OUTPUT:

    a polynomial in one variable and a tuple of vectors

    EXAMPLES::

        sage: from sage.matrix.generalized_jordan_form import min_orb
        sage: V = FreeModule(QQ, 4)
        sage: X = matrix(QQ, 4, 4, [0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1])
        sage: b = V.basis()
        sage: v = b[0] + 4 * b[1]
        sage: min_orb(X, v)
        (X^3 - 1, [(1, 4, 0, 0), (0, 1, 4, 0), (4, 0, 1, 0)])
        sage: min_orb(X, b[3])
        (X - 1, [(0, 0, 0, 1)])

        sage: v = vector(QQ,[1,4,3,0])
        sage: min_orb(X,v)
        (X^3 - 1, [(1, 4, 3, 0), (3, 1, 4, 0), (4, 3, 1, 0)])
    """
    V = X._row_ambient_module()
    liste = [v]
    last = v
    ring = V.base_ring()
    while True:
        last *= X
        liste.append(last)
        relation = V.linear_dependence(liste)
        if relation:
            break
    poly = PolynomialRing(ring, 'X')(list(relation[0]))
    return poly / poly.leading_coefficient(), liste[:-1]


def decompose(X):
    """
    Return a decomposition of the vector space `V` under the action of `X`.

    THIS IS NOT MATHEMATICALLY CHECKED YET ! BEWARE !

    INPUT:

    - `X` -- matrix acting on a vector space `V`

    OUTPUT:

    - a direct sum decomposition of `V` into `X`-invariant subspaces
      `V_1`, ..., `V_k`

    - irreducible polynomials `p_1`, ..., `p_k`

    - integers `e_1`, ..., `e_k`

    - vectors v_1, ..., v_k, with v_i either zero or generating V_i under X

    The minimal polynomial of `X` acting on `V_i` is `p_i^{e_i}`.

    EXAMPLES::

        sage: from sage.matrix.generalized_jordan_form import decompose
        sage: X = matrix(QQ,4,4,[0,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1])
        sage: decompose(X)
        ([Vector space of degree 4 and dimension 1 over Rational Field
          Basis matrix:
          [1 1 1 0], Vector space of degree 4 and dimension 1 over
          Rational Field
          Basis matrix:
          [0 0 0 1], Vector space of degree 4 and dimension 2 over
          Rational Field
          Basis matrix:
          [ 1  0 -1  0]
          [ 0  1 -1  0]],
         [X - 1, X - 1, X^2 + X + 1],
         [1, 1, 1],
         [(1, 1, 1, 0), (0, 0, 0, 1), (-1, 1, 0, 0)])
    """
    V = X._row_ambient_module()
    ring = V.base_ring()
    polynomes = PolynomialRing(ring, 'X')
    Gen = {}
    g = []
    L = []
    dimV = V.dimension()
    nullspace = V.span([])
    subV = V.span([])
    while subV.dimension() < dimV:
        for v in V.basis():
            if v not in subV:
                break
        m = min_orb(X, v)[0]
        p = polynomes.one()
        for i in range(len(L)):
            g_i = g[i]
            if m.gcd(g_i) != 1:
                f = g_i ** m.valuation(g_i)
                w = v * (m // f)(X)
                p = p * f
                W = V.span(min_orb(X, w)[1])
                # now compare W with known subspaces in L[i]
                sub_U = []
                notsub_U = []
                for U in L[i]:
                    if U.is_subspace(W):
                        sub_U.append(U)
                    else:
                        notsub_U.append(U)
                full_span = sum([U for U in L[i]], nullspace)
                if not W.intersection(full_span).dimension():
                    Gen[W] = w
                    L[i].append(W)
                elif W.is_subspace(full_span):
                    pass
                elif not W.intersection(sum(notsub_U, nullspace)).dimension():
                    Gen[W] = w
                    L[i] = notsub_U + [W]
                else:
                    # todo, find optimal set of U in L[i]
                    T = W + full_span
                    Gen[T] = 0
                    L[i] = [T]

        v *= p(X)
        m = m // p

        for qi, ei in m.factor():
            w = v * (m // (qi ** ei))(X)
            W = V.span(min_orb(X, w)[1])
            Gen[W] = w
            L.append([W])
            g.append(qi)

        subV = sum(U for Li in L for U in Li)

    V = []
    p = []
    e = []
    v = []
    for i, L_i in enumerate(L):
        g_i = g[i]
        for U in L_i:
            V.append(U)
            p.append(g_i)
            # HERE BELOW IS THE POINT THAT IS NOT CLEAR ENOUGH ! BEWARE !
            e.append(U.dimension() // g_i.degree())

            v.append(Gen[U])
    return V, p, e, v


def split_primary(W, X, f, q):
    """
    Split a primary subspace `W` into a direct sum of orbit spaces.

    INPUT:

    - `W` -- subspace of the ambient vector space `V`
    - `X` -- matrix acting on `V`
    - `f` -- irreducible polynomial in one variable
    - `q` -- integer

    OUTPUT:

    - `m_1`, ..., `m_q` -- integers
    - `w_{i,j}` -- vectors in `W`

    This finds a complete cyclic decomposition of `W` for `X`, assuming the
    subspace `W` is primary for `X`.

    EXAMPLES::

        sage: from sage.matrix.generalized_jordan_form import split_primary
        sage: V = FreeModule(QQ,4)
        sage: W = V
        sage: X = matrix(QQ,4,4,[1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,1])
        sage: f = polygen(QQ,'x') - 1
        sage: q = 3
        sage: split_primary(W, X, f, q)
        ({1: 1, 2: 0, 3: 1}, {(1, 1): (1, 0, 0, 0), (3, 1): (0, 0, 0, 1)})
    """
    V = X._row_ambient_module()
    N = f(X)
    WW = [(N ** i).left_kernel().intersection(W) for i in range(q + 1)]

    B = {i: [v for v in WW[i].basis() if v not in WW[i - 1]]
         for i in range(1, q + 1)}

    m = {}
    ww = {}
    d = f.degree()
    S = {q + 1: []}
    for i in range(q, 0, -1):
        S[i] = [v * N for v in S[i + 1]]
        m[i] = ZZ.zero()
        for w in B[i]:
            if w not in V.span(S[i]) + WW[i - 1]:
                S[i] += [w * X ** k for k in range(d)]
                m[i] += 1
                ww[i, m[i]] = w
    return m, ww


def N_matrix(ring, d):
    """
    Return an auxiliary square matrix.

    INPUT:

    - ``ring`` -- a ring

    - `d` -- an integer

    OUTPUT:

    a square matrix of size `d` with just one non-zero coefficient

    EXAMPLES::

        sage: from sage.matrix.generalized_jordan_form import N_matrix
        sage: N_matrix(QQ,4)
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [1 0 0 0]
    """
    m = matrix(ring, d, d, 0)
    m[d - 1, 0] = ring.one()
    return m


def J_matrix(p, e):
    """
    Return an auxiliary square matrix.

    INPUT:

    - p -- a polynomial in one variable
    - e -- an integer

    EXAMPLES::

        sage: from sage.matrix.generalized_jordan_form import J_matrix
        sage: x = polygen(ZZ,'x')
        sage: J_matrix(x**2-3*x+5,3)
        [ 0  1  0  0  0  0]
        [-5  3  1  0  0  0]
        [ 0  0  0  1  0  0]
        [ 0  0 -5  3  1  0]
        [ 0  0  0  0  0  1]
        [ 0  0  0  0 -5  3]
    """
    ring = p.base_ring()
    C = companion_matrix(p, format='bottom')
    N = N_matrix(ring, p.degree())
    step = [C, N] + [0] * (e - 1)
    step = step * (e - 1) + [C]
    return block_matrix(ring, e, e, step, subdivide=False)


def TJ_matrix(X, v, p, e):
    """
    Return an auxiliary rectangular matrix.

    INPUT:

    - X -- a square matrix
    - v -- a vector
    - p -- a polynomial in one variable
    - e -- an integer

    EXAMPLES::

        sage: from sage.matrix.generalized_jordan_form import TJ_matrix
        sage: x = polygen(QQ,x)
        sage: g = (x-4)*(x-3)
        sage: C = companion_matrix(x**4-1)
        sage: TJ_matrix(C, vector(QQ,[1,1,0,1]), g**2, 2)
        [   1    1    0    1]
        [   1    0    1    1]
        [   0    1    1    1]
        [   1    1    1    0]
        [ -37  204 -109   50]
        [ 204 -109   50  -37]
        [-109   50  -37  204]
        [  50  -37  204 -109]
    """
    d = p.degree()
    pX = p(X)
    step = [v * pX ** j * X ** i for j in range(e) for i in range(d)]
    return matrix(p.base_ring(), d * e, len(list(v)), step)
