from sage.rings.infinity import Infinity
from sage.rings.all import GF, ZZ, mod, IntegerModRing
from copy import copy
from sage.matrix.all import matrix
from sage.modules.all import vector
from sage.groups.all import GO,Sp

def Hensel_qf(G, F, a, b):
    r"""
    Lift ``F`` modulo ``p^n`` satisfying ``G == F * G * F.T``.

    INPUT:

     - ``G`` -- a block diagonal matrix of the form
    ``[G0*p^n0, G1*p^n1, ... , Gk*p^nk]``
    with integers `nk < .... < n1 < n0`
    and `Gi` unimodular and symmetric.
    - ``F`` -- invertible `p`-adic matrix
      such that `(G, G, F)` is `a`-adapted
    - ``a`` -- integer the starting precision
    - ``b```-- integer the target precision

    OUTPUT:

    - ``Fk`` -- the lift of `F` such that
    ``Z == F * G * F.T`` modulo `p^n` with `n = prec`

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import Hensel_qf, _min_val
        sage: R = Zp(3, type='fixed-mod', prec=7, print_mode='terse', show_prec=False)
        sage: G = matrix(R, 6, 6, [0, 243, 0, 0, 0, 0, 243, 0, 0, 0, 0, 0, 0, 0, 0, 27, 0, 0, 0, 0, 27, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])
        sage: F = matrix(R, 6, 6, [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 2, 1, 0, 0, 0, 2, -1, 6, 3, 0, 1, 1, 9, 3, 6, 1, 0])
        sage: G
        [  0 243   0   0   0   0]
        [243   0   0   0   0   0]
        [  0   0   0  27   0   0]
        [  0   0  27   0   0   0]
        [  0   0   0   0   0   1]
        [  0   0   0   0   1   0]
        sage: F
        [   0    1    0    0    0    0]
        [   1    0    0    0    0    0]
        [   2    1    0    1    0    0]
        [   1    2    1    0    0    0]
        [   2 2186    6    3    0    1]
        [   1    9    3    6    1    0]
        sage: Flift = Hensel_qf(G, F, 1, 6)
        sage: _min_val(Flift*G*Flift.T-G) >= 6
        True

        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False)
        sage: U = matrix(R, 2, [0, 1, 1 ,0])
        sage: V = matrix(R, 2, [2, 1, 1 ,2])
        sage: G = matrix.block_diagonal([2*U, 2*U, V])
        sage: F = matrix(R, 6, [1, 0, 0, 0, 0, 0,
        ....:                   1, 1, 1, 1, 0, 0,
        ....:                   1, 0, 1, 0, 0, 0,
        ....:                   1, 0, 0, 1, 0, 0,
        ....:                   1, 0, 0, 0, 1, 1,
        ....:                   0, 0, 0, 1, 0, 1])
        sage: Fl = Hensel_qf(G, F, 1, 6)
        sage: _min_val(Fl*G*Fl.T-G)
        10
    """
    # Input checks
    if F.determinant().valuation() != 0:
        raise ValueError("F must be invertible")
    if not (G.ncols() == F.ncols() and G.nrows() == F.nrows()):
        raise ValueError("G, F must have the same size")
    if not G.base_ring() == F.base_ring():
        raise ValueError("not the same basering")
    if G != G.T:
        raise ValueError("G must be symmetric")
    R = G.base_ring()
    n = R.precision_cap()
    if b > n:
        raise ValueError("Desired precision is higher than base ring precision")
    for k in range(0,G.ncols()-1):
        n1 = _min_val(G[k,:])
        n2 = _min_val(G[k+1,:])
        if n1 < n2:
            raise ValueError("block valuations must be descending")
    # this is not enough ... mhhh
    if _min_val(F*G*F.T-G) < a:
        raise ValueError("F must satisfy Z == F * G * F.T  modulo p^a.")

    # the real worker
    F = copy(F) # leave input unchanged
    F = _Hensel_qf(G, G, F, a, b) #works inplace
    return F

def _reverse_homogeneous_blocks(G):
    r"""

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _reverse_homogeneous_blocks
        sage: R = Zp(3,type='fixed-mod',prec=4,print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R, [3^3,3^3,3^1,1,1,1])
        sage: U = _reverse_homogeneous_blocks(G)
        sage: G
        [27  0  0  0  0  0]
        [ 0 27  0  0  0  0]
        [ 0  0  3  0  0  0]
        [ 0  0  0  1  0  0]
        [ 0  0  0  0  1  0]
        [ 0  0  0  0  0  1]
        sage: Gr = U.T*G*U
        sage: Gr
        [ 1  0  0  0  0  0]
        [ 0  1  0  0  0  0]
        [ 0  0  1  0  0  0]
        [ 0  0  0  3  0  0]
        [ 0  0  0  0 27  0]
        [ 0  0  0  0  0 27]
        sage: Ur = _reverse_homogeneous_blocks(Gr)
        sage: Ur.T*Gr*Ur == G
        True
    """
    from sage.combinat.permutation import Permutation
    ind, _ = _block_indices_vals(G)
    ind.append(G.ncols())
    perm = []
    for k in range(len(ind)-1):
        perm = range(ind[k] + 1, ind[k+1] + 1) + perm
    perm = Permutation(perm)
    perm = perm.to_matrix()
    return perm

def _last_block_index(G):
    r"""
    Return the starting index of the last modular block.

    INPUT:

    - ``G`` -- a `p`-adic block diagonal matrix
      with valuation of the blocks descending

    OUTPUT:

    - ``(index, valuation)`` -- with ``i`` the maximal index such that
      G[i:,i:] is modular and valuation is the minimum valuation of ``G[i:,i:]``

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _last_block_index
        sage: G = matrix(Zp(2), 3, [4,0,0, 0,2,1, 0,1,2])
        sage: _last_block_index(G)
        (1, 0, 2)
    """
    n = G.ncols()
    val = _min_val(G[n-1,:])
    val_current = val
    for k in range(2, n+1):
        val_current = _min_val(G[n-k,:])
        if val != val_current:
            return (n-k+1, val, val_current)
    return 0, val, val_current

def _block_indices_vals(G):
    r"""

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _block_indices_vals
        sage: G = matrix(Zp(2), 3, [4,0,0, 0,2,1, 0,1,2])
        sage: _block_indices_vals(G)
        ([0, 1], [2, 0])
    """
    indices = []
    valuations = []
    while G.ncols() != 0:
        i, val, _ = _last_block_index(G)
        indices.append(i)
        valuations.append(val)
        G = G[:i,:i]
    indices.reverse()
    valuations.reverse()
    return indices, valuations

def _min_val(M):
    r"""
    Return the minimum valuation of an entry of M

    Helper function for :meth:`Hensel` No input checks.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _min_val
        sage: G = matrix.diagonal(Zp(2), [2,6,8,20,14])
        sage: _min_val(G)
        1

    TESTS::

        sage: G = matrix(Zp(2),[])
        sage: _min_val(G)
        +Infinity
    """
    L = [x.valuation() for x in M.list()]
    if len(L) == 0:
        return Infinity
    else:
        return min(L)

def _Hensel_qf(Z, G, F, a, b):
    r"""
    The real worker for :meth:`Hensel`.

    INPUT:

    - ``Z`` -- symmetric `p`-adic `n \times n` matrix
    - ``G`` -- symmetric `p`-adic `n \times n` matrix
    - ``F`` -- symmetric `p`-adic `n \times n` matrix
    - ``a`` -- integer
    - ``b`` -- integer

    We require that the triple `(Z,G,F)` is `a`-adapted.

    OUTPUT:

    Fl such that ``(Z,G,Fl)`` is `b`-adapted.

    EXAMPLES::
    """
    i, s1, s2 = _last_block_index(G)
    s = s2 - s1
    if s == 0:
        assert i == 0
        s = b-a
    if G.base_ring().prime() == 2:
        _Hensel_qf_modular = _Hensel_qf_modular_even
    else:
        _Hensel_qf_modular = _Hensel_qf_modular_odd
    Zn = Z[i:,i:] - F[i:,:i]*G[:i,:i]*F[i:,:i].T
    G2 = G[i:,i:]
    F[i:,i:] = _Hensel_qf_modular(Zn, G[i:,i:], F[i:,i:], a, b)
    if i == 0:
        return F
    while a < b:
        F[:i,i:] = (Z[:i,i:] - F[:i,:i]*G[:i,:i]*F[i:,:i].T) * (G[i:,i:]*F[i:,i:].T).inverse()
        Zn = Z[:i,:i] - F[:i,i:]*G[i:,i:]*F[:i,i:].T
        F[:i,:i] = _Hensel_qf(Zn, G[:i,:i], F[:i,:i], a, a+s)
        F[:i,i:] = (Z[:i,i:] - F[:i,:i]*G[:i,:i]*F[i:,:i].T) * (G[i:,i:]*F[i:,i:].T).inverse()
        a = a + s
    return F

def _Hensel_qf_modular_odd(Z, G, F, a, b):
    r"""
    Helper function for :meth:`_Hensel_qf`.

    No input checks. Let `p` be an odd prime number.

    INPUT:

    - ``Z`` -- symmetric `p`-adic `n \times n` matrix
    - ``G`` -- symmetric `p`-adic `n \times n` matrix
    - ``F`` -- symmetric `p`-adic `n \times n` matrix
    - ``a`` -- integer
    - ``b`` -- integer

    We require that the triple `(Z,G,F)` is `a`-adapted.

    OUTPUT:

    `Fl` such that ``(Z,G,Fl)`` is `b`-adapted.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _Hensel_qf_modular_odd,_min_val
        sage: R = Zp(3, type='fixed-mod', print_mode='terse', show_prec=False)
        sage: Gr = GO(3, 3, e=1)
        sage: G = Gr.invariant_bilinear_form().change_ring(R)
        sage: Flist = [g.matrix().change_ring(R).T for g in Gr]
        sage: error_vals = []
        sage: for F in Flist:
        ....:     Flift = _Hensel_qf_modular_odd(G, G, F, 1, 4)
        ....:     error = Flift*G*Flift.T - G
        ....:     error_vals.append(_min_val(error))
        sage: all([e>=4 for e in error_vals])
        True
    """
    while a < b:
        Y = (Z - F*G*F.T) / 2
        F = F + Y*(G*F.T).inverse()
        a = 2 * a
    return F

def _Hensel_qf_modular_even(Z, G, F, a, b):
    r"""
    Helper function for :meth:`_Hensel_qf`.

    No input checks. Let `p` be an odd prime number.

    INPUT:

    - ``Z`` -- symmetric `p`-adic `n \times n` matrix
    - ``G`` -- symmetric `p`-adic `n \times n` matrix
    - ``F`` -- symmetric `p`-adic `n \times n` matrix
    - ``a`` -- integer
    - ``b`` -- integer

    We require that the triple `(Z,G,F)` is `a`-adapted.

    OUTPUT:

    - `Fl` such that ``(Z,G,Fl)`` is `b`-adapted
    - raises a ``ValueError`` if ``F`` cannot be lifted

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _Hensel_qf_modular_even, _min_val
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False)
        sage: G = matrix(R, 4, [0,1,0,0, 1,0,0,0, 0,0,1,0, 0,0,0,5])
        sage: G
        [0 1 0 0]
        [1 0 0 0]
        [0 0 1 0]
        [0 0 0 5]
        sage: F = matrix(R, 4, [1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0])
        sage: F
        [1 0 0 0]
        [1 1 1 1]
        [1 0 0 1]
        [1 0 1 0]
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 4)
        sage: Fl
        [   1    0    0    0]
        [ 605    1    3  411]
        [ 383    0 1210  597]
        [ 795    0  297  198]
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (6, 6)
        sage: F = matrix(R, 4, [0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1])
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 4)
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (+Infinity, +Infinity)

    ::

        sage: G = matrix(R, 4, [2,1,0,0, 1,2,0,0, 0,0,3,0, 0,0,0,7])
        sage: F = matrix(R, 4, [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0])
        sage: G
        [2 1 0 0]
        [1 2 0 0]
        [0 0 3 0]
        [0 0 0 7]
        sage: F
        [1 0 0 0]
        [1 1 0 0]
        [0 0 0 1]
        [0 0 1 0]
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 5)
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (8, 8)
        sage: Z = G + matrix.identity(R, 4)*2^2
        sage: Fl = _Hensel_qf_modular_even(Z, G, F, 1, 5)
        sage: Er = Fl*G*Fl.T - Z
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (6, 6)

    TESTS::

        sage: G = matrix(R, 3, [2,1,0, 1,2,0, 0,0,7])
        sage: F = matrix(R, 3, [0, 1, 0, 1, 1, 0, 0, 0, 1])
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 8)
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (+Infinity, +Infinity)
    """
    R = Z.base_ring()
    n = Z.ncols()
    v = _min_val(G)
    G = (G/2**v).change_ring(R)
    Z = (Z/2**v).change_ring(R)
    if a == 0:
        raise ValueError("a must be a non-zero integer")
    if a == 1:
        Y = Z - F*G*F.T
        X = _solve_X(Y/2, (Y/4).diagonal(), G.inverse().diagonal())
        X = 2 * X.change_ring(R)
        F = F + X*(G*F.T).inverse()
        a = 2
    while a < b:
        Y = Z - F*G*F.T
        D = matrix.diagonal(Y.diagonal())

        # L is the the strictly lower triangular part of Y
        L = copy(Y)
        for i in range(0,n):
            L[i,i:] = 0

        F = F + (L + D/2)*(G*F.T).inverse()
        a = 2*a - 1
    # confirm computation
    assert _min_val(Z-F*G*F.T) >= b
    assert _min_val(matrix((Z-F*G*F.T).diagonal())) >= b + 1
    return F

def _solve_X(Y, b, g, ker=False):
    r"""
    Solve a certain linear equation.

    This is a helper function for :meth:`_Hensel_qf_modular_even`.

    :MATH::

        Y = X + X^T

    :MATH::

        b_i = X_{ii} + \sum_{j=1}^n X_{ij}g_j \quad i \in \{1, \dots, n\}

    INPUT:

    - ``Y`` - converts to an `n \times n` matrix over `\FF_2`
    - ``b`` - converts to `\FF_2^n`
    - ``g`` - converts to `\FF_2^n`

    OUTPUT:

    - ``X`` - an `n \times n` matrix over `\FF_2` the solution of the
      linear equation above.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _solve_X
        sage: k = GF(2)
        sage: Y = matrix(k, 3, [0,0,1, 0,0,1, 1,1,0])
        sage: b = vector(k, [1, 0, 0])
        sage: g = vector(k, [0, 1, 0])
        sage: X = _solve_X(Y, b, g)
        sage: X
        [1 0 1]
        [0 0 1]
        [0 0 0]
    """
    k = GF(2)
    Y = matrix(k, Y)
    if Y != Y.T:
        raise ValueError("Y must be symmetric")
    b = vector(k, b)
    g = vector(k, g)
    n = Y.ncols()

    equations = []
    for i in range(0, n):
        R = matrix.zero(k, n, n)
        R[i,:] = g
        R[i,i] += 1
        eqn = R.list() + [b[i]]
        equations.append(eqn)
    # equations to be a symmetric matrix
    for i in range(0, n):
        for j in range(i+1, n):
            R = matrix.zero(n, n)
            R[i, j] = 1
            R[j, i] = 1
            eq = R.list() + [Y[i,j]]
            equations.append(eq)
    A = matrix(k, equations)
    c = A[:,-1]
    A = A[:,:-1]
    # A*Xcoeff == c
    Xcoeff = A.solve_right(c)
    X = matrix(k, n, n, Xcoeff.list())
    if ker:
        Ker = []
        for Xcoeff in A.right_kernel().basis():
            X = matrix(k, n, n, Xcoeff.list())
            Ker.append(X)
        return Ker
    # confirm the computation
    assert Y == X + X.T
    for i in range(n):
        assert b[i] == X[i,i] + sum([X[i,j]*g[j] for j in range(n)])
    return X

def _mod_p_to_a_kernel(G, a):
    r"""
    """
    n = G.ncols()
    R = G.base_ring()
    E = G.parent().one()
    p = G.base_ring().prime()
    ind, val = _block_indices_vals(G)

    par = []
    ind.append(n)
    for k in range(0,len(ind)-1):
        i = ind[k+1] - 1    # last index of block i
    # add a virtual even block at the end
    val.append(val[-1]-1)

    # block diagonal contribution
    gens = []
    for k in range(len(ind)-1):
        i1 = ind[k]
        i2 = ind[k+1]
        Gk_inv = (G[i1:i2,i1:i2]/p**val[k]).inverse().change_ring(R)
        ni = i2 - i1
        Ek = matrix.identity(R, i2 - i1)
        Zk = matrix.zero(R, i2 - i1)
        if p == 2 and a == 1:
            gensk = _solve_X(matrix.zero(R,ni,ni),matrix.zero(R,1,ni).list(),Gk_inv.diagonal(),ker=True)
        else:
            # basis for the space of anti-symmetric matrices
            gensk = []
            for i in range(ni):
                for j in range(i):
                    gk = copy(Zk)
                    gk[i,j] = 1
                    gk[j,i] = -1
                    gensk.append(gk)
        for h in gensk:
            g = copy(E)
            g[i1:i2,i1:i2] = Ek + p**a*h.change_ring(R)*Gk_inv
            gens.append(g)

    # generators below the block diagonal.
    for i in range(n):
        for j in range(i):
            g = copy(E)
            g[i,j] = p**a
            flag = True
            for k in range(len(ind)-1):
                if ind[k] <= i < ind[k+1] and ind[k] <= j < ind[k+1]:
                    flag = False
                    break
            if flag:
                gens.append(g)
    return gens

def _normal(G):
    r"""
    Return the transformation to normal form.

    INPUT:

    - ``G`` -- `p`-adic symmetric matrix.

    OUTPUT:

    - ``D`` -- the normal form
    - ``B`` -- a transformation matrix

    EXAMPLES::

        sage:
    """
    from sage.quadratic_forms.genera.normal_form import (_jordan_odd_adic,
                                                         _normalize,
                                                         _jordan_2_adic,
                                                         _two_adic_normal_forms)
    p = G.base_ring().prime()
    if p == 2:
        D1, B1 = _jordan_2_adic(G)
        D2, B2 = _normalize(D1)
        D3, B3 = _two_adic_normal_forms(D2)
        B = B3 * B2 * B1
        return D3, B
    else:
        D1, B1 = _jordan_odd_adic(G)
        D2, B2 = _normalize(D1)
        B = B2 * B1
        return D2, B

def _orthogonal_grp_gens_odd(G):
    r"""

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _orthogonal_grp_gens_odd
        sage: R = Zp(3,type='fixed-mod', prec=2, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R,[])
        sage: _orthogonal_grp_gens_odd(G)
        []
        sage: G = matrix.diagonal(R,[1])
        sage: _orthogonal_grp_gens_odd(G)
        [[-1]]
        sage: G = matrix.diagonal(R,[1,1])
        sage: _orthogonal_grp_gens_odd(G)
        [
        [ 3 -2]  [  6  -2]
        [-1  0], [  1 -12]
        ]
        sage: G = matrix.diagonal(R,[1,2])
        sage: _orthogonal_grp_gens_odd(G)
        [
        [2 0]  [-1  0]
        [0 2], [ 0  1]
        ]
    """
    from sage.quadratic_forms.genera.normal_form import p_adic_normal_form
    from sage.quadratic_forms.genera.normal_form import _min_nonsquare
    from sage.arith.misc import legendre_symbol
    R = G.base_ring()
    p = R.prime()
    ug = G.det()
    if ug.valuation() != 0:
        raise ValueError("G is not of scale 1.")
    r = G.ncols()
    if r == 0:
        return []
    if r == 1:
        return [matrix(R,[-1])]
    O = GO(r, p, e=1)
    uo = O.invariant_bilinear_form().det()
    if legendre_symbol(uo, p) != legendre_symbol(ug,p):
        if r % 2 == 0:
            # there are two inequivalent orthogonal_groups
            O = GO(r, p, e=-1)
        else:
            # There is only a single orthogonal group up to isomorphism since
            # O(G) = O(cG). But G and cG are not isomorphic if c is not a square
            c = ZZ(_min_nonsquare(p))
            G = c * G # destroys the normal form of G
    b = O.invariant_bilinear_form()
    b = b.change_ring(R)
    # compute an isomorphism of b and G
    bn, Ub = _normal(b)
    Dn, Ud = _normal(G)
    U = Ud.inverse() * Ub
    assert bn == Dn
    Uinv = U.adjoint()*U.det().inverse_of_unit()
    gens = [U * g.matrix().change_ring(R) * Uinv for g in O.gens()]
    for g in gens:
        err = g*G*g.T - G
        assert _min_val(err) >= 1
    return gens

def _orthogonal_gens_bilinear(G):
    r"""

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _orthogonal_gens_bilinear
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: _orthogonal_gens_bilinear(U)
        sage: _orthogonal_gens_bilinear(V)
        sage: _orthogonal_gens_bilinear(W0)
        sage: _orthogonal_gens_bilinear(matrix.block_diagonal([V,W1])
    """

    r = G.ncols()
    R = G.base_ring()
    def gens_normal_form(O):
        gens = []
        b = O.invariant_form().change_ring(R)
        _, U = _normal(b)
        gens = [g.matrix().change_ring(R) for g in O.gens()]
        Uinv = U.change_ring(GF(2)).inverse().change_ring(R)
        gens = [U * g * Uinv for g in gens]
        return gens
    # corner cases
    if r <= 1:
        gens = []
    elif r == 2 and mod(G[-1,-1],2) == 1:
        return [matrix(R,2,[0,1,1,0])]
    # odd cases
    elif r % 2 == 1:
        O = Sp(r - 1, 2)
        gens = gens_normal_form(O)
        E1 = matrix.identity(R,1)
        gens = [matrix.block_diagonal([g,E1]) for g in gens]
    elif G[-1,-1].valuation() == 0:
        O = Sp(r - 2, 2)
        gens = gens_normal_form(O)
        gens = [matrix.block_diagonal(g,matrix.identity(R,2)) for g in gens]
        E = matrix.identity(R,r)
        for a in (R**(r-2)).basis():
            g = copy(E)
            g[:-2,-2] = a
            g[:-2,-1] = a
            g[-2:,:-2] = g[:-2,-2:].T * G[:-2,:-2]
            gens.append(g)
        g = copy(E)
        g[-2:,-2:] = matrix(R,2,[0,1,1,0])
        gens.append(g)
    # even case
    else:
        O = Sp(r, 2)
        sp = O.invariant_form().change_ring(R)
        _, U = _normal(sp)
        Uinv = U.change_ring(GF(2)).inverse().change_ring(R)
        gens = [U * g.matrix().change_ring(R) * Uinv for g in O.gens()]
    # check that generators are isometries
    for g in gens:
        err = g*G*g.T-G
        assert _min_val(err) >= 1, (g.change_ring(GF(2)),err.change_ring(GF(2)))
    return gens

def _orthogonal_grp_quadratic(G):
    r"""

    INPUT:

    - ``G`` -- homogeneous `2`-adic matrix of scale `1` and in normal form representing a
      quadratic form on an $\FF_2$ vector space with values in `\Zmod{4}`.

    OUTPUT:

    - a list of matrices. Generators of the orthogonal group modulo `2`.

    TESTS::

        sage: from sage.groups.fqf_orthogonal.lift import _orthogonal_grp_quadratic
        sage: R = Zp(2,type='fixed-mod', prec=3, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R,2,[0,1,1,0])
        sage: V = matrix(R,2,[2,1,1,2])
        sage: W0 = matrix(R,2,[1,0,0,3])
        sage: W1 = matrix(R,2,[1,0,0,1])
        sage: G = matrix.block_diagonal([W0])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: gen
        []
        sage: G = matrix.block_diagonal([U, W0])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: gen
        [
        [0 1 0 0]  [1 0 1 1]  [1 0 0 0]
        [1 0 0 0]  [0 1 0 0]  [0 1 1 1]
        [0 0 1 0]  [0 7 1 0]  [7 0 1 0]
        [0 0 0 1], [0 1 0 1], [1 0 0 1]
        ]
        sage: G = matrix.block_diagonal([U, U, W0])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: G = matrix.block_diagonal([U, V, W0])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: G = matrix(R,1,[1])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: G = matrix.block_diagonal([W1])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: gen
        [
        [0 1]
        [1 0]
        ]
        sage: G = matrix.block_diagonal([U, W1])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: gen
        [
        [1 0 0 0]  [4 1 0 0]  [1 0 0 0]
        [1 1 3 3]  [1 4 0 0]  [0 1 0 0]
        [5 0 1 0]  [0 0 1 0]  [0 0 0 7]
        [3 0 0 1], [0 0 0 1], [0 0 1 2]
        ]
        sage: G = matrix.block_diagonal([U, U, W1])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: G = matrix.block_diagonal([U, V, W1])
        sage: gen = _orthogonal_grp_quadratic(G)

    TESTS:

        sage: _orthogonal_grp_quadratic(matrix([]))
        []
    """
    r = G.ncols()
    R = G.base_ring()
    # corner cases
    if r == 0:
        return []
    v = _min_val(G)
    G = (G/R.prime()**v).change_ring(R)
    if r <= 1:
        gens = []
    elif r == 2:
        if G[-1,-1].valuation() == 0:
            if mod(G[-1,-1] + G[-2,-2], 4) == 0:
                gens = []
            else:
                gens = [matrix(R, 2, [0, 1, 1, 0])]
        elif G[-1,-1].valuation() == 1:
            gens = [matrix(R, 2, [0, 1, 1, 0]),
                    matrix(R,2,[0, 1, 1, 1])]
        else:
            gens = [matrix(R, 2, [0, 1, 1, 0])]
    # normal cases
    elif r % 2 == 1:
        # the space of points of even square is preserved
        # so is its orthogonal complement -> invariant vector
        if G[-2,-2] == 0:
            e = 1
        else:
            e = -1
        O = GO(r - 1, 2, e)
        b = O.invariant_quadratic_form().change_ring(R)
        b = b + b.T
        _, U = _normal(b)
        gens = [g.matrix().change_ring(R) for g in O.gens()]
        Uinv = U.adjoint()
        gens = [U * g * Uinv for g in gens]
        E1 = matrix.identity(R,1)
        gens = [matrix.block_diagonal([g,E1]) for g in gens]
    # now r % 2 == 0
    elif G[-1,-1].valuation() == 0:
        # odd case
        if mod(G[-1,-1] + G[-2,-2], 4) == 0:
            gens = _orthogonal_grp_quadratic(G[:-2,:-2])
            gens = [_lift(G, g, 0) for g in gens]
            Id = matrix.identity(R, r - 2)
            gens += [_lift(G, Id, a) for a in (R**(r-2)).basis()]
        else:
            assert mod(G[-1,-1] + G[-2,-2], 4) == 2
            O = Sp(r - 2, 2)
            sp = O.invariant_form().change_ring(R)
            _, U = _normal(sp)
            Uinv = U.adjoint()*U.det().inverse_of_unit()
            gens = [U * g.matrix().change_ring(R) * Uinv for g in O.gens()]
            gens = [_lift(G, g, 0) for g in gens]
            gens += [_lift(G, matrix.identity(R,r-2), 1)]
    else:
        # even case
        if G[-1,-1] == 0:
            e = 1
        else:
            e = -1
        O = GO(r, 2, e)
        b = O.invariant_quadratic_form().change_ring(R)
        b = b + b.T
        _, U = _normal(b)
        Uinv = U.adjoint()*U.det().inverse_of_unit()
        gens = [U * g.matrix().change_ring(R) * Uinv for g in O.gens()]
    # check that generators are isometries
    for g in gens:
        err = g*G*g.T-G
        assert _min_val(err) >= 1, err.change_ring(IntegerModRing(4))
        assert _min_val(matrix.diagonal(err.diagonal())) >= 2, err.change_ring(IntegerModRing(4))
    return gens

def _lift(q, f, a):
    r"""

    INPUT:

    - ``q`` -- of scale `1` in homogeneous normal form
    - ``f`` -- the `n-2 \times n-2` matrix to be lifted
    - ``a`` -- ``0`` or ``1`` there are two possible lifts

    OUTPUT:

    - ``g`` -- the lift of ``f`` as determined by ``a``

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _lift
        sage: R = Zp(2,type='fixed-mod',prec=2,print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R,2,[0,1,1,0])
        sage: W0 = matrix(R,2,[1,0,0,3])
        sage: W1 = matrix(R,2,[1,0,0,1])
        sage: q0 = matrix.block_diagonal([U,W0])
        sage: g0 = matrix(R,2,[0,1,1,0])
        sage: g0l = _lift(q0,g0,vector([1,1]))
        sage: g0l
        [0 1 1 1]
        [1 0 1 1]
        [3 3 0 3]
        [1 1 1 2]

    The essential property of the lifts is that is preserves the bilinear form
    `\mod 2` and the quadratic `\mod 4`::

        sage: (g0l.T*q0*g0l - q0)
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 0]

    The parameter ``a`` is ``g0l1[-2,:-2]``::

        sage: _lift(q0,g0,vector([0,1]))
        [0 1 0 0]
        [1 0 1 1]
        [0 3 1 0]
        [0 1 0 1]

    In the second case one can lift any form preserving the bilinear form on the
    small part. This is the whole symplectic group::

        sage: q1 = matrix.block_diagonal([U,W1])
        sage: g1 = matrix(R,2,[1,1,1,0])
        sage: g1l = _lift(q1,g1,1)
        sage: (g1l.T*q1*g1l - q1)
        [0 0 2 0]
        [0 0 0 0]
        [2 0 0 2]
        [0 0 2 0]
    """
    if mod(q[-1,-2], 2) != 0:
        raise ValueError("The form must be odd.")
    # notation
    g = matrix.block_diagonal([f, matrix.identity(2)])
    b = g.parent()(1)
    b[-2,-1] = 1
    qb = b * q * b.T
    G = qb[:-2,:-2]
    fG = f * G
    fGinv = fG.adjoint() * fG.det().inverse_of_unit()

    if mod(q[-2,-2] + q[-1,-1], 4) == 2:
        g[-1, -2] = a
        g[:-2,-2] = vector(((G - f*G*f.T)/2).diagonal())
        g[-1,:-2] = (fGinv * g[:-2,-2]).T
    else:
        g[:-2,-2] = a
        g[-1,:-2] = (fGinv * g[:-2,-2]).T
        g[-1,-2] = (g[-1,:-2]*G*g[-1,:-2].T)[0,0].expansion(1)
    err = g*qb*g.T-qb
    # check that lifting succeeded
    assert _min_val(err) >= 1
    assert _min_val(matrix.diagonal(err.diagonal())) >= 2
    binv = b.adjoint() * b.det().inverse_of_unit()
    return binv * g * b

def _gens_mod_p(G):
    r"""
    Return the generators of the orthogonal groups of ``G`` modulo `p`.

    Let `V = \Zp^n` and `b: V \times V \rightarrow \Zp` be the bilinear form
    `b(x,y)= x^T G y`. This method computes generators of the image of
    the orthogonal group `O(V,b)` under

    ..MATH:

        O(V,b) \rightarrow GL(V/pV)

    INPUT::

        -``G`` -- gram matrix of a non-degenerate, symmetric, bilinear
          `p`-adic form.

    OUTPUT::

        - generators modulo `p`

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _gens_mod_p
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: R3 = Zp(3, type='fixed-mod', prec=6, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.block_diagonal([U,2*V])
        sage: _gens_mod_p(G)
        [
        [0 1 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]
        [1 0 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]
        [0 0 1 0]  [0 0 0 1]  [0 0 0 1]  [1 0 1 0]  [0 1 1 0]  [0 0 1 0]
        [0 0 0 1], [0 0 1 0], [0 0 1 1], [0 0 0 1], [0 0 0 1], [1 0 0 1],
        <BLANKLINE>
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 1 0 1]
        ]
        sage: G = matrix.block_diagonal([U,2*W1])
        sage: _gens_mod_p(G)
        [
        [0 1 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]
        [1 0 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]
        [0 0 1 0]  [0 0 0 1]  [1 0 1 0]  [0 1 1 0]  [0 0 1 0]  [0 0 1 0]
        [0 0 0 1], [0 0 1 0], [0 0 0 1], [0 0 0 1], [1 0 0 1], [0 1 0 1]
        ]

    TESTS::

        sage: R = Zp(3, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R,[])
        sage: gens = _gens_mod_p(G)
        sage: G = matrix.diagonal(R,[3*1])
        sage: gens = _gens_mod_p(G)
        sage: G = matrix.diagonal(R,[3*1,3*1])
        sage: gens = _gens_mod_p(G)
        sage: G = matrix.diagonal(R,[1,3*1,3*1,9,2*27])
        sage: gens = _gens_mod_p(G)
        """
    n = G.ncols()
    R = G.base_ring()
    E = G.parent().one()
    p = R.prime()
    indices, valuations = _block_indices_vals(G)
    indices.append(n)
    gens = []
    for k in range(len(indices)-1):
        i1 = indices[k]
        i2 = indices[k+1]
        Gi = (G[i1:i2,i1:i2]/ZZ(p)**valuations[k]).change_ring(R)
        gens_homog = _orthogonal_grp_gens_odd(Gi)
        for f in gens_homog:
            g = copy(E)
            g[i1:i2, i1:i2] = f
            gens.append(g)
    # generators below the block diagonal.
    for i in range(n):
        for j in range(i):
            g = copy(E)
            g[i,j] = 1
            flag = True
            for k in range(len(indices)-1):
                if indices[k] <= i < indices[k+1] and indices[k] <= j < indices[k+1]:
                    g[i,j] = 0
                    flag = False
                    break
            if flag:
                gens.append(g)
    return gens

def _gens_mod_2(G):
    r"""
    Return the generators of the orthogonal groups of ``G`` modulo `2`.

    Let `V = \FF_2^n` and `b: V \times V \rightarrow \FF_2` be the bilinear form
    `b(x,y)= x^T G y`. Compute generators of `O(V,b)`.
    INPUT::

    -``G`` -- gram matrix of a non-degenerate, symmetric, bilinear
     `2` form over `\FF_2` in normal form.

    OUTPUT::

    - a list of matrices -- the generators

    EXAMPLES::

        sage: R = Zp(2,type='fixed-mod',print_mode='terse',show_prec=False,prec=6)
        sage: U = matrix(R,2,[0,1,1,0])
        sage: V = matrix(R,2,[2,1,1,2])
        sage: W0 = matrix(R,2,[1,0,0,3])
        sage: W1 = matrix(R,1,[1])
        sage: W2 = matrix(R,2,[1,0,0,1])
        sage: G = matrix.block_diagonal([4*U,4*W0,2*U,V,W1])
        sage: _gens_mod_2(G)
    """
    n = G.ncols()
    R = G.base_ring()
    E = G.parent().one()
    p = G.base_ring().prime()
    ind0, val0 = _block_indices_vals(G)
    par0 = []
    ind0.append(n)
    for k in range(0,len(ind0)-1):
        i = ind0[k+1] - 1    # last index of block i
        pa = mod(G[i,i]//ZZ(2)**val0[k], 2)
        par0.append(pa)

    ind = []
    val = []
    par = []
    k = 0
    for v in range(val0[0]+2, val0[-1]-2,-1):
        try:
            k = val0.index(v)
            ind.append((ind0[k],ind0[k+1]))
            val.append(v)
            par.append(par0[k])
        except ValueError:
            ind.append((ind0[k+1],ind0[k+1]))
            val.append(v)
            par.append(0)
    val[-1] = 0

    gens = []
    for k in range(2,len(val)-1):
        if par[k+1] == 1:
            i1 = ind[k][0]
            i2 = ind[k][1]
            i3 = ind[k+1][1]
            Gk = (G[i1:i3,i1:i3]/ZZ(2)**val[k+1]).change_ring(R)
            gens_k = _gens_pair(Gk, i2-i1, on_second=False)
        elif par[k-1] == 1:
            i1 = ind[k-1][0]
            i2 = ind[k][0]
            i3 = ind[k][1]
            Gk = (G[i1:i3,i1:i3]/ZZ(2)**val[k]).change_ring(R)
            gens_k = _gens_pair(Gk, i2-i1, on_second=True)
        else:
            i1 = ind[k][0]
            i3 = ind[k][1]
            Gk = (G[i1:i3,i1:i3]/ZZ(2)**val[k]).change_ring(R)
            gens_k = _orthogonal_grp_quadratic(Gk)

        for h in gens_k:
            g = copy(E)
            g[i1:i3,i1:i3] = h
            gens.append(g)

    # a change in convention
    trafo = copy(E)
    for k in range(2, len(ind)-1):
        if par[k] == 1 and mod(ind[k][1]-ind[k][0],2) == 0:
           i = ind[k][1]
           trafo[i-2:i,i-2:i] = matrix(R,2,[1,1,0,1])
    trafoinv = trafo.inverse().change_ring(R)
    Gt = trafo * G * trafo.T

    # ker
    # row wise starting with the last row
    for k in range(len(ind)-2,2,-1):
        pa = par[k-2:k+1]
        i = ind[k][1]
        Gi = (Gt[:i,:i]/ZZ(2)**val[k]).change_ring(R)
        gensK = _ker_gens(Gi, ind[k-1][0], ind[k-1][1], pa)
        E = matrix.identity(R,n-i)
        gensK = [matrix.block_diagonal([g,E]) for g in gensK]
        gensK = [trafoinv * g * trafo for g in gensK]
        gens += gensK
    return gens

def _gens_pair(G, k, on_second):
    r"""
    """
    gen = []
    R = G.base_ring()
    n = G.ncols()
    G1 = G[:k,:k]  # 2^1 - modular
    G2 = G[k:,k:]  # 2^0 - modular
    E = G.parent()(1)
    if on_second:
        for f in _orthogonal_gens_bilinear(G2):
            a = vector(((f*G2*f.T-G2)/2).diagonal())
            g = copy(E)
            g[k:,k:] = f
            g[k:,k-1] = a
            gen.append(g)
    else:
        for f in _orthogonal_gens_bilinear((G1/2).change_ring(R)):
            a = vector((f*G1*f.T-G1).diagonal())/2
            g = copy(E)
            g[:k,:k] = f
            g[:k,-1] = a
            g[k:,:k] = - G2 * g[:k,k:].T * (G1*f.T).inverse()
            gen.append(g)
    return gen

def _ker_gens(G, i1, i2, parity):
    r"""

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _ker_gens
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [0, 1, 1, 1])
        sage: W1 = matrix(R, 1, [1])
        sage: W2 = matrix(R, 2, [2, 1, 1, 1])
        sage: G = matrix.block_diagonal([V*4,U*2,U])
        sage: gen = _ker_gens(G,2,4,[0,0,0])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([V*4,W1*2,W2])
        sage: gen = _ker_gens(G,2,3,[0,1,1])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([V*4,U*2,W1])
        sage: gen = _ker_gens(G,2,4,[0,0,1])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([V*4,W2*2,U])
        sage: gen = _ker_gens(G,2,4,[0,1,0])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([W1*4,U*2,U])
        sage: gen = _ker_gens(G,1,3,[1,0,0])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([W1*4,U*2,W0])
        sage: gen = _ker_gens(G,1,3,[1,0,1])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([W0*4,W0*2,U])
        sage: gen = _ker_gens(G,2,4,[1,1,0])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
        sage: G = matrix.block_diagonal([W1*4,U*2,W2*2,W1])
        sage: gen = _ker_gens(G,2,6,[1,1,1])
        sage: gen = [Hensel_qf(G,g,1,2) for g in gen]
        sage: len(gen)
    """
    n = G.nrows()
    E = G.parent()(1)
    gens = []
    e = n - 1
    if parity[2]==1 and mod(n - i2, 2)==0:
        e = n-2
    for i in range(i2,n):
        for j in range(i2):
            g = copy(E)
            if parity == [0,0,0] or parity == [1,0,0]:
                g[i,j] = 1
            elif parity == [0,0,1]:
                if (j < i1) or (i != e):
                    g[i,j] = 1
            elif parity == [0,1,0] or parity == [1,1,0]:
                if not (j == i2 - 1):
                    g[i,j] = 1
            elif parity == [0,1,1]:
                if not ((j == i2 - 1) or (i == e and j >= i1)):
                    g[i,j] = 1
            elif parity == [1,0,1]:
                if not (i == e and j == i1 - 1):
                    g[i,j] = 1
            elif parity == [1,1,1]:
                if not ((i == e and j == i1 - 1) or (j == i2 - 1)):
                    g[i,j] = 1
            if parity[0]==1 and parity[2]==1:
                # compensate
                G[e, i1-1] = (g[e,i1:i2]*G[i1:i2,i1:i2]*g[e,i1:i2].T)[0,0]//4
                # the second row depends on the third
                g[i1:i2,i1-1] = - (G[i1:i2,i1:i2]* g[i2:,i1:i2].T * G[i2:,i2:].inverse())[:,-1]/2
            if g[i,j] == 1:   # no need to append the identity
                gens.append(g)
    return gens

def _gens(G, b):
    r"""
    Return generators.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _gens
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: G = matrix.block_diagonal([2*U,V])
        sage: gens = _gens(G,2)
        sage: G = matrix.block_diagonal([2*U,W1])
        sage: gens = _gens(G,2)
        sage: G = matrix.block_diagonal([2*V,V])
        sage: gens = _gens(G,2)
        sage: G = matrix.diagonal(R,[2,1])
        sage: gens = _gens(G,2)
        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2,2,4,4])
        sage: aut = A.aut()
        sage: gens = [aut(g) for g in gens]
        sage: Oq = aut.subgroup(gens)
        sage: Oq.order()
        1152

    TESTS::

        sage: R = Zp(3, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R,[3*1])
        sage: gens = _gens(G,1)
        sage: G = matrix.diagonal(R,[3*1,3*1])
        sage: gens = _gens(G,2)
        sage: G = matrix.diagonal(R,[2*27,9,3*1,3*1,1])
        sage: gens = _gens(G,4)
    """
    k = 1
    gens = []
    while k <= b:
        gen = _mod_p_to_a_kernel(G, k)
        gen = [Hensel_qf(G, f, k+1, b+1) for f in gen]
        k *= 2
        gens += gen

    if G.base_ring().prime() == 2:
        gen = _gens_mod_2(G)
    else:
        gen = _gens_mod_p(G)
    gens += [Hensel_qf(G, g, 1, b) for g in gen]
    return gens
