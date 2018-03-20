from sage.rings.infinity import Infinity
from sage.rings.all import GF, ZZ, mod
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

        sage: from sage.groups.torsion_orthogonal.lift import Hensel_qf, _min_val
        sage: R = Qp(3, type='floating-point', prec=7, print_mode='terse')
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

        sage: R = Qp(2,type='floating-point', prec=10, print_mode='terse')
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
        9
    """
    # Input checks
    if F.determinant().valuation() != 0:
        raise ValueError("F must be invertible")
    if not (G.parent() == F.parent()):
        raise ValueError("G, F must have the same parent")
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
    F = _Hensel_qf(G, G, F, a, b)
    return F

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

        sage: from sage.groups.torsion_orthogonal.lift import _last_block_index
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

        sage: from sage.groups.torsion_orthogonal.lift import _block_indices_vals
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

        sage: from sage.groups.torsion_orthogonal.lift import _min_val
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
    F[i:,i:] = _Hensel_qf_modular(Zn, G[i:,i:], F[i:,i:], a, b)
    if i == 0:
        return F
    while a < b:
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

        sage: from sage.groups.torsion_orthogonal.lift import _Hensel_qf_modular_odd,_min_val
        sage: R = Qp(3, type='floating-point', print_mode='terse')
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

        sage: from sage.groups.torsion_orthogonal.lift import _Hensel_qf_modular_even, _min_val
        sage: R = Qp(2, type='floating-point', prec=10, print_mode='terse')
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
        [1021    1    1    1]
        [ 925    0 1210  597]
        [ 761    0  297  198]
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
        (5, 5)
        sage: Z = G + matrix.identity(R, 4)*2^2
        sage: Fl = _Hensel_qf_modular_even(Z, G, F, 1, 5)
        sage: Er = Fl*G*Fl.T - Z
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (5, 5)

    TESTS::

        sage: G = matrix(R, 3, [2,1,0, 1,2,0, 0,0,7])
        sage: F = matrix(R, 3, [0, 1, 0, 1, 1, 0, 0, 0, 1])
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 8)
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(matrix(Er.diagonal()))
        (9, 9)
    """
    R = Z.base_ring()
    n = Z.ncols()
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
    assert _min_val(Z-F*G*F.T)>=b # sanity check
    return F

def _solve_X(Y, b, g):
    r"""
    Solve a certain linear equation.

    This is a helper function for :meth:`_Hensel_qf_modular_even`.

    :MATH::

        Y = X + X^T

    :MATH::

        b_i = X_{ii} + \sum_{k=1}^n X_{ij}g_j \quad i \in \{1, \dots, n\}

    INPUT:

    - ``Y`` - converts to an `n \times n` matrix over `\FF_2`
    - ``b`` - converts to `\FF_2^n`
    - ``g`` - converts to `\FF_2^n`

    OUTPUT:

    - ``X`` - an `n \times n` matrix over `\FF_2` the solution of the
      linear equation above.

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _solve_X
        sage: k = GF(2)
        sage: Y = matrix(k, 3, [0,0,1, 1,0,1, 1,0,0])
        sage: b = vector(k, [1, 0, 0])
        sage: g = vector(k, [0, 1, 0])
        sage: _solve_X(Y, b, g)
        [1 0 0]
        [0 0 0]
        [0 0 0]
    """
    k = GF(2)
    Y = matrix(k, Y)
    b = vector(k, b)
    g = vector(k, g)
    n = Y.ncols()

    equations = []
    for i in range(0, n):
        R = matrix.zero(k, n, n)
        R[i,:] = g
        R[i,i] += 1
        eqn = R.list() + [b[i]]
        equations .append(eqn)
    # equations to be a symmetric matrix
    for i in range(0, n):
        for j in range(i+1, n):
            R = matrix.zero(n, n)
            R[i, j] = 1
            R[j, i] = 1
            eq = R.list() + [Y[i,j]]
            equations.append(eqn)
    A = matrix(k, equations)
    b = A[:,-1]
    A = A[:,:-1]
    # A*Xcoeff == b
    Xcoeff = A.solve_right(b)
    X = matrix(k, n, n, Xcoeff.list())
    return X

def _basis_kernel(R, N, s=0):
    r"""

    INPUT:

    - ``R`` -- a space of `n \times n` matrices
    - ``N`` -- an ascending list of integers ``0==N[0] < ... < N[-1]==n``
    - ``s`` -- an integer (default: ``0``)

    OUTPUT:

    - a list of matrices

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _basis_kernel
        sage: R = MatrixSpace(ZZ, 4, 4)
        sage: N = [0, 2, 4]
        sage: _basis_kernel(R, N)
        [
        [0 1 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]
        [1 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]
        [0 0 0 0]  [1 0 0 0]  [0 1 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 1]
        [0 0 0 0], [0 0 0 0], [0 0 0 0], [1 0 0 0], [0 1 0 0], [0 0 1 0]
        ]
        sage: _basis_kernel(R, N, 2)
        [
        [0 0 0 0]
        [0 0 0 0]
        [0 0 0 1]
        [0 0 1 0]
        ]
    """
    basis = []
    n = R.ncols()
    assert R.ncols() == R.nrows()
    for i in range(n):
        for j in range(s, i):
            b = copy(R.zero())
            for k in range(len(N)-1):
                if N[k] <= i < N[k+1] and N[k] <= j < N[k+1]:
                    b[j,i] = 1
                    break
            b[i,j] = 1
            basis.append(b)
    return basis

def _mod_p_kernel(G, b, lift=True):
    r"""
    Return generators of the kernel of O(A/p^bA) --> O(A/pA).

    INPUT:

    - ``G`` - a symmetric `p`-adic matrix
    - ``b`` - an integer

    OUTPUT:

    - a list of matrices

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _mod_p_kernel
        sage: R = Zp(3, type='fixed-mod', prec=4, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R, [3^2, 3^2, 1, 1])
        sage: ker_gens = _mod_p_kernel(G, 4, lift=False)
        sage: ker_gens
        [
        [1 3 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]
        [3 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]
        [0 0 1 0]  [3 0 1 0]  [0 3 1 0]  [0 0 1 0]  [0 0 1 0]  [0 0 1 3]
        [0 0 0 1], [0 0 0 1], [0 0 0 1], [3 0 0 1], [0 3 0 1], [0 0 3 1],
        <BLANKLINE>
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 9]
        [0 0 9 1]
        ]
    """
    n = G.ncols()
    R = G.parent()
    E = R.one()
    p = G.base_ring().prime()

    indices, valuations = _block_indices_vals(G)
    assert b >= valuations[0]
    valuations = [b - v for v in valuations]
    indices.append(n)
    k = 1
    gens = []
    while k < b:
        for i in range(len(valuations)):
            if k < valuations[i]:
                break
        s1 = indices[i]
        basis = _basis_kernel(R, indices, s1)
        gen = [E + p**k*F for F in basis]
        if lift:
            gen = [Hensel_qf(G, f, k, b) for f in gens]
        gens += gen
        k *= 2
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

def _gens_homogeneous(G):
    r"""
    Return generators of the orthogonal group modulo `p`

    INPUT:

    - ``G`` -- homogeneous `p`-adic matrix in normal form representing a
      quadratic form on an $\FF_p$ vector space with values in `\Zmod{p}`.

    OUTPUT:

    - a list of matrices. Generators of the orthogonal group modulo `p`.

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _gens_homogeneous
        sage: R = Zp(2,type='fixed-mod', prec=4, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R,2,[0,1,1,0])
        sage: V = matrix(R,2,[2,1,1,2])
        sage: W0 = matrix(R,2,[1,0,0,3])
        sage: W1 = matrix(R,2,[1,0,0,1])
        sage: G = matrix.block_diagonal([W0])
        sage: gen = _gens_homogeneous(G)
        sage: len(gen)
        0
        sage: G = 2 * matrix.block_diagonal([U, W0])
        sage: gen = _gens_homogeneous(G)
        sage: len(gen)
        2
        sage: G = 2 * matrix.block_diagonal([U, U, W0])
        sage: gen = _gens_homogeneous(G)
        sage: len(gen)
        6
        sage: G = 2 * matrix.block_diagonal([U, V, W0])
        sage: gen = _gens_homogeneous(G)
        sage: len(gen)
        6

    Some examples at odd primes::

        sage: R = Zp(3,type='fixed-mod', prec=6, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R,[1,3,3,9])
    """
    R = G.base_ring()
    p = R.prime()
    v = _min_val(G)
    Ghomogeneous = (G/p**v).change_ring(R)
    if p == 2:
        return _gens_homogeneous_even(Ghomogeneous)
    else:
        return _gens_homogeneous_odd(Ghomogeneous)

def _gens_homogeneous_odd(G):
    r"""

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _gens_homogeneous_odd
        sage: R = Zp(3,type='fixed-mod', prec=2, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R,[])
        sage: _gens_homogeneous_odd(G)
        []
        sage: G = matrix.diagonal(R,[1])
        sage: _gens_homogeneous_odd(G)
        [[-1]]
        sage: G = matrix.diagonal(R,[1,1])
        sage: _gens_homogeneous_odd(G)
        [
        [ 3 -2]  [-3 -2]
        [-1  0], [ 1 -3]
        ]
        sage: G = matrix.diagonal(R,[1,2])
        sage: _gens_homogeneous_odd(G)
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
            G = c * G
    b = O.invariant_bilinear_form()
    b = b.change_ring(R)
    # compute an isomorphism
    bn, U = _normal(b)
    assert bn == G
    Uinv = U.adjoint()*U.det().inverse_of_unit()
    gens = [U * g.matrix().change_ring(R) * Uinv for g in O.gens()]
    for g in gens:
        err = g*G*g.T - G
        assert _min_val(err) >= 1
    return gens

def _gens_homogeneous_even(G):
    r"""

    INPUT:

    - ``G`` -- homogeneous `2`-adic matrix of scale `1` and in normal form representing a
      quadratic form on an $\FF_2$ vector space with values in `\Zmod{4}`.

    OUTPUT:

    - a list of matrices. Generators of the orthogonal group modulo `2`.

    TESTS::

        sage: from sage.groups.torsion_orthogonal.lift import _gens_homogeneous_even
        sage: R = Zp(2,type='fixed-mod', prec=3, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R,2,[0,1,1,0])
        sage: V = matrix(R,2,[2,1,1,2])
        sage: W0 = matrix(R,2,[1,0,0,3])
        sage: W1 = matrix(R,2,[1,0,0,1])
        sage: G = matrix.block_diagonal([W0])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([U, W0])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([U, U, W0])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([U, V, W0])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix(R,1,[1])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([W1])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([U, W1])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([U, U, W1])
        sage: gen = _gens_homogeneous_even(G)
        sage: G = matrix.block_diagonal([U, V, W1])
        sage: gen = _gens_homogeneous_even(G)

    TESTS:

        sage: _gens_homogeneous_even(matrix([]))
        []
    """
    r = G.ncols()
    R = G.base_ring()
    if r <= 1:
        gens = []
    elif r == 2:
        if mod(G[0,0] + G[1,1], 4) == 0:
            gens = []
        else:
            gens = []
    elif r % 2 == 1:
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
        if mod(G[-1,-1] + G[-2,-2], 4) == 0:
            gens = _gens_homogeneous_even(G[:-2,:-2])
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
            gens += [_lift(G, g.parent().one(), 1)]
    else:
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

        sage: from sage.groups.torsion_orthogonal.lift import _lift
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
    b = copy(g.parent().one())
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

        sage: from sage.groups.torsion_orthogonal.lift import _gens_mod_p
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: R3 = Zp(3, type='fixed-mod', prec=6, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.block_diagonal([U,2*V])
        sage: _gens_mod_p(G)
        [
        [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]
        [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]
        [1 0 1 0]  [0 1 1 0]  [0 0 1 0]  [0 0 1 0]
        [0 0 0 1], [0 0 0 1], [1 0 0 1], [0 1 0 1]
        ]
        sage: G = matrix.block_diagonal([U,2*W1])
        sage: _gens_mod_p(G)
        [
        [1 0 0 0]  [1 0 0 0]  [1 0 0 0]  [1 0 0 0]
        [0 1 0 0]  [0 1 0 0]  [0 1 0 0]  [0 1 0 0]
        [1 0 1 0]  [0 1 1 0]  [0 0 1 0]  [0 0 1 0]
        [0 0 0 1], [0 0 0 1], [1 0 0 1], [0 1 0 1]
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
    R = G.parent()
    E = R.one()
    p = G.base_ring().prime()
    indices, valuations = _block_indices_vals(G)
    indices.append(n)
    gens = []
    for k in range(len(indices)-1):
        i1 = indices[k]
        i2 = indices[k+1]
        Gi = G[i1:i2,i1:i2]
        gens_homog = _gens_homogeneous(Gi)
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

def _gens(G, b):
    r"""
    Return generators.

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _gens
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: G = matrix.block_diagonal([2*U,V])
        sage: gens = _gens(G,2)
        sage: G = matrix.block_diagonal([2*U,W1])
        sage: gens = _gens(G,2)
        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2,2,4,4])
        sage: aut = A.aut()
        sage: gens = [aut(g) for g in gens]
        sage: Oq = aut.subgroup(gens)
        sage: Oq.order()
        32

    TESTS::

        sage: R = Zp(3, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R,[3*1])
        sage: gens = _gens(G,1)
        sage: G = matrix.diagonal(R,[3*1,3*1])
        sage: gens = _gens(G,2)
        sage: G = matrix.diagonal(R,[2*27,9,3*1,3*1,1])
        sage: gens = _gens(G,4)
    """
    gensK = _mod_p_kernel(G, b)
    gens = _gens_mod_p(G)
    gens = [Hensel_qf(G, g, 1, b) for g in gens]
    return gens + gensK
