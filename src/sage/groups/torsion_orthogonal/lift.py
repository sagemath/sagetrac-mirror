from sage.rings.infinity import Infinity
from sage.rings.all import GF
from copy import copy
from sage.matrix.all import Matrix
from sage.modules.all import vector

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
        sage: G = Matrix(R, 6, 6, [0, 243, 0, 0, 0, 0, 243, 0, 0, 0, 0, 0, 0, 0, 0, 27, 0, 0, 0, 0, 27, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])
        sage: F = Matrix(R, 6, 6, [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 2, 1, 0, 0, 0, 2, -1, 6, 3, 0, 1, 1, 9, 3, 6, 1, 0])
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
        sage: U = Matrix(R, 2, [0, 1, 1 ,0])
        sage: V = Matrix(R, 2, [2, 1, 1 ,2])
        sage: G = Matrix.block_diagonal([2*U, 2*U, V])
        sage: F = Matrix(R, 6, [1, 0, 0, 0, 0, 0,
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
        sage: G = Matrix(Zp(2), 3, [4,0,0, 0,2,1, 0,1,2])
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
        sage: G = Matrix(Zp(2), 3, [4,0,0, 0,2,1, 0,1,2])
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

        sage: G = Matrix(Zp(2),[])
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
        sage: G = Matrix(R, 4, [0,1,0,0, 1,0,0,0, 0,0,1,0, 0,0,0,5])
        sage: G
        [0 1 0 0]
        [1 0 0 0]
        [0 0 1 0]
        [0 0 0 5]
        sage: F = Matrix(R, 4, [1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0])
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
        sage: _min_val(Er), _min_val(Matrix(Er.diagonal()))
        (6, 6)
        sage: F = Matrix(R, 4, [0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1])
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 4)
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(Matrix(Er.diagonal()))
        (+Infinity, +Infinity)

    ::

        sage: G = Matrix(R, 4, [2,1,0,0, 1,2,0,0, 0,0,3,0, 0,0,0,7])
        sage: F = Matrix(R, 4, [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0])
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
        sage: _min_val(Er), _min_val(Matrix(Er.diagonal()))
        (5, 5)
        sage: Z = G + matrix.identity(R, 4)*2^2
        sage: Fl = _Hensel_qf_modular_even(Z, G, F, 1, 5)
        sage: Er = Fl*G*Fl.T - Z
        sage: _min_val(Er), _min_val(Matrix(Er.diagonal()))
        (5, 5)

    TESTS::

        sage: G = Matrix(R, 3, [2,1,0, 1,2,0, 0,0,7])
        sage: F = Matrix(R, 3, [0, 1, 0, 1, 1, 0, 0, 0, 1])
        sage: Fl = _Hensel_qf_modular_even(G, G, F, 1, 8)
        sage: Er = Fl*G*Fl.T - G
        sage: _min_val(Er), _min_val(Matrix(Er.diagonal()))
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
        D = Matrix.diagonal(Y.diagonal())

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
        sage: Y = Matrix(k, 3, [0,0,1, 1,0,1, 1,0,0])
        sage: b = vector(k, [1, 0, 0])
        sage: g = vector(k, [0, 1, 0])
        sage: _solve_X(Y, b, g)
        [1 0 0]
        [0 0 0]
        [0 0 0]
    """
    k = GF(2)
    Y = Matrix(k, Y)
    b = vector(k, b)
    g = vector(k, g)
    n = Y.ncols()

    equations = []
    for i in range(0, n):
        R = Matrix.zero(k, n, n)
        R[i,:] = g
        R[i,i] += 1
        eqn = R.list() + [b[i]]
        equations .append(eqn)
    # equations to be a symmetric matrix
    for i in range(0, n):
        for j in range(i+1, n):
            R = Matrix.zero(n, n)
            R[i, j] = 1
            R[j, i] = 1
            eq = R.list() + [Y[i,j]]
            equations.append(eqn)
    A = Matrix(k, equations)
    b = A[:,-1]
    A = A[:,:-1]
    # A*Xcoeff == b
    Xcoeff = A.solve_right(b)
    X = Matrix(k, n, n, Xcoeff.list())
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

def _mod_p_kernel(G, b):
    r"""
    Return generators of the kernel of O(A/p^bA) --> O(A/pA).

    INPUT:

    - ``G`` - a symmetric `p`-adic matrix
    - ``b`` - an integer

    OUTPUT:

    - a list of matrices

    EXAMPLES::

        sage: from sage.groups.torsion_orthogonal.lift import _mod_p_kernel
        sage: R = Qp(3, type='floating-point', prec=7, print_mode='terse')
        sage: G = Matrix(R, 6, 6, [0, 243, 0, 0, 0, 0, 243, 0, 0, 0, 0, 0, 0, 0, 0, 27, 0, 0, 0, 0, 27, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0])
        sage: F = Matrix(R, 6, 6, [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 1, 2, 1, 0, 0, 0, 2, -1, 6, 3, 0, 1, 1, 9, 3, 6, 1, 0])
        sage: _mod_p_kernel(G, 5)
        [
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]
        [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]
        [0 0 1 3 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]
        [0 0 3 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]
        [0 0 0 0 1 0]  [0 0 3 0 1 0]  [0 0 0 3 1 0]  [0 0 0 0 1 0]
        [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 3 0 0 1],
        <BLANKLINE>
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]
        [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]
        [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 9 0 0]  [0 0 1 0 0 0]
        [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 9 1 0 0]  [0 0 0 1 0 0]
        [0 0 0 0 1 0]  [0 0 0 0 1 3]  [0 0 0 0 1 0]  [0 0 9 0 1 0]
        [0 0 0 3 0 1], [0 0 0 0 3 1], [0 0 0 0 0 1], [0 0 0 0 0 1],
        <BLANKLINE>
        [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]
        [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]
        [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]
        [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]
        [0 0 0 9 1 0]  [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 0 1 9]
        [0 0 0 0 0 1], [0 0 9 0 0 1], [0 0 0 9 0 1], [0 0 0 0 9 1],
        <BLANKLINE>
        [ 1  0  0  0  0  0]
        [ 0  1  0  0  0  0]
        [ 0  0  1  0  0  0]
        [ 0  0  0  1  0  0]
        [ 0  0  0  0  1 81]
        [ 0  0  0  0 81  1]
        ]
    """
    n = G.ncols()
    R = G.parent()
    E = R.one()
    p = G.base_ring().prime()

    indices, valuations = _block_indices_vals(G)
    assert b >= valuations[0]
    indices.append(n)
    k = 1
    gens = []
    while k <= b:
        for i in range(len(valuations)):
            if k <= b - valuations[i]:
                break
        s1 = indices[i]
        basis = _basis_kernel(R, indices, s1)
        gen = [E + p**k*F for F in basis]
        gens += gen
        k *= 2
    return gens