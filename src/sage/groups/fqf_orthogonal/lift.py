r"""
Lifting of quadratic matrix equations over the `p`-adics.

<Paragraph description>

EXAMPLES::

    sage:

AUTHORS:

- Simon Brandhorst (2018-05-15): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.infinity import Infinity
from sage.rings.all import GF
from copy import copy
from sage.matrix.all import matrix
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
    if F.ncols() == 0:
        return F
    # the real worker
    F = copy(F) # leave input unchanged
    F = _Hensel_qf(G, G, F, a, b) #works inplace
    return F

def _reverse_homogeneous_blocks(G):
    r"""
    Reverses the order of the homogeneous blocks of ``G``.

    Input:

    - ``G`` -- a block diagonal matrix over the `p`-adics

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
    Return a list of indices and a list of valuation of the homogeneous blocks.

    Input:

    - ``G`` -- a symmetric `p`-adic block diagonal matrix with modular blocks
      which have descending valuations

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
    Return the minimum valuation of an entry of ``M``

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
    The real worker for :meth:`Hensel_qf`.

    No input checks.

    INPUT:

    - ``Z`` -- symmetric `p`-adic `n \times n` matrix
    - ``G`` -- symmetric `p`-adic `n \times n` matrix
    - ``F`` -- symmetric `p`-adic `n \times n` matrix
    - ``a`` -- integer
    - ``b`` -- integer

    We require that the triple `(Z, G, F)` is `a`-adapted.

    OUTPUT:

    a matrix ``Fl`` such that `(Z, G, Fl)` is `b`-adapted.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal.lift import _Hensel_qf
        sage: R = Zp(3, type='fixed-mod', prec=5, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R, [3^2,1,1])
        sage: Z = G + matrix(R, 3, [0,3^2,0, 3^2,0,0, 0,0,3])
        sage: F = matrix(R,3, [1,0,0, 0,0,1, 0,1,0])
        sage: Z - F*G*F.T
        [0 9 0]
        [9 0 0]
        [0 0 3]
        sage: Fn = _Hensel_qf(Z, G, F, 1, 4)#
        sage: Z - F*G*F.T
        [  0   0   0]
        [  0   0   0]
        [  0   0 -81]
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
    Gn = G[i:,i:]
    F[i:,i:] = _Hensel_qf_modular(Zn, Gn, F[i:,i:], a, b)
    K = (G[i:,i:]*F[i:,i:].T).inverse()
    if i == 0:
        return F
    while a < b:
        # an extra line that recomputes the upper block diagonal
        # if the input really is adapted this should not be necessary
        # but in any case it does not hurt
        F[:i,i:] = (Z[:i,i:] - F[:i,:i]*G[:i,:i]*F[i:,:i].T) * K
        Zn = Z[:i,:i] - F[:i,i:]*G[i:,i:]*F[:i,i:].T
        F[:i,:i] = _Hensel_qf(Zn, G[:i,:i], F[:i,:i], a, a+s)
        F[:i,i:] = (Z[:i,i:] - F[:i,:i]*G[:i,:i]*F[i:,:i].T) * K
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

    Deals with the case that `G` is modular and `p=2`.

    INPUT:

    - ``Z`` -- symmetric `p`-adic `n \times n` matrix
    - ``G`` -- symmetric `p`-adic `n \times n` matrix
    - ``F`` -- symmetric `p`-adic `n \times n` matrix
    - ``a`` -- integer
    - ``b`` -- integer

    We require that the triple `(Z, G, F)` is `a`-adapted.

    OUTPUT:

    - `Fl` such that ``(Z, G, Fl)`` is `b`-adapted
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
    n = Z.ncols()
    if a == 0:
        raise ValueError("a must be a non-zero integer")
    if a == 1:
        R = Z.base_ring()
        v = _min_val(G)
        G = (G/2**v).change_ring(R)
        Z = (Z/2**v).change_ring(R)
        Y = Z - F*G*F.T
        X = _solve_X(Y/2, (Y/4).diagonal(), G.inverse().diagonal())
        X = 2 * X.change_ring(R)
        F = F + X*(G*F.T).inverse()
        a = 2
    while a < b:
        Y = Z - F*G*F.T
        for i in range(0,n):
            Y[i,i+1:] = 0
            Y[i,i] = Y[i,i]//2
        F = F + Y*(G*F.T).inverse()
        a = 2*a - 1
    # confirm computation
    # assert _min_val(Z-F*G*F.T) >= b
    # assert _min_val(matrix((Z-F*G*F.T).diagonal())) >= b + 1
    return F

def _solve_X(Y, b, g, ker=False):
    r"""
    Solve a certain linear equation mod `2`.

    This is a helper function for :meth:`_Hensel_qf_modular_even`.

    :MATH::

        Y = X + X^T

    :MATH::

        b_i = X_{ii} + \sum_{j=1}^n X_{ij}g_j \quad i \in \{1, \dots, n\}

    INPUT:

    - ``Y`` -- converts to an `n \times n` matrix over `\FF_2`
    - ``b`` -- converts to `\FF_2^n`
    - ``g`` -- converts to `\FF_2^n`

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
