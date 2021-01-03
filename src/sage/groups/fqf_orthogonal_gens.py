r"""
Generators of the orthogonal group of a torsion quadratic form.

<Paragraph description>

EXAMPLES::

<Lots and lots of examples>

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
from sage.rings.all import GF, ZZ, mod, IntegerModRing, Zp
from copy import copy
from sage.matrix.all import matrix
from sage.modules.all import vector
from sage.groups.all import GO, Sp
from sage.groups.fqf_orthogonal_lift import (Hensel_qf,
                                             _block_indices_vals,
                                             _min_val,
                                             _solve_X)
from sage.groups.matrix_gps.linear import GL
from sage.matrix.matrix_space import MatrixSpace
from copy import copy
from sage.groups.fqf_orthogonal import _isom_fqf

def _mod_p_to_a_kernel(G, a):
    r"""
    Return the kernel of the orthogonal group under mod `p` reduction.

    Let `L` be a `p`-adic lattice. For a positive integer `a` define
    `O(L/p^aL)` as the image of `O(L)` in `GL(L/p^aL)`.
    This function returns generators of the kernel of `O(L/p^aL)`--> `O(L/pL)`.

    INPUT:

    - ``G`` -- a non-singular, symmetric, integral, block diagonal
      `p`-adic matrix with modular blocks of descending valuation
    - ``a``-- an integer

    OUTPUT:

    - a list of `p`-adic matrices

    EXAMPLES::

        sage:
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
            gensk = _solve_X(matrix.zero(R,ni,ni), matrix.zero(R,1,ni).list(),Gk_inv.diagonal(), ker=True)
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

def _gens_normal_form(O, R):
    r"""
    Return generators for the orthogonal group `O` with normal invariant form.

    INPUT:

    - ``O`` -- an orthogonal or sympletic group
    - ``R`` -- a `p`-adic ring

    OUTPUT:

    - a list of matrices over `R`

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _gens_normal_form
        sage: R = Zp(3,type='fixed-mod', prec=2, print_mode='terse', show_prec=False, print_pos=False)
        sage: O = GO(4,3, e=1)
        sage: _gens_normal_form(O, R)
        [
        [2 2 2 1]  [2 1 1 2]
        [2 1 1 1]  [1 1 2 2]
        [2 2 1 2]  [2 1 2 1]
        [1 2 1 1], [1 1 1 1]
        ]
    """
    p = R.prime()
    gens = []
    try:
        b = O.invariant_quadratic_form().change_ring(R)
        b = b + b.T
    except AttributeError:
        b = O.invariant_form().change_ring(R)
    _, U = _normal(b)
    U = U.change_ring(GF(p))
    Uinv = U.inverse()
    gens = [U * g.matrix() * Uinv for g in O.gens()]
    gens = [g.change_ring(R) for g in gens]
    return gens

def _orthogonal_grp_gens_odd(G):
    r"""
    Return generators of the orthogonal group modulo `p` for odd `p`.

    `G` is the gram matrix of the bilinear form.

    INPUT:

    - ``G`` -- a unimodular, symmetric, integral `p`-adic matrix; `p` odd

    OUTPUT:

    - a list of `p`-adic matrices well defined modulo `p`

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _orthogonal_grp_gens_odd
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
        [0 1]  [0 1]
        [2 0], [1 0]
        ]
        sage: G = matrix.diagonal(R,[1,2])
        sage: _orthogonal_grp_gens_odd(G)
        [
        [2 0]  [2 0]
        [0 2], [0 1]
        ]

    TESTS:

        sage: _orthogonal_grp_gens_odd(matrix(R,[]))
        []
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
    Ud = G.parent().one()
    if legendre_symbol(uo, p) != legendre_symbol(ug, p):
        if r % 2 == 0:
            # there are two inequivalent orthogonal_groups
            O = GO(r, p, e=-1)
        else:
            # There is only a single orthogonal group up to isomorphism since
            # O(G) = O(cG). But G and cG are not isomorphic if c is not a square
            c = ZZ(_min_nonsquare(p))
            G = c * G # destroys the normal form of G
            _, Ud = _normal(G) # restore it
    gens = _gens_normal_form(O, R)
    Uinv = Ud.adjugate()*Ud.det().inverse_of_unit()
    gens = [Uinv * g * Ud for g in gens]
    # check that generators are indeed isometries
    for g in gens:
        err = g*G*g.T - G
        assert _min_val(err) >= 1
    return gens

def _orthogonal_gens_bilinear(G):
    r"""
    Return generators of the orthogonal group of the bilinear form modulo `2`.

    Here `G` is the gram matrix of the bilinear form.

    INPUT:

    - ``G`` -- a modular, symmetric `2`-adic matrix in homogeneous normal form

    OUTPUT:

    - a list of `2`-adic matrices well defined modulo `2`

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _orthogonal_gens_bilinear
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: _orthogonal_gens_bilinear(U) == _orthogonal_gens_bilinear(V)
        True
        sage: _orthogonal_gens_bilinear(W0)
        [
        [0 1]
        [1 0]
        ]
        sage: _orthogonal_gens_bilinear(matrix.block_diagonal([V,W1]))
        [
        [1 0|0 0]  [0 1|0 0]
        [1 1|0 0]  [1 0|0 0]  [1 0 1 1]  [1 0 0 0]  [1 0 0 0]
        [---+---]  [---+---]  [0 1 0 0]  [0 1 1 1]  [0 1 0 0]
        [0 0|1 0]  [0 0|1 0]  [2 1 1 0]  [1 2 1 0]  [0 0 0 1]
        [0 0|0 1], [0 0|0 1], [2 1 0 1], [1 2 0 1], [0 0 1 0]
        ]

    TESTS:

        sage: _orthogonal_gens_bilinear(matrix([]))
        []
    """
    r = G.ncols()
    R = G.base_ring()
    # corner cases
    if r <= 1:
        gens = []
    elif r == 2 and mod(G[-1,-1],2) == 1:
        return [matrix(R,2,[0,1,1,0])]
    # odd cases
    elif r % 2 == 1:
        # the space of points of 0 square is non degenerate and preserved
        # there the bilinear form is alternating and the group is a
        # sympletic group
        # the orthogonal complement has rank one and is invariant as well
        O = Sp(r - 1, 2)
        gens = _gens_normal_form(O, R)
        E1 = matrix.identity(R,1)
        gens = [matrix.block_diagonal([g, E1]) for g in gens]
    elif G[-1,-1].valuation() == 0:
        # the space of points of 0 square is degenerate with kernel of dim 1
        # the quotient by the kernel is again sympletic
        # but we obtain several possible lifts
        O = Sp(r - 2, 2)
        gens = _gens_normal_form(O, R)
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
    # just a symplectic group nothing special
    else:
        O = Sp(r, 2)
        gens = _gens_normal_form(O, R)
    # check that generators are isometries
    for g in gens:
        err = g*G*g.T-G
        assert _min_val(err) >= 1, (g.change_ring(GF(2)),err.change_ring(GF(2)))
    return gens

def _orthogonal_grp_quadratic(G):
    r"""
    Return generators of the orthogonal group of the quadratic form modulo `4`.

    Here `G`represents a quadratic form on an $\FF_2$ vector space with
    values in `\Zmod{4}` given by `q(x) = xGx^T \mod 4`.

    INPUT:

    - ``G`` -- a homogeneous, symmetric `2`-adic matrix in normal form

    OUTPUT:

    - a list of matrices. Generators of the orthogonal group modulo `2`.

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _orthogonal_grp_quadratic
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
        [1 0 0 0]  [0 1 0 0]  [1 0 0 0]
        [1 1 3 3]  [1 0 0 0]  [0 1 0 0]
        [5 0 1 0]  [0 0 1 0]  [0 0 0 7]
        [3 0 0 1], [0 0 0 1], [0 0 1 2]
        ]
        sage: G = matrix.block_diagonal([U, U, W1])
        sage: gen = _orthogonal_grp_quadratic(G)
        sage: G = matrix.block_diagonal([U, V, W1])
        sage: gen = _orthogonal_grp_quadratic(G)

    TESTS::

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
    # usual cases
    elif r % 2 == 1:
        # an odd case
        # the space of points of even square is preserved
        # there the quadratic form is classical
        # so is its orthogonal complement -> we get an invariant vector
        if G[-2,-2] == 0:
            e = 1
        else:
            e = -1
        O = GO(r - 1, 2, e)
        gens = _gens_normal_form(O, R)
        # the invariant vector is the last one
        # continue the isometry as the identity there
        E1 = matrix.identity(R,1)
        gens = [matrix.block_diagonal([g,E1]) for g in gens]
    # now r % 2 == 0
    elif G[-1,-1].valuation() == 0:
        # odd case
        # the space of points of even square is preserved
        # but it is degenerate
        # modulo the degenerate part the form is
        if mod(G[-1,-1] + G[-2,-2], 4) == 0:
            # the degenerate part v has q(v) = 0
            # thus q is preserved an we get
            # a classical form
            gens = _orthogonal_grp_quadratic(G[:-2,:-2])
            gens = [_lift(G, g, 0) for g in gens]
            Id = matrix.identity(R, r - 2)
            gens += [_lift(G, Id, a) for a in (R**(r-2)).basis()]
        else:
            # here the degenerate part v has q(v) = 2
            # thus only the bilinear form descends to the quotient
            # and we get a sympletic group
            assert mod(G[-1,-1] + G[-2,-2], 4) == 2
            O = Sp(r - 2, 2)
            gens = _gens_normal_form(O, R)
            gens = [_lift(G, g, 0) for g in gens]
            gens += [_lift(G, matrix.identity(R,r-2), 1)]
    else:
        # even case
        # the groups are classical we just need to bring them to the right basis
        if G[-1,-1] == 0:
            e = 1
        else:
            e = -1
        O = GO(r, 2, e)
        gens = _gens_normal_form(O, R)
    # check that generators are isometries
    for g in gens:
        err = g*G*g.T-G
        assert _min_val(err) >= 1, err.change_ring(IntegerModRing(4))
        assert _min_val(matrix.diagonal(err.diagonal())) >= 2, err.change_ring(IntegerModRing(4))
    return gens

def _lift(q, f, a):
    r"""
    Return a list of the lifts of f.

    INPUT:

    - ``q`` -- odd of scale `1` in homogeneous normal form
    - ``f`` -- the `n-2 \times n-2` matrix to be lifted
    - ``a`` -- ``0``, ``1`` or a row matrix depending on the case

    OUTPUT:

    - ``g`` -- the lift of ``f`` as determined by ``a``

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _lift
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
    fGinv = fG.adjugate() * fG.det().inverse_of_unit()

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
    binv = b.adjugate() * b.det().inverse_of_unit()
    return binv * g * b

def _gens_mod_p(G):
    r"""
    Return generators of the orthogonal groups of ``G`` modulo `p`.

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

        sage: from sage.groups.fqf_orthogonal_gens import _gens_mod_p
        sage: R = Zp(3, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: G = matrix.diagonal(R, [3*1, 3*1])
        sage: _gens_mod_p(G)
        [
        [0 1]  [0 1]
        [2 0], [1 0]
        ]
        sage: G = matrix.diagonal(R, [1, 3, 3, 9, 2*27])
        sage: _gens_mod_p(G)
        [
        [-1  0  0  0  0]  [1 0 0 0 0]  [1 0 0 0 0]  [ 1  0  0  0  0]
        [ 0  1  0  0  0]  [0 0 1 0 0]  [0 0 1 0 0]  [ 0  1  0  0  0]
        [ 0  0  1  0  0]  [0 2 0 0 0]  [0 1 0 0 0]  [ 0  0  1  0  0]
        [ 0  0  0  1  0]  [0 0 0 1 0]  [0 0 0 1 0]  [ 0  0  0 -1  0]
        [ 0  0  0  0  1], [0 0 0 0 1], [0 0 0 0 1], [ 0  0  0  0  1],
        <BLANKLINE>
        [ 1  0  0  0  0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]
        [ 0  1  0  0  0]  [1 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]
        [ 0  0  1  0  0]  [0 0 1 0 0]  [1 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]
        [ 0  0  0  1  0]  [0 0 0 1 0]  [0 0 0 1 0]  [1 0 0 1 0]  [0 1 0 1 0]
        [ 0  0  0  0 -1], [0 0 0 0 1], [0 0 0 0 1], [0 0 0 0 1], [0 0 0 0 1],
        <BLANKLINE>
        [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]  [1 0 0 0 0]
        [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]  [0 1 0 0 0]
        [0 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]  [0 0 1 0 0]
        [0 0 1 1 0]  [0 0 0 1 0]  [0 0 0 1 0]  [0 0 0 1 0]  [0 0 0 1 0]
        [0 0 0 0 1], [1 0 0 0 1], [0 1 0 0 1], [0 0 1 0 1], [0 0 0 1 1]
        ]

    TESTS::

        sage: G = matrix.diagonal(R, [])
        sage: _gens_mod_p(G)
        []
        sage: G = matrix.diagonal(R, [3*1])
        sage: _gens_mod_p(G)
        [[-1]]
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


        sage: from sage.groups.fqf_orthogonal_gens import _gens_mod_2
        sage: R = Zp(2,type='fixed-mod',print_mode='terse',show_prec=False,prec=6)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 1, [1])
        sage: W2 = matrix(R, 2, [1, 0, 0, 1])
        sage: G = matrix.block_diagonal([4*U, 4*W0, 2*U, V, W1])
        sage: gen = _gens_mod_2(G)
        sage: len(gen)
        30
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

        sage: from sage.groups.fqf_orthogonal_gens import _ker_gens
        sage: from sage.groups.fqf_orthogonal_lift import Hensel_qf
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [0, 1, 1, 1])
        sage: W1 = matrix(R, 1, [1])
        sage: W2 = matrix(R, 2, [2, 1, 1, 1])
        sage: G = matrix.block_diagonal([V*4, U*2, U])
        sage: gen = _ker_gens(G, 2, 4, [0, 0, 0])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        8
        sage: G = matrix.block_diagonal([V*4, W1*2, W2])
        sage: gen = _ker_gens(G, 2, 3, [0, 1, 1])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        4
        sage: G = matrix.block_diagonal([V*4, U*2, W1])
        sage: gen = _ker_gens(G, 2, 4, [0, 0, 1])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        2
        sage: G = matrix.block_diagonal([V*4, W2*2, U])
        sage: gen = _ker_gens(G, 2, 4, [0, 1, 0])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        6
        sage: G = matrix.block_diagonal([W1*4, U*2, U])
        sage: gen = _ker_gens(G, 1, 3, [1,0,0])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        6
        sage: G = matrix.block_diagonal([W1*4, U*2, W0])
        sage: gen = _ker_gens(G, 1, 3, [1,0,1])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        5
        sage: G = matrix.block_diagonal([W0*4, W0*2, U])
        sage: gen = _ker_gens(G, 2, 4, [1,1,0])
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
        sage: len(gen)
        6
        sage: G = matrix.block_diagonal([W1*4, U*2, W2*2, W1])
        sage: gen = _ker_gens(G, 2, 6, [1,1,1])
        sage: gen
        []
        sage: gen = [Hensel_qf(G, g, 1, 2) for g in gen]
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
                g[e, i1-1] = (g[e,i1:i2]*G[i1:i2,i1:i2]*g[e,i1:i2].T)[0,0]//4
                # the second row depends on the third
                g[i1:i2,i1-1] = - (G[i1:i2,i1:i2]* g[i2:,i1:i2].T * G[i2:,i2:].inverse())[:,-1]/2
            if g[i,j] == 1:   # no need to append the identity
                gens.append(g)
    return gens

def _gens(G, b):
    r"""
    Return generators.

    INPUT:

    - ``G`` -- a non-singular, integral, symmetric `p`-adic matrix in
      descending normal form
    - ``b`` -- a positive integer

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _gens
        sage: R = Zp(2, type='fixed-mod', prec=10, print_mode='terse', show_prec=False, print_pos=False)
        sage: U = matrix(R, 2, [0, 1, 1, 0])
        sage: V = matrix(R, 2, [2, 1, 1, 2])
        sage: W0 = matrix(R, 2, [1, 0, 0, 3])
        sage: W1 = matrix(R, 2, [1, 0, 0, 1])
        sage: G = matrix.block_diagonal([2*U, V])
        sage: gens = _gens(G, 2)
        sage: G = matrix.block_diagonal([2*U, W1])
        sage: gens = _gens(G, 2)
        sage: G = matrix.diagonal(R, [2, 1])
        sage: gens = _gens(G, 2)
        sage: G = matrix.block_diagonal([2*V, V])
        sage: gens = _gens(G, 2)
        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2, 2, 4, 4])
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

def _compute_gens(T, deg=True):
    r"""
    Return generators of the orthogonal group of ``T``.

    INPUT:

    - ``T`` -- torsion orthogonal module in normal form.
    - ``deg`` -- avoids an infinite recursion

    OUTPUT:

    - a list of matrices -- the generators

    EXAMPLES::

        sage: from sage.groups.fqf_orthogonal_gens import _compute_gens
        sage: T = TorsionQuadraticForm(matrix.diagonal([2/3,2/3]))
        sage: _compute_gens(T)
        [
        [  1 726]  [0 1]  [0 1]
        [  3   1], [2 0], [1 0]
        ]
    """
    T = T.normal_form()
    # corner case
    invs = T.invariants()
    if len(invs) == 0:
        return []
    # a well behaved degenerate case
    if deg and T.is_degenerate() and T.invariants()[-1].is_prime():
        from sage.modules.torsion_quadratic_module import _normalize
        p = T.invariants()[-1]
        n = len(T.invariants())
        N = _normalize(T)
        q = N.gram_matrix_quadratic()
        k = len([i for i in range(n) if q[:,i]==0])
        r = n - k
        Idk = matrix.identity(k)
        Idr = matrix.identity(r)
        NR = N.submodule_with_gens(N.gens()[k:])
        gensNR = [NR._to_smith()*g*NR._to_gens() for g in _compute_gens(NR,deg=False)]
        if k > 0:
            gens = [matrix.block_diagonal([Idk,g]) for g in gensNR]
            gens += [matrix.block_diagonal([g.matrix(),Idr]) for g in GL(k,p).gens()]
        else:
            gens = gensNR
        if k>0 and r>0:
            h = matrix.identity(n)
            for g in MatrixSpace(GF(p),r,k).basis():
                h[k:,:k] = g
                gens.append(h)
                gens.append(copy(h))
        return [N._to_gens() * g * N._to_smith() for g in gens]
    elif deg and T.is_degenerate():
        return _isom_fqf(T)

    # normal form gens for the different primes
    blocks = []
    gensT_orders = [t.order() for t in T.gens()]
    n = len(T.gens())
    P = T.annihilator().gen().prime_divisors()
    for p in P:
        indices = []
        for k in range(len(gensT_orders)):
            if mod(gensT_orders[k], p) == 0:
                indices.append(k)
        blocks.append([p, tuple(indices)])

    to_smith = T.to_smith()
    to_normal = T.to_gens()

    # compute generators of the orthogonal groups
    gens = []
    for bl in blocks:
        # compute the generators of the p-primary part
        # the whole group is the direct product of the p-primary parts
        p = bl[0]
        indices = bl[1]
        q_p = T.gram_matrix_quadratic()[indices, indices]
        b = invs[-1].valuation(p)
        G_p = q_p * p**b
        if p != 2:
            # make sure each homogeneous block of G_p stays in normal form
            G_p = G_p/2

        R = Zp(p, type='fixed-mod', prec=b+5, show_prec=False)
        G_p = G_p.change_ring(R)

        # the generators in matrix form
        gens_mat = _gens(G_p, b)

        # extend as identity on the orthogonal complement
        E1 = matrix.identity(indices[0])
        E2 = matrix.identity(n - indices[-1] - 1)
        for g in gens_mat:
            g = g.change_ring(ZZ)
            g = matrix.block_diagonal([E1,g,E2])
            g = to_normal * g * to_smith
            gens.append(g)
    return gens
