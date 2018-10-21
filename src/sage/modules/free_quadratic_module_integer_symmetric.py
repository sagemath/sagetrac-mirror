r"""
Integral lattices

An integral lattice is a finitely generated free abelian group
`L \cong \ZZ^r` equipped with a non-degenerate, symmetric bilinear
form `L \times L \colon \rightarrow \ZZ`.

Here, lattices have an ambient quadratic space `\QQ^n` and
a distinguished basis.

EXAMPLES::

    sage: M = Matrix(ZZ, [[0,1], [1,0]])
    sage: IntegralLattice(M)
    Lattice of degree 2 and rank 2 over Integer Ring
    Basis matrix:
    [1 0]
    [0 1]
    Inner product matrix:
    [0 1]
    [1 0]

A lattice can be defined by an inner product matrix of the
ambient space and a basis::

    sage: G = matrix.identity(3)
    sage: basis = [[1,-1,0], [0,1,-1]]
    sage: L = IntegralLattice(G, basis)
    sage: L
    Lattice of degree 3 and rank 2 over Integer Ring
    Basis matrix:
    [ 1 -1  0]
    [ 0  1 -1]
    Inner product matrix:
    [1 0 0]
    [0 1 0]
    [0 0 1]

    sage: L.gram_matrix()
    [ 2 -1]
    [-1  2]

AUTHORS:

- Simon Brandhorst (2017-09): First created
- Paolo Menegatti (2018-03): Added IntegralLatticeDirectSum, IntegralLatticeGluing
"""

#*****************************************************************************
#       Copyright (C) 2017 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.modules.free_quadratic_module import FreeQuadraticModule_submodule_with_basis_pid, FreeQuadraticModule
from sage.matrix.constructor import matrix
from sage.structure.element import is_Matrix
from sage.arith.misc import gcd
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.misc.cachefunc import cached_method

###############################################################################
#
# Constructor functions
#
###############################################################################

def IntegralLattice(data, basis=None):
    r"""
    Return the integral lattice spanned by ``basis`` in the ambient space.

    A lattice is a finitely generated free abelian group `L \cong \ZZ^r`
    equipped with a non-degenerate, symmetric bilinear form
    `L \times L \colon \rightarrow \ZZ`. Here, lattices have an
    ambient quadratic space `\QQ^n` and a distinguished basis.

    INPUT:

    The input is a descriptor of the lattice and a (optional) basis.
    - ``data`` -- can be one of the following:

      * a symmetric matrix over the rationals -- the inner product matrix
      * an integer -- the dimension for a euclidian lattice
      * a symmetric Cartan type or anything recognized by
        :class:`CartanMatrix` (see also
        :mod:`Cartan types <sage.combinat.root_system.cartan_type>`)
        -- for a root lattice
      * the string ``"U"`` or ``"H"`` -- for hyperbolic lattices

    - ``basis`` -- (optional) a matrix whose rows form a basis of the
      lattice,  or a list of module elements forming a basis

    OUTPUT:

    A lattice in the ambient space defined by the inner_product_matrix.
    Unless specified, the basis of the lattice is the standard basis.

    EXAMPLES::

        sage: H5 = Matrix(ZZ, 2, [2,1,1,-2])
        sage: IntegralLattice(H5)
        Lattice of degree 2 and rank 2 over Integer Ring
        Basis matrix:
        [1 0]
        [0 1]
        Inner product matrix:
        [ 2  1]
        [ 1 -2]

    A basis can be specified too::

        sage: IntegralLattice(H5, Matrix([1,1]))
        Lattice of degree 2 and rank 1 over Integer Ring
        Basis matrix:
        [1 1]
        Inner product matrix:
        [ 2  1]
        [ 1 -2]

    We can define a Euclidian lattice just by its dimension::

        sage: IntegralLattice(3)
        Lattice of degree 3 and rank 3 over Integer Ring
        Basis matrix:
        [1 0 0]
        [0 1 0]
        [0 0 1]
        Inner product matrix:
        [1 0 0]
        [0 1 0]
        [0 0 1]

    Here is an example of the `A_2` root lattice in Euclidian space::

        sage: basis = Matrix([[1,-1,0], [0,1,-1]])
        sage: A2 = IntegralLattice(3, basis)
        sage: A2
        Lattice of degree 3 and rank 2 over Integer Ring
        Basis matrix:
        [ 1 -1  0]
        [ 0  1 -1]
        Inner product matrix:
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: A2.gram_matrix()
        [ 2 -1]
        [-1  2]

    We use ``"U"`` or ``"H"`` for defining a hyperbolic lattice::

        sage: L1 = IntegralLattice("U")
        sage: L1
        Lattice of degree 2 and rank 2 over Integer Ring
        Basis matrix:
        [1 0]
        [0 1]
        Inner product matrix:
        [0 1]
        [1 0]
        sage: L1 == IntegralLattice("H")
        True

    We can construct root lattices by specifying their type
    (see :mod:`Cartan types <sage.combinat.root_system.cartan_type>`
    and :class:`CartanMatrix`)::

        sage: IntegralLattice(["E", 7])
        Lattice of degree 7 and rank 7 over Integer Ring
        Basis matrix:
        [1 0 0 0 0 0 0]
        [0 1 0 0 0 0 0]
        [0 0 1 0 0 0 0]
        [0 0 0 1 0 0 0]
        [0 0 0 0 1 0 0]
        [0 0 0 0 0 1 0]
        [0 0 0 0 0 0 1]
        Inner product matrix:
        [ 2  0 -1  0  0  0  0]
        [ 0  2  0 -1  0  0  0]
        [-1  0  2 -1  0  0  0]
        [ 0 -1 -1  2 -1  0  0]
        [ 0  0  0 -1  2 -1  0]
        [ 0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0 -1  2]
        sage: IntegralLattice(["A", 2])
        Lattice of degree 2 and rank 2 over Integer Ring
        Basis matrix:
        [1 0]
        [0 1]
        Inner product matrix:
        [ 2 -1]
        [-1  2]
        sage: IntegralLattice("D3")
        Lattice of degree 3 and rank 3 over Integer Ring
        Basis matrix:
        [1 0 0]
        [0 1 0]
        [0 0 1]
        Inner product matrix:
        [ 2 -1 -1]
        [-1  2  0]
        [-1  0  2]
        sage: IntegralLattice(["D", 4])
        Lattice of degree 4 and rank 4 over Integer Ring
        Basis matrix:
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        Inner product matrix:
        [ 2 -1  0  0]
        [-1  2 -1 -1]
        [ 0 -1  2  0]
        [ 0 -1  0  2]

    We can specify a basis as well::

        sage: G = Matrix(ZZ, 2, [0,1,1,0])
        sage: B = [vector([1,1])]
        sage: IntegralLattice(G, basis=B)
        Lattice of degree 2 and rank 1 over Integer Ring
        Basis matrix:
        [1 1]
        Inner product matrix:
        [0 1]
        [1 0]
        sage: IntegralLattice(["A", 3], [[1,1,1]])
        Lattice of degree 3 and rank 1 over Integer Ring
        Basis matrix:
        [1 1 1]
        Inner product matrix:
        [ 2 -1  0]
        [-1  2 -1]
        [ 0 -1  2]
        sage: IntegralLattice(4, [[1,1,1,1]])
        Lattice of degree 4 and rank 1 over Integer Ring
        Basis matrix:
        [1 1 1 1]
        Inner product matrix:
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        sage: IntegralLattice("A2", [[1,1]])
        Lattice of degree 2 and rank 1 over Integer Ring
        Basis matrix:
        [1 1]
        Inner product matrix:
        [ 2 -1]
        [-1  2]

    TESTS::

        sage: IntegralLattice(["A", 1, 1])
        Traceback (most recent call last):
        ...
        ValueError: lattices must be nondegenerate; use FreeQuadraticModule instead
        sage: IntegralLattice(["D", 3, 1])
        Traceback (most recent call last):
        ...
        ValueError: lattices must be nondegenerate; use FreeQuadraticModule instead
    """
    if is_Matrix(data):
        inner_product_matrix = data
    elif isinstance(data, Integer):
        inner_product_matrix = matrix.identity(ZZ, data)
    elif data == "U" or data == "H":
        inner_product_matrix = matrix([[0,1],[1,0]])
    else:
        inner_product_matrix = CartanMatrix(data)
    if basis is None:
        basis = matrix.identity(ZZ, inner_product_matrix.ncols())
    if inner_product_matrix != inner_product_matrix.transpose():
        raise ValueError("the inner product matrix must be symmetric\n%s"
                         % inner_product_matrix)

    A = FreeQuadraticModule(ZZ,
                            inner_product_matrix.ncols(),
                            inner_product_matrix=inner_product_matrix)
    return FreeQuadraticModule_integer_symmetric(ambient=A,
                                                 basis=basis,
                                                 inner_product_matrix=A.inner_product_matrix(),
                                                 already_echelonized=False)

def IntegralLatticeDirectSum(Lattices, return_embeddings=False):
    r"""
    Return the direct sum of the lattices contained in the list ``Lattices``.

    INPUT:

    - ``Lattices`` -- a list of lattices ``[L_1,...,L_n]``
    - ``return_embeddings`` -- (default: ``False``) a boolean

    OUTPUT:

    The direct sum of the `L_i` if ``return_embeddings`` is ``False`` or
    the tuple ``[L, phi]`` where `L` is the direct sum of `L_i` and ``phi``
    is the list of embeddings from `L_i` to `L`.

    EXAMPLES::

        sage: from sage.modules.free_quadratic_module_integer_symmetric import IntegralLatticeDirectSum
        sage: L1 = IntegralLattice("D4")
        sage: L2 = IntegralLattice("A3", [[1, 1, 2]])
        sage: L3 = IntegralLattice("A4", [[0, 1, 1, 2], [1, 2, 3, 1]])
        sage: Lattices = [L1, L2, L3]
        sage: IntegralLatticeDirectSum([L1, L2, L3])
        Lattice of degree 11 and rank 7 over Integer Ring
        Basis matrix:
        [1 0 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 1 1 2 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 1 2]
        [0 0 0 0 0 0 0 1 2 3 1]
        Inner product matrix:
        [ 2 -1  0  0  0  0  0  0  0  0  0]
        [-1  2 -1 -1  0  0  0  0  0  0  0]
        [ 0 -1  2  0  0  0  0  0  0  0  0]
        [ 0 -1  0  2  0  0  0  0  0  0  0]
        [ 0  0  0  0  2 -1  0  0  0  0  0]
        [ 0  0  0  0 -1  2 -1  0  0  0  0]
        [ 0  0  0  0  0 -1  2  0  0  0  0]
        [ 0  0  0  0  0  0  0  2 -1  0  0]
        [ 0  0  0  0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0  0  0  0 -1  2]
        sage: [L, phi] = IntegralLatticeDirectSum([L1, L2, L3], True)
        sage: LL3 = L.sublattice(phi[2].image().basis_matrix())
        sage: L3.discriminant() == LL3.discriminant()
        True
        sage: x = L3([1, 2, 3, 1])
        sage: phi[2](x).inner_product(phi[2](x)) == x.inner_product(x)
        True

    TESTS::

        sage: IntegralLatticeDirectSum([IntegralLattice("D4")])
        Lattice of degree 4 and rank 4 over Integer Ring
        Basis matrix:
        [1 0 0 0]
        [0 1 0 0]
        [0 0 1 0]
        [0 0 0 1]
        Inner product matrix:
        [ 2 -1  0  0]
        [-1  2 -1 -1]
        [ 0 -1  2  0]
        [ 0 -1  0  2]

        sage: L1 = IntegralLattice(2 * matrix.identity(2), [[1/2, 1/2]])
        sage: L2 = IntegralLattice("A3", [[1, 1, 2]])
        sage: [L, phi] = IntegralLatticeDirectSum([L1, L2], True)
        sage: L
        Lattice of degree 5 and rank 2 over Integer Ring
        Basis matrix:
        [1/2 1/2   0   0   0]
        [  0   0   1   1   2]
        Inner product matrix:
        [ 2  0  0  0  0]
        [ 0  2  0  0  0]
        [ 0  0  2 -1  0]
        [ 0  0 -1  2 -1]
        [ 0  0  0 -1  2]
    """
    for L in Lattices:
        if not isinstance(L, FreeQuadraticModule_integer_symmetric):
            raise ValueError("Lattices must be a list of lattices")
    N = len(Lattices)
    dims = [L_i.dimension() for L_i in Lattices]
    degrees = [L_i.degree() for L_i in Lattices]
    degree_tot = sum(degrees)
    sum_degree = [sum(degrees[:i]) for i in range(N+1)]
    inner_product_list = [copy(L_i.inner_product_matrix()) for L_i in Lattices]
    IM = matrix.block_diagonal(inner_product_list)
    ambient = FreeQuadraticModule(ZZ,
                                  degree_tot,
                                  inner_product_matrix=IM)
    basis = [matrix.block(1, 3, [matrix.zero(dims[i], sum_degree[i]),
                                 Lattices[i].basis_matrix(),
                                 matrix.zero(dims[i], sum_degree[-1] - sum_degree[i+1])
                                ])  for i in range(N)]
    basis_matrix = matrix.block(N, 1, basis)
    ipm = ambient.inner_product_matrix()
    direct_sum = FreeQuadraticModule_integer_symmetric(ambient=ambient,
                                                       basis=basis_matrix,
                                                       inner_product_matrix=ipm,
                                                       already_echelonized=False)
    if not return_embeddings:
        return direct_sum
    sum_dims = [sum(dims[:i]) for i in range(N+1)]
    phi = [Lattices[i].hom(direct_sum.basis()[sum_dims[i]:sum_dims[i+1]])
           for i in range(N)]
    return [direct_sum, phi]

def IntegralLatticeGluing(Lattices, glue, return_embeddings=False):
    r"""
    Return an overlattice of the direct sum as defined by ``glue``.

    INPUT:

    - ``Lattices`` -- a list of lattices `[L_1,...,L_n]`
    - ``glue`` -- a list where the elements are lists in the form `[g_1,...,g_n]`;
      here `g_i` is an element of the discriminant group of `L_i`and
      the overlattice is spanned by the additional ``[sum(g) for g in glue]``
    - ``return_embeddings`` -- (default: ``False``) a boolean


    OUTPUT:

    The glued lattice given by `L_i` and ``glue`` if ``return_embeddings``
    is ``False`` or the tuple ``[L, phi]`` where `L` is the glued lattice and
    ``phi`` the list of embeddings from `L_i` to `L`

    EXAMPLES:

    A single lattice can be glued. This is the same as taking an overlattice::

        sage: from sage.modules.free_quadratic_module_integer_symmetric import IntegralLatticeGluing
        sage: L1 = IntegralLattice(matrix([[4]]))
        sage: g1 = L1.discriminant_group().gens()[0]
        sage: glue = [[2 * g1]]
        sage: L = IntegralLatticeGluing([L1], glue)
        sage: L
        Lattice of degree 1 and rank 1 over Integer Ring
        Basis matrix:
        [1/2]
        Inner product matrix:
        [4]
        sage: L.gram_matrix()
        [1]
        sage: IntegralLatticeGluing([L1], glue, return_embeddings=True)
        [Lattice of degree 1 and rank 1 over Integer Ring
         Basis matrix:
         [1/2]
         Inner product matrix:
         [4], [Free module morphism defined by the matrix
          [2]
          Domain: Lattice of degree 1 and rank 1 over Integer Ring
          Basis matrix:
          [1]
          Inner product matrix:
          [4]
          Codomain: Lattice of degree 1 and rank 1 over Integer Ring
          Basis matrix:
          [1/2]
          Inner product matrix:
          [4]]]

        sage: L1 = IntegralLattice([[2]])
        sage: L2 = IntegralLattice([[2]])
        sage: AL1 = L1.discriminant_group()
        sage: AL2 = L2.discriminant_group()
        sage: AL1
        Finite quadratic module over Integer Ring with invariants (2,)
        Gram matrix of the quadratic form with values in Q/2Z:
        [1/2]
        sage: g1 = L1.discriminant_group().gens()[0]
        sage: g2 = L2.discriminant_group().gens()[0]
        sage: glue = [[g1, g2]]
        sage: IntegralLatticeGluing([L1, L2], glue)
        Lattice of degree 2 and rank 2 over Integer Ring
        Basis matrix:
        [1/2 1/2]
        [  0   1]
        Inner product matrix:
        [2 0]
        [0 2]

        sage: L1 = IntegralLattice("A4")
        sage: L2 = IntegralLattice("A4")
        sage: g1 = L1.discriminant_group().gens()[0]
        sage: g2 = L2.discriminant_group().gens()[0]
        sage: glue = [[g1, 2 * g2]]
        sage: [V, phi] = IntegralLatticeGluing([L1, L2], glue, True)
        sage: V
        Lattice of degree 8 and rank 8 over Integer Ring
        Basis matrix:
        [1/5 2/5 3/5 4/5 2/5 4/5 1/5 3/5]
        [  0   1   0   0   0   0   0   0]
        [  0   0   1   0   0   0   0   0]
        [  0   0   0   1   0   0   0   0]
        [  0   0   0   0   1   0   0   0]
        [  0   0   0   0   0   1   0   0]
        [  0   0   0   0   0   0   1   0]
        [  0   0   0   0   0   0   0   1]
        Inner product matrix:
        [ 2 -1  0  0  0  0  0  0]
        [-1  2 -1  0  0  0  0  0]
        [ 0 -1  2 -1  0  0  0  0]
        [ 0  0 -1  2  0  0  0  0]
        [ 0  0  0  0  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]
        sage: V.sublattice(phi[0].image().basis_matrix())
        Lattice of degree 8 and rank 4 over Integer Ring
        Basis matrix:
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]
        Inner product matrix:
        [ 2 -1  0  0  0  0  0  0]
        [-1  2 -1  0  0  0  0  0]
        [ 0 -1  2 -1  0  0  0  0]
        [ 0  0 -1  2  0  0  0  0]
        [ 0  0  0  0  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]

    Different gluings can be composed::

        sage: D4 = IntegralLattice("D4")
        sage: D4.discriminant_group()
        Finite quadratic module over Integer Ring with invariants (2, 2)
        Gram matrix of the quadratic form with values in Q/2Z:
        [  1 1/2]
        [1/2   1]
        sage: L2 = IntegralLattice(2 * matrix.identity(2))
        sage: L2.discriminant_group()
        Finite quadratic module over Integer Ring with invariants (2, 2)
        Gram matrix of the quadratic form with values in Q/2Z:
        [1/2   0]
        [  0 1/2]
        sage: g1 = D4.discriminant_group().gens()[0]
        sage: g2 = L2.discriminant_group().gens()[0] + L2.discriminant_group().gens()[1]
        sage: D6, phi = IntegralLatticeGluing([D4, L2], [[g1, g2]], True)
        sage: AD6 = D6.discriminant_group()
        sage: AD6.normal_form()
        Finite quadratic module over Integer Ring with invariants (2, 2)
        Gram matrix of the quadratic form with values in Q/2Z:
        [3/2   0]
        [  0 3/2]
        sage: f1, g1 = AD6.normal_form().gens()
        sage: f2, g2 = L2.discriminant_group().gens()
        sage: E8, psi = IntegralLatticeGluing([D6, L2], [[f1, f2], [g1, g2]], True)
        sage: D4embed = E8.sublattice(psi[0](phi[0].image()).basis_matrix())
        sage: x = D4([1, 0, 0, 0])
        sage: psi[0](phi[0](x)).inner_product(psi[0](phi[0](x)))==x.inner_product(x)
        True
        sage: D4embed
        Lattice of degree 8 and rank 4 over Integer Ring
        Basis matrix:
        [1 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0 0]
        Inner product matrix:
        [ 2 -1  0  0  0  0  0  0]
        [-1  2 -1 -1  0  0  0  0]
        [ 0 -1  2  0  0  0  0  0]
        [ 0 -1  0  2  0  0  0  0]
        [ 0  0  0  0  2  0  0  0]
        [ 0  0  0  0  0  2  0  0]
        [ 0  0  0  0  0  0  2  0]
        [ 0  0  0  0  0  0  0  2]

    The input may be a list of three or more lattices::

        sage: A7 = IntegralLattice("A7")
        sage: D5 = IntegralLattice("D5")
        sage: gA7 = A7.discriminant_group().gens()[0]
        sage: gD5 = D5.discriminant_group().gens()[0]
        sage: [L, phi] = IntegralLatticeGluing([A7, A7, D5, D5],
        ....:                          [[gA7, gA7, gD5, 2 * gD5],
        ....:                          [gA7, 7 * gA7, 2 * gD5, gD5]], True)
        sage: L.determinant()
        1
        sage: B = phi[0].matrix()
        sage: B*L.gram_matrix()*B.transpose()==A7.gram_matrix()
        True

    The gluing takes place in the direct sum of the respective ambient spaces::

        sage: L1 = IntegralLattice("D4", [[1, 1, 0, 0], [0, 1, 1, 0]])
        sage: L2 = IntegralLattice("E6", [[0, 2, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1]])
        sage: [f1, f2] = L1.discriminant_group().gens()
        sage: [g1, g2] = L2.discriminant_group().gens()
        sage: [L, phi] = IntegralLatticeGluing([L1, L2], [[f1, g1], [f2, 2 * g2]], True)
        sage: phi[0]
        Free module morphism defined by the matrix
        [ 2  2 -1 -2]
        [ 0  2  0 -1]
        Domain: Lattice of degree 4 and rank 2 over Integer Ring
        Basis matrix:
        [1 1 0 0]
        [0 1 1 0]
        Inner product matrix:
        [ 2 -1  0  0]
        [-1  2 -1 -1]
        [ 0 -1  2  0]
        [ 0 -1  0  2]
        Codomain: Lattice of degree 10 and rank 4 over Integer Ring
        Basis matrix:
        [ 1/2    0 -1/2    0    0  1/2    0    0  1/2  1/2]
        [   0  1/2  1/2    0    0    0    0    0  1/2  1/2]
        [   0    0    0    0    0    1    0    0    0    0]
        [   0    0    0    0    0    0    0    0    1    1]
        Inner product matrix:
        [ 2 -1  0  0  0  0  0  0  0  0]
        [-1  2 -1 -1  0  0  0  0  0  0]
        [ 0 -1  2  0  0  0  0  0  0  0]
        [ 0 -1  0  2  0  0  0  0  0  0]
        [ 0  0  0  0  2  0 -1  0  0  0]
        [ 0  0  0  0  0  2  0 -1  0  0]
        [ 0  0  0  0 -1  0  2 -1  0  0]
        [ 0  0  0  0  0 -1 -1  2 -1  0]
        [ 0  0  0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0  0  0 -1  2]
        sage: B = phi[0].matrix()
        sage: B * L.gram_matrix() * B.transpose()==L1.gram_matrix()
        True
    """
    [direct_sum, phi] = IntegralLatticeDirectSum(Lattices, return_embeddings=True)
    N = len(Lattices)
    for g in glue:
        if not len(g)==N:
            raise ValueError("the lengths of the lists do not match")
    for i in range(N):
        ALi = Lattices[i].discriminant_group()
        for g in glue:
            try:
                ALi(g[i])
            except:
                raise ValueError("the gluing vectors must be in the"
                                 "corresponding discriminant groups")
    generators = [sum(phi[i](g[i].lift()*g[i].order())/g[i].order()
                      for i in range(N))
                  for g in glue]
    glued_lattice = direct_sum.overlattice(generators)
    if not return_embeddings:
        return glued_lattice
    HomSpaces = [Lattices[i].Hom(glued_lattice) for i in range(N)]
    G = glued_lattice.basis_matrix().solve_left(direct_sum.basis_matrix())
    f = [HomSpaces[i](phi[i].matrix() * G) for i in range(N)]
    return [glued_lattice, f]

###############################################################################
#
# Base class for Lattices
#
###############################################################################

class FreeQuadraticModule_integer_symmetric(FreeQuadraticModule_submodule_with_basis_pid):
    r"""
    This class represents non-degenerate, integral,
    symmetric free quadratic `\ZZ`-modules.

    INPUT:

    - ``ambient`` -- an ambient free quadratic module
    - ``basis`` -- a list of elements of ambient or a matrix
    - ``inner_product_matrix`` -- a symmetric matrix over the rationals

    EXAMPLES::

        sage: IntegralLattice("U",basis=[vector([1,1])])
        Lattice of degree 2 and rank 1 over Integer Ring
        Basis matrix:
        [1 1]
        Inner product matrix:
        [0 1]
        [1 0]
    """
    def __init__(self, ambient, basis, inner_product_matrix,
                 check=True, already_echelonized=False):
        r"""
        Create the integral lattice spanned by ``basis`` in the ambient space.

        TESTS::

            sage: L = IntegralLattice("U")
            sage: TestSuite(L).run()
        """
        FreeQuadraticModule_submodule_with_basis_pid.__init__(
                                        self,
                                        ambient,
                                        basis,
                                        inner_product_matrix,
                                        check=check,
                                        already_echelonized=already_echelonized)
        if self.determinant() == 0:
            raise ValueError("lattices must be nondegenerate; "
                            "use FreeQuadraticModule instead")
        if self.gram_matrix().base_ring() is not ZZ:
            if self.gram_matrix().denominator() != 1:
                raise ValueError("lattices must be integral; "
                            "use FreeQuadraticModule instead")

    def _mul_(self, other, switch_sides=False):
        r"""
        Multiplication of the basis by ``other``.

        EXAMPLES::

            sage: M = Matrix(ZZ,2,[1,2,2,-1])
            sage: L = IntegralLattice(M)
            sage: 2 * L
            Lattice of degree 2 and rank 2 over Integer Ring
            Basis matrix:
            [2 0]
            [0 2]
            Inner product matrix:
            [ 1  2]
            [ 2 -1]
            sage: L * matrix(ZZ,2,[1,2,3,4])
            Lattice of degree 2 and rank 2 over Integer Ring
            Basis matrix:
            [1 2]
            [3 4]
            Inner product matrix:
            [ 1  2]
            [ 2 -1]
        """
        B = self.basis_matrix()
        B = other * B if switch_sides else B * other
        # check whether it is integral
        if other in ZZ or other.denominator()==1:
            return self.sublattice(B.rows())
        else:
            return self.span(B.rows())

    def _repr_(self):
        r"""
        The print representation of this lattice.

        EXAMPLES::

            sage: A2 = IntegralLattice("A2")
            sage: A2
            Lattice of degree 2 and rank 2 over Integer Ring
            Basis matrix:
            [1 0]
            [0 1]
            Inner product matrix:
            [ 2 -1]
            [-1  2]
        """
        if self.is_sparse():
            s = "Sparse lattice of degree %s and rank %s over %s\n"%(
                self.degree(), self.rank(), self.base_ring()) + \
                "Basis matrix:\n%s\n" % self.basis_matrix() + \
                "Inner product matrix:\n%s" % self.inner_product_matrix()
        else:
            s = "Lattice of degree %s and rank %s over %s\n"%(
                self.degree(), self.rank(), self.base_ring()) + \
                "Basis matrix:\n%s\n" % self.basis_matrix() + \
                "Inner product matrix:\n%s" % self.inner_product_matrix()
        return s

    @cached_method
    def is_even(self):
        r"""
        Return whether the diagonal entries of the Gram matrix are even.

        EXAMPLES::

            sage: G = Matrix(ZZ,2,2,[-1,1,1,2])
            sage: L = IntegralLattice(G)
            sage: L.is_even()
            False
            sage: L = IntegralLattice("A2")
            sage: L.is_even()
            True
        """
        return all(d % 2 == 0 for d in self.gram_matrix().diagonal())

    @cached_method
    def dual_lattice(self):
        r"""
        Return the dual lattice as a :class:`FreeQuadraticModule`

        Let `L` be a lattice. Its dual lattice is

        .. MATH::

            L^\vee = \{x \in L \otimes \QQ :  (x, l) \in \ZZ \; \forall l \in L \}.

        EXAMPLES::

            sage: L = IntegralLattice("A2")
            sage: Ldual=L.dual_lattice()
            sage: Ldual
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/3 2/3]
            [  0   1]

        Since our lattices are always integral, a lattice is contained in its dual::

            sage: L.is_submodule(Ldual)
            True
        """
        return self.span(self.gram_matrix().inverse()*self.basis_matrix())

    @cached_method
    def discriminant_group(self, s=0):
        r"""
        Return the discriminant group `L^\vee / L` of this lattice.

        INPUT:

        - ``s`` -- an integer (default: 0)

        OUTPUT:

        The `s` primary part of the discriminant group.
        If `s=0`, returns the whole discriminant group.

        EXAMPLES::

            sage: L = IntegralLattice(Matrix(ZZ,2,2,[2,1,1,-2])*2)
            sage: L.discriminant_group()
            Finite quadratic module over Integer Ring with invariants (2, 10)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2 9/5]
            sage: L.discriminant_group(2)
            Finite quadratic module over Integer Ring with invariants (2, 2)
            Gram matrix of the quadratic form with values in Q/2Z:
            [  1 1/2]
            [1/2   1]
            sage: L.discriminant_group(5)
            Finite quadratic module over Integer Ring with invariants (5,)
            Gram matrix of the quadratic form with values in Q/2Z:
            [6/5]

        TESTS::

            sage: L = IntegralLattice("H")
            sage: L.discriminant_group()
            Finite quadratic module over Integer Ring with invariants ()
            Gram matrix of the quadratic form with values in Q/2Z:
            []
        """
        from sage.modules.torsion_quadratic_module import TorsionQuadraticModule
        D = TorsionQuadraticModule(self.dual_lattice(), self)
        d = D.annihilator().gen()
        a = d.prime_to_m_part(s)
        Dp_gens = [a*g for g in D.gens()]
        return D.submodule(Dp_gens)

    def signature(self):
        r"""
        Return the signature of this lattice, which is defined as
        the difference between the number of positive eigenvalues and
        the number of negative eigenvalues in the Gram matrix.

        EXAMPLES::

            sage: U = IntegralLattice("U")
            sage: U.signature()
            0
        """
        sig = self.signature_pair()
        return sig[0] - sig[1]

    @cached_method
    def signature_pair(self):
        r"""
        Return the signature tuple `(n_+,n_-)` of this lattice.

        Here `n_+` (resp. `n_-`) is the number of positive (resp. negative)
        eigenvalues of the Gram matrix.

        EXAMPLES::


            sage: A2 = IntegralLattice("A2")
            sage: A2.signature_pair()
            (2, 0)
        """
        from sage.quadratic_forms.quadratic_form import QuadraticForm
        return QuadraticForm(QQ, self.gram_matrix()).signature_vector()[:2]

    def direct_sum(self, M):
        r"""
        Return the direct sum of this lattice with ``M``.

        INPUT:

        - ``M`` -- a module over `\ZZ`

        EXAMPLES::

            sage: A = IntegralLattice(1)
            sage: A.direct_sum(A)
            Lattice of degree 2 and rank 2 over Integer Ring
            Basis matrix:
            [1 0]
            [0 1]
            Inner product matrix:
            [1 0]
            [0 1]
        """
        IM = matrix.block_diagonal([self.inner_product_matrix(),
                                    M.inner_product_matrix()])
        ambient = FreeQuadraticModule(ZZ,
                                      self.degree() + M.degree(), IM)
        smzero = matrix.zero(self.rank(), M.degree())
        mszero = matrix.zero(M.rank(), self.degree())
        basis = self.basis_matrix().augment(smzero).stack(
                            mszero.augment(M.basis_matrix()))
        ipm = ambient.inner_product_matrix()
        return FreeQuadraticModule_integer_symmetric(ambient=ambient,
                                                     basis=basis,
                                                     inner_product_matrix=ipm,
                                                     already_echelonized=False)

    def is_primitive(self, M):
        r"""
        Return whether ``M`` is a primitive submodule of this lattice.

        A `\ZZ`-submodule ``M`` of a `\ZZ`-module ``L`` is called primitive if
        the quotient ``L/M`` is torsion free.

        INPUT:

        - ``M`` -- a submodule of this lattice

        EXAMPLES::

            sage: U = IntegralLattice("U")
            sage: L1 = U.span([vector([1,1])])
            sage: L2 = U.span([vector([1,-1])])
            sage: U.is_primitive(L1)
            True
            sage: U.is_primitive(L2)
            True
            sage: U.is_primitive(L1+L2)
            False

        We can also compute the index::

            sage: (L1+L2).index_in(U)
            2
        """
        return (gcd((self/M).invariants()) == 0)

    def orthogonal_complement(self, M):
        r"""
        Return the orthogonal complement of ``M`` in this lattice.

        INPUT:

        - ``M`` -- a module in the same ambient space or
          a list of elements of the ambient space

        EXAMPLES::

            sage: H5 = Matrix(ZZ,2,[2,1,1,-2])
            sage: L = IntegralLattice(H5)
            sage: S = L.span([vector([1,1])])
            sage: L.orthogonal_complement(S)
            Lattice of degree 2 and rank 1 over Integer Ring
            Basis matrix:
            [1 3]
            Inner product matrix:
            [ 2  1]
            [ 1 -2]

            sage: L = IntegralLattice(2)
            sage: L.orthogonal_complement([vector(ZZ,[1,0])])
            Lattice of degree 2 and rank 1 over Integer Ring
            Basis matrix:
            [0 1]
            Inner product matrix:
            [1 0]
            [0 1]
        """
        from sage.modules.free_module import FreeModule_generic
        if not isinstance(M,FreeModule_generic):
            M = self.span(M)
        elif M.ambient_vector_space() != self.ambient_vector_space():
            raise ValueError("M must have the same "
                             "ambient vector space as this lattice")

        K = (self.inner_product_matrix() * M.basis_matrix().transpose()).kernel()
        K = self.span(K.basis())
        K = K.base_extend(QQ)
        return self.sublattice(self.intersection(K).basis())

    def sublattice(self, basis):
        r"""
        Return the sublattice spanned by ``basis``.

        INPUT:

        - ``basis`` -- A list of elements of this lattice.

        EXAMPLES::

            sage: U = IntegralLattice("U")
            sage: S = U.sublattice([vector([1,1])])
            sage: S
            Lattice of degree 2 and rank 1 over Integer Ring
            Basis matrix:
            [1 1]
            Inner product matrix:
            [0 1]
            [1 0]
            sage: U.sublattice([vector([1,-1])/2])
            Traceback (most recent call last):
            ...
            ValueError: lattices must be integral; use FreeQuadraticModule instead
            sage: S.sublattice([vector([1,-1])])
            Traceback (most recent call last):
            ...
            ValueError: the basis (= [(1, -1)]) does not span a submodule
        """
        M = FreeQuadraticModule_integer_symmetric(
            ambient=self.ambient_module(), basis=basis,
            inner_product_matrix=self.inner_product_matrix(),
            already_echelonized=False)
        if not M.is_submodule(self):
            raise ValueError("the basis (= %s) does not span "
                             "a submodule" % basis)
        return M

    def overlattice(self, gens):
        r"""
        Return the lattice spanned by this lattice and ``gens``.

        INPUT:

        - ``gens`` -- a list of elements or a rational matrix

        EXAMPLES::

            sage: L = IntegralLattice(Matrix(ZZ,2,2,[2,0,0,2]))
            sage: M = L.overlattice([vector([1,1])/2])
            sage: M.gram_matrix()
            [1 1]
            [1 2]
        """
        basis = (self + self.span(gens)).basis()
        return FreeQuadraticModule_integer_symmetric(
            ambient=self.ambient_module(), basis=basis,
            inner_product_matrix=self.inner_product_matrix(),
            already_echelonized=False)

    def orthogonal_group(self, gens=None, is_finite=None):
        """
        Return the orthogonal group of this lattice as a matrix group.

        The elements are isometries of the ambient vector space
        which preserve this lattice. They are represented by
        matrices with respect to the standard basis.

        INPUT:

        - ``gens`` -- a list of matrices (default:``None``)
        - ``is_finite`` -- bool (default: ``None``) If set to ``True``,
          then the group is placed in the category of finite groups. Sage does not check this.

        OUTPUT:

        The matrix group generated by ``gens``.
        If ``gens`` is not specified, then generators of the full
        orthogonal group of this lattice are computed. They are
        continued as the identity on the orthogonal complement of
        the lattice in its ambient space. Currently, we can only
        compute the orthogonal group for positive definite lattices.

        EXAMPLES::

            sage: A4 = IntegralLattice("A4")
            sage: Aut = A4.orthogonal_group()
            sage: Aut
            Group of isometries with 5 generators (
            [-1  0  0  0]  [0 0 0 1]  [-1 -1 -1  0]  [ 1  0  0  0]  [ 1  0  0  0]
            [ 0 -1  0  0]  [0 0 1 0]  [ 0  0  0 -1]  [-1 -1 -1 -1]  [ 0  1  0  0]
            [ 0  0 -1  0]  [0 1 0 0]  [ 0  0  1  1]  [ 0  0  0  1]  [ 0  0  1  1]
            [ 0  0  0 -1], [1 0 0 0], [ 0  1  0  0], [ 0  0  1  0], [ 0  0  0 -1]
            )

        The group acts from the right on the lattice and its discriminant group::

            sage: x = A4.an_element()
            sage: g = Aut.an_element()
            sage: g
            [ 1  1  1  0]
            [ 0  0 -1  0]
            [ 0  0  1  1]
            [ 0 -1 -1 -1]
            sage: x*g
            (1, 1, 1, 0)
            sage: (x*g).parent()==A4
            True
            sage: (g*x).parent()
            Vector space of dimension 4 over Rational Field
            sage: y = A4.discriminant_group().an_element()
            sage: y*g
            (1)

        If the group is finite we can compute the usual things::

            sage: Aut.order()
            240
            sage: conj = Aut.conjugacy_classes_representatives()
            sage: len(conj)
            14
            sage: Aut.structure_description()   # optional - database_gap
            'C2 x S5'

        The lattice can live in a larger ambient space::

            sage: A2 = IntegralLattice(matrix.identity(3),Matrix(ZZ,2,3,[1,-1,0,0,1,-1]))
            sage: A2.orthogonal_group()
            Group of isometries with 3 generators (
            [-1/3  2/3  2/3]  [ 2/3  2/3 -1/3]  [1 0 0]
            [ 2/3 -1/3  2/3]  [ 2/3 -1/3  2/3]  [0 0 1]
            [ 2/3  2/3 -1/3], [-1/3  2/3  2/3], [0 1 0]
            )

        It can be negative definite as well::

            sage: A2m = IntegralLattice(-Matrix(ZZ,2,[2,1,1,2]))
            sage: G = A2m.orthogonal_group()
            sage: G.order()
            12

        If the lattice is indefinite, sage does not know how to compute generators.
        Can you teach it?::

            sage: U = IntegralLattice(Matrix(ZZ,2,[0,1,1,0]))
            sage: U.orthogonal_group()
            Traceback (most recent call last):
            ...
            NotImplementedError: currently, we can only compute generators for orthogonal groups over definite lattices.

        But we can define subgroups::

            sage: S = IntegralLattice(Matrix(ZZ,2,[2, 3, 3, 2]))
            sage: f = Matrix(ZZ,2,[0,1,-1,3])
            sage: S.orthogonal_group([f])
            Group of isometries with 1 generator (
            [ 0  1]
            [-1  3]
            )

        TESTS:

        We can handle the trivial group::

            sage: S = IntegralLattice(Matrix(ZZ,2,[2, 3, 3, 2]))
            sage: S.orthogonal_group([])
            Group of isometries with 1 generator (
            [1 0]
            [0 1]
            )
        """
        from sage.categories.groups import Groups
        from sage.groups.matrix_gps.isometries import GroupOfIsometries
        sig = self.signature_pair()
        if gens is None:
            gens = []
            if sig[1]==0 or sig[0]==0: #definite
                from sage.quadratic_forms.quadratic_form import QuadraticForm
                is_finite = True
                # Compute transformation matrix to the ambient module.
                L = self.overlattice(self.ambient_module().gens())
                Orthogonal = L.orthogonal_complement(self)
                B = self.basis_matrix().stack(Orthogonal.basis_matrix())
                if sig[0] == 0: #negative definite
                    q = QuadraticForm(ZZ, -2*self.gram_matrix())
                else:    # positve definite
                    q = QuadraticForm(ZZ, 2*self.gram_matrix())
                identity = matrix.identity(Orthogonal.rank())
                for g in q.automorphism_group().gens():
                    g = g.matrix().T
                    # We continue g as identity on the orthogonal complement.
                    g = matrix.block_diagonal([g, identity])
                    g = B.inverse()*g*B
                    gens.append(g)
            else: #indefinite
                raise NotImplementedError(
                    "currently, we can only compute generators "
                    "for orthogonal groups over definite lattices.")
        deg = self.degree()
        base = self.ambient_vector_space().base_ring()
        inv_bil = self.inner_product_matrix()
        if is_finite:
            cat = Groups().Finite()
        else:
            cat = Groups()
        D = self.discriminant_group()
        G = GroupOfIsometries(deg,
                              base,
                              gens,
                              inv_bil,
                              category=cat,
                              invariant_submodule=self,
                              invariant_quotient_module=D)
        return G

    automorphisms=orthogonal_group

    def genus(self):
        r"""
        Return the genus of this lattice.

        EXAMPLES::

            sage: L = IntegralLattice("U")
            sage: L.genus()
            Genus of
            [0 1]
            [1 0]
            Signature:  (1, 1)
            Genus symbol at 2:    1^2
        """
        from sage.quadratic_forms.genera.genus import Genus
        return Genus(self.gram_matrix())

    def tensor_product(self, other, discard_basis=False):
        r"""
        Return the tensor product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- an integral lattice
        - ``discard_basis`` -- a boolean (default: ``False``). If ``True``, then the lattice
                            returned is equipped with the standard basis.

        EXAMPLES::

            sage: L = IntegralLattice("D3", [[1,-1,0], [0,1,-1]])
            sage: L1 = L.tensor_product(L)
            sage: L2 = L.tensor_product(L, True)
            sage: L1
            Lattice of degree 9 and rank 4 over Integer Ring
            Basis matrix:
            [ 1 -1  0 -1  1  0  0  0  0]
            [ 0  1 -1  0 -1  1  0  0  0]
            [ 0  0  0  1 -1  0 -1  1  0]
            [ 0  0  0  0  1 -1  0 -1  1]
            Inner product matrix:
            [ 4 -2 -2 -2  1  1 -2  1  1]
            [-2  4  0  1 -2  0  1 -2  0]
            [-2  0  4  1  0 -2  1  0 -2]
            [-2  1  1  4 -2 -2  0  0  0]
            [ 1 -2  0 -2  4  0  0  0  0]
            [ 1  0 -2 -2  0  4  0  0  0]
            [-2  1  1  0  0  0  4 -2 -2]
            [ 1 -2  0  0  0  0 -2  4  0]
            [ 1  0 -2  0  0  0 -2  0  4]
            sage: L1.gram_matrix()
            [ 36 -12 -12   4]
            [-12  24   4  -8]
            [-12   4  24  -8]
            [  4  -8  -8  16]
            sage: L2
            Lattice of degree 4 and rank 4 over Integer Ring
            Basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Inner product matrix:
            [ 36 -12 -12   4]
            [-12  24   4  -8]
            [-12   4  24  -8]
            [  4  -8  -8  16]
        """
        if not isinstance(other, FreeQuadraticModule_integer_symmetric):
            raise ValueError("other (=%s) must be an integral lattice" % other)
        if discard_basis:
            gram_matrix = self.gram_matrix().tensor_product(other.gram_matrix())
            return IntegralLattice(gram_matrix)
        else:
            inner_product_matrix = self.inner_product_matrix().tensor_product(other.inner_product_matrix())
            basis_matrix = self.basis_matrix().tensor_product(other.basis_matrix())
            n = self.degree()
            m = other.degree()
            ambient = FreeQuadraticModule(self.base_ring(), m * n, inner_product_matrix)
            return FreeQuadraticModule_integer_symmetric(ambient=ambient,
                                                         basis=basis_matrix,
                                                         inner_product_matrix=ambient.inner_product_matrix())

    def twist(self, s, discard_basis=False):
        r"""
        Return the lattice with inner product matrix scaled by ``s``.

        INPUT:

        - ``s`` -- a nonzero integer
        - ``discard_basis`` -- a boolean (default: ``False``).
          If ``True``, then the lattice returned is equipped
          with the standard basis.

        EXAMPLES::

            sage: L = IntegralLattice("A4")
            sage: L.twist(3)
            Lattice of degree 4 and rank 4 over Integer Ring
            Basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            Inner product matrix:
            [ 6 -3  0  0]
            [-3  6 -3  0]
            [ 0 -3  6 -3]
            [ 0  0 -3  6]
            sage: L = IntegralLattice(3,[[2,1,0],[0,1,1]])
            sage: L
            Lattice of degree 3 and rank 2 over Integer Ring
            Basis matrix:
            [2 1 0]
            [0 1 1]
            Inner product matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: L.twist(1)
            Lattice of degree 3 and rank 2 over Integer Ring
            Basis matrix:
            [2 1 0]
            [0 1 1]
            Inner product matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: L.twist(1, True)
            Lattice of degree 2 and rank 2 over Integer Ring
            Basis matrix:
            [1 0]
            [0 1]
            Inner product matrix:
            [5 1]
            [1 2]
        """
        try:
            s = self.base_ring()(s)
        except:
            ValueError("the scaling factor must be an element of the base ring.")
        if (s==0):
            raise ValueError("the scaling factor must be non zero")
        if discard_basis:
            return IntegralLattice(s * self.gram_matrix())
        else:
            n = self.degree()
            inner_product_matrix = s * self.inner_product_matrix()
            ambient = FreeQuadraticModule(self.base_ring(), n, inner_product_matrix)
            return FreeQuadraticModule_integer_symmetric(ambient=ambient, basis=self.basis(), inner_product_matrix=inner_product_matrix)

