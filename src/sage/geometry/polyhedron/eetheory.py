r"""
Module for Equivariant Ehrhart Theory computations.

This module is centered around the computation of the equivariant `H^*` function
for a polytope `P` which is invariant under a linear group action using the
function :func:`Hstar_function`.
The equivariant Ehrhart series is the formal power series
`\sum_{m \geq 0} \chi_{mP}t^m`, where `\chi_{mP}` is
the permutation character of the group action on the lattice points in the
`m`th dilate of `P`. The equivariant Ehrhart series has the rational generating
function
`\sum_{m \geq 0} \chi_{mP}t^m = \frac{H^*(t)} {(1-t)det(Id-\rho(t)}.`
For details, see [1]_.
There are also functions to test whether the `H^*` series is polynomial,
:func:`is_polynomial`, and effective, :func:`is_effective`.

One can also use this module to compute fixed subpolytopes of a polytope under
the action of a group with the functions :func:`fixed_subpolytopes` and
:func:`fixed_subpolytope`.`
Finally, the module contains the function :func:`match_perms_to_mats` which
returns a dictionary between elements of the group written as permutations of a
polytope's vertices and written as a matrices.

REFERENCES:

    - [1] Alan Stapledon. Equivariant Ehrhart Theory. Advances in Mathematics 226 (2011), no. 4, 3622-3654
    - [2] Federico Ardila, Mariel Supina, and Andr\ ́{e}s R. Vindas-Mel\ ́{e}ndez, The equivariant ehrhart theory of the permutahedron, 2020.

EXAMPLES:

A first example follows Prop 6.1 of [1]_  which states that for a simplex,
`H^*` is given by the permutation character of the group on the box points at
each height. The one dimensional simplex [(-1,1),(1,1)] has two box points,
one at (0,0) and the other at (0,1). The group Z/2Z acts linearly on the simplex
by reflection across the y-axis. The box points are fixed under this action,
so we expect an `H^*` polynomial of `\chi_triv + \chi_triv t`::

    sage: simplex = Polyhedron(vertices = [[-1,1],[1,1]], backend = 'normaliz') # optional - pynormaliz
    sage: G = simplex.restricted_automorphism_group(output = 'permutation')     # optional - pynormaliz
    sage: G.order()   # optional - pynormaliz
    2
    sage: Hstar_function(simplex,G) # optional - pynormaliz
    chi_1*t + chi_1
    sage: G.character_table()       # optional - pynormaliz
    [ 1 -1]
    [ 1  1]

A second example is the action of the symmetric group on the 2-dimensional
permutahedron in 3-dimensional space, given by permuting the three basis
vectors. As shown in [2]_, the corresponding Hstar function is polynomial and
effective, while for all higher dimensional permutahedra under the symmetric
group action, it is not. The following computation agrees with the
`H^*`-polynomial recovered in [2]_::

    sage: p2 = polytopes.permutahedron(3, backend = 'normaliz')      # optional - pynormaliz
    sage: G = p2.restricted_automorphism_group(output='permutation') # optional - pynormaliz
    sage: H = G.subgroup(gens=[G.gens()[1],G.gens()[2]])             # optional - pynormaliz
    sage: H.order()                                                  # optional - pynormaliz
    6
    sage: [Hstar, Hlin] = [Hstar_function(p2,H), Hstar_function(p2,H, output = 'Hstar_as_lin_comb')]; Hstar # optional - pynormaliz
    chi_0*t^2 + (chi_0 + chi_1 + chi_2)*t + chi_0
    sage: H.character_table()        # optional - pynormaliz
    [ 1  1  1]
    [ 1 -1  1]
    [ 2  0 -1]
    sage: is_polynomial(Hstar)       # optional - pynormaliz
    True
    sage: is_effective(Hstar,Hlin)   # optional - pynormaliz
    True

AUTHORS:

- Sophia Elia (2021): Initial version
- Jean-Philippe Labbé (2021): Initial version
"""

##############################################################################
#     Copyright (C) 2021 Sophia Elia         <sophiae56 at math.fu-berlin.de>
#                   2021 Jean-Philippe Labbe <labbe at math.fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################


def Hstar_function(polytope, acting_group=None, output=None):
    r"""
    Return `H^*` as a rational function in `t` with coefficients in
    the ring of class functions of the ``acting_group`'
    of ``polytope``.

    Here, `\H^*(t) = \sum\_{m} \chi\_{m\text{``polytope``}}t^m det(Id-\rho(t))`.
    The irreducible characters of ``acting_group`` form an orthonormal basis
    for the ring (algebra) of class functions with values in `\mathbb C`.
    The coefficients of `H^*(t)` are expressed in this basis.

    INPUT:

    - ``polytope`` -- polyhedron object. A full dimensional lattice polytope
      with backend='normaliz'.

    - ``acting_group`` -- (default=None) a permgroup object. A subgroup of
      `polytope`'s `restricted_automorphism_group` output as a permutation.
       If ``None``, it is set to the full `restricted_automorphism_group`
       of `polytope`. The acting group should always use output = 'permutation'.

    - ``output`` -- string. an output option. The allowed values are:

        * ``None`` (default): returns the rational function `H^*(t)`. `H^*` is
          a rational function in `t` with coefficients in the ring of
          class functions.
        * ``'e_series_list'``: string. Returns a list of the ehrhart_series
          for the fixed_subpolytopes of each conjugacy class representative.
        * ``'determinant_vec'``: string. Returns a list of the determinants
          of `Id-\rho*t` for each conjugacy class representative.
        * ``'Hstar_as_lin_comb'``: string. Returns a vector of the coefficients
          of the irreducible representations in the expression of `H^*`.
        * ``'prod_det_es'``: string. Returns a vector of the product of
          determinants and the Ehrhart series.

    OUTPUT:

    The default output is the rational function `H^*`. `H^*` is a rational
    function in `t` with coefficients in the ring of class functions.
    There are several output options to see the intermediary outputs of the
    function.


    EXAMPLES:

    The `H^*`-polynomial of the standard `d-1` dimensional simplex
    `S = conv(e_1, \dots, e_d)` under its ``restricted_automorphism_group``
    is equal to 1 = `\chi_{trivial}` (Prop 6.1 [1]_).
    Here is the computation for the 3-dimensional standard simplex::

        sage: S = polytopes.simplex(3, backend = 'normaliz'); S              # optional - pynormaliz
        A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
        sage: G = S.restricted_automorphism_group(output = 'permutation'); G # optional - pynormaliz
        Permutation Group with generators [(2,3), (1,2), (0,1)]
        sage: len(G)                                                         # optional - pynormaliz
        24
        sage: phi = Hstar_function(S,G); phi                                 # optional - pynormaliz
        chi_4
        sage: G.character_table()                                            # optional - pynormaliz
        [ 1 -1  1  1 -1]
        [ 3 -1  0 -1  1]
        [ 2  0 -1  2  0]
        [ 3  1  0 -1 -1]
        [ 1  1  1  1  1]

    The next example is Example 7.6 in [1]_, and shows that `H^*` is not always
    a polynomial. Let P be the polytope with vertices
    `\pm(0,0,1),\pm(1,0,1), \pm(0,1,1), \pm(1,1,1)` and let
    G = `\mathbb Z / 2\mathbb Z` act on P as follows::

        sage: P = Polyhedron(vertices=[[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[0,1,1],
        ....: [0,-1,-1],[1,1,1],[-1,-1,-1]],backend='normaliz')           # optional - pynormaliz
        sage: K = P.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
        sage: G = K.subgroup(gens = [K[6]]); G                            # optional - pynormaliz
        Subgroup generated by [(0,2)(1,3)(4,6)(5,7)] of (Permutation Group with generators [(2,4)(3,5), (1,2)(5,6), (0,1)(2,3)(4,5)(6,7), (0,7)(1,3)(2,5)(4,6)])
        sage: conj_reps = G.conjugacy_classes_representatives()           # optional - pynormaliz
        sage: Dict = match_perms_to_mats(P,conj_reps, acting_group = G)   # optional - pynormaliz
        sage: list(Dict.keys())[0]                                        # optional - pynormaliz
        (0,2)(1,3)(4,6)(5,7)
        sage: list(Dict.values())[0]                                      # optional - pynormaliz
        [-1  0  1  0]
        [ 0  1  0  0]
        [ 0  0  1  0]
        [ 0  0  0  1]
        sage: len(G)                                                      # optional - pynormaliz
        2
        sage: G.character_table()                                         # optional - pynormaliz
        [ 1  1]
        [ 1 -1]

    Then we calculate the rational function `H^*(t)`::

        sage: Hst = Hstar_function(P,G); Hst     # optional - pynormaliz
        (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3 + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)

    To see the exact as written in [1]_, we can format it as
    ``'Hstar_as_lin_comb'``. The first coordinate is the coefficient of the
    trivial character; the second is the coefficient of the sign character::

        sage: lin = Hstar_function(P,G,output = 'Hstar_as_lin_comb'); lin  # optional - pynormaliz
        ((t^4 + 3*t^3 + 8*t^2 + 3*t + 1)/(t + 1), (3*t^3 + 2*t^2 + 3*t)/(t + 1))
    """
    # Setting the group
    G_perm = polytope.restricted_automorphism_group(output='permutation')

    if acting_group is not None:
        if not acting_group.is_subgroup(G_perm):
            raise TypeError("The 'acting_group' should be a subgroup of the 'restricted_automorphism_group'.")
        G_perm = acting_group
    # Create the Gap group one time only (each creation has different conj reps)
    G_perm_gap = G_perm._libgap_()

    # Fixing the conjugacy classes representatives once and for all
    cls = G_perm_gap.ConjugacyClasses()
    L = [cl.Representative() for cl in cls]
    conj_classes = [ConjugacyClassGAP(G_perm, G_perm.element_class(rep, G_perm, check=False)) for rep in L]
    conj_reps = [cl[0] for cl in conj_classes]
    # print(conj_reps)

    # Creating the Character Table
    n_classes = len(conj_reps)
    irrG_perm_gap = G_perm_gap.Irr()
    ct = [[irrG_perm_gap[i, j] for j in range(n_classes)] for i in range(n_classes)]
    from sage.rings.all import CyclotomicField
    e = irrG_perm_gap.Flat().Conductor()
    K = CyclotomicField(e)
    ct = [[K(x) for x in v] for v in ct]
    # Finally return the result as a matrix.
    from sage.matrix.all import MatrixSpace
    MS = MatrixSpace(K, n_classes)
    char_initial = MS(ct)

    # A check on whether the character table has permuted columns
    tbl = G_perm_gap.CharacterTable()
    perm = tbl.IdentificationOfConjugacyClasses()
    ident_perm = [i for i in range(1, 1 + n_classes)]
    if perm == ident_perm:
        pass
    else:
        raise ValueError("The conjugacy classes don't match with the character table")

    # Create fixed subpolytopes and their Ehrhart series
    group_dict = match_perms_to_mats(polytope, conj_reps, acting_group)
    fix_polys = fixed_subpolytopes(polytope, conj_reps)
    list_es = [fix_polys[g].ehrhart_series() for g in conj_reps]

    # get the list of the denominators det([Id - rho (t)])
    Ring = PolynomialRing(QQbar, 't')
    det_vector = list()
    dim = group_dict[G_perm.gens()[0]].dimensions()[0]
    t = Ring.gens()[0]
    ts_matrix = t * identity_matrix(Ring, dim)
    identity = identity_matrix(Ring, dim)

    # create a flag to fix the determinant if polytope isn't full dimensional
    flag = False
    if not polytope.is_full_dimensional():
        flag = True
        codim = polytope.ambient_dim() - polytope.dim()

    for perm in conj_reps:
        mat = group_dict[perm]
        mat = mat.change_ring(Ring)
        new_matrix = identity - mat*ts_matrix
        if flag:
            det = (1-t)**-codim*(new_matrix.determinant())
        else:
            det = new_matrix.determinant()
        det_vector.append(det)

    FF = Ring.fraction_field()
    initial_result = vector(FF, [a*b for a, b in zip(det_vector, list_es)])
    Char = char_initial.change_ring(FF)
    new_result = Char.solve_left(initial_result)

    new_new_result = _express_phi_as_polynomial_in_t(new_result)
    if output is None:
        return(new_new_result)
    elif output is 'e_series_list':
        return(list_es)
    elif output is 'determinant_vec':
        return(det_vector)
    elif output is 'Hstar_as_lin_comb':
        return(new_result)
    elif output is 'prod_det_es':
        return(initial_result)


def fixed_subpolytope(polytope, vertex_permutation):
    r"""
    Return the fixed subpolytope of ``polytope`` by the cyclic action of
    ``vertex_permutation``.

    INPUT:

    - ``polytope`` -- polyhedron object; a compact lattice polytope

    - ``vertex_permutation`` -- permutation; a permutation of the vertices of
      ``polytope``.

    OUTPUT:

    A subpolytope of ``polytope``.

    .. NOTE::

        The vertex_permutation is obtained as a permutation of the vertices
        represented as a permutation. For example, vertex_permutation =
        polytope.restricted_automorphism_group(output='permutation').

        Requiring a lattice polytope as opposed to a rational polytope as
        input is purely conventional.

    EXAMPLES:

    The fixed subpolytopes of the cube can be obtained as follows::

        sage: Cube = polytopes.cube(backend = 'normaliz')                   # optional - pynormaliz
        sage: AG = Cube.restricted_automorphism_group(output='permutation') # optional - pynormaliz
        sage: reprs = AG.conjugacy_classes_representatives()                # optional - pynormaliz

    The fixed subpolytope of the identity element of the group is the entire
    cube::
        sage: reprs[0]                                                      # optional - pynormaliz
        ()
        sage: fixed_subpolytope(Cube,reprs[0])                              # optional - pynormaliz
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8
        vertices
        sage: _.vertices()                                                  # optional - pynormaliz
        (A vertex at (-1, -1, -1),
         A vertex at (-1, -1, 1),
         A vertex at (-1, 1, -1),
         A vertex at (-1, 1, 1),
         A vertex at (1, -1, -1),
         A vertex at (1, -1, 1),
         A vertex at (1, 1, -1),
         A vertex at (1, 1, 1))

    You can obtain non-trivial examples::

        sage: fsp1 = fixed_subpolytope(Cube,reprs[8]);fsp1 # optional - pynormaliz
        A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex
        sage: fsp1.vertices()                              # optional - pynormaliz
        (A vertex at (0, 0, 0),)
        sage: fsp2 = fixed_subpolytope(Cube,reprs[3]);fsp2 # optional - pynormaliz
        A 2-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
        sage: fsp2.vertices()                              # optional - pynormaliz
        (A vertex at (-1, -1, 0),
        A vertex at (-1, 1, 0),
        A vertex at (1, -1, 0),
        A vertex at (1, 1, 0))

    The next example shows that fixed_subpolytope still works for rational polytopes::

       sage: P = Polyhedron(vertices = [[0,0],[3/2,0],[3/2,3/2],[0,3/2]], backend ='normaliz') # optional - pynormaliz
       sage: P.vertices() # optional - pynormaliz
       (A vertex at (0, 0),
        A vertex at (0, 3/2),
        A vertex at (3/2, 0),
        A vertex at (3/2, 3/2))
       sage: G = P.restricted_automorphism_group(output = 'permutation');G  # optional - pynormaliz
       Permutation Group with generators [(1,2), (0,1)(2,3), (0,3)]
       sage: len(G) # optional - pynormaliz
       8
       sage: G[2] # optional - pynormaliz
       (0,1)(2,3)
       sage: fixed_set = fixed_subpolytope(P,G[2]); fixed_set # optional - pynormaliz
       A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices
       sage: fixed_set.vertices() # optional - pynormaliz
       (A vertex at (0, 3/4), A vertex at (3/2, 3/4))
    """
    orbits = Set([Set(i) for i in vertex_permutation.cycle_tuples(singletons=True)])

    # If its the identity, returns the polytope
    if not orbits:
        return polytope

    # Make an index shift flag
    shift = True
    if 0 in vertex_permutation.domain():
        shift = False

    vertices = []
    for orbit in orbits:
        size = len(orbit)
        if shift:
            # in this case, the indices in the orbit are 1 more than the index in the V
            s = sum([(polytope.Vrepresentation()[i-1]).vector() for i in orbit])
        else:
            s = sum([(polytope.Vrepresentation()[i]).vector() for i in orbit])
        orbit_barycenter = (1/QQ(size)) * s
        vertices += [orbit_barycenter]

    return Polyhedron(vertices=vertices, backend='normaliz')


def fixed_subpolytopes(polytope, conj_class_reps):
    r"""
    Return the fixed subpolytope of an element in each conjugacy class of a
    given subgroup of the automorphism group of ``polytope``.

    INPUT:

    - ``polytope`` -- polyhedron object. a lattice polytope.

    - ``conj_class_reps`` -- a list of representatives of the conjugacy classes
      of the subgroup of the ``restricted_automorphism_group`` of the polytope.
      Each element is written as a permutation of the vertices of the polytope.

    OUTPUT:

    A dictionary with conj_class_reps as keys and the fixed subpolytopes
    as values.

    NOTE:

    Two elements in the same conjugacy class fix lattice-isomorphic
    subpolytopes.

    EXAMPLES:

    Here is an example for the square::

        sage: p = polytopes.hypercube(2, backend = 'normaliz'); p               # optional - pynormaliz
        A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
        sage: aut_p = p.restricted_automorphism_group(output = 'permutation')   # optional - pynormaliz
        sage: aut_p.order()                                                     # optional - pynormaliz
        8
        sage: conj_list = aut_p.conjugacy_classes_representatives(); conj_list  # optional - pynormaliz
        [(), (1,2), (0,1)(2,3), (0,1,3,2), (0,3)(1,2)]
        sage: fixed_subpolytopes(p,conj_list)                                   # optional - pynormaliz
        {(): A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices,
        (1,2): A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices,
        (0,1)(2,3): A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices,
        (0,1,3,2): A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex,
        (0,3)(1,2): A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex}
    """
    fixed_subpolytopes = {}

    for element in conj_class_reps:
        fixed_subpoly = fixed_subpolytope(polytope, element)
        fixed_subpolytopes[element] = fixed_subpoly
    return fixed_subpolytopes


def _express_phi_as_polynomial_in_t(initial_phi):
    r"""
    Rewrite the vector representing `H^*(t)` given as a linear combination of
    the irreducible representations as a rational function in `t`.

    The function `H^*` is first calculated as a linear combination of the irreducible
    representations of ``polytope.restricted_automorphism_group``. The convention is
    to express `H^*` as a rational function in `t` with coefficients in the ring of
    class functions.

    INPUT:

    - ``result`` -- a vector of rational functions in `t`.

    OUTPUT:

    A rational function in `t` with coefficients in the ring of class functions
    of ``polytope.restricted_automorphism_group()``.

    EXAMPLES:

    The expression of `H^*` as a polynomial in `t` for a 3-dimensional simplex
    is computed as follows::

        sage: simplex = Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],backend='normaliz') # optional - pynormaliz
        sage: phi = Hstar_function(simplex); phi # optional - pynormaliz
        chi_4

    The polynomial is `\chi_4 \cdot t^0`. We can see which irreducible
    representation `\chi_4` corresponds to by looking at the character table::

        sage: G = simplex.restricted_automorphism_group(output='permutation') # optional - pynormaliz
        sage: char = G.character_table();char # optional - pynormaliz
        [ 1 -1  1  1 -1]
        [ 3 -1  0 -1  1]
        [ 2  0 -1  2  0]
        [ 3  1  0 -1 -1]
        [ 1  1  1  1  1]

    Thus `\chi_4` corresponds to the trivial representation of the group, and
    for every element in the group, it evaluates to 1.

    As another example, we can look at `H^*(t)` for the `\pm 1` square::

        sage: square = Polyhedron(vertices = [[1,1],[-1,1],[-1,-1],[1,-1]], backend ='normaliz') # optional - pynormaliz
        sage: phi = Hstar_function(square) ; phi       # optional - pynormaliz
        chi_0*t^2 + (2*chi_0 + chi_2 + chi_3 + chi_4)*t + chi_0

    Plugging in the values from the first column of the character table below
    yields the `h^*`-polynomial of the square, `t^2+6t+1`::

        sage: G = square.restricted_automorphism_group(output='permutation') # optional - pynormaliz
        sage: G.character_table() # optional - pynormaliz
        [ 1  1  1  1  1]
        [ 1 -1 -1  1  1]
        [ 1 -1  1 -1  1]
        [ 1  1 -1 -1  1]
        [ 2  0  0  0 -2]
    """
    chi_vars = var(','.join('chi_{}'.format(i) for i in range(len(initial_phi))))
    Chi_ring = PolynomialRing(QQbar, chi_vars)
    virtual_ring = PolynomialRing(Chi_ring, initial_phi.base_ring().gens())
    fraction_virtual_ring = virtual_ring.fraction_field()
    new_result = initial_phi.change_ring(fraction_virtual_ring)*vector(fraction_virtual_ring, chi_vars)
    return new_result


def match_perms_to_mats(polytope, conj_class_reps, acting_group=None, additional_elts=None):
    r"""
    An element of ``polytope``'s ``restricted_autormorphism_group`` may
    be represented either as a permutation of the vertices of ``polytope``
    or as a matrix. This function returns a dictionary with permutations
    as keys and matrices as values.

    When ``additional_elts`` is ``None``, the dictionary is returned for the
    generators and conjugacy classes representatives in conj_class_reps of the
    ``restricted_automorphism_group`` or the ``acting_group``.
    When ``additional_elts`` is not ``None``, each element in
    ``additional_elts`` also becomes a key.

    INPUT:

    - ``polytope`` -- polyhedron object. a lattice polytope.

    - ``conj_class_reps`` -- list. A list of representatives of the conjugacy
      classes of the ``restricted_automorphism_group``.

    - ``acting_group`` -- a subgroup of the ``polytope``'s
      ``restricted_automorphism_group``.

    - ``additional_elts`` -- list (default=None). a subset of the
      ``restricted_automorphism_group`` of ``polytope`` expressed as
      permutations.

    OUTPUT:

    A dictionary between elements of ``the restricted_automorphism_group``
    or ``acting_group`` expressed as permutations (keys) and matrices (values).

    EXAMPLES:

    This example shows the dictionary between permutations and matrices
    for the generators of the ``restricted_automorphism_group`` of the
    `\pm 1` 2-dimensional square::

        sage: square = Polyhedron(vertices=[[1,1],[-1,1],[-1,-1],[1,-1]], backend='normaliz') # optional - pynormaliz
        sage: aut_square = square.restricted_automorphism_group(output = 'permutation')       # optional - pynormaliz
        sage: conj_reps = aut_square.conjugacy_classes_representatives()                      # optional - pynormaliz
        sage: gens_dict = match_perms_to_mats(square,conj_reps); gens_dict                    # optional - pynormaliz
        {(): [1 0 0]
         [0 1 0]
         [0 0 1],
         (1,2): [0 1 0]
         [1 0 0]
         [0 0 1],
         (0,1)(2,3): [ 1  0  0]
         [ 0 -1  0]
         [ 0  0  1],
         (0,1,3,2): [ 0  1  0]
         [-1  0  0]
         [ 0  0  1],
         (0,3): [ 0 -1  0]
         [-1  0  0]
         [ 0  0  1],
         (0,3)(1,2): [-1  0  0]
         [ 0 -1  0]
         [ 0  0  1]}
        sage: square.vertices() # optional - pynormaliz
        (A vertex at (-1, -1),
        A vertex at (-1, 1),
        A vertex at (1, -1),
        A vertex at (1, 1))

    This example tests the functionality for extra elements::

        sage: C = polytopes.cross_polytope(2)
        sage: G = C.restricted_automorphism_group(output = 'permutation')
        sage: conj_reps = G.conjugacy_classes_representatives()
        sage: add_elt = [G[6]]; add_elt
        [(0,2,3,1)]
        sage: match_perms_to_mats(C,conj_reps,additional_elts = add_elt)
        {(): [1 0 0]
         [0 1 0]
         [0 0 1],
         (1,2): [ 1  0  0]
         [ 0 -1  0]
         [ 0  0  1],
         (0,1)(2,3): [0 1 0]
         [1 0 0]
         [0 0 1],
         (0,1,3,2): [ 0 -1  0]
         [ 1  0  0]
         [ 0  0  1],
         (0,2,3,1): [ 0  1  0]
         [-1  0  0]
         [ 0  0  1],
         (0,3): [-1  0  0]
         [ 0  1  0]
         [ 0  0  1],
         (0,3)(1,2): [-1  0  0]
         [ 0 -1  0]
         [ 0  0  1]}
    """
    V = [v.homogeneous_vector() for v in polytope.Vrepresentation()]
    Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()
    Vplus = list(matrix(V) * Qplus)
    W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

    G = polytope.restricted_automorphism_group(output='permutation')
    if acting_group is not None:
        G = acting_group

    group_dict = {}

    for perm in G.gens():
        group_dict[perm] = _match_perm(perm, V, Vplus, W)

    for perm in conj_class_reps:
        group_dict[perm] = _match_perm(perm, V, Vplus, W)

    if additional_elts is not None:
        for perm in additional_elts:
            group_dict[perm] = _match_perm(perm, V, Vplus, W)
    return group_dict


def _match_perm(permutation, V, Vplus, W):
    r"""
    Return the matrix representation of a permutation in the
    ``restricted_autormorphism_group`` of ``polytope``.

    INPUT:

    - ``polytope`` -- polyhedron object. A lattice polytope.

    - ``V`` -- list. a list of vectors from the ``match_perms_to_mats`` function.

    - ``Vplus`` -- list. from the ``match_perms_to_mats`` function.

    - ``W`` -- matrix. from the ``match_perms_to_mats`` function.

    OUTPUT:

    A matrix that acts in the same way on the polytope as the ``permutation``.

    EXAMPLES:

    This example shows a reflection across the x-axis::

        sage: cross = polytopes.cross_polytope(2, backend = 'normaliz') # optional - pynormaliz
        sage: cross.vertices() # optional - pynormaliz
        (A vertex at (-1, 0),
         A vertex at (0, -1),
         A vertex at (0, 1),
         A vertex at (1, 0))
        sage: G = cross.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
        sage: flip = G.gens()[0]; flip   # optional - pynormaliz
        (1,2)
        sage: V = [v.homogeneous_vector() for v in cross.Vrepresentation()]   # optional - pynormaliz
        sage: Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()    # optional - pynormaliz
        sage: Vplus = list(matrix(V) * Qplus)   # optional - pynormaliz
        sage: W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V))) # optional - pynormaliz
        sage: _match_perm(flip, V, Vplus, W)   # optional - pynormaliz
        [ 1  0  0]
        [ 0 -1  0]
        [ 0  0  1]
    """
    A = sum(V[permutation(i)].column() * Vplus[i].row() for i in range(len(V)))
    return A + W


def is_polynomial(Hstar):
    r"""
    Checks if the equivariant `H^*`-series is a polynomial.

    INPUT:

    - ``Hstar`` -- a rational function in `t` with coefficients in the ring of
                    class functions.

    OUTPUT:

    Boolean. Whether the `H^*` series is a polynomial.

    EXAMPLES:

    Example of a simplex under symmetric group action::

        sage: S = polytopes.simplex(3); S
        A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
        sage: G = S.restricted_automorphism_group(output = 'permutation');
        sage: len(G)
        24
        sage: Hstar_simplex = Hstar_function(S,G); Hstar_simplex # optional - pynormaliz
        chi_4
        sage: is_polynomial(Hstar_simplex) # optional - pynormaliz
        True

    The following is example 7.6 in Stapledon::

        sage: P = Polyhedron(vertices=[[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[0,1,1],
        ....: [0,-1,-1],[1,1,1],[-1,-1,-1]],backend='normaliz')           # optional - pynormaliz
        sage: G = P.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
        sage: H = G.subgroup(gens = [G[6]])                               # optional - pynormaliz
        sage: Hstar = Hstar_function(P,H); Hstar                          # optional - pynormaliz
        (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3 + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)
        sage: is_polynomial(Hstar)                  # optional - pynormaliz
        False
    """
    flag = True
    if Hstar.denominator() == 1:
        pass
    else:
        flag = False
    return flag


def is_effective(Hstar, Hstar_as_lin_comb):
    r"""
    Check if the `H^*` series is effective.

    The `H^*` series is effective if it is a polynomial and the coefficients
    of each $t^i$ are effective representations. The coefficients of each
    irreducible representation must be non-negative. The `H^*` series must be
    polynomial in order for it to be effective.

    INPUT:

    - ``Hstar`` -- a rational function in `t` with coefficients in the ring of
                   class functions.

    OUTPUT:

    Boolean. Whether the `H^*` series is effective.

    EXAMPLES:

    The `H^*` series of the two-dimensional permutahedron under the action
    of the symmetric group is effective::

        sage: p2 = polytopes.permutahedron(3, backend = 'normaliz')      # optional - pynormaliz
        sage: G = p2.restricted_automorphism_group(output='permutation') # optional - pynormaliz
        sage: H = G.subgroup(gens=[G.gens()[1],G.gens()[2]])             # optional - pynormaliz
        sage: H.order()                                                  # optional - pynormaliz
        6
        sage: [Hstar, Hlin] = [Hstar_function(p2,H), Hstar_function(p2,H, output = 'Hstar_as_lin_comb')] # optional - pynormaliz
        sage: is_effective(Hstar,Hlin)   # optional - pynormaliz
        True

    The `H^*` series must be a polynomial in order to be effective. If it is
    not polynomial, a value error is returned::

        sage: P = Polyhedron(vertices=[[0,0,1],[0,0,-1],[1,0,1],[-1,0,-1],[0,1,1],
        ....: [0,-1,-1],[1,1,1],[-1,-1,-1]],backend='normaliz')           # optional - pynormaliz
        sage: G = P.restricted_automorphism_group(output = 'permutation') # optional - pynormaliz
        sage: H = G.subgroup(gens = [G[6]])                               # optional - pynormaliz
        sage: Hstar = Hstar_function(P,H); Hstar                          # optional - pynormaliz
        (chi_0*t^4 + (3*chi_0 + 3*chi_1)*t^3 + (8*chi_0 + 2*chi_1)*t^2 + (3*chi_0 + 3*chi_1)*t + chi_0)/(t + 1)
        sage: Hstar_lin = Hstar_function(P,H, output = 'Hstar_as_lin_comb') # optional - pynormaliz
        sage: is_effective(Hstar, Hstar_lin)  # optional - pynormaliz
        Traceback (most recent call last):
        ...
        ValueError: The Hstar vector must be polynomial
    """
    if not is_polynomial(Hstar):
        raise ValueError("The Hstar vector must be polynomial")
    flag = True
    for irrep in range(len(Hstar_as_lin_comb)):
        coeffs = Hstar_as_lin_comb[irrep].numerator().coefficients()
        for i in coeffs:
            if i < 0:
                flag = False
    return flag
