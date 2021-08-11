r"""
The Varchenko-Gelfand ring and the Moseley ring of a central hyperplane
arrangement.

EXAMPLES:

(The examples are written for the Varchenko-Gelfand ring below, but it works
similarly for the Moseley ring.)

We construct the Varchenko-Gelfand ring of the reflection arrangement of
type `A_2`::

    sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
    sage: A = hyperplane_arrangements.braid(3)
    sage: VG = VarchenkoGelfandRing(QQ, A)
    sage: VG
    Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over Rational Field

In order to work with the elements of this ring, we first need to pick
a basis. There are currently three bases implemented.

- *NBC Basis*. This is the basis consisting of the monomials corresponding to
  the NBC sets of the underlying matroid::

    sage: N = VG.nbc_basis()
    sage: N
    Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over the Rational Field on the NBC basis

  Often, we will want to loop over the elements of the basis. For example,
  we loop through the basis elements and convert each to polynomials
  in the underlying polynomial ring::

    sage: for b in N.basis():
    ....:     print("{} --> {}".format(b, b.to_polynomial()))
    N[] --> 1
    N[0] --> x0
    N[1] --> x1
    N[2] --> x2
    N[0, 1] --> x0*x2 + x1*x2 - x2
    N[0, 2] --> x0*x2

  The elements N[0], N[1], N[2] correspond to the monomials x0, x1, x2 in the
  ring. As a conventient shorthand for defining an element, you can do the
  following::

    sage: N[0,2]
    N[0, 2]

  If the input is not an NBC set, then the corresponding element will be
  convert to an element in the NBC basis::

    sage: N[0,1,2]              # this is x0 * x1 * x2 in the ring
    N[0, 1]
    sage: N[0] * N[1] * N[2]    # this is x0 * x1 * x2 in the ring
    N[0, 1]
    sage: N[1] * N[2]           # this is x1 * x2 in the ring, which is not part of the NBC basis
    N[0, 1] + N[2] - N[0, 2]
    sage: N.one()               # writing N[] doesn't work
    N[]

- *Covector basis*. This is the basis corresponding to the functions that take
  value 1 on a specific chamber and 0 on all other chambers (via the map that
  sends the variable `x_i` to the function that takes value
  1 on all chambers in the positive half-space of the `i`-th hyperplane and
  -1 on all chambers in the negative half-space of the `i`-th hyperplane.

  Since the chambers are easily indexed by covectors, this basis is indexed by
  covectors::

    sage: C = VG.covector_basis()
    sage: for b in C.basis():
    ....:     print("{} --> {}".format(b, b.to_polynomial()))
    C[1, 1, 1] --> x0*x2 + x1*x2 - x2
    C[1, -1, 1] --> -x1*x2 + x2
    C[1, -1, -1] --> -x0*x2 + x0
    C[-1, 1, 1] --> -x0*x2 + x2
    C[-1, 1, -1] --> -x1*x2 + x1
    C[-1, -1, -1] --> x0*x2 + x1*x2 - x0 - x1 - x2 + 1

- *Change of basis*. You can also easily change between representations in the
  different bases. For example, to convert an element `a` to the NBC basis, use
  `N(a)`::

    sage: N(C[1, -1, 1])
    -N[0, 1] + N[0, 2]

  Here we expand all the elements of the covector basis in the NBC basis::

    sage: for b in C.basis():
    ....:     print("{} --> {}".format(b, N(b)))
    C[1, 1, 1] --> N[0, 1]
    C[1, -1, 1] --> -N[0, 1] + N[0, 2]
    C[1, -1, -1] --> N[0] - N[0, 2]
    C[-1, 1, 1] --> N[2] - N[0, 2]
    C[-1, 1, -1] --> -N[0, 1] - N[2] + N[1] + N[0, 2]
    C[-1, -1, -1] --> N[0, 1] - N[1] + N[] - N[0]

  And the other way:

    sage: for b in N.basis():
    ....:     print("{} --> {}".format(b, C(b)))
    N[] --> C[-1, -1, -1] + C[-1, 1, -1] + C[-1, 1, 1] + C[1, -1, -1] + C[1, -1, 1] + C[1, 1, 1]
    N[0] --> C[1, -1, -1] + C[1, -1, 1] + C[1, 1, 1]
    N[1] --> C[-1, 1, -1] + C[-1, 1, 1] + C[1, 1, 1]
    N[2] --> C[-1, 1, 1] + C[1, -1, 1] + C[1, 1, 1]
    N[0, 1] --> C[1, 1, 1]
    N[0, 2] --> C[1, -1, 1] + C[1, 1, 1]

  This looks good since the identity element is ``N[]``, which is the sum of
  all the indicator functions ranging over all the chambers of the arrangement.

- *Normal Basis*. This basis is used for internally storing the underlying
  polynomial ring. The basis itself is what is computed from Sage and depends
  on various factors, including the monomial order on the polynomial ring.

  It is used to have a specific realization of this ring, and all other bases
  are defined in relation to this basis.

    sage: x = VG.normal_basis()
    sage: for b in x.basis():
    ....:     print("{} --> {}".format(b, b.to_polynomial()))
    X[1, 2] --> x1*x2
    X[0, 2] --> x0*x2
    X[2] --> x2
    X[1] --> x1
    X[0] --> x0
    X[] --> 1

- *New bases*. If you want to define a new basis, then you can define it in
  relation to any of the above bases. The simplest way is to define a method
  called ``to_polynomial_on_basis`` that converts a key of a basis element to
  a polynomial in the underlying polynomial ring. See the class ``NBC`` for an
  example; in particular, ``to_polynomial_on_basis`` is defined for the NBC
  basis by mapping an NBC set `S` to the product of `x_i` for `i` in `S` (for
  an explicit example, see the method below in which `S` is called ``nbc``).


EXAMPLES (MoseleyRing):

These are all the examples from above, but with the MoseleyRing in place of the
VarchenkoGelfandRing (and without the commentary)::

    sage: from sage.algebras.varchenko_gelfand import MoseleyRing
    sage: A = hyperplane_arrangements.braid(3)
    sage: MR = MoseleyRing(QQ, A)
    sage: MR
    Moseley ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over Rational Field

    sage: N = MR.nbc_basis()
    sage: for b in N.basis():
    ....:     print("{} --> {}".format(b, b.to_polynomial()))
    N[] --> 1
    N[0] --> x0
    N[1] --> x1
    N[2] --> x2
    N[0, 1] --> x0*x2 + x1*x2
    N[0, 2] --> x0*x2

    sage: N[0,2]
    N[0, 2]
    sage: N[0,1,2]              # this is x0 * x1 * x2 in the ring
    0
    sage: N[0] * N[1] * N[2]    # this is x0 * x1 * x2 in the ring
    0
    sage: N[1] * N[2]           # this is x1 * x2 in the ring, which is not part of the NBC basis
    N[0, 1] - N[0, 2]
    sage: N.one()               # writing N[] doesn't work
    N[]

    sage: C = MR.covector_basis()
    sage: for b in C.basis():
    ....:     print("{} --> {}".format(b, b.to_polynomial()))
    C[1, 1, 1] --> x0*x2 + x1*x2
    C[1, -1, 1] --> -x1*x2 + x2
    C[1, -1, -1] --> -x0*x2 + x0
    C[-1, 1, 1] --> -x0*x2 + x2
    C[-1, 1, -1] --> -x1*x2 + x1
    C[-1, -1, -1] --> x0*x2 + x1*x2 - x0 - x1 + 1

    sage: for b in C.basis():
    ....:     print("{} --> {}".format(b, N(b)))
    C[1, 1, 1] --> N[0, 1]
    C[1, -1, 1] --> -N[0, 1] + N[0, 2] + N[2]
    C[1, -1, -1] --> N[0] - N[0, 2]
    C[-1, 1, 1] --> N[2] - N[0, 2]
    C[-1, 1, -1] --> -N[0, 1] + N[0, 2] + N[1]
    C[-1, -1, -1] --> N[0, 1] - N[0] - N[1] + N[]

    sage: for b in N.basis():
    ....:     print("{} --> {}".format(b, C(b)))
    N[] --> C[-1, -1, -1] + C[-1, 1, -1] + C[1, -1, -1]
    N[0] --> -1/2*C[-1, 1, 1] + C[1, -1, -1] + 1/2*C[1, -1, 1] + 1/2*C[1, 1, 1]
    N[1] --> C[-1, 1, -1] + 1/2*C[-1, 1, 1] - 1/2*C[1, -1, 1] + 1/2*C[1, 1, 1]
    N[2] --> 1/2*C[-1, 1, 1] + 1/2*C[1, -1, 1] + 1/2*C[1, 1, 1]
    N[0, 1] --> C[1, 1, 1]
    N[0, 2] --> -1/2*C[-1, 1, 1] + 1/2*C[1, -1, 1] + 1/2*C[1, 1, 1]

    sage: x = MR.normal_basis()
    sage: for b in x.basis():
    ....:     print("{} --> {}".format(b, b.to_polynomial()))
    X[1, 2] --> x1*x2
    X[0, 2] --> x0*x2
    X[2] --> x2
    X[1] --> x1
    X[0] --> x0
    X[] --> 1

AUTHORS:

- Franco Saliola, Galen Dorpalen-Barry (2021): initial version

"""

# ****************************************************************************
#       Copyright (C) 2021 Franco Saliola <saliola.franco at uqam.ca>
#                          Galen Dorpalen-Barry
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.parent import Parent
from sage.misc.bindable_class import BindableClass
from sage.categories.all import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.structure.unique_representation import UniqueRepresentation
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.sets.set import Set
from sage.categories.realizations import Category_realization_of_parent
from sage.categories.rings import Rings
from sage.categories.fields import Fields

class VarchenkoGelfandRing(UniqueRepresentation, Parent):
    r"""
    The Varchenko-Gelfand ring of a central hyperplane arrangement.
    """
    def __init__(self, base_ring, arrangement, homogenized=False):
        """
        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: A = hyperplane_arrangements.braid(3); A
            Arrangement <t1 - t2 | t0 - t1 | t0 - t2>
            sage: VG = VarchenkoGelfandRing(QQ, A); VG
            Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over Rational Field
        """
        # test that the arrangement is central
        if not arrangement.is_central():
            raise ValueError("the hyperplane arrangement must be central")

        # initiate the parent object by specifying the category and that we
        # will be defining multiple bases (WithRealizations)
        if not (base_ring in Fields() or base_ring in Rings()):
            raise ValueError(f"{base_ring=} must be a Ring")
        self._base_ring = base_ring
        Parent.__init__(self, category=AlgebrasWithBasis(base_ring).WithRealizations())

        # save data for later use
        self._arrangement = arrangement
        self._homogenized = homogenized
        self._matroid = self._arrangement.matroid()

        # register coercions (how to convert between the bases)
        self._register_coercions()

    def _register_coercions(self):
        N = self.nbc_basis()
        C = self.covector_basis()
        X = self.normal_basis()
        N.module_morphism(on_basis=lambda I : C.from_polynomial(N.to_polynomial_on_basis(I)), codomain=C).register_as_coercion()
        N.module_morphism(on_basis=lambda I : X.from_polynomial(N.to_polynomial_on_basis(I)), codomain=X).register_as_coercion()
        C.module_morphism(on_basis=lambda I : C.from_polynomial(C.to_polynomial_on_basis(I)), codomain=C).register_as_coercion()
        C.module_morphism(on_basis=lambda I : X.from_polynomial(C.to_polynomial_on_basis(I)), codomain=X).register_as_coercion()
        X.module_morphism(on_basis=lambda I : N.from_polynomial(X.to_polynomial_on_basis(I)), codomain=N).register_as_coercion()
        X.module_morphism(on_basis=lambda I : C.from_polynomial(X.to_polynomial_on_basis(I)), codomain=C).register_as_coercion()

    def _repr_(self):
        r"""
        String representation of this object.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: VG = VarchenkoGelfandRing(QQ, A)
            sage: VG._repr_()
            'Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over Rational Field'
        """
        return "Varchenko-Gelfand ring of {} over {}".format(self.hyperplane_arrangement(), self.base_ring())

    def a_realization(self):
        r"""
        Returns the default realization.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: VG = VarchenkoGelfandRing(QQ, A)
            sage: VG.a_realization()
            Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over the Rational Field on the Covector basis

        """
        return self.covector_basis()

    def base_ring(self):
        r"""

        """

        return self._base_ring

    def hyperplane_arrangement(self):
        r"""
        The hyperplane arrangement defining this Varchenko-Gelfand ring.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: VG = VarchenkoGelfandRing(QQ, A)
            sage: VG.hyperplane_arrangement()
            Arrangement <t1 - t2 | t0 - t1 | t0 - t2>
        """
        return self._arrangement

    @cached_method
    def underlying_polynomial_ring(self):
        r"""
        The Varchenko-Gelfand ring as a polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: VG = VarchenkoGelfandRing(QQ, A)
            sage: VG.underlying_polynomial_ring()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field by the ideal (x0^2 - x0, x1^2 - x1, x2^2 - x2, -x0*x1 + x0*x2 + x1*x2 - x2)

        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        A = self._arrangement
        n = A.n_hyperplanes()

        var_string = ["x{}".format(i) for i in range(n)]
        if self._homogenized:
            var_string.append("z0")

        R = PolynomialRing(self.base_ring(), var_string, len(var_string))
        if self._homogenized:
            x = R.gens()[:-1]
            z0 = R.gens()[-1]
        else:
            x = R.gens()
            z0 = R.one()

        # Build the ideal that we will quotient by.
        # First, there is one generator of the form xi^2 - xi for each variable xi.
        ideal_gens = [xi * xi - xi * z0 for xi in x]
        if self._homogenized:
            ideal_gens.append(z0 * z0 - 1)

        # Second, there is one generator for each empty intersection of the form
        #   `\cap_{i \in I} H_i^+ \cap \cap_{j \in J} H_j^- = \emptyset`,
        # which is given by`\prod_{I}
        #   `\prod_{i \in I} x_i \prod_{j \in J} (x_j - 1) - \prod_{i \in I} (x_i - 1) \prod_{j \in J} x_j`.
        # It suffices to take one such generator for each circuit of the associated
        # matroid.

        for circuit in self._matroid.circuits():
            circuit = sorted(circuit)
            m = matrix([A[i].normal() for i in circuit])
            for v in m.left_kernel().basis():
                I = [circuit[i] for (i, vi) in enumerate(v) if vi > 0]
                J = [circuit[j] for (j, vj) in enumerate(v) if vj < 0]
                gen = R.prod(x[i] for i in I) * R.prod(x[j] - z0 for j in J) - R.prod(x[i] - z0 for i in I) * R.prod(x[j] for j in J)
                ideal_gens.append(gen)

        I = R.ideal(ideal_gens)
        Q = R.quotient(I, names=R.variable_names())
        return Q

    @cached_method
    def change_of_basis_matrix(self, basis0, basis1):
        r"""
        The change of basis matrix from ``basis0`` to ``basis1``.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: VG = VarchenkoGelfandRing(QQ, A)
            sage: C = VG.covector_basis()
            sage: N = VG.nbc_basis()
            sage: X = VG.normal_basis()
            sage: VG.change_of_basis_matrix(X, X)
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            sage: VG.change_of_basis_matrix(C, X)
            [ 1  1 -1  0  0  0]
            [-1  0  1  0  0  0]
            [ 0 -1  0  0  1  0]
            [ 0 -1  1  0  0  0]
            [-1  0  0  1  0  0]
            [ 1  1 -1 -1 -1  1]
            sage: VG.change_of_basis_matrix(N, X)
            [ 0  0  0  0  0  1]
            [ 0  0  0  0  1  0]
            [ 0  0  0  1  0  0]
            [ 0  0  1  0  0  0]
            [ 1  1 -1  0  0  0]
            [ 0  1  0  0  0  0]

        """
        return matrix([basis1.from_polynomial(basis0.to_polynomial_on_basis(key)).to_vector()
                                                        for key in basis0.basis().keys()])

    def covectors_in_cone(self, cone_mask):
        r"""
        The covectors of the regions that lie in a cone.

        The cone is specified by prescribing which entries of the sign vectors
        must be positive or negative.

        INPUT:

        - ``cone_mask`` -- list with entries in `\{1, 0, -1\}`, where the entry
          in position `i` specifies whether we should consider the positive
          (+1), the negative (-1), or neither (0) halfspace associated with the
          `i`-th hyperplane.

        OUTPUT:

        - ``cone_covectors`` -- list of covectors

        *WARNING* When constructing a hyperplane arrangement, the order in
        which the hyperplanes are specified is not the order with which the
        hyperplanes are stored internally, or the order with which the sign
        vectors are computed. Instead, they are sorted according to their
        (implicitly defined) normal vector.

        EXAMPLES:

        Construct a hyperplane arrangement::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H(-x + z, z, x - 5*y, x - y, x + y - 2*z, x + y)
            sage: VG = VarchenkoGelfandRing(QQ, A)

        To compute the covectors of the regions lying in the cone defined by
        the halfspaces `z > 0`, `x - y > 0`, `x + y > 0`, identify the
        positions of the hyperplanes in ``A.hyperplanes()``::

            sage: A.hyperplanes()
            (Hyperplane -x + 0*y + z + 0,
             Hyperplane 0*x + 0*y + z + 0,
             Hyperplane x - 5*y + 0*z + 0,
             Hyperplane x - y + 0*z + 0,
             Hyperplane x + y - 2*z + 0,
             Hyperplane x + y + 0*z + 0)

        So, we are interested in the positive halfspaces of the hyperplanes in
        positions 1, 3, 5, so we set ``cone_mask=(0,1,0,1,0,1)``:

            sage: cone_covectors = VG.covectors_in_cone((0,1,0,1,0,1))
            sage: cone_covectors
            [(1, 1, 1, 1, -1, 1),
             (1, 1, -1, 1, -1, 1),
             (-1, 1, 1, 1, 1, 1),
             (-1, 1, 1, 1, -1, 1),
             (-1, 1, -1, 1, 1, 1),
             (-1, 1, -1, 1, -1, 1)]

        Here are the polynomials corresponding to heaviside function for each
        of these regions::

            sage: C = VG.covector_basis()
            sage: heavisides = [C.to_polynomial_on_basis(covector) for covector in cone_covectors]
            sage: heavisides
            [x0*x2*x5,
             -x0*x2*x5 - x0*x4*x5 - x3*x4*x5 + x0*x5 + x3*x5 + x4*x5 - x5,
             x1*x2*x5 + x2*x4*x5 - x2*x5,
             -x0*x2*x5 - x2*x4*x5 + x2*x5,
             -x1*x2*x5 + x1*x3*x5 - x2*x4*x5 + x3*x4*x5 + x2*x5 - x3*x5,
             x0*x2*x5 + x0*x4*x5 + x2*x4*x5 - x0*x5 - x2*x5 - x4*x5 + x5]

        Here are the same elements, but expressed in terms of the covector basis::

            sage: [C.from_polynomial(heaviside) for heaviside in heavisides]
            [C[1, 1, 1, 1, -1, 1],
             C[1, 1, -1, 1, -1, 1],
             C[-1, 1, 1, 1, 1, 1],
             C[-1, 1, 1, 1, -1, 1],
             C[-1, 1, -1, 1, 1, 1],
             C[-1, 1, -1, 1, -1, 1]]

        Here are the same elements, but expressed in terms of the nbc basis::

            sage: N = VG.nbc_basis()
            sage: [N.from_polynomial(heaviside) for heaviside in heavisides]
            [N[0, 2, 5],
             N[0, 3] - N[0, 2, 3],
             N[0, 1] - N[2] - N[0, 1, 2] - N[0, 1, 4] - N[1] + N[0, 2] - N[4] - N[0, 2, 4] + N[] + N[1, 4] + N[0, 4] + N[1, 2] - N[0] + N[2, 4],
             N[1, 5] + N[2] + N[0, 1, 4] - N[0, 2] + N[4] + N[0, 2, 4] - N[] - N[1, 4] - N[0, 4] + N[0] - N[0, 1, 5] - N[2, 4],
             N[2] + N[0, 1, 2] + N[1, 3] - N[0, 2] + N[4] + N[0, 2, 4] - N[] - N[0, 4] - N[1, 2] + N[0] - N[0, 1, 3] - N[2, 4],
             -N[2] + N[0, 2] - N[4] - N[0, 2, 4] + N[] + N[0, 4] - N[0] + N[2, 4]]

        Test that the ideal generated by the heaviside functions is equal to
        the ideal generated by the product of the hyperplanes defining the
        cone::

            sage: R = VG.underlying_polynomial_ring()
            sage: x = R.gens()
            sage: R.ideal(heavisides) == R.ideal(x[1] * x[3] * x[5])
            True

        The product `x_1 x_3 x_4` defines a heaviside function which is `0` on
        regions outside the cone (this can be seen by plotting the arrangement
        in the hyperplane `z = 1`). Hence, this element should belong to the
        above ideal. Let's vertify this::

            sage: x[1] * x[3] * x[4] in R.ideal(heavisides)
            True

        """
        A = self.hyperplane_arrangement()
        cone_covectors = []
        for region in A.regions():
            covector = A.sign_vector(region.representative_point())
            if all(cone_mask[i] / covector[i] >= 0 for i in range(len(covector))):
                cone_covectors.append(tuple(covector))
        return cone_covectors

    def rref_basis_of_cone_ideal(self, cone_mask):
        r"""

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = x | x - y | x - z | y | y - z
            sage: VG = VarchenkoGelfandRing(QQ, A)

            sage: VG.rref_basis_of_cone_ideal((1, 0, 1, 0, 0))
            [N[0, 2], N[0, 4] - N[0, 1, 4], N[0, 1, 2]]

            sage: VG.rref_basis_of_cone_ideal((0, 1, 0, 0, 0))
            [N[1],
             N[0, 1],
             N[1, 2],
             N[1, 3],
             N[1, 4],
             -N[0, 3, 4] + N[0, 4],
             N[0, 1, 2],
             N[0, 1, 3],
             N[0, 1, 4]]

        """
        cone_covectors = self.covectors_in_cone(cone_mask)
        covec_to_poly = self.covector_basis().to_polynomial_on_basis
        heavisides = [covec_to_poly(covector) for covector in cone_covectors]
        N = self.nbc_basis()
        heavisides_in_nbc_basis = [N.from_polynomial(heaviside) for heaviside in heavisides]
        m = []
        for element in heavisides_in_nbc_basis:
            m.append(element.to_vector())
        m = matrix(m)
        return [N.from_vector(row) for row in m.rref()]

    class Bases(Category_realization_of_parent):
        r"""
        The category of the realizations. This contains code that will be
        inherited by all the different bases.
        """
        def super_categories(self):
            category = AlgebrasWithBasis(self.base().base_ring())
            return [self.base().Realizations(),
                    category.Realizations().WithBasis()]

        class ParentMethods:
            def underlying_polynomial_ring(self):
                return self.realization_of().underlying_polynomial_ring()

            def hyperplane_arrangement(self):
                return self.realization_of().hyperplane_arrangement()

            def _repr_(self):
                return "Varchenko-Gelfand ring of {} over the {} on the {} basis".format(self.hyperplane_arrangement(), self.base_ring(), self._realization_name())

            @cached_method
            def product_on_basis(self, x, y):
                r"""
                Product of two basis elements.

                This is achieved by converting the elements to polynomials in
                the underlying polynomial ring and then expanding the result in
                the basis.

                EXAMPLES::

                    sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                    sage: A = hyperplane_arrangements.braid(3)
                    sage: VG = VarchenkoGelfandRing(QQ, A)
                    sage: X = VG.normal_basis()
                    sage: X[0, 1] * X[1]
                    X[0, 2] - X[2] + X[1, 2]
                    sage: X[0] * X[2]
                    X[0, 2]

                    sage: A = hyperplane_arrangements.braid(3)
                    sage: VG = VarchenkoGelfandRing(QQ, A)
                    sage: N = VG.nbc_basis()
                    sage: poly1 = N[0,2].to_polynomial(); poly1
                    x0*x2
                    sage: poly2 = N[1].to_polynomial(); poly2
                    x1
                    sage: poly1 * poly2
                    x0*x2 + x1*x2 - x2
                    sage: N.product_on_basis(Set([0,2]), Set([1]))
                    N[0, 1]
                    sage: _.to_polynomial()
                    x0*x2 + x1*x2 - x2

                A bigger example for sanity checking::

                    sage: for b0 in N.basis():
                    ....:     for b1 in N.basis():
                    ....:         print("{} * {} = {}".format(b0, b1, b0 * b1))
                    N[] * N[] = N[]
                    N[] * N[0] = N[0]
                    N[] * N[1] = N[1]
                    N[] * N[2] = N[2]
                    N[] * N[0, 1] = N[0, 1]
                    N[] * N[0, 2] = N[0, 2]
                    N[0] * N[] = N[0]
                    N[0] * N[0] = N[0]
                    N[0] * N[1] = N[0, 1]
                    N[0] * N[2] = N[0, 2]
                    N[0] * N[0, 1] = N[0, 1]
                    N[0] * N[0, 2] = N[0, 2]
                    N[1] * N[] = N[1]
                    N[1] * N[0] = N[0, 1]
                    N[1] * N[1] = N[1]
                    N[1] * N[2] = N[0, 1] + N[2] - N[0, 2]
                    N[1] * N[0, 1] = N[0, 1]
                    N[1] * N[0, 2] = N[0, 1]
                    N[2] * N[] = N[2]
                    N[2] * N[0] = N[0, 2]
                    N[2] * N[1] = N[0, 1] + N[2] - N[0, 2]
                    N[2] * N[2] = N[2]
                    N[2] * N[0, 1] = N[0, 1]
                    N[2] * N[0, 2] = N[0, 2]
                    N[0, 1] * N[] = N[0, 1]
                    N[0, 1] * N[0] = N[0, 1]
                    N[0, 1] * N[1] = N[0, 1]
                    N[0, 1] * N[2] = N[0, 1]
                    N[0, 1] * N[0, 1] = N[0, 1]
                    N[0, 1] * N[0, 2] = N[0, 1]
                    N[0, 2] * N[] = N[0, 2]
                    N[0, 2] * N[0] = N[0, 2]
                    N[0, 2] * N[1] = N[0, 1]
                    N[0, 2] * N[2] = N[0, 2]
                    N[0, 2] * N[0, 1] = N[0, 1]
                    N[0, 2] * N[0, 2] = N[0, 2]

                """
                poly1 = self.to_polynomial_on_basis(x)
                poly2 = self.to_polynomial_on_basis(y)
                return self.from_polynomial(poly1 * poly2)

            def to_polynomial(self, element):
                r"""
                EXAMPLES::

                    sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                    sage: A = hyperplane_arrangements.braid(3)
                    sage: VG = VarchenkoGelfandRing(QQ, A)
                    sage: X = VG.normal_basis()
                    sage: X.to_polynomial(X[0, 2] - X[0])
                    x0*x2 - x0

                    sage: N = VG.nbc_basis()
                    sage: N.to_polynomial(N[0, 2] - N[0])
                    x0*x2 - x0

                    sage: C = VG.normal_basis()
                    sage: C.to_polynomial(C[1, -1, 1] - C[1, 1, 1])
                    x1*x2 - x1

                """
                R = self.underlying_polynomial_ring()
                return R.sum(coeff * self.to_polynomial_on_basis(monom) for (monom, coeff) in element)

            def from_polynomial(self, polynomial):
                r"""
                Express a polynomial in the Varchenko-Gelfand ring as a linear
                combination of the elements of the nbc basis.

                EXAMPLES::

                    sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                    sage: A = hyperplane_arrangements.braid(3)
                    sage: VG = VarchenkoGelfandRing(QQ, A)

                Convert the elements of the nbc basis to polynomials::

                    sage: N = VG.nbc_basis()
                    sage: for b in N.basis():
                    ....:     poly = b.to_polynomial()
                    ....:     print("{} = {}".format(N.from_polynomial(poly), poly))
                    N[] = 1
                    N[0] = x0
                    N[1] = x1
                    N[2] = x2
                    N[0, 1] = x0*x2 + x1*x2 - x2
                    N[0, 2] = x0*x2

                Convert the elements of the normal basis to polynomials::

                    sage: X = VG.normal_basis()
                    sage: for b in X.basis():
                    ....:     poly = b.to_polynomial()
                    ....:     print("{} = {}".format(X.from_polynomial(poly), poly))
                    X[1, 2] = x1*x2
                    X[0, 2] = x0*x2
                    X[2] = x2
                    X[1] = x1
                    X[0] = x0
                    X[] = 1

                Convert the elements of the nbc basis to elements of the normal
                basis::

                    sage: for b in N.basis():
                    ....:     print("{} = {}".format(b, X.from_polynomial(b.to_polynomial())))
                    N[] = X[]
                    N[0] = X[0]
                    N[1] = X[1]
                    N[2] = X[2]
                    N[0, 1] = -X[2] + X[0, 2] + X[1, 2]
                    N[0, 2] = X[0, 2]

                and conversely::

                    sage: for b in X.basis():
                    ....:     print("{} = {}".format(b, N.from_polynomial(b.to_polynomial())))
                    X[1, 2] = N[0, 1] - N[0, 2] + N[2]
                    X[0, 2] = N[0, 2]
                    X[2] = N[2]
                    X[1] = N[1]
                    X[0] = N[0]
                    X[] = N[]

                """
                X = self.realization_of().normal_basis()
                M = self.realization_of().change_of_basis_matrix(self, X)
                v = X.from_polynomial(polynomial).to_vector()
                return self.from_vector(M.solve_left(v))

        class ElementMethods:
            def to_polynomial(self):
                r"""
                The polynomial corresponding to this element.

                EXAMPLES::

                    sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                    sage: A = hyperplane_arrangements.braid(3)
                    sage: VG = VarchenkoGelfandRing(QQ, A)
                    sage: X = VG.normal_basis()
                    sage: X[0,2].to_polynomial()
                    x0*x2
                    sage: X[0,1,2].to_polynomial()
                    x0*x2 + x1*x2 - x2

                    sage: N = VG.nbc_basis()
                    sage: N[0, 1].to_polynomial()
                    x0*x2 + x1*x2 - x2

                    sage: C = VG.covector_basis()
                    sage: C[1, -1, 1].to_polynomial()
                    -x1*x2 + x2

                """
                return self.parent().to_polynomial(self)

    class Normal(CombinatorialFreeModule, BindableClass):
        def __init__(self, VG):
            r"""
            Return a basis of the Varchenko-Gelfand polynomial ring consisting
            of monomials.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: M = VG.normal_basis()
                sage: M
                Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over the Rational Field on the Normal basis

            """
            PR = VG.underlying_polynomial_ring()
            normal_basis = PR.defining_ideal().normal_basis()
            basis_keys = [Set([i for (i, j) in enumerate(monom.exponents()[0]) if j]) for monom in normal_basis]
            CombinatorialFreeModule.__init__(self,
                                             VG.base_ring(),
                                             basis_keys,
                                             category=VG.Bases(),
                                             prefix='X')


        def __getitem__(self, indices):
            """
            This method implements a conventient shorthand for defining an
            element of this basis.

            Given a list of indices corresponding to the monomials `i_1, i_2,
            ...`, we compute the product of the generators of the polynomial
            ring `x_{i_1} x_{i_2} ...`, and express the result in the normal
            basis of the defining ideal.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: X = VG.normal_basis()
                sage: X[0,2]
                X[0, 2]

            Note that not all monomials are part of the normal monomial basis;
            explicitly, the monomial `x_0 x_1 x_2` does not appear in the basis
            for the above Varchenko-Gelfand::

                sage: R = VG.underlying_polynomial_ring()
                sage: R.defining_ideal().normal_basis()
                [x1*x2, x0*x2, x2, x1, x0, 1]

            However, it still does define an element of this ring,
            obtained by expressing this element in the normal basis::

                sage: X[0,1,2]
                -X[2] + X[0, 2] + X[1, 2]
                sage: x = R.gens()
                sage: x[0] * x[1] * x[2]
                x0*x2 + x1*x2 - x2

            """
            if isinstance(indices, (int, Integer)):
                indices = [indices]
            R = self.underlying_polynomial_ring()
            x = R.gens()
            poly = R.prod(x[i] for i in indices)
            return self.from_polynomial(poly)

        def _repr_term(self, m):
            r"""
            The string representation of a single basis element.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: X = VarchenkoGelfandRing(QQ, A).normal_basis()
                sage: X._repr_term([0,1])
                'X[0, 1]'
            """
            return self.prefix() + '[' + ", ".join(str(i) for i in sorted(m)) + ']'

        @cached_method
        def one_basis(self):
            r"""
            The multiplicative identity element of the algebra.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: X = VG.normal_basis()
                sage: X.one()
                X[]
                sage: X.one_basis()
                {}

            """
            return Set([])


        def to_polynomial_on_basis(self, indices):
            r"""
            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: X = VG.normal_basis()
                sage: X.to_polynomial_on_basis([0,2])
                x0*x2
            """
            R = self.underlying_polynomial_ring()
            x = R.gens()
            return R.prod(x[i] for i in indices)

        def from_polynomial(self, poly):
            r"""
            Create an element from a polynomial in the underlying polynomial
            ring.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: R = VG.underlying_polynomial_ring()
                sage: x = R.gens()
                sage: x[0] * x[1]
                x0*x2 + x1*x2 - x2

                sage: X = VG.normal_basis()
                sage: X.from_polynomial(x[0] * x[1])
                -X[2] + X[0, 2] + X[1, 2]

            """
            poly = poly.lift()
            d = {}
            for monom in poly.monomials():
                coeff = poly.monomial_coefficient(monom)
                monom = Set([i for (i, j) in enumerate(monom.exponents()[0]) if j])
                d[monom] = coeff

            return self._from_dict(d)

        def _polynomial_to_normal_basis_vector(self, polynomial):
            r"""
            Express a polynomial in the Varchenko-Gelfand ring as a linear
            combination of the elements of the monomial basis.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: X = VG.normal_basis()
                sage: [X._polynomial_to_normal_basis_vector(b.to_polynomial()) for b in X.basis()]
                [(1, 0, 0, 0, 0, 0),
                 (0, 1, 0, 0, 0, 0),
                 (0, 0, 1, 0, 0, 0),
                 (0, 0, 0, 1, 0, 0),
                 (0, 0, 0, 0, 1, 0),
                 (0, 0, 0, 0, 0, 1)]
                sage: N = VG.nbc_basis()
                sage: [X._polynomial_to_normal_basis_vector(b.to_polynomial()) for b in N.basis()]
                [(0, 0, 0, 0, 0, 1),
                 (0, 0, 0, 0, 1, 0),
                 (0, 0, 0, 1, 0, 0),
                 (0, 0, 1, 0, 0, 0),
                 (1, 1, -1, 0, 0, 0),
                 (0, 1, 0, 0, 0, 0)]
            """
            poly = polynomial.lift()
            return vector([poly.monomial_coefficient(b.to_polynomial().lift()) for b in self.basis()])

    normal_basis = Normal

    class NBC(CombinatorialFreeModule, BindableClass):
        def __init__(self, VG):
            r"""
            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: VG.nbc_basis()
                Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over the Rational Field on the NBC basis

            """
            basis_keys = [Set(nbc) for nbc in VG._matroid.no_broken_circuits_sets()]
            CombinatorialFreeModule.__init__(self,
                                             VG.base_ring(),
                                             basis_keys,
                                             category=VG.Bases(),
                                             prefix='N')

        def __getitem__(self, indices):
            """
            This method implements a conventient shorthand for defining an
            element of this basis.

            Given a list of indices corresponding to the monomials `i_1, i_2,
            ...`, we compute the product of the generators of the polynomial
            ring `x_{i_1} x_{i_2} ...`, and express the result in the normal
            basis of the defining ideal.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: N = VG.nbc_basis()
                sage: N[0,2]
                N[0, 2]

            If the input is not an NBC set, then the corresponding element will
            be convert to an element in the NBC basis::

                sage: N[0,1,2]
                N[0, 1]
                sage: N[0] * N[1] * N[2]
                N[0, 1]

            """
            if isinstance(indices, (int, Integer)):
                indices = [indices]
            R = self.underlying_polynomial_ring()
            x = R.gens()
            poly = R.prod(x[i] for i in indices)
            return self.from_polynomial(poly)


        def _repr_term(self, m):
            r"""
            The string representation of a single basis element.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: N = VarchenkoGelfandRing(QQ, A).nbc_basis()
                sage: N._repr_term([0,1])
                'N[0, 1]'

            """
            return self.prefix() + '[' + ", ".join(str(i) for i in sorted(m)) + ']'

        @cached_method
        def one(self):
            r"""
            The multiplicative identity element of the algebra.
            """
            return self.monomial(Set([]))

        def to_polynomial_on_basis(self, nbc):
            r"""
            The polynomial associated with the no broken circuit set ``nbc``.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: N = VG.nbc_basis()
                sage: for nbc in A.matroid().no_broken_circuits_sets():
                ....:     poly = N.to_polynomial_on_basis(nbc)
                ....:     print("{:>6}: {}".format(str(Set(nbc)), poly))
                    {}: 1
                   {0}: x0
                   {1}: x1
                   {2}: x2
                {0, 1}: x0*x2 + x1*x2 - x2
                {0, 2}: x0*x2
            """
            PR = self.underlying_polynomial_ring()
            x = PR.gens()
            return PR.prod(x[i] for i in nbc)

    nbc_basis = NBC

    class Covector(CombinatorialFreeModule, BindableClass):
        def __init__(self, VG):
            r"""
            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: VG.covector_basis()
                Varchenko-Gelfand ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over the Rational Field on the Covector basis

            """
            A = VG.hyperplane_arrangement()
            self._region_to_covector = dict((region, tuple(A.sign_vector(region.representative_point()))) for region in A.regions())
            self._covector_to_region = dict((value, key) for (key, value) in self._region_to_covector.items())
            basis_keys = [self._region_to_covector[region] for region in A.regions()]
            CombinatorialFreeModule.__init__(self,
                                             VG.base_ring(),
                                             basis_keys=basis_keys,
                                             category=VG.Bases(),
                                             prefix='C')

        def __getitem__(self, covector):
            """
            This method implements a conventient shorthand for defining an
            element of this basis.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: C = VG.covector_basis()
                sage: C[-1, 1, 1]
                C[-1, 1, 1]

            TESTS::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: C = VG.covector_basis()
                sage: C[-1, -1, 1]
                Traceback (most recent call last):
                ...
                ValueError: input is not a covector of a region of the hyperplane arrangement

            """
            if covector not in self.basis().keys():
                raise ValueError("input is not a covector of a region of the hyperplane arrangement")
            return self.monomial(covector)

        def _repr_term(self, m):
            r"""
            The string representation of a single basis element.

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: C = VarchenkoGelfandRing(QQ, A).covector_basis()
                sage: C._repr_term([-1,1,-1])
                'C[-1, 1, -1]'

            """
            return self.prefix() + '[' + ", ".join(str(i) for i in m) + ']'

        @cached_method
        def one(self):
            r"""
            The multiplicative identity element of the algebra.
            """
            return NotImplemented

        def to_polynomial_on_basis(self, covector):
            r"""
            The polynomial associated with a chamber (where the chamber is
            expressed via its covector).

            EXAMPLES::

                sage: from sage.algebras.varchenko_gelfand import VarchenkoGelfandRing
                sage: A = hyperplane_arrangements.braid(3)
                sage: VG = VarchenkoGelfandRing(QQ, A)
                sage: C = VG.covector_basis()
                sage: for covector in C.basis().keys():
                ....:     poly = C.to_polynomial_on_basis(covector)
                ....:     print("{:>12}: {}".format(str(covector), poly))
                   (1, 1, 1): x0*x2 + x1*x2 - x2
                  (1, -1, 1): -x1*x2 + x2
                 (1, -1, -1): -x0*x2 + x0
                  (-1, 1, 1): -x0*x2 + x2
                 (-1, 1, -1): -x1*x2 + x1
                (-1, -1, -1): x0*x2 + x1*x2 - x0 - x1 - x2 + 1

                sage: N = VG.nbc_basis()
                sage: for covector in sorted(C.basis().keys()):
                ....:     poly = C.to_polynomial_on_basis(covector)
                ....:     print("{:>12}: {}".format(str(covector), N.from_polynomial(poly)))
                (-1, -1, -1): N[0, 1] - N[0] - N[1] + N[]
                 (-1, 1, -1): -N[0, 1] + N[0, 2] + N[1] - N[2]
                  (-1, 1, 1): -N[0, 2] + N[2]
                 (1, -1, -1): - N[0, 2] + N[0]
                  (1, -1, 1): -N[0, 1] + N[0, 2]
                   (1, 1, 1): N[0, 1]
            """
            if 0 in covector:
                raise(ValueError, "covector cannot contain 0")

            PR = self.underlying_polynomial_ring()
            x = PR.gens()

            region = self._covector_to_region[covector]
            A = self.hyperplane_arrangement()
            H = A.parent()
            hyperplanes = A.hyperplanes()

            polynomial = PR.one()
            for ieq in region.Hrepresentation():
                wall = H(ieq).hyperplanes()[0]
                if wall in hyperplanes:
                    i = hyperplanes.index(wall)
                    polynomial *= x[i]
                else:
                    i = hyperplanes.index(-wall)
                    polynomial *= (1-x[i])

            return polynomial

    covector_basis = Covector

class MoseleyRing(VarchenkoGelfandRing):
    r"""
    The Moseley ring of a central hyperplane arrangement.
    """
    def __init__(self, base_ring, arrangement):
        """
        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import MoseleyRing
            sage: A = hyperplane_arrangements.braid(3); A
            Arrangement <t1 - t2 | t0 - t1 | t0 - t2>
            sage: MR = MoseleyRing(QQ, A); MR
            Moseley ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over Rational Field
        """

        # test that the arrangement is central
        if not arrangement.is_central():
            raise ValueError("the hyperplane arrangement must be central")

        # initiate the parent object by specifying the category and that we
        # will be defining multiple bases (WithRealizations)
        if not (base_ring in Fields() or base_ring in Rings()):
            raise ValueError(f"{base_ring=} must be a Ring")
        self._base_ring = base_ring
        Parent.__init__(self, category=AlgebrasWithBasis(base_ring).WithRealizations())

        # save data for later use
        self._arrangement = arrangement
        self._matroid = self._arrangement.matroid()

        # register coercions (how to convert between the bases)
        self._register_coercions()

    def _repr_(self):
        r"""
        String representation of this object.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import MoseleyRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: MoseleyRing(QQ, A)._repr_()
            'Moseley ring of Arrangement <t1 - t2 | t0 - t1 | t0 - t2> over Rational Field'
        """
        return "Moseley ring of {} over {}".format(self.hyperplane_arrangement(), self.base_ring())

    @cached_method
    def underlying_polynomial_ring(self):
        r"""
        The Varchenko-Gelfand ring as a polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import MoseleyRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: MoseleyRing(QQ, A).underlying_polynomial_ring()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2 over Rational Field by the ideal (x0^2, x1^2, x2^2, x0*x1 - x0*x2 - x1*x2)

        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        A = self._arrangement
        n = A.n_hyperplanes()

        var_string = ["x{}".format(i) for i in range(n)]

        R = PolynomialRing(self.base_ring(), var_string, len(var_string))
        x = R.gens()

        # Build the ideal that we will quotient by.
        # First, there is one generator of the form xi^2 - xi for each variable xi.
        ideal_gens = [xi * xi for xi in x]

        # Second, there is one generator for each empty intersection of the form
        #   `\cap_{i \in I} H_i^+ \cap \cap_{j \in J} H_j^- = \emptyset`,
        # which is given by`\prod_{I}
        #   `\prod_{i \in I} x_i \prod_{j \in J} (x_j - 1) - \prod_{i \in I} (x_i - 1) \prod_{j \in J} x_j`.
        # It suffices to take one such generator for each circuit of the associated
        # matroid.
        # Unlike the usual VG ring, we take the leading term of the difference.

        for circuit in self._matroid.circuits():
            circuit = sorted(circuit)
            m = matrix([A[i].normal() for i in circuit])
            for v in m.left_kernel().basis():
                I = [circuit[i] for (i, vi) in enumerate(v) if vi > 0]
                J = [circuit[j] for (j, vj) in enumerate(v) if vj < 0]
                gen = R.sum(R.prod(x[i] for i in I) * R.prod(x[j] for j in J if j != k) for k in J) \
                    - R.sum(R.prod(x[i] for i in I if i != k) * R.prod(x[j] for j in J) for k in I)
                ideal_gens.append(gen)

        I = R.ideal(ideal_gens)
        Q = R.quotient(I, names=R.variable_names())
        return Q

    def hilbert_series(self):
        r"""
        The Hilbert series of the polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import MoseleyRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: MR = MoseleyRing(QQ, A)
            sage: MR.hilbert_series()
            2*t^2 + 3*t + 1

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H(x, y, z, x+y, x+z, y+z, x-y, x-z, y-z)
            sage: MR = MoseleyRing(QQ, A)
            sage: MR.hilbert_series()
            15*t^3 + 23*t^2 + 9*t + 1

        """
        I = self.underlying_polynomial_ring().defining_ideal()
        return I.hilbert_series()

    def hilbert_polynomial(self):
        r"""
        The Hilbert series of the polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import MoseleyRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: MR = MoseleyRing(QQ, A)
            sage: MR.hilbert_polynomial()
            0

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H(x, y, z, x+y, x+z, y+z, x-y, x-z, y-z)
            sage: MR = MoseleyRing(QQ, A)
            sage: MR.hilbert_polynomial()
            0

        """
        I = self.underlying_polynomial_ring().defining_ideal()
        return I.hilbert_polynomial()

    def hilbert_numerator(self):
        r"""
        The Hilbert series of the polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.varchenko_gelfand import MoseleyRing
            sage: A = hyperplane_arrangements.braid(3)
            sage: MR = MoseleyRing(QQ, A)
            sage: MR.hilbert_numerator()
            -2*t^5 + 3*t^4 + 2*t^3 - 4*t^2 + 1

            sage: H.<x,y,z> = HyperplaneArrangements(QQ)
            sage: A = H(x, y, z, x+y, x+z, y+z, x-y, x-z, y-z)
            sage: MR = MoseleyRing(QQ, A)
            sage: MR.hilbert_numerator()
            -15*t^12 + 112*t^11 - 342*t^10 + 512*t^9 - 273*t^8 - 288*t^7 + 588*t^6 - 384*t^5 + 63*t^4 + 48*t^3 - 22*t^2 + 1

        """
        I = self.underlying_polynomial_ring().defining_ideal()
        return I.hilbert_numerator()
