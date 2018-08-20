r"""
Finite Dimensional Graded Lie Algebras With Basis

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.graded_modules import GradedModulesCategory

class FiniteDimensionalGradedLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    r"""
    Category of finite dimensional graded Lie algebras with a basis.

    A grading of a Lie algebra `\mathfrak{g}` is a direct sum decomposition
    `\mathfrak{g} = \bigoplus_{i} V_i` such that `[V_i,V_j] \subset V_{i+j}`.

    EXAMPLES::

        sage: C = LieAlgebras(ZZ).WithBasis().FiniteDimensional().Graded(); C
        Category of finite dimensional graded lie algebras with basis over Integer Ring
        sage: C.super_categories()
        [Category of graded lie algebras with basis over Integer Ring,
         Category of finite dimensional lie algebras with basis over Integer Ring]

        sage: C is LieAlgebras(ZZ).WithBasis().FiniteDimensional().Graded()
        True

    TESTS::

        sage: C = LieAlgebras(QQ).FiniteDimensional().WithBasis().Graded()
        sage: TestSuite(C).run()
    """
    class ParentMethods:
        def _test_grading(self, **options):
            r"""
            Tests that the Lie bracket respects the grading.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES::

                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: C = LieAlgebras(QQ).WithBasis().Graded()
                sage: C = C.FiniteDimensional().Stratified().Nilpotent()
                sage: L = NilpotentLieAlgebra(QQ, {('x','y'): {'z': 1}},
                ....:                         category=C)
                sage: L._test_grading()
                sage: L = NilpotentLieAlgebra(QQ, {('x','y'): {'x': 1}},
                ....:                         category=C)
                sage: L._test_grading()
                Traceback (most recent call last):
                ...
                AssertionError: Lie bracket [x, y] has degree 1, not degree 2

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)

            from sage.misc.misc import some_tuples
            for X,Y in some_tuples(self.basis(), 2, tester._max_runs):
                i = X.degree()
                j = Y.degree()
                Z = self.bracket(X, Y)
                if Z == 0:
                    continue
                tester.assertEquals(Z.degree(), i + j,
                    msg="Lie bracket [%s, %s] has degree %d, not degree %d " %
                        (X, Y, Z.degree(), i + j))
                tester.assertTrue(
                    Z.to_vector() in self.homogeneous_component_as_submodule(i + j),
                    msg="Lie bracket [%s, %s] is not in the "
                        "homogeneous component of degree %d" % (X, Y, i + j))

        @cached_method
        def homogeneous_component_as_submodule(self, d):
            r"""
            Return the ``d``-th homogeneous component of ``self``
            as a submodule.

            EXAMPLES::

                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: C = LieAlgebras(QQ).WithBasis().Graded()
                sage: C = C.FiniteDimensional().Stratified().Nilpotent()
                sage: L = NilpotentLieAlgebra(QQ, {('x','y'): {'z': 1}},
                ....:                         category=C)
                sage: L.homogeneous_component_as_submodule(2)
                Sparse vector space of degree 3 and dimension 1 over Rational Field
                Basis matrix:
                [0 0 1]
            """
            B = self.homogeneous_component_basis(d)
            return self.module().submodule([X.to_vector() for X in B])

    class Stratified(CategoryWithAxiom_over_base_ring):
        r"""
        Category of finite dimensional stratified Lie algebras with a basis.

        A stratified Lie algebra is a graded Lie algebra that is generated
        as a Lie algebra by its homogeneous component of degree 1.

        TESTS::

            sage: C = LieAlgebras(QQ).Graded().FiniteDimensional().WithBasis().Stratified()
            sage: TestSuite(C).run()
        """
        class ParentMethods:
            def _test_generated_by_degree_one(self, **options):
                r"""
                Tests that the Lie algebra is generated by the homogeneous component
                of degree one.

                INPUT:

                - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

                EXAMPLES::

                    sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                    sage: C = LieAlgebras(QQ).WithBasis().Graded()
                    sage: C = C.FiniteDimensional().Stratified().Nilpotent()
                    sage: sc = {('x','y'): {'z': 1}}
                    sage: L.<x,y,z> = NilpotentLieAlgebra(QQ, sc, category=C)
                    sage: L._test_generated_by_degree_one()

                We modify ``L`` with a different grading::

                    sage: L._basis_degrees = {x: 1, y: 2, z: 3}
                    sage: L.homogeneous_component_as_submodule.clear_cache()
                    sage: L._test_generated_by_degree_one()
                    Traceback (most recent call last):
                    ...
                    AssertionError: [x] does not generate Nilpotent Lie algebra on
                    3 generators (x, y, z) over Rational Field

                See the documentation for :class:`TestSuite` for more information.
                """
                tester = self._tester(**options)

                V1 = self.homogeneous_component_as_submodule(1)
                B1 = V1.basis()
                m = self.module()

                V = V1
                d = 0
                i = 0
                while V.dimension() > d:
                    if i > tester._max_runs:
                        return
                    B = V.basis()
                    d = V.dimension()
                    V = m.submodule(B + [self.bracket(X, Y).to_vector()
                                         for X in B1 for Y in B])

                tester.assertEqual(V, m,
                    msg="%s does not generate %s" % ([self(X) for X in B1], self))

            def degree_on_basis(self, m):
                r"""
                Return the degree of the basis element indexed by ``m``.

                If the degrees of the basis elements are not defined,
                they will be computed. By assumption the stratification
                `V_1 \oplus \cdots \oplus V_s` of ``self`` is such that each
                component `V_k` is spanned by some subset of the basis.

                The degree of a basis element `X` is therefore the largest
                index `k`such that `X \in V_k \oplus \cdots \oplus V_s`. The
                space  `V_k \oplus \cdots \oplus V_s` is by assumption the
                `k`-th term of the lower central series.

                EXAMPLES::

                    sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                    sage: C = LieAlgebras(QQ).WithBasis().Graded()
                    sage: C = C.FiniteDimensional().Stratified().Nilpotent()
                    sage: sc = {('X','Y'): {'Z': 1}}
                    sage: L.<X,Y,Z> = NilpotentLieAlgebra(QQ, sc, category=C)
                    sage: X.degree()
                    1
                    sage: Y.degree()
                    1
                    sage: L[X, Y]
                    Z
                    sage: Z.degree()
                    2
                """
                if not hasattr(self, '_basis_degrees'):
                    lcs = self.lower_central_series(submodule=True)
                    self._basis_degrees = {}

                    for k in reversed(range(len(lcs) - 1)):
                        for X in self.basis():
                            if X in self._basis_degrees:
                                continue
                            if X.to_vector() in lcs[k]:
                                self._basis_degrees[X] = k + 1

                m = self.basis()[m]
                return self._basis_degrees[m]

