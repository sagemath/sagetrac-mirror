r"""
Super Hopf algebras with basis
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory

class SuperHopfAlgebrasWithBasis(SuperModulesCategory):
    """
    The category of super Hopf algebras with a distinguished basis.

    EXAMPLES::

        sage: C = HopfAlgebras(ZZ).WithBasis().Super(); C
        Category of super hopf algebras with basis over Integer Ring
        sage: sorted(C.super_categories(), key=str)
        [Category of graded hopf algebras with basis over Integer Ring,
         Category of super algebras with basis over Integer Ring,
         Category of super coalgebras with basis over Integer Ring,
         Category of super hopf algebras over Integer Ring]

    TESTS::

        sage: TestSuite(C).run()
    """
    def extra_super_categories(self):
        return [self.base_category().Graded()]

    class ParentMethods:

        def _test_antipode(self, **options):
            r"""
            Test the antipode.

            An *antipode* `S` of a (super) Hopf algebra is a linear
            endomorphism of the Hopf algebra that satisfies the
            following conditions (see :wikipedia:`HopfAlgebra`).

            - If `\mu` and `\Delta` denote the product and coproduct of the
              Hopf algebra, respectively, then `S` satisfies

              .. MATH::

                  \mu \circ (S \tensor 1) \circ \Delta = unit \circ counit
                  \mu \circ (1 \tensor S) \circ \Delta = unit \circ counit

            - `S` is an *anti*-homomorphism:

               .. MATH::

                   S(ab) = (-1)^{\deg a \deg b} S(b) S(a)

               for homogeneous `a` and `b`.

            These properties are tested on :meth:`some_elements`.

            TESTS::

                sage: A = SteenrodAlgebra(7)
                sage: A._test_antipode()
            """
            tester = self._tester(**options)

            S = self.antipode

            IS = lambda x: self.sum(c * self.monomial(t1) * S(self.monomial(t2))
                                for ((t1, t2), c) in x.coproduct())

            SI = lambda x: self.sum(c * S(self.monomial(t1)) * self.monomial(t2)
                                for ((t1, t2), c) in x.coproduct())

            sign = lambda x, y: (-1)**(x.degree() * y.degree())

            for x in tester.some_elements():

                # antipode is an anti-homomorphism
                for y in tester.some_elements():
                    tester.assertTrue(S(x) * S(y) == sign(x, y) * S(y * x))

                # mu * (S # I) * delta == counit * unit
                tester.assertTrue(SI(x) == self.counit(x) * self.one())

                # mu * (I # S) * delta == counit * unit
                tester.assertTrue(IS(x) == self.counit(x) * self.one())

