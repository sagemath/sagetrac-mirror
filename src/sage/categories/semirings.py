r"""
Semirngs
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas Borie <nicolas.borie@math.u-psud.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from magmas_and_additive_magmas import MagmasAndAdditiveMagmas
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
#from sage.misc.lazy_attribute import lazy_attribute

class Semirings(CategoryWithAxiom):
    """
    The category of semirings.

    A semiring `(S,+,*)` is similar to a ring, but without the
    requirement that each element must have an additive inverse. In
    other words, it is a combination of a commutative additive monoid
    `(S,+)` and a multiplicative monoid `(S,*)`, where `*` distributes
    over `+`.

    .. SEEALSO: :wikipedia:`Semiring`

    EXAMPLES::

        sage: Semirings()
        Category of semirings
        sage: Semirings().super_categories()
        [Category of associative additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of monoids]

        sage: sorted(Semirings().axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveUnital', 'Associative', 'Distributive', 'Unital']

        sage: Semirings() is (CommutativeAdditiveMonoids() & Monoids()).Distributive()
        True

        sage: Semirings().AdditiveInverse()
        Category of rings


    TESTS::

        sage: TestSuite(Semirings()).run()
    """
    _base_category_class_and_axiom = (MagmasAndAdditiveMagmas.Distributive.AdditiveAssociative.AdditiveCommutative.AdditiveUnital.Associative, "Unital")

    class Differential(CategoryWithAxiom):
        class ParentMethods:
            @abstract_method
            def differential(self):
                """
                Return the differential of this semiring as a morphism.
                """

            def _test_differential(self, **options):
                """
                Verify that the differential of ``self`` satisfies the
                Leibniz rule.

                INPUT:

                - ``options`` -- any keyword arguments accepted
                  by :meth:`_tester`
                """
                tester = self._tester(**options)
                D = self.differential
                S = self.some_elements()
                for x in S:
                    for y in S:
                        tester.assertEquals(D(x*y), D(x)*y + x*D(y))

        class ElementMethods:
            def differential(self):
                """
                Return the differential of ``self``.
                """
                return self.parent().differential(self)

    class SubcategoryMethods:
        @cached_method
        def Differential(self):
            r"""
            Return the full subcategory of the objects of ``self`` equipped
            with a differential.

            A differential semiring `S` is a semiring equipped with at least
            one *derivation*, that is a linear morphism `\partial` satisfying
            the Leibniz rule:

            .. MATH:: \partial (xy) = (\partial x) y + x (\partial y)

            for all `x,y \in S`.

            .. SEEALSO:: :wikipedia:`Differential_algebra`

            EXAMPLES::

                sage: Semirings().Differential()
                Category of semigroups

            TESTS::

                sage: TestSuite(Semirings().Differential()).run()
                sage: Rings().Associative.__module__
                'sage.categories.semirings'
            """
            return self._with_axiom('Differential')

