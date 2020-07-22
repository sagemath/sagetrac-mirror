r"""
Additive monoids
"""
#*****************************************************************************
#  Copyright (C) 2013-2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_import import LazyImport
from sage.categories.category_with_axiom import CategoryWithAxiom_singleton
from sage.categories.additive_semigroups import AdditiveSemigroups
from sage.categories.homsets import HomsetsCategory
from .normed_additive_or_multiplicative_monoids import (NormedMonoidsCategory,
                                                        NormedAdditiveOrMultiplicativeMonoids)

class AdditiveMonoids(CategoryWithAxiom_singleton):
    """
    The category of additive monoids.

    An *additive monoid* is a unital :class:`additive semigroup
    <sage.categories.additive_semigroups.AdditiveSemigroups>`, that
    is a set endowed with a binary operation `+` which is associative
    and admits a zero (see :wikipedia:`Monoid`).

    EXAMPLES::

        sage: from sage.categories.additive_monoids import AdditiveMonoids
        sage: C = AdditiveMonoids(); C
        Category of additive monoids
        sage: C.super_categories()
        [Category of additive unital additive magmas, Category of additive semigroups]
        sage: sorted(C.axioms())
        ['AdditiveAssociative', 'AdditiveUnital']
        sage: from sage.categories.additive_semigroups import AdditiveSemigroups
        sage: C is AdditiveSemigroups().AdditiveUnital()
        True

    TESTS::

        sage: C.Algebras(QQ).is_subcategory(AlgebrasWithBasis(QQ))
        True
        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (AdditiveSemigroups, "AdditiveUnital")

    AdditiveCommutative = LazyImport('sage.categories.commutative_additive_monoids', 'CommutativeAdditiveMonoids', at_startup=True)
    AdditiveInverse = LazyImport('sage.categories.additive_groups', 'AdditiveGroups', at_startup=True)

    class ParentMethods:
        def sum(self, args):
            r"""
            Return the sum of the elements in ``args``, as an element
            of ``self``.

            INPUT:

            - ``args`` -- a list (or iterable) of elements of ``self``

            EXAMPLES::

                sage: S = CommutativeAdditiveMonoids().example()
                sage: (a,b,c,d) = S.additive_semigroup_generators()
                sage: S.sum((a,b,a,c,a,b))
                3*a + 2*b + c
                sage: S.sum(())
                0
                sage: S.sum(()).parent() == S
                True
            """
            return sum(args, self.zero())

    class Homsets(HomsetsCategory):

        def extra_super_categories(self):
            """
            Implement the fact that a homset between two monoids is
            associative.

            EXAMPLES::

                sage: from sage.categories.additive_monoids import AdditiveMonoids
                sage: AdditiveMonoids().Homsets().extra_super_categories()
                [Category of additive semigroups]
                sage: AdditiveMonoids().Homsets().super_categories()
                [Category of homsets of additive unital additive magmas, Category of additive monoids]

            .. TODO::

                This could be deduced from
                :meth:`AdditiveSemigroups.Homsets.extra_super_categories`.
                See comment in :meth:`Objects.SubcategoryMethods.Homsets`.
            """
            return [AdditiveSemigroups()]

    class SubcategoryMethods:
        def Normed(self):
            r"""
            Return the subcategory of the normed objects of ``self``.
            """
            return NormedMonoidsCategory.category_of(self)

    class Normed(NormedMonoidsCategory):
        def extra_super_categories(self):
            return [NormedAdditiveOrMultiplicativeMonoids()]

        class ParentMethods:

            def _test_norm_positive(self, **options):
                tester = self._tester(**options)
                S = tester.some_elements()
                o = self.zero()
                norm = self.norm_function()
                # Test additive neutral
                tester.assertEqual(norm(o), 0)
                # Test positivity
                for a in S:
                    d = norm(a)
                    if a != o:
                        tester.assertGreater(d, 0)
                    else:
                        tester.assertEqual(d, 0)

            def _test_norm_subadditive(self, **options):
                r"""
                Test that this normed additive monoid has a properly
                implemented norm.

                INPUT:

                - ``options`` -- any keyword arguments accepted
                  by :meth:`_tester`

                EXAMPLES:

                The "0-norm" is not a norm::

                """
                tester = self._tester(**options)
                S = tester.some_elements()
                o = self.zero()
                norm = self.norm_function()
                # Test triangle inequality
                for a in S:
                    for b in S:
                        tester.assertLessEqual(norm(a + b), norm(a) + norm(b))

