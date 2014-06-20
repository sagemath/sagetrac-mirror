r"""
Nonassociative Rngs
"""
#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.misc.lazy_import import LazyImport
from magmas_and_additive_magmas import MagmasAndAdditiveMagmas

class NonassociativeRngs(CategoryWithAxiom):
    """
    The category of nonassocitiave rngs.

    A nonassociative rng `(R,+,*)` is a rng (a ring without requiring
    a multiplicative identity), but without the requirement that each
    elements be associative. In other words, it is a combination of an
    abelian group `(R, +)` and a multiplicative magma `(R, *)`,
    where `*` distributes over `+`.

    .. SEEALSO: :wikipedia:`Nonassociative_ring`

    EXAMPLES::

        sage: from sage.categories.nonassociative_rngs import NonassociativeRngs
        sage: NonassociativeRngs()
        Category of nonassociative rngs
        sage: NonassociativeRngs().super_categories()
        [Category of additive commutative additive associative additive unital distributive magmas and additive magmas,
         Category of commutative additive groups]

        sage: sorted(NonassociativeRngs().axioms())
        ['AdditiveAssociative', 'AdditiveCommutative', 'AdditiveInverse',
         'AdditiveUnital', 'Distributive']

        sage: NonassociativeRngs() is (CommutativeAdditiveGroups() & Magmas()).Distributive()
        True

        sage: NonassociativeRngs().Associative()
        Category of rngs

    TESTS::

        sage: from sage.categories.nonassociative_rngs import NonassociativeRngs
        sage: TestSuite(NonassociativeRngs()).run()
    """
    _base_category_class_and_axiom = (MagmasAndAdditiveMagmas.Distributive.AdditiveAssociative.AdditiveCommutative.AdditiveUnital, "AdditiveInverse")

    Associative = LazyImport('sage.categories.rngs', 'Rngs', at_startup=True)

    class Unital(CategoryWithAxiom):
        """
        Category of nonassociative rings.

        This is a nonassociative rng with a multiplicative identity.
        """
        def _repr_object_names(self):
            """
            EXAMPLES::

                sage: from sage.categories.nonassociative_rngs import NonassociativeRngs
                sage: NonassociativeRngs().Unital() # indirect doctest
                Category of nonassociative rings
            """
            return "nonassociative rings"

