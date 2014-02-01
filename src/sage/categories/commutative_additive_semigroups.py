r"""
Commutative additive semigroups
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.additive_semigroups import AdditiveSemigroups

class CommutativeAdditiveSemigroups(CategoryWithAxiom):
    """
    The category of additive abelian semigroups, i.e. sets with an
    associative and abelian operation +.

    EXAMPLES::

        sage: C = CommutativeAdditiveSemigroups(); C
        Category of commutative additive semigroups
        sage: C.example()
        An example of a commutative monoid: the free commutative monoid generated by ('a', 'b', 'c', 'd')

        sage: sorted(C.super_categories(), key=str)
        [Category of additive commutative additive magmas,
         Category of additive semigroups]
        sage: sorted(C.axioms())
        ['AdditiveAssociative', 'AdditiveCommutative']
        sage: C is AdditiveMagmas().AdditiveAssociative().AdditiveCommutative()
        True

    .. NOTE::

        This category is currently empty and only serves as a place
        holder to make ``C.example()`` work.

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = [AdditiveSemigroups, "AdditiveCommutative"]
