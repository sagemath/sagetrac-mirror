r"""
Integral domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_contains_method_by_parent_class
from sage.categories.domains import Domains

class IntegralDomains(CategoryWithAxiom):
    """
    The category of integral domains

    An integral domain is commutative ring with no zero divisors, or
    equivalently a commutative domain.

    EXAMPLES::

        sage: C = IntegralDomains(); C
        Category of integral domains
        sage: sorted(C.super_categories(), key=str)
        [Category of commutative rings, Category of domains]
        sage: C is Domains().Commutative()
        True
        sage: C is Rings().Commutative().NoZeroDivisors()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (Domains, "Commutative")

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in IntegralDomains()
            True
            sage: QQ in IntegralDomains()
            True
            sage: ZZ in IntegralDomains()
            True
            sage: IntegerModRing(4) in IntegralDomains()
            False
            sage: IntegerModRing(5) in IntegralDomains()
            True

        This implementation will not be needed anymore once every
        field in Sage will be properly declared in the category
        :class:`IntegralDomains`().
        """
        try:
            return self._contains_helper(x) or x.is_integral_domain()
        except Exception:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of integral domains.

        This helper just tests whether the given object's category
        is already known to be a sub-category of the category of
        integral domains. There are, however, rings that are initialised
        as plain commutative rings and found out to be integral domains
        only afterwards. Hence, this helper alone is not enough
        for a proper containment test.

        TESTS::

            sage: R = Zmod(7)
            sage: R.category()
            Join of Category of finite commutative rings
             and Category of subquotients of monoids
             and Category of quotients of semigroups
            sage: ID = IntegralDomains()
            sage: ID._contains_helper(R)
            False
            sage: R in ID  # This changes the category!
            True
            sage: ID._contains_helper(R)
            True
        """
        return Category_contains_method_by_parent_class(cls())

    class ParentMethods:
        def is_integral_domain(self):
            """
            Return True, since this in an object of the category of integral domains.

            EXAMPLES::

                sage: QQ.is_integral_domain()
                True
                sage: Parent(QQ,category=IntegralDomains()).is_integral_domain()
                True

            """
            return True

    class ElementMethods:
        pass
