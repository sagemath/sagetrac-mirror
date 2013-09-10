r"""
Integral domains

AUTHORS:

- Teresa Gomez-Diaz (2008): initial version

- Julian Rueth (2013-09-10): added category refinement similar to the one done by Fields

"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton, Category_contains_method_by_parent_class
from sage.misc.cachefunc import cached_method
from sage.categories.commutative_rings import CommutativeRings
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.domains import Domains
import sage.rings.ring

class IntegralDomains(Category_singleton):
    """
    The category of integral domains, i.e., non-trivial commutative rings with
    no zero divisors

    EXAMPLES::

        sage: IntegralDomains()
        Category of integral domains
        sage: IntegralDomains().super_categories()
        [Category of commutative rings, Category of domains]

    TESTS::

        sage: TestSuite(IntegralDomains()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: IntegralDomains().super_categories()
            [Category of commutative rings, Category of domains]
        """
        return [CommutativeRings(), Domains()]

    def __contains__(self, x):
        """
        Return whether ``x`` is in the category of integral domains.

        EXAMPLES::

            sage: GF(4, "a") in IntegralDomains()
            True
            sage: QQ in IntegralDomains()
            True
            sage: ZZ in IntegralDomains()
            True
            sage: IntegerModRing(4) in IntegralDomains()
            False

        Test that :trac:`15183` has been fixed::

            sage: IntegerModRing(2) in IntegralDomains()
            True

        """
        try:
            return self._contains_helper(x) or sage.rings.ring._is_IntegralDomain(x)
        except StandardError:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of integral domains.

        This helper just tests whether the given object's category is already
        known to be a sub-category of the category of integral domains. There
        are, however, rings that are initialised as plain commutative rings and
        found out to be integral domains only afterwards. Hence, this helper
        alone is not enough for a proper containment test.

        TESTS::

            sage: P.<x> = QQ[]
            sage: Q = P.quotient(x^2+2)
            sage: Q.category()
            Join of Category of commutative algebras over Rational Field and Category of subquotients of monoids and Category of quotients of semigroups
            sage: ID = IntegralDomains()
            sage: ID._contains_helper(Q)
            False
            sage: Q in ID  # This changes the category!
            True
            sage: ID._contains_helper(Q)
            True

        """
        return Category_contains_method_by_parent_class(cls())

    class ParentMethods:
        def is_integral_domain(self):
            """
            Return True, since this in an object of the category of integral domains.

            EXAMPLES::

                sage: Parent(QQ,category=IntegralDomains()).is_integral_domain()
                True

            """
            return True

    class ElementMethods:
        pass
