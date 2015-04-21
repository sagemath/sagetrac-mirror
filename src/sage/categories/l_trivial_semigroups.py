r"""
L-Trivial Semigroups
"""
#*****************************************************************************
#  Copyright (C) 2009-2010 Florent Hivert <florent.hivert at univ-rouen.fr>
#                2009-2015 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from magmas import Magmas
from semigroups import Semigroups

class LTrivialSemigroups(CategoryWithAxiom):
    def extra_super_categories(self):
        r"""
        Implement the fact that a `L`-trivial semigroup is `H`-trivial.

        EXAMPLES::

            sage: Semigroups().LTrivial().extra_super_categories()
            [Category of htrivial semigroups]
        """
        return [Semigroups().HTrivial()]

    def RTrivial_extra_super_categories(self):
        r"""
        Implement the fact that an `L`-trivial and `R`-trivial monoid
        is `J`-trivial.

        EXAMPLES::

            sage: Semigroups().LTrivial().RTrivial_extra_super_categories()
            [Category of jtrivial magmas]

        TESTS::

            sage: Semigroups().LTrivial().RTrivial() is Semigroups().JTrivial()
            True
        """
        return [Magmas().JTrivial()]

    def Commutative_extra_super_categories(self):
        r"""
        Implement the fact that a commutative `R`-trivial semigroup is `J`-trivial.

        EXAMPLES::

            sage: Semigroups().LTrivial().Commutative_extra_super_categories()
            [Category of jtrivial semigroups]

        TESTS::

            sage: Semigroups().LTrivial().Commutative() is Semigroups().JTrivial().Commutative()
            True
        """
        return [self.JTrivial()]
