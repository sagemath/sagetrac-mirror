r"""
Super Algebras
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.super_modules import SuperModulesCategory
from sage.categories.algebras import Algebras
from sage.categories.modules import Modules
from sage.misc.lazy_import import LazyImport

class SuperAlgebras(SuperModulesCategory):
    """
    The category of super algebras.

    An `R`-*super algebra* is an `R`-super module `A` endowed with an
    `R`-algebra structure satisfying

    .. MATH::

        A_0 A_0 \subseteq A_0, \qquad
        A_0 A_1 \subseteq A_1, \qquad
        A_1 A_0 \subseteq A_1, \qquad
        A_1 A_1 \subseteq A_0

    and `1 \in A_0`.

    EXAMPLES::

        sage: Algebras(ZZ).Super()
        Category of super algebras over Integer Ring

    TESTS::

        sage: TestSuite(Algebras(ZZ).Super()).run()
    """
    def extra_super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(ZZ).Super().super_categories() # indirect doctest
            [Category of graded algebras over Integer Ring,
             Category of super modules over Integer Ring]
        """
        return [self.base_category().Graded()]

