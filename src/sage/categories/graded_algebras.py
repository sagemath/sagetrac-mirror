r"""
Graded Algebras
"""
#*****************************************************************************
#  Copyright (C) 2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.graded_modules import GradedModulesCategory
from sage.misc.cachefunc import cached_method

class GradedAlgebras(GradedModulesCategory):
    """
    The category of graded algebras

    EXAMPLES::

        sage: GradedAlgebras(ZZ)
        Category of graded algebras over Integer Ring
        sage: GradedAlgebras(ZZ).super_categories()
        [Category of algebras over Integer Ring,
         Category of graded modules over Integer Ring]

    TESTS::

        sage: TestSuite(GradedAlgebras(ZZ)).run()
    """

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class SubcategoryMethods:
        @cached_method
        def Differential(self):
            r"""
            Return the full subcategory of the objects of ``self`` equipped
            with a differential.

            A differential graded algebra `A` is a `\ZZ / 2 \ZZ` graded
            algebra equipped with at least one *derivation*, that is a
            linear morphism `\partial` satisfying the *graded* Leibniz rule:

            .. MATH::

                \partial (xy) = (\partial x) y + (-1)^{\deg(x)} x (\partial y)

            for all `x,y \in A`.

            .. SEEALSO:: :wikipedia:`Differential_algebra`

            EXAMPLES::

                sage: Algebras(QQ).Graded().Differential()
                Category of graded differential algebras

            TESTS::

                sage: TestSuite(Algebra(QQ).Graded().Differential()).run()
                sage: Rings().Associative.__module__
                'sage.categories.graded_algebras'
            """
            return self._with_axiom('Differential')

