# -*- coding: utf-8 -*-
from chasage.all import *
r"""
Quadri algebras

A quadri-algebra _[AguiarLoday] is a vector space `Q` together with four
operations `\searrow`, `\nearrow`, `\nwarrow` and `\swarrow` (called
respectively south east, north east, north west and south west products):
`Q \otimes Q \to Q` satisfying the axioms below. In order to state them,
consider the following operations

MATH::

    \begin{align}
        x \succ y  &:= x \nearrow y + x \searrow y\\
        x \prec y  &:= x \nwarrow y + x \swarrow y\\
        x \vee  y  &:= x \searrow y + x \swarrow y\\
        x \wedge y &:= x \nearrow y + x \nwarrow y
    \end{align}

called respectively east, west, south and north product; and

MATH::

    \begin{align}
    x \star y &:= x \searrow y + x \nearrow y + x \nwarrow y + x \swarrow y\\
              &= x \succ y + x \prec y = x \vee y + x \wedge y
    \end{align}

The axioms are

MATH::

    \begin{array}{ccc}
        (x \nwarrow y) \nwarrow z = x \nwarrow (y \star z) &
        (x \nwarrow y) \nearrow z = x \nearrow (y \prec z) &
        (x \wedge y) \nearrow z   = x \nearrow (y \succ z)\\
        (x \swarrow y) \nwarrow z = x \swarrow (y \wedge z) &
        (x \searrow y) \nwarrow z = x \searrow (y \nwarrow z)&
        (x \vee y) \nearrow z     = x \searrow (y \nearrow z)\\
        (x \prec y) \swarrow z    = x \swarrow (y \vee z) &
        (x \succ y) \swarrow z    = x \searrow (y \swarrow z) &
        (x \star y) \searrow z    = x \searrow (y \searrow z)
    \end{array}

References:

.. [Loday] Quadri-algebras,
    Marcelo Aguiar and
    Jean-Louis Loday

AUTHORS:

 - Jean-Baptiste Priez (first version)
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import LazyImport
from sage.categories.dual import DualObjectsCategory
from sage.structure.sage_object import have_same_parent
from sage.categories.realizations import RealizationsCategory
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.misc.lazy_attribute import lazy_attribute

from sage.categories.category_with_axiom import \
    CategoryWithAxiom_over_base_ring
from sage.categories.category_types import Category_over_base_ring

from sage.categories.algebras_over_operads import AlgebrasOverOperad


class QuadriAlgebras(Category_over_base_ring):
    """
    The category of quadri-algebras over a given base ring.

    A quadri-algebra _[AguiarLoday] is a vector space `Q` together with four
    operations `\searrow`, `\nearrow`, `\nwarrow` and `\swarrow` (called
    respectively south east, north east, north west and south west products):
    `Q \otimes Q \to Q` satisfying the axioms below. In order to state them,
    consider the following operations

    MATH::

        \begin{align}
            x \succ y  &:= x \nearrow y + x \searrow y\\
            x \prec y  &:= x \nwarrow y + x \swarrow y\\
            x \vee  y  &:= x \searrow y + x \swarrow y\\
            x \wedge y &:= x \nearrow y + x \nwarrow y
        \end{align}

    called respectively east, west, south and north products; and

    MATH::

        \begin{align}
        x \star y &:= x \searrow y + x \nearrow y + x \nwarrow y + x \swarrow y\\
                  &= x \succ y + x \prec y = x \vee y + x \wedge y
        \end{align}

    The axioms are

    MATH::

        \begin{array}{ccc}
            (x \nwarrow y) \nwarrow z = x \nwarrow (y \star z) &
            (x \nwarrow y) \nearrow z = x \nearrow (y \prec z) &
            (x \wedge y) \nearrow z   = x \nearrow (y \succ z)\\
            (x \swarrow y) \nwarrow z = x \swarrow (y \wedge z) &
            (x \searrow y) \nwarrow z = x \searrow (y \nwarrow z)&
            (x \vee y) \nearrow z     = x \searrow (y \nearrow z)\\
            (x \prec y) \swarrow z    = x \swarrow (y \vee z) &
            (x \succ y) \swarrow z    = x \searrow (y \swarrow z) &
            (x \star y) \searrow z    = x \searrow (y \searrow z)
        \end{array}

    EXAMPLES::

        sage: from chasage.combinat.cha.categories.dendriform_algebras import \
        ....:     DendriformAlgebras
        sage: DendriformAlgebras(ZZ)
        Category of dendriform algebras over Integer Ring

    TESTS::

        sage: TestSuite(DendriformAlgebras(ZZ)).run()
    """

    def super_categories(self):
        """
        TESTS::

            sage: from chasage.combinat.cha.categories.dendriform_algebras \
            ....:     import DendriformAlgebras
            sage: DendriformAlgebras(ZZ).super_categories()
            [Category of algebras over Integer Ring]
        """
        R = self.base_ring()
        return [AlgebrasOverOperad(R)]

    class ParentMethods:

        def _test_quadri_algebras(self, **options):
            r"""
            Run generic tests on the methods quadri-algebras.

            See also: :class:`TestSuite`.

            TESTS::

                sage: FQSym(QQ)._test_quadri_algebras()
            """
            pass

        # east, west, south and north products
        def east_product(self, x, y):
            r"""
            The east product of two elements defined by
            `x \succ y  := x \nearrow y + x \searrow y`.
            """
            return x.ne_product(y) + x.se_product(y)

        def west_product(self, x, y):
            r"""
            The west product of two elements defined by
            `x \prec y  := x \nwarrow y + x \swarrow y`.
            """
            return x.nw_product(y) + x.sw_product(y)

        def south_product(self, x, y):
            r"""
            The south product of two elements defined by
            `x \vee  y  := x \searrow y + x \swarrow y`.
            """
            return x.se_product(y) + x.sw_product(y)

        def north_product(self, x, y):
            r"""
            The north product of two elements defined by
            `x \wedge y := x \nearrow y + x \nwarrow y`.
            """
            return x.ne_product(y) + x.nw_product(y)

        def ne_product(self, x, y):
            pass

        def nw_product(self, x, y):
            pass

        def se_product(self, x, y):
            pass

        def sw_product(self, x, y):
            pass

    class ElementMethods:

        pass

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def east_product(self, x, y):
                r"""
                The east product of two elements defined by
                `x \succ y  := x \nearrow y + x \searrow y`.
                """
                return x.ne_product(y) + x.se_product(y)

            def west_product(self, x, y):
                r"""
                The west product of two elements defined by
                `x \prec y  := x \nwarrow y + x \swarrow y`.
                """
                return x.nw_product(y) + x.sw_product(y)

            def south_product(self, x, y):
                r"""
                The south product of two elements defined by
                `x \vee  y  := x \searrow y + x \swarrow y`.
                """
                return x.se_product(y) + x.sw_product(y)

            def north_product(self, x, y):
                r"""
                The north product of two elements defined by
                `x \wedge y := x \nearrow y + x \nwarrow y`.
                """
                return x.ne_product(y) + x.nw_product(y)

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            pass

from chasage.combinat.cha.categories.tools import expand_binary_operator
expand_binary_operator(QuadriAlgebras, "ne_product")
expand_binary_operator(QuadriAlgebras, "nw_product")
expand_binary_operator(QuadriAlgebras, "se_product")
expand_binary_operator(QuadriAlgebras, "sw_product")
