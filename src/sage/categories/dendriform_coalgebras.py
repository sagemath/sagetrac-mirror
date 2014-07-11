# -*- coding: utf-8 -*-
r"""
Dendriform coalgebras

A dendriform coalgebra `C` over a field `K` is a non-counitary coassociative
coalgebra, such that the coproduct `\delta` can be split into two part
`\delta_\prec` and `\delta_\succ` (called left and right coproducts) with good
compatibilities.
That structure is the dual of a dendriform algebra.
(see :mod:`sage.combinat.hopf_algebras.categories.dendriform_algebras`)

AUTHORS:

 - Jean-Baptiste Priez
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category_types import Category_over_base_ring
from sage.categories.modules import Modules
from sage.misc.abstract_method import abstract_method
from sage.categories.dual import DualObjectsCategory
from sage.categories.realizations import RealizationsCategory
from sage.categories.all import tensor
from sage.categories.with_realizations import WithRealizationsCategory

from sage.categories.category_with_axiom import \
    CategoryWithAxiom_over_base_ring
from sage.misc.cachefunc import cached_method


class DendriformCoalgebras(Category_over_base_ring):
    """
    The category of dendriform coalgebras over a given base ring.

    An dendriform coalgebra over a ring `R` is an coalgebra with `\delta` is a
    coproduct defined by `\delta := \delta_\prec + \delta_\succ + \bar\delta`
    three operators with `\bar\delta(x) = 1\otimes x + x\otimes 1`
    such that respect following properties, for `x,y,z` three basis element:

    MATH::

        \begin{align}
            ( \delta_\prec \otimes Id ) \circ \delta_\prec( x ) &=
            ( Id \otimes \bar\delta )\circ \delta_\prec(x)\\
            ( \delta_\succ \otimes Id ) \circ \delta_\prec( x ) &=
            ( Id \otimes \delta_\prec )\circ \delta_\succ(x)\\
            ( \bar\delta \otimes Id ) \circ \delta_\succ( y ) &=
            ( Id \otimes \delta_\succ )\circ \delta_\succ(x)\\
        \end{align}

    EXAMPLES::

        sage: from sage.combinat.hopf_algebras.categories.dendriform_coalgebras \
        ....:     import DendriformCoalgebras
        sage: DendriformCoalgebras(ZZ)
        Category of dendriform coalgebras over Integer Ring

    TESTS::
        sage: from sage.combinat.hopf_algebras.categories.dendriform_coalgebras \
        ....:     import DendriformCoalgebras
        sage: TestSuite(DendriformCoalgebras(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        TESTS::

            sage: from sage.categories.dendriform_algebras import DendriformAlgebras
            sage: DendriformAlgebras(QQ)
            Category of dendriform algebras over Rational Field
        """
        R = self.base_ring()
        return [Modules(R)]

    class ParentMethods:

        def _test_codendriform(self, **options):
            r"""
            Run generic tests on the methods :meth:`.left_coproduct` and
            `.right_coproduct'.

            See also: :class:`TestSuite`.

            TESTS::

                sage: # FQSym(QQ)._test_codendriform()
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            F = S[0].parent()
            FxFxF = tensor((F, F, F))
            delta = lambda x: F.m_coproduct(x)  # - tensor((self.one(), x)) \
            #                                    - tensor((x, self.one()))
            dprec = lambda x: x.left_coproduct()
            dsucc = lambda x: x.right_coproduct()
            mlm = lambda f, x: x.apply_multilinear_morphism(f, FxFxF)
            for e in S:
                # ( \delta_\prec \otimes Id ) \circ \delta_\prec (c) &=
                # ( Id \otimes \bar\delta ) \circ \delta_\prec (c)
                tester.assertEqual(
                    mlm(lambda l, r: tensor((dprec(l), r)), dprec(e)),
                    mlm(lambda l, r: tensor((l, delta(r))), dprec(e))
                )
                # ( \delta_\succ \otimes Id ) \circ \delta_\prec(c) &=
                # ( Id \otimes \delta_\prec ) \circ \delta_\succ(c)
                tester.assertEqual(
                    mlm(lambda l, r: tensor((dsucc(l), r)), dprec(e)),
                    mlm(lambda l, r: tensor((l, dprec(r))), dsucc(e))
                )
                # ( \bar\delta \otimes Id ) \circ \delta_\succ(c) &=
                # ( Id \otimes \delta_\succ ) \circ \delta_\succ(c)
                tester.assertEqual(
                    mlm(lambda l, r: tensor((delta(l), r)), dsucc(e)),
                    mlm(lambda l, r: tensor((l, dsucc(r))), dsucc(e))
                )

        @abstract_method
        def left_coproduct(self, elem):
            """
            The left coproduct (`\delta_\prec`) of an element.
            """

        @abstract_method
        def right_coproduct(self, elem):
            """
            The right coproduct (`\delta_\prec`) of an element.
            """

        def m_coproduct(self, elem):
            """
            By default a dendriform coalgebra is a non-counital coalgebra with 
            coproduct `\tilde{\delta}` is defined for any `x` by:

            MATH::

                \tilde{\delta}(x) := \delta_\prec(x) + \delta_\succ(x)
            """
            return self.right_coproduct(elem) + self.left_coproduct(elem)

    class ElementMethods:
        pass
#         def m_coproduct(self):
#             return self.parent().m_coproduct(self)

    class WithRealizations(WithRealizationsCategory):

        class ParentMethods:
            pass

    class Realizations(RealizationsCategory):

        class ParentMethods:

            def m_coproduct(self, x):
                """
                (See :meth:`DendriformAlgebras.ParentMethods.product`.)
                """
                return self.left_coproduct(x) + self.right_coproduct(x)

    class WithBasis(CategoryWithAxiom_over_base_ring):

        class ParentMethods:

            pass

    class DualObjects(DualObjectsCategory):

        def extra_super_categories(self):
            r"""
            Returns the dual category

            EXAMPLES:

            The category of algebras over the Rational Field is dual
            to the category of coalgebras over the same field::

                sage: C = Algebras(QQ)
                sage: C.dual()
                Category of duals of algebras over Rational Field
                sage: C.dual().extra_super_categories()
                [Category of coalgebras over Rational Field]
            """
            from sage.categories.dendriform_algebras import DendriformAlgebras
            return [DendriformAlgebras(self.base_category().base_ring())]

from sage.categories.operators_tools import expand_unary_operator

expand_unary_operator(DendriformCoalgebras, "left_coproduct")
expand_unary_operator(DendriformCoalgebras, "right_coproduct")
