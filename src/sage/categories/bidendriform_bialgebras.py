# -*- coding: utf-8 -*-
r"""
Bidendriform bialgebras

A bidendriform bialgebra is both a dendriform bialgebra and a codendriform,
with some compatibilities. It is a family `(A, \ll, \gg, \delta_\prec,
\delta_\succ)` such that:

    1 - `(A, \ll, \gg)` is a dendriform algebra,
    2 - `(A, \delta_\prec, \delta_\succ)` is a dendriform algebra,
    3 - the following compatibilities are satisfied: for all `a, b \in A`,

        MATH::

        \begin{align}
            \delta_\succ(a \gg b) &= a' b_\succ' \otimes a'' \gg b_\succ''
                + a' \otimes a'' \gg b
                + b_\succ' \otimes a \gg b_\succ''
                + a b_\succ' \otimes b_\succ''
                + a \otimes b\\
            \delta_\succ(a \ll b) &= a' b_\succ' \otimes a'' \ll b_\succ''
                + a' \otimes a'' \ll b
                + b_\succ' \otimes a \ll b_\succ''\\
            \delta_\prec(a \gg b) &= a' b_\prec' \otimes a'' \gg b_\prec''
                + ab_\prec' \otimes b_\prec''
                + b_\prec \otimes a \gg b_\prec''\\
            # \delta_\prec(a \ll b) &= a'b_\prec' \otimes a'' \ll b_\prec''
                + a'b \otimes a''
                + b_\prec' \otimes a \ll b_\prec''
                + b \otimes a
        \end{align}

AUTHORS:

 - Jean-Baptiste Priez
"""
#*****************************************************************************
#  Copyright (C) 2013      Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************
from sage.categories.category import Category
from sage.categories.category_types import Category_over_base_ring
from sage.categories.tensor import tensor

from sage.categories.dendriform_algebras import \
    DendriformAlgebras
from sage.categories.dendriform_coalgebras import \
    DendriformCoalgebras
from sage.categories.algebras import Algebras


class BidendriformBialgebras(Category_over_base_ring):
    """
    The category of bidendriform bialgebras over a given base ring.

    A bidendriform bialgebra over a ring `R` is an bialgebra which is a
    dendriform algebra and a dendriform coalgebra such that respect these
    following conditions:


    TODO: should `R` be a commutative ring?

    EXAMPLES::

        sage: from sage.combinat.hopf_algebras.categories.\
        ....:     bidendriform_bialgebras import BidendriformBialgebras
        sage: BidendriformBialgebras(QQ)
        Category of bidendriform bialgebras over Rational Field

    TESTS::

        sage: from sage.combinat.hopf_algebras.categories.\
        ....:     bidendriform_bialgebras import BidendriformBialgebras
        sage: TestSuite(BidendriformBialgebras(ZZ)).run()
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.hopf_algebras.categories.\
            ....:     bidendriform_bialgebras import BidendriformBialgebras
            sage: BidendriformBialgebras(ZZ).super_categories()
            [Category of dendriform algebras over Integer Ring,
             Category of dendriform coalgebras over Integer Ring]
        """
        R = self.base_ring()
        return Category.join([DendriformAlgebras(R), DendriformCoalgebras(R)], as_list=True)

    class ParentMethods:

        def _test_bidendriform(self, **options):
            r"""
            Run generic tests on the methods :meth:`.left_coproduct` and
            :meth:`.right_coproduct' and :meth:`.left_product` and
            :meth:`.right_product'.

            See also: :class:`TestSuite`.

            TESTS::

                sage: FQSym(QQ)._test_bidendriform()
            """
            return
            #FIXME: replace tests
            import itertools
            tester = self._tester(**options)
            S = tester.some_elements()
            F = S[0].parent()
            FxF = F.tensor_square()
            delta = lambda x: F.m_coproduct(x)  # - tensor((self.one(), x)) \
            #                                    - tensor((x, self.one()))
            dprec = lambda x: x.left_coproduct()
            dsucc = lambda x: x.right_coproduct()
            mprod = lambda x, y: F.m_product(x, y)
            if Algebras(self.base_ring()).Unital() in self.categories():
                minus = lambda x: x[self.one().support()[0]] * self.one()
            else: minus = lambda x: x
            for a in S:
                a = minus(a)
                for b in S:
                    b = minus(b)
                    # \delta_\succ(a \gg b) =
                    #     a' b_\succ' \otimes a'' \gg b_\succ''
                    #   + a' \otimes a'' \gg b
                    #   + b_\succ' \otimes a \gg b_\succ''
                    #   + a b_\succ' \otimes b_\succ''
                    #   + a \otimes b
                    tester.assertEqual(
                        # \delta_\succ(a \gg b)
                        dsucc(a >> b),
                        # a' b_\succ' \otimes a'' \gg b_\succ''
                        FxF.sum(ca * cb * tensor([mprod(F(al), F(bl)), F(ar) >> F(br)])
                            for (((al, ar), ca), ((bl, br), cb)) in \
                            itertools.product(delta(a), dsucc(b))
                        ) +
                        # a' \otimes a'' \gg b
                        FxF.sum(ca * tensor([F(al), F(ar) >> b])
                            for ((al, ar), ca) in delta(a)
                        ) +
                        # b_\succ' \otimes a \gg b_\succ''
                        FxF.sum(cb * tensor([F(bl), a >> F(br)])
                            for ((bl, br), cb) in dsucc(b)
                        ) +
                        # a b_\succ' \otimes b_\succ''
                        FxF.sum(cb * tensor([mprod(a, F(bl)), F(br)])
                            for ((bl, br), cb) in dsucc(b)
                        ) +
                        # a \otimes b
                        tensor([a, b]),
                        "First condition (see :mod:sage.categories.bidendriform_bialgebras)"
                    )
                    # \delta_\succ(a \ll b) =
                    #       a' b_\succ' \otimes a'' \ll b_\succ''
                    #     + a' \otimes a'' \ll b
                    #     + b_\succ' \otimes a \ll b_\succ''
                    tester.assertEqual(
                        # \delta_\succ(a \ll b)
                        dsucc(a << b),
                        # a' b_\succ' \otimes a'' \ll b_\succ''
                        FxF.sum(ac * bc * tensor([mprod(F(al), F(bl)), F(ar) << F(br)])
                            for (((al, ar), ac), ((bl, br), bc)) in \
                            itertools.product(delta(a), dsucc(b))
                        ) +
                        # a' \otimes a'' \ll b
                        FxF.sum(ac * tensor([F(al), F(ar) << b])
                            for ((al, ar), ac) in delta(a)
                        ) +
                        # b_\succ' \otimes a \ll b_\succ''
                        FxF.sum(bc * tensor([F(bl), a << F(br)])
                            for ((bl, br), bc) in dsucc(b)),
                        "Second condition (see :mod:sage.categories.bidendriform_bialgebras)"
                    )
                    # \delta_\prec(a \gg b) =
                    #     a' b_\prec' \otimes a'' \gg b_\prec''
                    #   + ab_\prec' \otimes b_\prec''
                    #   + b_\prec \otimes a \gg b_\prec''
                    tester.assertEqual(
                        # \delta_\prec(a \gg b)
                        dprec(a >> b),
                        # a' b_\prec' \otimes a'' \gg b_\prec''
                        FxF.sum(ac * bc * tensor([mprod(F(al), F(bl)), F(ar) >> F(br)])
                            for (((al, ar), ac), ((bl, br), bc)) in \
                            itertools.product(delta(a), dprec(b))
                        ) +
                        # ab_\prec' \otimes b_\prec''
                        FxF.sum(bc * tensor([mprod(a, F(bl)), F(br)])
                            for ((bl, br), bc) in dprec(b)
                        ) +
                        # b_\prec \otimes a \gg b_\prec''
                        FxF.sum(bc * tensor([F(bl), a >> F(br)])
                            for ((bl, br), bc) in dprec(b)),
                        "Third condition (see :mod:sage.categories.bidendriform_bialgebras)"
                    )
                    # \delta_\prec(a \ll b) =
                    #     a'b_\prec' \otimes a'' \ll b_\prec''
                    #   + a'b \otimes a''
                    #   + b_\prec' \otimes a \ll b_\prec''
                    #   + b \otimes a
                    tester.assertEqual(
                        dprec(a << b),
                        # a'b_\prec' \otimes a'' \ll b_\prec''
                        FxF.sum(ac * bc * tensor([mprod(F(al), F(bl)), F(ar) << F(br)])
                            for (((al, ar), ac), ((bl, br), bc)) in \
                            itertools.product(delta(a), dprec(b))
                        ) +
                        # a'b \otimes a''
                        FxF.sum(ac * tensor([mprod(F(al), b), F(ar)])
                            for ((al, ar), ac) in delta(a)
                        ) +
                        # b_\prec' \otimes a \ll b_\prec''
                        FxF.sum(bc * tensor([F(bl), a << F(br)])
                            for ((bl, br), bc) in dprec(b)
                        ) +
                        # b \otimes a
                        tensor([b, a]),
                        "Fourth condition (see :mod:sage.categories.bidendriform_bialgebras)"
                    )
        class ElementMethods:
            pass