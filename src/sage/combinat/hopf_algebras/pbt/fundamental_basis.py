# -*- coding: utf-8 -*-
r"""
The fundamental basis of PBT Hopf algebra.
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.pbt import PlanarBinaryTreeFunctions
from tree_shuffle_product import TreeShuffleProduct


class Fundamental(PlanarBinaryTreeFunctions.Bases.Base):
    """
    This is the fundamental basis of ``PBT``.

    .. see: [HNT05]_ §3.7 A sylvester description of PBT.

    EXAMPLES::

        sage: P = PBT(QQ).P(); P
        The combinatorial Hopf algebra of Planar Binary Trees Functions over the Rational Field on the Fundamental basis

    .. MATH::

        \mathbf{P}_{T} =
        \sum_{\sigma\in \mathfrak{S}; BST(\sigma) = T} \mathbf{F}_\sigma

    where `(\mathbf{F}_\sigma)` is the fundamental basis of ``FQSym``.

    The product of `(\mathbf{P}_T)_{T}` is described by

    .. MATH::

        \mathbf{P}_{T'} \\times \mathbf{P}_{T''} =
        \sum_{T\in T' \bullet T''} \mathbf{F}_{\gamma}

    where `T' \bullet T''` is the shuffle of trees:

    .. MATH::

        T' \bullet T'' := T_1' \wedge (T_2' \bullet T'') +
        (T' \bullet T_1'') \wedge T_2''

    with `T' := T_1' \wedge T_2'` and `T'' := T_1'' \wedge T_2''`.

    EXAMPLES::

        P  *  P    = P      + P     + P      + P        + P      + P
         o       o    o        o         o        o          o           o
          \     /      \        \         \      /          /           /
           o   o        o        o         o    o          o          _o_
                \        \      /         /      \          \        /   \
                 o        o    o         o        o          o      o     o
                         /      \       / \        \        / \      \
                        o        o     o   o        o      o   o      o
                         \        \                   \
                          o        o                   o

    And the coproduct is described by

    .. MATH::

        \Delta(\mathbf{P}_{T}) =
        \sum_{i\in [n]} \mathbb{F}_{\sigma_{\mid i[}} \otimes
        \mathbb{F}_{\sigma_{\mid [i}}

    where `[n] := [1,..,n]`, `i[` is the sub-interval of `[n]` defined
    by `[1,..,i-1]` and `[i` the sub-interval defined by `[i,..,n]`.

    EXAMPLES::

        sage: ascii_art(P[3,1,2,4].coproduct())
        1 # P      + P  # P    + P  # P      + P    # P    + P    # P    + P      # P
                 o    o      o    o        o    o        o      o      o      o      o
                /           /             /      \      /      /      /      / \
               o           o             o        o    o      o      o      o   o
              / \           \           /
             o   o           o         o
        <BLANKLINE>
         + P      # 1
                o
               /
              o
             / \
            o   o

    (See [HNT05]_, [LR98]_ or [LR02]_).

    TESTS::

        sage: P = PBT(QQ).P()
        sage: # TestSuite(P).run()
    """
    _prefix = "P"

    def build_morphisms(self):
        self._P_to_F_FQSym()

    def _P_to_F_FQSym(self):
        """
        TEST::

            sage: P = PBT(QQ).P()
            sage: F = FQSym(QQ).F()
            sage: F(P[3,1,2])
            F[1, 3, 2] + F[3, 1, 2]
            sage: F(P[3,1,2,4])
            F[1, 3, 2, 4] + F[3, 1, 2, 4]

        TESTS::

            sage: P = PBT(QQ, True).P()
            sage: F = FQSym(QQ).F()
            sage: F(P[3,1,2])
            F[3, 2, 1]
            sage: F(P[3,1,2,4])
            F[3, 2, 1, 4] + F[3, 2, 4, 1] + F[3, 4, 2, 1]
        """
        from sage.combinat.hopf_algebras.all import FQSym

        # fundamental basis of PBT to fundamental basis of FQSym
        F = FQSym(self.base()).F()
        self.module_morphism(
            on_basis=lambda bt: F.sum_of_monomials(
                self._get_sylvester_class(bt)
            ), codomain=F
        ).register_as_coercion()

    def dual_basis(self):
        return self.realization_of().Q()

    def product_on_basis(self, bt1, bt2):
        r"""
        The product of the P basis : shuffle of tree

        TESTS::

            sage: P = PBT(QQ).P()
            sage: ascii_art(P[4,2,1,3] * P[3,1,2])
            P            + P              + P              + P              +
               _o_            __o___           ___o___             __o___
              /   \          /      \         /       \           /      \
             o     o        o       _o_      o         o        _o_       o
              \     \        \     /   \      \       / \      /   \
               o     o        o   o     o      o     o   o    o     o
                    / \            \                /          \     \
                   o   o            o              o            o     o
            <BLANKLINE>
            <BLANKLINE>
            <BLANKLINE>
            P              + P
                   __o__             o_
                  /     \           /  \
               __o__     o         o    o
              /     \             /
             o       o          _o_
              \     /          /   \
               o   o          o     o
                               \
                                o

        """
        return self.sum_of_monomials(TreeShuffleProduct(bt1, bt2))

    def _extract_sub_tree_from_root(self, bt):
        B = self.basis().keys()
        if bt.is_empty():
            yield ([], bt)
            return
        yield ([bt], B())
        for (l0, b0) in self._extract_sub_tree_from_root(bt[0]):
            for (l1, b1) in self._extract_sub_tree_from_root(bt[1]):
                yield (l0 + l1, B([b0, b1]))

    def coproduct_on_basis(self, bt):
        r"""
        The coproduct of the P basis : deconcatenation of tree

        TESTS::

            sage: P = PBT(QQ).P()
            sage: P[[]].coproduct()
            P[] # P[]
            sage: P[1].coproduct()
            P[] # P[1] + P[1] # P[]
            sage: P[1,3,2].coproduct()
            P[] # P[1, 3, 2] + P[1] # P[2, 1] + P[1] # P[1, 2] + P[2, 1] # P[1] + P[1, 2] # P[1] + P[1, 3, 2] # P[]

        EXAMPLES::

            sage: P = PBT(QQ).P()
            sage: ascii_art(P[4,2,1,3].coproduct())
            1 # P        + P  # P      + P  # P    + P    # P    + P    # P    + P      #
                   _o_      o      o      o      o    o      o      o        o    o
                  /   \           / \           /      \      \      \      /      \
                 o     o         o   o         o        o      o      o    o        o
                  \                             \                                    \
                   o                             o                                    o
            <BLANKLINE>
            P  + P    # P  + P    # P    + P    # P  + P        # 1
             o    o      o      o      o      o    o      _o_
                   \           /      /      /           /   \
                    o         o      o      o           o     o
                   /                         \           \
                  o                           o           o
        """
        from sage.categories.tensor import tensor

        return self.tensor_square().sum(
            map(lambda (li, bt): tensor(
                [self.prod(map(self.monomial, li)),
                 self.monomial(bt)]),
                self._extract_sub_tree_from_root(bt)))

    def left_product_on_basis(self, t1, t2):
        """
        TESTS::

            sage: P = PBT(QQ).P()
            sage: P[1] << P[1]
            P[2, 1]
            sage: P[1] << P[[]]
            P[1]
            sage: P[[]] << P[1]
            0
            sage: b1 = BinaryTrees(4).random_element()
            sage: b2 = BinaryTrees(7).random_element()
            sage: (P(b1)<<P(b2)) + (P(b1)>>P(b2)) == P(b1) * P(b2)
            True
        """
        if t1.is_empty():
            return self.zero()
        if t2.is_empty():
            return self.monomial(t1)

        B = self.basis().keys()
        return self.sum_of_monomials(
            B([t1[0], t]) for t in TreeShuffleProduct(t1[1], t2)
        )

    def right_product_on_basis(self, t1, t2):
        """
        TESTS::

            sage: P = PBT(QQ).P()
            sage: P[1] >> P[1]
            P[1, 2]
            sage: P[1] >> P[[]]
            0
            sage: P[[]] >> P[1]
            P[1]
            sage: b1 = BinaryTrees(4).random_element()
            sage: b2 = BinaryTrees(7).random_element()
            sage: (P(b1)<<P(b2)) + (P(b1)>>P(b2)) == P(b1) * P(b2)
            True
        """
        if t1.is_empty():
            return self.monomial(t2)
        if t2.is_empty():
            return self.zero()

        B = self.basis().keys()
        return self.sum_of_monomials(
            B([t, t2[1]]) for t in TreeShuffleProduct(t1, t2[0])
        )

    def _extract_max(self, bt):
        """
        TESTS::

            sage: P = PBT(QQ).P()
            sage: bt = BinaryTree([[None, []],[[], None]]); ascii_art([bt])
            [   __o__   ]
            [  /     \  ]
            [ o       o ]
            [  \     /  ]
            [   o   o   ]
            sage: ascii_art(P._extract_max(bt))
            (   o    o )
            (  /    /  )
            ( o    o   )
            (  \       )
            (   o,     )
        """
        if bt.is_empty():
            return (bt, bt)
        B = self.basis().keys()
        if bt[1].is_empty():
            return (B(), bt)
        (bl, br) = self._extract_max(bt[1])
        return (B([bt[0], bl]), br)

    def left_coproduct_on_basis(self, bt):
        """
        TESTS::

            sage: P = PBT(QQ).P()
            sage: ascii_art([P[3,1,2].left_coproduct()])
            [ P  # P    + P    # P  ]
            [  o      o    o      o ]
            [        /      \       ]
            [       o        o      ]
        """
        from sage.categories.tensor import tensor
        if bt.is_empty() or all(sub_bt.is_empty() for sub_bt in bt):
            return self.base_ring().zero()
        (bt_min, bt_max) = self._extract_max(bt)
        return self.tensor_square().sum_of_terms(map(
            lambda ((bl, br), coeff): ((bl / bt_max, br), coeff),
            self.monomial(bt_min).coproduct(). \
                monomial_coefficients().iteritems()
        )) - tensor([self.monomial(bt), self.one()])

    def right_coproduct_on_basis(self, bt):
        """
        TESTS::

            sage: P = PBT(QQ).P()
            sage: ascii_art([P[3,1,2].right_coproduct()])
            [ P  # P    + P    # P  ]
            [  o    o        o    o ]
            [        \      /       ]
            [         o    o        ]
            sage: ascii_art([P[3,1,2].coproduct() - \
            ....:     (P[3,1,2].right_coproduct() + \
            ....:     P[3,1,2].left_coproduct())])
            [ 1 # P      + P      # 1 ]
            [        o        o       ]
            [       / \      / \      ]
            [      o   o    o   o     ]
        """
        from sage.categories.tensor import tensor
        if bt.is_empty() or all(sub_bt.is_empty() for sub_bt in bt):
            return self.base_ring().zero()
        (bt_min, bt_max) = self._extract_max(bt)
        return self.tensor_square().sum_of_terms(map(
            lambda ((bl, br), coeff): ((bl, br / bt_max), coeff),
            self.monomial(bt_min).coproduct(). \
                monomial_coefficients().iteritems()
        )) - tensor([self.one(), self.monomial(bt)])
