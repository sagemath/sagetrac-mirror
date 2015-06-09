# -*- coding: utf-8 -*-
r"""
The fundamental basis of WQSym Hopf algebra.

S-basis of WQSym
"""
#*****************************************************************************
#       Copyright (C) 2014 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.hopf_algebras.wqsym import WordQuasiSymmetricFunctions
from sage.combinat.packed_word import to_pack
from sage.combinat.shuffle import ShuffleProduct


class Fundamental(WordQuasiSymmetricFunctions.Bases.Base):
    """
    This is the fundamental basis of ``WQSym``:

    EXAMPLES::

        sage: S = WQSym(QQ).S();S
        The combinatorial Hopf algebra of Word Quasi-Symmetric Functions over the Rational Field on the Fundamental basis


    The product of `(\mathbb{S}_\sigma)_{\sigma\in\mathfsf{PW}}` is
    described by the shifted shuffle on packed words
    (..see: :meth:`sage.combinat.packed_word.PackedWord.shifted_shuffle`)

    EXAMPLES::

        sage: S[1,2,1] * S[2,1]
        S[1, 2, 1, 4, 3] + S[1, 2, 4, 1, 3] + S[1, 2, 4, 3, 1] + S[1, 4, 2, 1, 3] + S[1, 4, 2, 3, 1] + S[1, 4, 3, 2, 1] + S[4, 1, 2, 1, 3] + S[4, 1, 2, 3, 1] + S[4, 1, 3, 2, 1] + S[4, 3, 1, 2, 1]

    And the coproduct is described by the deconcatenation and
    packization.

    EXAMPLES::

        sage: S[3,1,2,1].coproduct()
        S[] # S[3, 1, 2, 1] + S[1] # S[1, 2, 1] + S[2, 1] # S[2, 1] + S[3, 1, 2] # S[1] + S[3, 1, 2, 1] # S[]
    """
    _prefix = "S"

    def build_morphisms(self):
        self._S_WQS_to_F_FQS()

    def _S_WQS_to_F_FQS(self):
        """
        TESTS::

            sage: F = FQSym(QQ).F()
            sage: S = WQSym(QQ).S()
            sage: F(S[1,1,1,1])
            F[1, 2, 3, 4]
            sage: p1 = S.monomial(PackedWords(5).random_element())
            sage: p2 = S.monomial(PackedWords(4).random_element())
            sage: F(p1 * p2) == F(p1) * F(p2)
            True
            sage: p3 = S.monomial(PackedWords(8).random_element())
            sage: F.tensor_square()(p3.coproduct()) == F(p3).coproduct()
            True
        """
        #################################################
        # S-basis of WQSym to F-basis of FQSym
        from sage.combinat.hopf_algebras.all import FQSym
        from sage.combinat.permutation import to_standard
        F = FQSym(self.base()).F()
        self.module_morphism(
            on_basis=lambda pw: F.monomial(to_standard(pw)),
            codomain=F
        ).register_as_coercion()
        #################################################

    def dual_basis(self):
        return self.realization_of().M()

    def product_on_basis(self, sigma, mu):
        """
        The product of `(\mathbb{S}_\sigma)_{\sigma\in\mathfsf{PW}}` is
        described by the shifted shuffle on packed words
        (..see: :meth:`sage.combinat.packed_word.PackedWord.shifted_shuffle`).

        TESTS::

            sage: S = WQSym(QQ).S()
            sage: S[1,1]* S[1,2]
            S[1, 1, 2, 3] + S[1, 2, 1, 3] + S[1, 2, 3, 1] + S[2, 1, 1, 3] + S[2, 1, 3, 1] + S[2, 3, 1, 1]
            sage: S[[]] * S[1]
            S[1]
            sage: S[1] * S[[]]
            S[1]
        """
        return self.sum_of_monomials(sigma.shifted_shuffle(mu))

    def coproduct_on_basis(self, sigma):
        """
        The coproduct is described by the deconcatenation and
        packization.

        TESTS::

            sage: S = WQSym(QQ).S()
            sage: S[1].coproduct()
            S[] # S[1] + S[1] # S[]
            sage: S[1,1].coproduct()
            S[] # S[1, 1] + S[1] # S[1] + S[1, 1] # S[]
            sage: S[2,1,2,1,3].coproduct()
            S[] # S[2, 1, 2, 1, 3] + S[1] # S[1, 2, 1, 3] + S[2, 1] # S[2, 1, 3] + S[2, 1, 2] # S[1, 2] + S[2, 1, 2, 1] # S[1] + S[2, 1, 2, 1, 3] # S[]

        """
        from sage.combinat.packed_word import to_pack
        return self.tensor_square().sum_of_monomials(
            (to_pack(sigma[:i]), to_pack(sigma[i:]))
                for i in range(len(sigma) + 1)
        )

    def left_product_on_basis(self, sigma, mu):
        """
        Left dendriform product defines as the shifted
        shuffle of ``sigma`` and ``mu`` such that
        for all `\tau \in \sigma \Cup \mu` one has
        `\tau_{m+n}` (for `m` and `n` respectively the size of
        `\sigma` and `\mu`) comes from `\sigma`.

        TESTS::

            sage: S = WQSym(QQ).S()
            sage: S.left_product(S[1], S[1])
            S[2, 1]
            sage: S[1] << S[1]
            S[2, 1]
            sage: S[1,2] << S[1]
            S[1, 3, 2] + S[3, 1, 2]
            sage: (S[1,2] >> S[3,1,1,2]) + (S[1,2] << S[3,1,1,2]) == S[1,2] * S[3,1,1,2]
            True
        """
        if len(sigma) < 1:
            return self(self.base_ring().zero())

        return self.sum_of_monomials(map(
            lambda gamma: self.basis().keys()(list(gamma) + [sigma[-1]]),
            ShuffleProduct(sigma[:-1], [l + len(sigma) for l in mu], list)
        ))

    def right_product_on_basis(self, sigma, mu):
        """
        Right dendriform product defines as the shifted
        shuffle of ``sigma`` and ``mu`` such that
        for all `\tau \in \sigma \Cup \mu` one has
        `\tau_{m+n}` (for `m` and `n` respectively the size of
        `\sigma` and `\mu`) comes from `\mu`.

        TESTS::

            sage: S = WQSym(QQ).S()
            sage: S.right_product(S[1], S[1])
            S[1, 2]
            sage: S[1,2] >> S[1]
            S[1, 2, 3]
        """
        if len(mu) < 1:
            return self(self.base_ring().zero())

        return self.sum_of_monomials(map(
            lambda gamma: self.basis().keys()(
                list(gamma) + [mu[-1] + len(sigma)]),
            ShuffleProduct(sigma, [l + len(sigma) for l in mu[:-1]], list)
        ))

    def left_coproduct_on_basis(self, pw):
        """
        Right dendriform coproduct defines as the deconcatenation-packization
        of ``sigma``such that the (last) maximum of `\sigma` is "left" on the left
        of the tensor.

        TESTS::

            sage: S = WQSym(QQ).S()
            sage: S[1,1,2,2].left_coproduct()
            0
            sage: S[1,1,2,2,1].left_coproduct()
            S[1, 1, 2, 2] # S[1]
            sage: F = FQSym(QQ).F()
            sage: FxF = F.tensor_square()
            sage: FxF(S[1,1,2,2,1].left_coproduct()) == F(S[1,1,2,2,1]).left_coproduct()
            True
        """
        if pw.size() < 2:
            return self(self.base_ring().zero())

        l = list(pw); l.reverse()
        id_max = pw.size() - l.index(max(pw)) - 1
        return self.tensor_square().sum_of_monomials(
            [(to_pack(pw[:i+1]), to_pack(pw[i+1:]))
            for i in range(id_max, pw.size()-1)]
        )

    def right_coproduct_on_basis(self, pw):
        """
        Right dendriform coproduct defines as the deconcatenation-packization
        of ``sigma``such that the (last) maximum of `\sigma` is "left" on the right
        of the tensor.

        TESTS::

            sage: S = WQSym(QQ).S()
            sage: S[1,1,2,2].right_coproduct()
            S[1] # S[1, 2, 2] + S[1, 1] # S[1, 1] + S[1, 1, 2] # S[1]
            sage: S[1,1,2,2,1].right_coproduct()
            S[1] # S[1, 2, 2, 1] + S[1, 1] # S[2, 2, 1] + S[1, 1, 2] # S[2, 1]
            sage: F = FQSym(QQ).F()
            sage: FxF = F.tensor_square()
            sage: FxF(S[1,1,2,2,1].right_coproduct()) == F(S[1,1,2,2,1]).right_coproduct()
            True

        """
        if pw.size() < 2:
            return self(self.base_ring().zero())

        l = list(pw); l.reverse()
        id_max = pw.size() - l.index(max(pw)) - 1

        return self.tensor_square().sum_of_monomials(
            [(to_pack(pw[:i+1]), to_pack(pw[i+1:]))
            for i in range(min(id_max, pw.size()))]
        )