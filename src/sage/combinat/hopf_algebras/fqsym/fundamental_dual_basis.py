# -*- coding: utf-8 -*-
r"""
The fundamental dual basis of FQSym Hopf algebra.

G-basis of FQSym
"""
#*****************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import itertools
from sage.misc.misc_c import prod
from sage.combinat.ncsf_qsym.ncsf import NonCommutativeSymmetricFunctions
from sage.combinat.permutation import to_standard, Permutations
from sage.combinat.hopf_algebras.fqsym import FreeQuasiSymmetricFunctions

class FundamentalDual(FreeQuasiSymmetricFunctions.Bases.Base):
    """
    The `(\mathbb{G}_{\sigma})_{\sigma\in \mathfrak{G}}`-basis is the dual
    basis of `(\mathbb{F}_{\sigma})_{\mathfrak{G}}`-basis::

    EXAMPLES::

        sage: G = FQSym(QQ).G();G
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the FundamentalDual basis

    So FQSym is auto dual then we define this basis as following:

    .. MATH::

        \mathbb{G}_\sigma := \mathbb{F}_{\sigma^{-1}} = \sum_{w\in
        \mathfrak{A}^*; std(w) = \sigma} w

    where `\sigma^{-1}` is the inverse of `\sigma`.

    In this basis, the product is the convolution product:

    .. MATH::

        \mathbb{G}_\sigma \star \mathbb{G}_\mu :=
        \sum_{\gamma := u\cdot v \in \mathfrak G;
            std(u) = \sigma;
            std(v) = \mu} \mathbb{G}_\gamma\;.

    EXAMPLES::

        sage: G[1,2] * G[2,1]
        G[1, 2, 4, 3] + G[1, 3, 4, 2] + G[1, 4, 3, 2] + G[2, 3, 4, 1] + G[2, 4, 3, 1] + G[3, 4, 2, 1]

    And the coproduct describes one 'un-shuffling' function:

    .. MATH::

        \Delta (\mathbb{G}_\sigma) =
        1 \otimes \mathbb{G}_\sigma + \sum_{i\in [n]}
            \mathbb{G}_{\sigma_{\mid i[}} \otimes
            \mathbb{G}_{std(\sigma_{\mid [i})}\;.

    EXAMPLES::

        sage: G[3,1,2,4].coproduct()
        G[] # G[3, 1, 2, 4] + G[1] # G[2, 1, 3] + G[1, 2] # G[1, 2] + G[3, 1, 2] # G[1] + G[3, 1, 2, 4] # G[]

    (See [MalReut]_ and [NCSF-VI]_.)

     TESTS::

        sage: G = FQSym(QQ).G(); G
        The combinatorial Hopf algebra of Free Quasi-Symmetric Functions over the Rational Field on the FundamentalDual basis
        sage: TestSuite(G).run()
    """
    _prefix = "G"

    def dual_basis(self):
        return self.realization_of().F()

    def build_morphisms(self):
        self.morph_F_G()

        ### NCSF morphisms ###
        self.morph_E_NCSF_to_G()
        self.morph_S_NCSF_to_G()
        self.morph_R_NCSF_to_G()

    def morph_R_NCSF_to_G(self):
        """
        Morphism from the non commutative ribbon Schur functions to G-basis

        TESTS::

            sage: G = FQSym(QQ).G()
            sage: R = NonCommutativeSymmetricFunctions(QQ).R()
            sage: G(R[2,1,1])
            G[1, 4, 3, 2] + G[2, 4, 3, 1] + G[3, 4, 2, 1]
            sage: G(R[2,1])
            G[1, 3, 2] + G[2, 3, 1]
            sage: G(R[1,2])
            G[2, 1, 3] + G[3, 1, 2]
        """
        R = NonCommutativeSymmetricFunctions(self.base()).R()
        R.module_morphism(
            on_basis=lambda ribbon: \
                self.sum_of_monomials(Permutations(
                    descents=(ribbon.descents(), ribbon.size())
                )),
            codomain=self
        ).register_as_coercion()

    def morph_S_NCSF_to_G(self):
        """
        Morphism from the non commutative complete functions to G-basis.

        TESTS::

            sage: G = FQSym(QQ).G()
            sage: S = NonCommutativeSymmetricFunctions( QQ ).S()
            sage: G(S[1])
            G[1]
            sage: G(S[4])
            G[1, 2, 3, 4]
            sage: G(S[1,1])
            G[1, 2] + G[2, 1]
            sage: G(S[1,2])
            G[1, 2, 3] + G[2, 1, 3] + G[3, 1, 2]
        """
        S = NonCommutativeSymmetricFunctions(self.base()).S()
        S.module_morphism(
            on_basis=lambda compo:
                prod([
                    self.monomial(Permutations(n).identity())
                    for n in compo
                ]),
            codomain=self
        ).register_as_coercion()

    def morph_E_NCSF_to_G(self):
        """
        Morphism from the non commutative elementary functions to G-basis

        TESTS::

            sage: G = FQSym(QQ).G()
            sage: E = NonCommutativeSymmetricFunctions( QQ ).Elementary()
            sage: G(E[4])
            G[4, 3, 2, 1]
            sage: G(E[3,1])
            G[3, 2, 1, 4] + G[4, 2, 1, 3] + G[4, 3, 1, 2] + G[4, 3, 2, 1]
        """
        E = NonCommutativeSymmetricFunctions(self.base()).Elementary()
        E.module_morphism(
            on_basis=lambda compo:
                prod([
                    self.monomial(Permutations(n).last())
                    for n in compo
                ]),
            codomain=self
        ).register_as_coercion()

    def morph_F_G(self):
        """
        TESTS::

            sage: G = FQSym(QQ).G()
            sage: F = FQSym(QQ).F()
            sage: F(G[1,2,3])
            F[1, 2, 3]
            sage: G(F[1,2,3])
            G[1, 2, 3]
            sage: G(F[3,1,2])
            G[2, 3, 1]
            sage: F(G[3,1,2])
            F[2, 3, 1]
        """
        F = self.realization_of().F()

        # G -> F
        self._module_morphism(
            lambda sigma: F(sigma.inverse()),
            codomain=F
        ).register_as_coercion()
        # F -> G
        F._module_morphism(
            lambda sigma: self(sigma.inverse()),
            codomain=self
        ).register_as_coercion()

    def product_on_basis(self, sigma, mu):
        """
        TESTS::

            sage: G = FQSym(QQ).G()
            sage: p = Permutations(6).random_element()
            sage: G(p) * G[[]] == G[[]] * G(p) == G(p)
            True
            sage: G[1]**2 == G.sum(G(sigma) for sigma in Permutations(2))
            True
            sage: G[1]**3 == G.sum(G(sigma) for sigma in Permutations(3))
            True
            sage: G[1] * G[2, 1, 3]
            G[1, 3, 2, 4] + G[2, 3, 1, 4] + G[3, 2, 1, 4] + G[4, 2, 1, 3]
            sage: G[1, 2] * G[2, 1]
            G[1, 2, 4, 3] + G[1, 3, 4, 2] + G[1, 4, 3, 2] + G[2, 3, 4, 1] + G[2, 4, 3, 1] + G[3, 4, 2, 1]
        """
        return self.sum_of_monomials(itertools.imap(
            lambda tau: tau.inverse(),
            sigma.inverse().shifted_shuffle(mu.inverse())
        ))

    def coproduct_on_basis(self, sigma):
        """
        TESTS::

            sage: G = FQSym(ZZ).G()
            sage: G[[]].coproduct()
            G[] # G[]
            sage: G[1].coproduct()
            G[] # G[1] + G[1] # G[]
            sage: G[1, 2].coproduct()
            G[] # G[1, 2] + G[1] # G[1] + G[1, 2] # G[]
            sage: ( G[1, 2] - G[2, 1] ).coproduct()
            G[] # G[1, 2] - G[] # G[2, 1] + G[1, 2] # G[] - G[2, 1] # G[]
        """
        Keys = self.basis().keys()

        def restrict(k):
            left = []
            right = []
            for sigi in sigma:
                if sigi > k:
                    right.append(sigi)
                else:
                    left.append(sigi)
            return (Keys(left), to_standard(right))

        return self.tensor_square().sum_of_monomials(
            restrict(i) for i in range(len(sigma) + 1)
        )

    def diese_linear_operator_on_basis(self, k, sigma):
        """
        TESTS::

            sage: G = FQSym(QQ).G()
            sage: G[1,3,2].diese_product(G[2,3,1])
            G[1, 4, 3, 5, 2] + G[1, 5, 3, 4, 2] + G[2, 4, 3, 5, 1] + G[2, 5, 3, 4, 1]

        """
        if sigma[k] == sigma[k - 1] + 1:
            return self(to_standard(sigma[:k] + sigma[k + 1:]))
        else:
            return self.zero()

    def expand_to_polynomial_on_basis(self, sigma, k):
        """
        TESTS::

            sage: G = FQSym(QQ).G()
            sage: G[3,1,2].expand_to_polynomial(3)
            a2*a1^2 + a3*a1^2 + a3*a1*a2 + a3*a2^2
            sage: G[3,1,2].expand_to_polynomial(4)
            a2*a1^2 + a3*a1^2 + a3*a1*a2 + a3*a2^2 + a4*a1^2 + a4*a1*a2 + a4*a1*a3 + a4*a2^2 + a4*a2*a3 + a4*a3^2
            sage: G[4,2,3,1].expand_to_polynomial(3)
            a3*a2^2*a1
            sage: (G[1].expand_to_polynomial(2))^2 == (G[1]^2).expand_to_polynomial(2)
            True
            sage: G[4,1,3,2].expand_to_polynomial(3)
            a3*a1*a2*a1
            sage: G[4,1,3,2].expand_to_polynomial(2)
            0
            sage: G[4,1,3,2].expand_to_polynomial(4)
            a3*a1*a2*a1 + a4*a1*a2*a1 + a4*a1*a3*a1 + a4*a1*a3*a2 + a4*a2*a3*a2
            sage: (G[4,1,3,2]*G[1,2]).expand_to_polynomial(4) == G[4,1,3,2].expand_to_polynomial(4) * G[1,2].expand_to_polynomial(4)
            True
            sage: F = FQSym(QQ).F()
            sage: F[4,1,3,2].expand_to_polynomial(4)
            a2*a3*a2*a1 + a2*a4*a2*a1 + a2*a4*a3*a1 + a3*a4*a3*a1 + a3*a4*a3*a2
            sage: F[4,1,3,2].expand_to_polynomial(4) == G(F[4,1,3,2]).expand_to_polynomial(4)
            True
        """
        from sage.algebras.free_algebra import FreeAlgebra
        FA = FreeAlgebra(self.base_ring(), ["a%d"%i for i in range(1, k+1)])

        def gen_increasing_words(l, descents):
            if l == 0: yield []
            elif l == 1:
                for i in range(k):
                    yield [i]
            else:
                if len(descents) > 0 and descents[-1]==l-1:
                    shift = 1
                    descents.pop()
                else: shift = 0
                for w in gen_increasing_words(l-1, descents):
                    for i in range(w[-1] + shift, k):
                        yield w + [i]

        inv = sigma.inverse()
        l = len(sigma)
        descents = [p+1 for p in inv.descents()]
        sigma_action = lambda w: prod([FA.gens()[w[sigma[i]-1]] for i in range(l)])

        return sum(itertools.imap(sigma_action, gen_increasing_words(l, descents)))
