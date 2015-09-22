from sage.combinat.permutation import Permutations, to_standard
from sage.combinat.hopf_algebras import GenericGradedConnectedHopfAlgebra, \
    register_as_realization


class SimpleFQSym(GenericGradedConnectedHopfAlgebra):
    """
    TESTS::

        sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
        sage: TestSuite(SimpleFQSym(QQ)).run()

    """

    def _repr_(self):
        """
        TESTS::

            sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
            sage: SimpleFQSym(QQ)
            FQSym

        """
        return "FQSym"

    _default_basis_indices_ = Permutations()

class FBasis(SimpleFQSym._Basis):
    """
    TESTS::

        sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
        sage: F = SimpleFQSym(QQ).F()
        sage: TestSuite(F).run()

    """

    _prefix_ = "F"

    def product_on_basis(self, sigma, mu):
        """

        TESTS::

            sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
            sage: F = SimpleFQSym(QQ).F()
            sage: sage: F.monomial(Permutation([1,2])) * F.monomial(Permutation([2,1]))
            F[1, 2, 4, 3] + F[1, 4, 2, 3] + F[1, 4, 3, 2] + F[4, 1, 2, 3] + F[4, 1, 3, 2] + F[4, 3, 1, 2]
        """
        return self.sum_of_monomials(sigma.shifted_shuffle(mu))

    def coproduct_on_basis(self, sigma):
        """
        TESTS::

            sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
            sage: F = SimpleFQSym(QQ).F()
            sage: F.monomial(Permutation([3,2,1,4])).coproduct()
            F[] # F[3, 2, 1, 4] + F[1] # F[2, 1, 3] + F[2, 1] # F[1, 2] + F[3, 2, 1] # F[1] + F[3, 2, 1, 4] # F[]

        """
        return self.tensor_square().sum_of_monomials(
            (to_standard(sigma[:i]), to_standard(sigma[i:]))
            for i in range(sigma.grade() + 1)
        )

class GBasis(SimpleFQSym._Basis):
    """
    TESTS::

        sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
        sage: G = SimpleFQSym(QQ).G()
        sage: TestSuite(G).run()


    """

    _prefix_ = "G"

    def _morphisms_(self):
        """
        TESTS::

            sage: from sage.combinat.hopf_algebras.examples.fqsym import SimpleFQSym
            sage: fqs = SimpleFQSym(QQ)
            sage: F, G = fqs.F(), fqs.G()
            sage: F(G.monomial(Permutation([3,1,2])))
            F[2, 3, 1]
            sage: G(F.monomial(Permutation([3,1,2])))
            G[2, 3, 1]

            sage: G.monomial(Permutation([3,1,2])) * G.monomial(Permutation([1]))
            G[3, 1, 2, 4] + G[4, 1, 2, 3] + G[4, 1, 3, 2] + G[4, 2, 3, 1]

            sage: G.monomial(Permutation([3,1,2])).coproduct()
            G[] # G[3, 1, 2] + G[1] # G[2, 1] + G[1, 2] # G[1] + G[3, 1, 2] # G[]

        """
        F = self.realization_of().F()

        F.module_morphism(
            on_basis=lambda sigma: self.monomial(sigma.inverse()),
            codomain=self
        ).register_as_coercion()

        self.module_morphism(
            on_basis=lambda sigma: F.monomial(sigma.inverse()),
            codomain=F
        ).register_as_coercion()


###############################################################
# IMPORTANT:: The following lines register the bases in FQSym #
###############################################################
register_as_realization(SimpleFQSym, FBasis, "F")             #
register_as_realization(SimpleFQSym, GBasis, "G")             #
###############################################################