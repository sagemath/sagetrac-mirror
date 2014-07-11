from sage.combinat.permutation import Permutations, to_standard
from sage.combinat.hopf_algebras import GenericGradedConnexeHopfAlgebra, \
    register_as_realization


class SimpleFQSym(GenericGradedConnexeHopfAlgebra):

    def _repr_(self):
        return "FQSym"

    _default_basis_indices_ = Permutations()

class FBasis(SimpleFQSym._Basis):

    _prefix_ = "F"

    def product_on_basis(self, sigma, mu):
        return self.sum_of_monomials(sigma.shifted_shuffle(mu))

    def coproduct_on_basis(self, sigma):
        return self.tensor_square().sum_of_monomials(
            (to_standard(sigma[:i]), to_standard(sigma[i:]))
            for i in range(sigma.grade() + 1)
        )

class GBasis(SimpleFQSym._Basis):

    _prefix_ = "G"

    def _morphisms_(self):
        F = self.realization_of().F()

        F.module_morphism(
            on_basis=lambda sigma: self(sigma.inverse()),
            codomain=self
        ).register_as_coercion()

        self.module_morphism(
            on_basis=lambda sigma: F(sigma.inverse()),
            codomain=F
        ).register_as_coercion()


###############################################################
### should be automatically made by inheritage ###            #
#### that could be do in a __metaclass__ but ...              #
#### the function *dynamiclass* use by category fails!!!      #
register_as_realization(SimpleFQSym, FBasis, "F")             #
register_as_realization(SimpleFQSym, GBasis, "G")             #
###############################################################