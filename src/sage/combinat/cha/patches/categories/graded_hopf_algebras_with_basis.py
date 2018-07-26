"""
This patch implement a default implementation of the antipode for Graded (and
connexe) Hopf algebras.
"""
from sage.combinat.cha.patches.monkey_patching import MonkeyPatch
from sage.misc.lazy_import import LazyImport

from sage.categories.graded_hopf_algebras_with_basis import \
        GradedHopfAlgebrasWithBasis
from sage.misc.lazy_attribute import lazy_attribute
from sage.categories.realizations import RealizationsCategory


class _(MonkeyPatch, GradedHopfAlgebrasWithBasis):

    Realizations = LazyImport(
        "sage.combinat.cha.patches.categories.graded_hopf_algebras_with_basis",
        "Realizations")


class Realizations(RealizationsCategory):

    class ParentMethods:
        def antipode_by_coercion(self, left, right):
            r"""
            """
            for R in self.realization_of().realizations():
                self_to_R = R.coerce_map_from(self)
                R_to_self = self.coerce_map_from(R)
                if self_to_R is not None and R_to_self is not None and \
                   R.antipode != R.antipode_by_coercion:
                    return self(R(left).antipode())
            return NotImplementedError


class _(MonkeyPatch, GradedHopfAlgebrasWithBasis.ParentMethods):

    def antipode_on_basis(self, indice):
        """
        MATH::

            S(x) := -\sum_{x^L \neq x} S(x^L) \times x^R

        in general or `x` if `\mid x \mid = 0`.

        TESTS::

            sage: F = FQSym(QQ).F()
            sage: F.antipode_on_basis(Permutation([1]))
            -F[1]

        """
        if indice.size() == 0:
            return self.one()
        else:
            from sage.categories.tensor import tensor
            S = self.antipode_on_basis
            x__S_Id = tensor([self, self]).module_morphism(
                lambda (a, b): S(a) * self.monomial(b),
                codomain=self)
            return -x__S_Id(
                self.monomial(indice).coproduct()
                - tensor([self(indice), self.one()])
            )

    @lazy_attribute
    def antipode(self):
        r"""
        """
        if self.antipode_on_basis is NotImplemented:
            return self.antipode_by_coercion
        else:
            return self.module_morphism(
                self.antipode_on_basis,
                codomain=self
            )


class _(MonkeyPatch, GradedHopfAlgebrasWithBasis.ElementMethods):

    def antipode(self):
        return self.parent().antipode(self)
