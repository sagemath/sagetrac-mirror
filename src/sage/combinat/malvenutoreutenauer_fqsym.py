from sage.categories.rings import Rings
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.words.word import Word
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.realizations import Category_realization_of_parent
from sage.misc.bindable_class import BindableClass

class MalvenutoReutenauer(UniqueRepresentation, Parent):
    def __init__(self, R):
        assert(R in Rings())
        self._base = R
        from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
        Parent.__init__(self, category = GradedHopfAlgebras(R).WithRealizations())

    def _repr_(self): # could be taken care of by the category
        return "The Malvenuto-Reutenaur Hopf Algebra of Permutations "+\
               "(FQSym) over the %s"%self.base_ring()

    def a_realization(self):
        return self.Fundamental()

    class Bases(Category_realization_of_parent):
        # this is where all the code for all bases goes
        def super_categories(self):
            r"""
            """
            R = self.base().base_ring()
            from sage.categories.graded_modules import GradedModules
            from sage.categories.graded_hopf_algebras import GradedHopfAlgebras
            from sage.categories.graded_hopf_algebras_with_basis import \
              GradedHopfAlgebrasWithBasis
            # NB: self.base() is the CHA
            return [self.base().Realizations(),
                    GradedModules(R).Realizations(),
                    GradedHopfAlgebras(R).Realizations(),
                    GradedHopfAlgebrasWithBasis(R)]

        class ParentMethods:
            def _repr_(self):
                return "%s in the %s basis" % (self.realization_of(), self._realization_name())
            def __getitem__(self, p, *rest):
                # This method implements the abuses of notations e.g. F[1,2,3]
                if not isinstance(p, Permutation):
                    p = self._indices(p)
                assert len(rest)==0
                return self.monomial(p)
            def one_basis(self):
                return Permutations()([])
            def degree_on_basis(self, p):
                return p.size()
            def counit_on_basis(self, p):
                if p != []:
                    return self.base_ring().zero()
                else:
                    return self.base_ring().one()

            def to_quasisymmetric_function_on_basis(self, p):
                F = self.realization_of().a_realization()
                return F.to_quasisymmetric_function(F(self(p)))

            def to_quasisymmetric_function(self, elt):
                F = self.realization_of().a_realization()
                return F.to_quasisymmetric_function(F(elt))

            def duality_pairing_by_coercion(self, right, left):
                F = self.realization_of().a_realization()
                return F.duality_pairing(right, left)

        class ElementMethods:
            def to_quasisymmetric_function(self):
                return self.parent().to_quasisymmetric_function(self)
            def duality_pairing(self, right):
                return self.parent().duality_pairing_by_coercion(self, right)

    class Fundamental(CombinatorialFreeModule, BindableClass):
        def __init__(self, my_CHA):
            CombinatorialSet = Permutations()
            CombinatorialFreeModule.__init__(self, my_CHA.base_ring(),
                                             CombinatorialSet,
                                             prefix='Fu', bracket=False,
                                             category=my_CHA.Bases())
        def _realization_name(self):
            return "Fundamental"
        def product_on_basis(self, p, q):
            return sum(self(r) for r in p.shifted_shuffle(q))
        def coproduct_on_basis(self, p):
            T=self.tensor_square()
            st = lambda pp: Word(pp).standard_permutation()
            return T.sum_of_monomials((st(p[:i]),st(p[i:])) for i in range(len(p)+1))
        def antipode_on_basis(self, p):
            if p==self.one_basis():
                return self.one()
            S = self.antipode_on_basis
            comul = self.coproduct_on_basis
            return -self.sum(c*self(pp)*S(qq) for ((pp,qq),c) in comul(p) if pp!=self.one_basis())

        def to_quasisymmetric_function(self, elt):
            FQS = QuasiSymmetricFunctions(self.base_ring()).F()
            return FQS.sum(c*FQS(p.descents_composition()) for (p,c) in elt)
        def duality_pairing(self, left, right):
            if left.parent()!=self:
                left = self(left)
            if right.parent()!=self:
                right= self(right)
            return sum(c*right.coefficient(p.inverse()) for (p,c) in left)

        class Element(CombinatorialFreeModule.Element):
            pass

    class Monomial(CombinatorialFreeModule, BindableClass):
        def __init__(self, my_CHA):
            CombinatorialSet = Permutations()
            CombinatorialFreeModule.__init__(self, my_CHA.base_ring(),
                                             CombinatorialSet,
                                             prefix='Mu', bracket=False,
                                             category=my_CHA.Bases())
            F = self.realization_of().a_realization()
            from sage.categories.graded_hopf_algebras_with_basis \
              import GradedHopfAlgebrasWithBasis
            categ = GradedHopfAlgebrasWithBasis(self.base_ring())
            phi = F.module_morphism(self._from_Fundamental_on_basis,
                                    codomain=self, triangular="lower",
                                    unitriangular=True, category=categ)
            self.register_coercion(phi)
            F.register_coercion(~phi)
        def _realization_name(self):
            return "Monomial"
        def _from_Fundamental_on_basis(self, p):
            return self.sum(self.monomial(q) for q in
              p.permutohedron_greater(side='left'))

        class Element(CombinatorialFreeModule.Element):
            pass

    class FundamentalDual(CombinatorialFreeModule, BindableClass):
        def __init__(self, my_CHA):
            CombinatorialSet = Permutations()
            CombinatorialFreeModule.__init__(self, my_CHA.base_ring(),
                                             CombinatorialSet,
                                             prefix='Fd', bracket=False,
                                             category=my_CHA.Bases())
            F = self.realization_of().a_realization()
            from sage.categories.graded_hopf_algebras_with_basis \
              import GradedHopfAlgebrasWithBasis
            categ = GradedHopfAlgebrasWithBasis(self.base_ring())
            F.module_morphism(self._from_Fundamental_on_basis,
                                codomain=self, category=categ).register_as_coercion()
            self.module_morphism(self._to_Fundamental_on_basis,
                                codomain=F, category=categ).register_as_coercion()
        def _realization_name(self):
            return "Fundamental dual"
        def _from_Fundamental_on_basis(self, p):
            return self.monomial(p.inverse())
        def _to_Fundamental_on_basis(self, p):
            F = self.realization_of().a_realization()
            return F.monomial(p.inverse())

        class Element(CombinatorialFreeModule.Element):
            pass

    class MonomialDual(CombinatorialFreeModule, BindableClass):
        def __init__(self, my_CHA):
            CombinatorialSet = Permutations()
            CombinatorialFreeModule.__init__(self, my_CHA.base_ring(),
                                             CombinatorialSet,
                                             prefix='Md', bracket=False,
                                             category=my_CHA.Bases())
            Fd = self.realization_of().FundamentalDual()
            from sage.categories.graded_hopf_algebras_with_basis \
              import GradedHopfAlgebrasWithBasis
            categ = GradedHopfAlgebrasWithBasis(self.base_ring())
            phi = self.module_morphism(self._to_Fundamental_dual_on_basis,
                                    codomain=Fd, triangular="upper",
                                    unitriangular=True, category=categ)
            Fd.register_coercion(phi)
            self.register_coercion(~phi)
            
        def _realization_name(self):
            return "Monomial dual"
        def _to_Fundamental_dual_on_basis(self, p):
            Fd = self.realization_of().FundamentalDual()
            return Fd.sum(Fd(q) for q in p.permutohedron_smaller(side='left'))

        class Element(CombinatorialFreeModule.Element):
            pass

    Fu = Fundamental
    Mu = Monomial
    Fd = FundamentalDual
    Md = MonomialDual

FQSym = MalvenutoReutenauer
