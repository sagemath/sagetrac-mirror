from sage.algebras.symmetric_algebra import SymmetricAlgebra
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.sets.family import Family

class TwistedWeylGroupAlgebraElement(CombinatorialFreeModule.Element):

    def _mul_(self, other):
        par = self.parent()
        ans = par.zero()
        for (w, c) in self:
            ans += par.sum_of_terms([(w*w1,c*par.weyl_auto(w)(c1)) for (w1, c1) in other])
        return ans

class TwistedWeylGroupAlgebra(CombinatorialFreeModule):
    r"""
    EXAMPLES::

        sage: from sage.algebras.twisted_weyl_group_algebra import TwistedWeylGroupAlgebra
        sage: A = TwistedWeylGroupAlgebra(['A',2])
        sage: A.one()
        W[1]
        sage: W = A.weyl_group()
        sage: w = W.an_element(); w
        s1*s2
        sage: S = A.base_ring()
        sage: s = S.an_element(); s
        a1
        sage: A.weyl_auto(w)(s)
        a2
        sage: A(s)*A(w)
        a1*W[s1*s2]
        sage: A(w)*A(s)
        a2*W[s1*s2]
        sage: (A(s)*A(w)+A.one())**2
        a1*a2*W[s2*s1] + 2*a1*W[s1*s2] + W[1]
    """
    @staticmethod
    def __classcall__(cls, cartan_type):
        cartan_type = CartanType(cartan_type)
        return super(TwistedWeylGroupAlgebra, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        self._cartan_type = cartan_type
        if cartan_type.is_affine():
            self._lattice = cartan_type.root_system().root_lattice(extended=True)
        else:
            self._lattice = cartan_type.root_system().root_lattice()
        self._base_ring = SymmetricAlgebra(self._lattice, prefix="a")
        self._weyl_group = self._lattice.weyl_group(prefix="s")

        # set up the action of the weyl group on the lattice and its symmetric algebra

        CombinatorialFreeModule.__init__(self, self._base_ring, self._weyl_group, category=AlgebrasWithBasis(self._base_ring), prefix="W")
        
        self._S_to_self = SetMorphism(Hom(self._base_ring, self, Sets()), lambda x: self.sum_of_terms([(self.one_basis(),x)], distinct=True))
        self._S_to_self.register_as_coercion()
        self._W_to_self = SetMorphism(Hom(self._weyl_group, self, Sets()), lambda w: self.monomial(w))
        self._W_to_self.register_as_coercion()

    def weyl_group(self):
        return self._weyl_group

    def base_ring(self):
        return self._base_ring

    @cached_method
    def weyl_auto(self, w):
        r"""
        This requires that the Weyl group has been defined as the Weyl group of the lattice.
        """
        alg = self.base_ring()
        return alg.algebra_morphism(Family(alg.keys(), lambda k: alg.from_module_map()(w.action(alg.base_module().basis()[k]))))

    @cached_method
    def one_basis(self):
        return self.weyl_group().one()

    Element = TwistedWeylGroupAlgebraElement

