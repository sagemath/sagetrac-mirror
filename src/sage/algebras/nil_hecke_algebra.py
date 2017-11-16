from sage.algebras.symmetric_algebra import SymmetricAlgebra
from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets
from sage.sets.family import Family
from sage.categories.rings import Rings

class NilHeckeAlgebraElement(CombinatorialFreeModule.Element):

    def _mul_(self, other):
        r"""
        Multiply ``self`` on the right by ``other``.

        EXAMPLES::

            sage: from sage.algebras.nil_hecke_algebra import NilHeckeAlgebra
            sage: H = NilHeckeAlgebra(['A',2])
            sage: H.A(1)**2
            0
            sage: H.A(1)*H.A(2)
            A[s1*s2]
            sage: S = H.base_ring()
            sage: S.simple_root(1)**2 * H.A(2)
            a1^2*A[s2]
            sage: S.simple_root(1)**2*H.A([2,1])
            a1^2*A[s2*s1]
        """
        par = self.parent()
        ans = par.zero()
        for (w, c) in self:
            ans += c * par.divided_difference_weyl(other, w)
        return ans

class NilHeckeAlgebra(CombinatorialFreeModule):
    r"""
    EXAMPLES::

        sage: from sage.algebras.nil_hecke_algebra import NilHeckeAlgebra
        sage: H = NilHeckeAlgebra(['A',2])
        sage: H.one()
        A[1]
        sage: S = H.base_ring()
        sage: s = S.an_element(); s
        a1
        sage: H(s)
        a1*A[1]
        sage: W = H.weyl_group()
        sage: w = W.an_element(); w
        s1*s2
        sage: H.A(w)
        A[s1*s2]
        sage: H.A(2)*H.A(2)
        0
        sage: H.A(1)*(H(s**2))
        a1^2*A[s1]
    """
    @staticmethod
    def __classcall__(cls, cartan_type):
        cartan_type = CartanType(cartan_type)
        return super(NilHeckeAlgebra, cls).__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        self._base_ring = EquivariantCohomologyPoint(cartan_type, lattice='root')

        CombinatorialFreeModule.__init__(self, self._base_ring, self.weyl_group(), category=ModulesWithBasis(self._base_ring)&Rings(), prefix="A")

        self._S_to_self = SetMorphism(Hom(self._base_ring, self, Sets()), lambda x: self.sum_of_terms([(self.one_basis(),x)], distinct=True))
        self._S_to_self.register_as_coercion()

    def base_ring(self):
        return self._base_ring

    def cartan_type(self):
        return self._base_ring.cartan_type()

    def weyl_group(self):
        return self._base_ring.weyl_group()

    def divided_difference_operator(self, f, i):
        r"""
        Apply the `i`-th divided difference operator to the element `f`.

        EXAMPLES::

            sage: from sage.algebras.nil_hecke_algebra import NilHeckeAlgebra
            sage: M = NilHeckeAlgebra(CartanType(['A',2]))
            sage: a = M.base_ring().algebra_generators(); a
            Finite family {1: a1, 2: a2}
            sage: s = M.weyl_group().simple_reflections()
            sage: f = a[1]*M.A(2) +  a[1]*a[2]*M.A([2,1])
            sage: M.divided_difference_operator(f, 2)
            (2*a1+a2)*A[s2*s1] + (-1)*A[s2]
            sage: M.divided_difference_operator(f, 1)
            (a1+2*a2)*A[s2*s1] + (-a1^2-a1*a2)*A[s1*s2*s1] + (-a1)*A[s1*s2] + 2*A[s2]
        """
        terms = []
        S = self.base_ring()
        for (w, c) in f:
            terms += [(w, S.divided_difference(c, i))]
            if not w.has_descent(i, side='left'):
                terms += [(w.apply_simple_reflection(i, side='left'), S.simple_reflection_automorphism(i)(c))]
        return self.sum_of_terms(terms)

    def divided_difference_weyl(self, f, w):
        r"""
        EXAMPLES::

            sage: from sage.algebras.nil_hecke_algebra import NilHeckeAlgebra
            sage: M = NilHeckeAlgebra(CartanType(['A',2]))
            sage: s = M.weyl_group().simple_reflections()
            sage: a = M.base_ring().algebra_generators()
            sage: f = a[1]*M.A(2) +  a[1]*a[2]*M.A([2,1])
            sage: M.divided_difference_weyl(f, s[1]*s[2])
            3*A[s2*s1] + (-a1+a2)*A[s1*s2*s1] + (-1)*A[s1*s2]
        """
        rw = w.reduced_word()
        for i in reversed(rw):
            f = self.divided_difference_operator(f,i)
        return f

    @cached_method
    def one_basis(self):
        return self.weyl_group().one()

    def A(self, i):
        r"""
        The basis element indexed by `i`.

        INPUT:

        - `i` -- Either a Weyl group element, or an index for a simple reflection,
          or a reduced word. If the word is not reduced the corresponding element is zero.

        EXAMPLES::

            sage: from sage.algebras.nil_hecke_algebra import NilHeckeAlgebra
            sage: H = NilHeckeAlgebra(['A',2])
            sage: H.A(2)
            A[s2]
            sage: H.A([1,2])
            A[s1*s2]
            sage: H.A([1,1])
            0
            sage: H.A([])
            A[1]
            sage: H.A(H.weyl_group().an_element())
            A[s1*s2]
        """
        W = self.weyl_group()
        if i in self.cartan_type().index_set():
            return self.monomial(W.simple_reflection(i))
        if i in W:
            return self.monomial(i)
        w = W.from_reduced_word(i)
        if w.length() == len(i):
            return self.monomial(W.from_reduced_word(i))        
        else:
            return self.zero()

    def weyl_element(self, i):
        r"""
        The image of a Weyl group element in ``self``.
 
        INPUT:

        - `i` -- Either a Weyl group element, or an index for a simple reflection,
          or a reduced word. If the word is not reduced, unlike the behavior of :meth:`A`,
          the product of the given simple reflections is used.

        EXAMPLES::

            sage: from sage.algebras.nil_hecke_algebra import NilHeckeAlgebra
            sage: H = NilHeckeAlgebra(['A',2])
            sage: H.weyl_element(1)
            (-a1)*A[s1] + A[1]
            sage: H.weyl_element([2,1])
            (a1*a2+a2^2)*A[s2*s1] + (-a1-a2)*A[s1] + A[1] + (-a2)*A[s2]
            sage: H.weyl_element([1,1])
            A[1]
            sage: H.weyl_element([])
            A[1]
            sage: H.weyl_element(H.weyl_group().simple_reflection(2))
            A[1] + (-a2)*A[s2]
        """
        W = self.weyl_group()
        if i in self.cartan_type().index_set():
            i = W.simple_reflection(i)
        elif i in W:
            pass
        else:
            i = W.from_reduced_word(i)
        return self._weyl_element(i)

    @cached_method
    def _weyl_element(self, w):
        W = self.weyl_group()
        if w == W.one():
            return self.one()
        i = w.first_descent(side='right')
        return self._weyl_element(w.apply_simple_reflection(i, side='right')) * (self.one() - self.sum_of_terms([(W.simple_reflection(i),self.base_ring().simple_root(i))]))

    Element = NilHeckeAlgebraElement

