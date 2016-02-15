"""
Smash product of algebras

Author: Mark Shimozono
"""
#*****************************************************************************
#  Copyright (C) 2014 Mark Shimozono <mshimo at math.vt.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import partial
from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.all import AlgebrasWithBasis, ModulesWithBasis
from sage.categories.category import Category
from sage.categories.morphism import SetMorphism
from sage.categories.tensor import tensor
from sage.categories.homset import Hom
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModule_Tensor, \
CombinatorialFreeModule_TensorGrouped, CartesianProductWithFlattening, CartesianProductWithUnflattening

class SmashProductAlgebraElement(CombinatorialFreeModule.Element):
    r"""
    Element class for :class:`SmashProductAlgebra`.
    """
    def to_opposite(self):
        r"""
        Map ``self`` to the smash product with factors in opposite order.

        This uses coercion.

        In the following example, `A` is the polynomial ring in one variable (realized by a free monoid),
        which represents the fundamental weight in the weight lattice of `SL_2` and `B` is the 
        Weyl group of `SL_2`, which is the symmetric group on two letters. The reflection acts on the
        fundamental weight by negation.

        EXAMPLES::
        
            sage: F = FreeMonoid(names=['a'])
            sage: A = F.algebra(ZZ); A.rename("A")
            sage: A._print_options['prefix'] = 'A'
            sage: W = WeylGroup("A1",prefix="s")
            sage: B = W.algebra(ZZ); B.rename("B")
            sage: cat = ModulesWithBasis(ZZ)
            sage: AB = tensor([A,B],category=cat)
            sage: BA = tensor([B,A],category=cat)
            sage: def twist_func((b,a)):
            ...       if b != b.parent().one() and mod(a.to_word().length(),2) == 1:
            ...           return -AB.monomial((a,b))
            ...       return AB.monomial((a,b))
            sage: def untwist_func((b,a)):
            ...       if a != a.parent().one() and mod(b.to_word().length(),2) == 1:
            ...           return -BA.monomial((a,b))
            ...       return BA.monomial((a,b))
            sage: ASB = SmashProductAlgebra(A, B, twist_on_basis=twist_func)
            sage: BSA = SmashProductAlgebra(B, A, untwist_func)
            sage: ASB.register_opposite(BSA)
            sage: ab = ASB.an_element(); ab
            3*A[a] B[s1] + A[a]
            sage: BSA = ASB.opposite()
            sage: ba = BSA(ab); ba
            -3*B[s1] A[a] + A[a]
            sage: ASB(ba) == ab
            True
            sage: ASB == BSA.opposite()
            True
            sage: ba == ab.to_opposite()
            True
            sage: A._print_options['prefix'] = 'B'

        """
        return self.parent().opposite()(self)

class SmashProductAlgebra(CombinatorialFreeModule_TensorGrouped):
    r"""
    Smash products of algebras.

    INPUT:

        - `A` and `B` -- algebras with basis over the same commutative ring
        - A specification of the product. There are several options for this; see below
        - ``category`` -- optional (default: None) category for resulting algebra. This should be a subcategory of
          :class:`AlgebrasWithBasis(R).TensorProducts()` where `R` is the base ring.
        - ``suppress_ones`` -- optional (default: None) The default printing of tensor monomials
          will suppress tensor factors of the identity.

    Returns the smash product `S(A,B,twist)`. As a module it is `A` tensor `B`.
    The product is determined by a module homomorphism ``twist`` from `B` tensor `A`
    to `A` tensor `B` as described below.

    Particular rings that are realized as smash products are the cohomological and K-theoretic nilHecke rings
    of Kostant and Kumar, the affine Hecke algebra, and the double affine Hecke algebra. These all come from actions
    of one algebra on another. Another class of examples comes from quasi-triangular Hopf algebras and the R-matrix.

    ..RUBRIC::

    Let `R` be the common base ring of `A` and `B`. All tensor products are taken over `R`.
    A smash product `(A,B,twist)` is an `R`-algebra structure imposed on `A` tensor `B` such that
    `A` and `B` are subalgebras in the natural way, whose multiplication is defined
    via an `R`-module homomorphism ``twist`` from `B` tensor `A` to `A` tensor `B` 
    as follows. Let `I_A` be the identity map on `A` and `m_A` be the multiplication map from
    `A` tensor `A` to `A`. Form the map from `A` tensor `B` tensor `A` tensor `B` to
    `A` tensor `B` by the composition

    ..MATH::

        m = (m_A tensor m_B) \circ (I_A tensor twist tensor I_B)

    The reference [CIMZ]_ gives necessary and sufficient conditions on ``twist`` for `m` to
    define an `R`-algebra structure on `A` tensor `B`.

    The homomorphism ``twist`` can be specified in several ways. The first way is to provide
    a function ``twist_on_basis`` which has as input an ordered pair `(b,a)` from the index sets of `B` and `A`,
    with values in `A \otimes B`.

        sage: W = WeylGroup(['A',2],prefix="s")
        sage: A = W.algebra(ZZ)
        sage: AA = tensor([A,A],category=ModulesWithBasis(ZZ).TensorProducts())
        sage: twist_func = lambda (b,a): AA.monomial((b*a*b.inverse(),b))
        sage: B = SmashProductAlgebra(A,A,twist_on_basis=twist_func, suppress_ones=False)
        sage: b = B.from_direct_product((A.an_element(), A.one())); b
        2*B[s1*s2*s1] # B[1] + 4*B[s1*s2] # B[1] + B[1] # B[1]
        sage: b*b
        4*B[s1*s2*s1] # B[1] + 8*B[s1*s2] # B[1] + 16*B[s2*s1] # B[1] + 8*B[s1] # B[1] + 8*B[s2] # B[1] + 5*B[1] # B[1]

    The second method is to provide the actual morphism ``twist_morphism``::

        sage: twist = AA.module_morphism(on_basis=twist_func, codomain=AA)
        sage: C = SmashProductAlgebra(A,A,twist_morphism=twist, suppress_ones=False)
        sage: c = C.from_direct_product((A.an_element(), A.one())); c
        2*B[s1*s2*s1] # B[1] + 4*B[s1*s2] # B[1] + B[1] # B[1]
        sage: c*c
        4*B[s1*s2*s1] # B[1] + 8*B[s1*s2] # B[1] + 16*B[s2*s1] # B[1] + 8*B[s1] # B[1] + 8*B[s2] # B[1] + 5*B[1] # B[1]

    The third method is to give (on the bases) a left action of `B` on `A \otimes B`.
    This function takes three arguments `bb, a, b`::

        sage: left_action_func = lambda bb,a,b: AA.monomial((bb*a*bb.inverse(), bb*b))
        sage: D = SmashProductAlgebra(A,A,left_action=left_action_func, suppress_ones=False)
        sage: d = D.from_direct_product((A.an_element(), A.one())); d
        2*B[s1*s2*s1] # B[1] + 4*B[s1*s2] # B[1] + B[1] # B[1]
        sage: d*d
        4*B[s1*s2*s1] # B[1] + 8*B[s1*s2] # B[1] + 16*B[s2*s1] # B[1] + 8*B[s1] # B[1] + 8*B[s2] # B[1] + 5*B[1] # B[1]

    The product is then determined by letting `A` act on the left factor by left multiplication.

    The final method is to give (on the bases) a right action of `A` on `A \otimes B`::
    This function takes three arguments `aa, a, b`::

        sage: right_action_func = lambda aa, a, b: AA.monomial((a*aa, aa.inverse()*b*aa))
        sage: E = SmashProductAlgebra(A,A,right_action=right_action_func, suppress_ones=False)
        sage: e = E.from_direct_product((A.an_element(), A.one())); d
        2*B[s1*s2*s1] # B[1] + 4*B[s1*s2] # B[1] + B[1] # B[1]
        sage: e*e
        4*B[s1*s2*s1] # B[1] + 8*B[s1*s2] # B[1] + 16*B[s2*s1] # B[1] + 8*B[s1] # B[1] + 8*B[s2] # B[1] + 5*B[1] # B[1]

    The product is then determined by letting `B` act on the right factor by right multiplication.

    REFERENCE:

        .. [CIMZ] Caenepeel, S.; Ion, Bogdan; Militaru, G.; Zhu, Shenglin. The factorization problem
          and the smash biproduct of algebras and coalgebras. Algebr. Represent. Theory 3 (2000),
          no. 1, 19--42. 

    The tensor product of algebras is implemented as a grouped tensor product via the concrete class
    :class:`sage.combinat.free_module.CombinatorialFreeModule_TensorGrouped` while tensor products of modules
    are automatically flattened using :class:`sage.combinat.free_module.CombinatorialFreeModule_Tensor`.
    In reality the tensor product of algebras is also a flat tensor product but the class remembers the grouping
    information. The index set of the basis is given by the "flattened" tuples. 

    In the following example both `A` and `B` are the group algebra of the symmetric group `S_3` over the integers,
    and `B` has a left action on `A` induced by conjugation of symmetric group elements.

    EXAMPLES::

        sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
        sage: r = W.from_reduced_word
        sage: A = W.algebra(ZZ); A.rename("A")
        sage: AA = tensor([A,A], category=ModulesWithBasis(ZZ))
        sage: twist = AA.module_morphism(on_basis=lambda (b,a): AA.monomial((b*a*b.inverse(),b)),codomain=AA)
        sage: ba = AA.monomial((r([1]),r([2]))); ba
        B[s1] # B[s2]
        sage: ab = twist(ba); ab
        B[s1*s2*s1] # B[s1]
        sage: C = SmashProductAlgebra(A, A, twist_morphism=twist, suppress_ones=False); C
        Smash product of A and A
        sage: c = A.monomial(A.basis().keys().an_element()); c
        B[s1*s2]
        sage: d = A.an_element(); d
        2*B[s1*s2*s1] + 4*B[s1*s2] + B[1]
        sage: cd = C.from_direct_product((c,d)); cd
        2*B[s1*s2] # B[s1*s2*s1] + 4*B[s1*s2] # B[s1*s2] + B[s1*s2] # B[1]
        sage: cd * cd
        2*B[s2*s1] # B[s1*s2*s1] + 8*B[s2*s1] # B[s1*s2] + 16*B[s2*s1] # B[s2*s1] + 8*B[s2*s1] # B[s2] + B[s2*s1] # B[1] + 2*B[1] # B[s1*s2*s1] + 8*B[1] # B[s1] + 4*B[1] # B[1]

    ..RUBRIC:: Opposite order smash product

    If the map ``twist`` has an inverse ``untwist``, then ``untwist`` is an algebra isomorphism from
    `S(A,B,twist)` to `S(B,A,untwist)`. To set up coercions between these two algebras, first create both algebras
    (call them ``ASB`` and ``BSA`` for short), and then invoke `ASB.register_opposite(BSA)`.

    ..IMPLEMENTATION:: The smash product is the flattened tensor product of `A` and `B`. However it knows its
    tensor factors. Coercion is immediately allowed between this underlying flattened tensor product module and the smash product,
    with no additional definitions. Intentionally, the coercion registration is made entirely separate from the creation of the
    algebras, for convenience of use by algebras with realizations which carry additional category information and use both
    forms of the smash product.

    """

    def __init__(self, A, B, twist_on_basis=None, twist_morphism=None, left_action=None, right_action=None, category=None, suppress_ones=None):
        """
        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: B = A
            sage: BA = tensor([B,A],category=ModulesWithBasis(ZZ))
            sage: AB = BA
            sage: def t_map((b,a)):
            ...       return AB.monomial((a,b))
            sage: C = SmashProductAlgebra(A,B,twist_on_basis=t_map); C # indirect doctest
            Smash product of A and A
        """
        R = A.base_ring()
        module_category = ModulesWithBasis(R)
        algebra_category = AlgebrasWithBasis(R)
        tensor_category = module_category.TensorProducts()
        if R != B.base_ring():
            raise TypeError("%s and %s must have the same base rings: \n %s \n %s\n"%(A,B,A.base_ring(),B.base_ring()))
        if A not in algebra_category or B not in algebra_category:
            raise TypeError("Tensor factors should be AlgebrasWithBasis over %s"%(R))

        if suppress_ones is None:
            self._suppress_ones = True
        elif suppress_ones not in (True,False):
            raise ValueError("%s should be True or False"%suppress_ones)
        else:
            self._suppress_ones = suppress_ones

        default_category = AlgebrasWithBasis(R).TensorProducts()
        if not category:
            category = default_category
        elif not category.is_subcategory(default_category):
            raise TypeError("Category %s is not a subcategory of %s"%(category, default_category))
        self._category = category

        AB = tensor([A,B], category=tensor_category)
        BA = tensor([B,A], category=tensor_category)

        CombinatorialFreeModule_TensorGrouped.__init__(self, (A,B), category=category)

        if left_action is not None:
            twist_on_basis = lambda (c,a): left_action(c, a, B.one_basis())
        elif right_action is not None:
            twist_on_basis = lambda (b,d): right_action(d, A.one_basis(), b)
        if twist_on_basis is not None:
            key_splitter = tensor([B,A]).index_to_indices()
            twist_on_basis_func = lambda x: twist_on_basis(key_splitter(x))
            self._twist = BA.module_morphism(on_basis=twist_on_basis_func, codomain=AB, category=module_category)
        elif twist_morphism is not None:
            if twist_morphism.domain() != BA:
                raise TypeError("Codomain of twist (%s) should be %s"%(twist.codomain(), BA))
            if twist_morphism.codomain() != AB:
                raise TypeError("Codomain of twist (%s) should be %s"%(twist.codomain(), AB))
            if not twist_morphism.category().is_subcategory(module_category.Homsets()):
                raise TypeError("twist should be a module morphism")
            self._twist = twist_morphism
        else:
            raise TypeError("Twist is improperly specified")

        self._has_opposite = False # default: no "opposite" smash product has been registered
        # Define the product morphism using twist
        self._ItwistI = tensor([A._identity_map(), self._twist, B._identity_map()], category=tensor_category)
        self._mAmB = tensor([A._product_morphism(),B._product_morphism()], category=tensor_category)
        self._product_morphism_map = SetMorphism(Hom(self._ItwistI.domain(), self._mAmB.codomain(), category=tensor_category), self._mAmB * self._ItwistI)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ)
            sage: B = A
            sage: BA = tensor([B,A],category=ModulesWithBasis(ZZ))
            sage: AB = BA
            sage: def t_map((b,a)):
            ...       return AB.monomial((a,b))
            sage: C = SmashProductAlgebra(A,B,twist_on_basis=t_map); C
            Smash product of A and A

            Question: Should this name include twist?

        """
        return "Smash product of %s and %s"%(self.factor(0), self.factor(1))

    def _repr_term(self, term):
        r"""
        A string for a term.
        """
        # remember the grouping on the tuple for printing
        term = self.index_to_indices()(term)
        if self._suppress_ones:
            if term[0] == self.factor(0).one_basis():
                if term[1] == self.factor(1).one_basis():
                    return "1"
                return self.factor(1)._repr_term(term[1])
            if term[1] == self.factor(1).one_basis():
                return self.factor(0)._repr_term(term[0])
            symb = " "
        else:
            from sage.categories.tensor import tensor
            symb = self._print_options['tensor_symbol']
            if symb is None:
                symb = tensor.symbol
        return symb.join(algebra._repr_term(t) for (algebra, t) in zip(self.factors(), term))

    def twist(self):
        r"""
        The map that defines the smash product structure.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: AA = tensor([A,A],category=ModulesWithBasis(ZZ))
            sage: def t_map((b,a)):
            ...       return AA.monomial((a,b))
            sage: C = SmashProductAlgebra(A,A,twist_on_basis=t_map)
            sage: C.twist()
            Generic endomorphism of A # A

        """
        return self._twist

    def _product_morphism(self):
        r"""
        The multiplication map on the smash product.

        This is a module homomorphism from the twofold tensor product of ``self``, to ``self``.
        """
        return self._product_morphism_map

    def product_on_basis(self, p1, p2):
        r"""
        The product of basis elements indexed by the grouped tuples `p1` and `p2`.

        ..IMPLEMENTATION::

        The tensor map form of multiplication takes the flattened tuples as input and output.
        In the following example, `C` is a smash product and `C2` is the direct product of `C` with itself
        with componentwise product defined in :meth:`sage.categories.algebras_with_basis.AlgebrasWithBasis.TensorProducts`.

        EXAMPLES::

            sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: AA = tensor([A,A], category=ModulesWithBasis(ZZ))
            sage: twist = AA.module_morphism(on_basis=lambda (b,a): AA.monomial((b*a*b.inverse(),b)),codomain=AA)
            sage: C = SmashProductAlgebra(A, A, twist_morphism=twist,suppress_ones=False)
            sage: p1 = (r([1]),r([1,2]))
            sage: p2 = (r([1,2,1]),r([2]))
            sage: C.product_on_basis(p1,p2)
            B[1] # B[s1]
            sage: C2 = tensor([C,C])
            sage: C2.product_on_basis(C2.indices_to_index()(p1,p1), C2.indices_to_index()(p2,p2))
            B[1] # B[s1] # B[1] # B[s1]
        """
        mult = self._product_morphism()
        return self(mult(mult.domain().monomial(p1+p2)))

    def _the_coercion_map(self, other_algebra_module, twist, x):
        return self(twist(other_algebra_module(x)))

    def register_opposite(self, other_algebra):
        if not isinstance(other_algebra, SmashProductAlgebra):
            raise TypeError, "%s is not a smash product"
        if not other_algebra.factor(0) == self.factor(1) or not other_algebra.factor(1) == self.factor(0):
            raise TypeError, "Factors are not opposite"
        twist = self.twist()
        untwist = other_algebra.twist()
        ABmod = twist.codomain()
        BAmod = twist.domain()
        if untwist.domain() != ABmod or untwist.codomain() != BAmod:
            raise TypeError, "Twists are not inverse"
        self._opposite_product = other_algebra
        other_algebra._opposite_product = self
        self._has_opposite = True
        other_algebra._has_opposite = True

        # set up coercions between ASB and BSA

        module_category=ModulesWithBasis(self.base_ring())

        other_to_self = SetMorphism(Hom(other_algebra, self, category=module_category), partial(self._the_coercion_map, BAmod, twist))
        other_to_self.register_as_coercion()

        self_to_other = SetMorphism(Hom(self, other_algebra, category=module_category), partial(other_algebra._the_coercion_map, ABmod, untwist))
        self_to_other.register_as_coercion()

    def opposite(self):
        r"""
        The smash product based on the other ordering of tensor factors.

        For a doctest, see :meth:`SmashProductAlgebraElement.to_opposite`
        """
        if not self._has_opposite:
            raise NotImplementedError, "Opposite smash product is not defined"
        return self._opposite_product

SmashProductAlgebra.Element = SmashProductAlgebraElement

