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

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.all import AlgebrasWithBasis, ModulesWithBasis
from sage.categories.category import Category
from sage.categories.morphism import SetMorphism
from sage.categories.tensor import tensor
from sage.categories.homset import Hom
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModule_Tensor, \
CartesianProductWithFlattening, CartesianProductWithUnflattening

"""
It allows two ways to construct smash products: a single smash product
A#B or an isomorphic pair of them A#B and B#A.

The single smash product construction is straightforward. 

The isomorphic pair is trickier.
In addition to having each of A#B and B#A we want them
to know about each other via :meth:`opposite`()
and even to construct coercions between them.

Both A#B and B#A are constructed in the same kind of way but
with different input data.

The single smash product case (when ``untwist`` is None) is straightforward.

For the pair, first A#B is constructed.
To make B#A the class :class:`SmashProductAlgebra` is invoked recursively,
but with the single smash product option (and different parameters).
It returns with the instance B#A.
Now that both A#B and B#A have been constructed as full-fledged :class:`SmashProductAlgebra` instances,
the two instances are told about each other and their mutual coercions are
established. Then the job is done.

"""

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
        
            sage: F = FreeMonoid(n=1, names=['a'])
            sage: A = F.algebra(ZZ); A.rename("A")
            sage: A._print_options['prefix'] = 'A'
            sage: W = WeylGroup("A1",prefix="s")
            sage: B = W.algebra(ZZ); B.rename("B")
            sage: cat = ModulesWithBasis(ZZ)
            sage: AB = tensor([A,B],category=cat)
            sage: BA = tensor([B,A],category=cat)
            sage: def twist_func(x):
            ...       if x[0] != x[0].parent().one() and mod(x[1].to_word().length(),2) == 1:
            ...           return -AB.monomial((x[1],x[0]))
            ...       return AB.monomial((x[1],x[0]))
            sage: twist = BA.module_morphism(on_basis=twist_func, codomain=AB, category=cat)
            sage: def untwist_func(x):
            ...       if x[1] != x[1].parent().one() and mod(x[0].to_word().length(),2) == 1:
            ...           return -BA.monomial((x[1],x[0]))
            ...       return BA.monomial((x[1],x[0]))
            sage: untwist = AB.module_morphism(on_basis=untwist_func, codomain=BA, category=cat)

            sage: ASB = SmashProductAlgebra(A, B, twist, untwist) # finally!
            sage: ab = ASB.an_element(); ab
            2*A[a] # B[s1] + 2*A[a] # B[1]
            sage: BSA = ASB.opposite()
            sage: ba = BSA(ab); ba
            -2*B[s1] # A[a] + 2*B[1] # A[a]            
            sage: ASB(ba) == ab
            True
            sage: ASB == BSA.opposite()
            True
            sage: ba == ab.to_opposite()
            True
            sage: A._print_options['prefix'] = 'B'

        """
        return self.parent().opposite()(self)

class SmashProductAlgebra(CombinatorialFreeModule_Tensor):
    r"""
    Smash products of algebras.

    INPUT:

        - ``A`` and ``B`` -- algebras with basis over the same commutative ring
        - ``twist`` -- a module homomorphism from `B` tensor `A` to `A` tensor `B`
        - ``untwist`` -- optional (default: None); the inverse of ``twist``

    Returns the smash product `S(A,B,twist)`. As a module it is `A` tensor `B`.
    The product is determined by a module homomorphism ``twist`` from `B` tensor `A`
    to `A` tensor `B` as described below.
    If ``untwist`` is not None, it is assumed to be the inverse map of ``twist``.
    In this case ``untwist`` is an algebra isomorphism from `S(A,B,twist)` to `S(B,A,untwist)`.
    Both `S(A,B,twist)` and `S(B,A,untwist)` are created, and coercions are registered between them,
    and each is given a way to access the other via :meth:`opposite`.

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

    REFERENCE:

        .. [CIMZ] Caenepeel, S.; Ion, Bogdan; Militaru, G.; Zhu, Shenglin. The factorization problem
          and the smash biproduct of algebras and coalgebras. Algebr. Represent. Theory 3 (2000),
          no. 1, 19--42. 

    If `B` has a left action on `A` by algebra endomorphisms, defining `twist(b tensor a) = (b . a) tensor b`
    makes `m` into an algebra multiplication map. Similarly if `A` has a right action on `B` by algebra
    endomorphisms, then one may use `twist(b tensor a) = a tensor (b . a)`.

    In the following example both `A` and `B` are the group algebra of the symmetric group `S_3` over the integers,
    and `B` has a left action on `A` induced by conjugation of symmetric group elements.

    EXAMPLES::

        sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
        sage: r = W.from_reduced_word
        sage: A = W.algebra(ZZ); A.rename("A")
        sage: AA = tensor([A,A], category=ModulesWithBasis(ZZ))
        sage: twist = AA.module_morphism(on_basis=lambda x: AA.monomial((x[0]*x[1]*x[0]**(-1),x[0])),codomain=AA)
        sage: ba = AA.monomial((r([1]),r([2]))); ba
        B[s1] # B[s2]
        sage: ab = twist(ba); ab
        B[s1*s2*s1] # B[s1]
        sage: C = SmashProductAlgebra(A, A, twist); C
        Smash product of A and A
        sage: c = A.monomial(A.basis().keys().an_element()); c
        B[s1*s2]
        sage: d = A.an_element(); d
        B[s1*s2*s1] + 3*B[s1*s2] + 3*B[s2*s1]
        sage: cd = C.from_direct_product(c,d); cd
        B[s1*s2] # B[s1*s2*s1] + 3*B[s1*s2] # B[s1*s2] + 3*B[s1*s2] # B[s2*s1]
        sage: cd * cd
        9*B[s2*s1] # B[s1*s2] + 9*B[s2*s1] # B[s2*s1] + 3*B[s2*s1] # B[s1] + 3*B[s2*s1] # B[s2] + 18*B[s2*s1] # B[1] + 3*B[1] # B[s1] + 3*B[1] # B[s2] + B[1] # B[1]

    Particular rings that are realized as smash products are the cohomological and K-theoretic nilHecke rings
    of Kostant and Kumar, the affine Hecke algebra, and the double affine Hecke algebra. These all come from actions
    of one algebra on another. Another class of examples comes from quasi-triangular Hopf algebras and the R-matrix.

    Some implementation details: The smash product is the flattened tensor product of `A` and `B`.
    Coercion is allowed between this flattened tensor product and the smash product.
    """

    @staticmethod
    def __classcall__(cls, A, B, twist, untwist=None):
        # type checking
        R = A.base_ring()
        module_category = ModulesWithBasis(R)
        algebra_category = AlgebrasWithBasis(R)
        tensor_category = module_category.TensorProducts()
        if R != B.base_ring():
            raise TypeError, "%s and %s must have the same base ring"%(A,B)
        if A not in algebra_category or B not in algebra_category:
            raise TypeError, "Tensor factors should be AlgebrasWithBasis over %s"%(R)
        if not twist.category().is_subcategory(module_category.hom_category()):
            raise TypeError, "twist should be a module morphism"
        AB = tensor([A,B], category=module_category)
        BA = tensor([B,A], category=module_category)
        if twist.codomain() != AB or twist.domain() != BA:
            raise TypeError, "Domain or codomain of twist is incorrect"

        if untwist is not None:
            if not untwist.category().is_subcategory(module_category.hom_category()):
                raise TypeError, "untwist should be a module morphism"
            if untwist.domain() != AB or untwist.codomain() != BA:
                raise TypeError, "Domain or codomain of untwist is incorrect"

        ASB = super(SmashProductAlgebra, cls).__classcall__(cls, A, B, twist)

        # the single smash product case
        if untwist is None:
            return ASB

        BSA = SmashProductAlgebra(B, A, untwist)
        
        ASB._opposite_product = BSA
        BSA._opposite_product = ASB
        ASB._has_opposite = True
        BSA._has_opposite = True

        # set up coercions
        ABmod = twist.codomain()
        BAmod = twist.domain()
#        def coerce_to_ABmod(x):
#            return ABmod(x)
#        def coerce_to_BAmod(x):
#            return BAmod(x)
#        def coerce_to_ASB(x):
#            return ASB(x)
#        def coerce_to_BSA(x):
#            return BSA(x)
#
#        ASB_to_ABmod = SetMorphism(Hom(ASB, ABmod, category=module_category), coerce_to_ABmod)
#        ABmod_to_ASB = SetMorphism(Hom(ABmod, ASB, category=module_category), coerce_to_ASB)
#        BSA_to_BAmod = SetMorphism(Hom(BSA, BAmod, category=module_category), coerce_to_BAmod)
#        BAmod_to_BSA = SetMorphism(Hom(BAmod, BSA, category=module_category), coerce_to_BSA)
#
#        ASB_to_BSA = SetMorphism(Hom(ASB, BSA, category=module_category), BAmod_to_BSA * untwist * ASB_to_ABmod)
#        BSA_to_ASB = SetMorphism(Hom(BSA, ASB, category=module_category), ABmod_to_ASB * twist * BSA_to_BAmod)

        def exalted_untwist(x):
            return BSA(untwist(ABmod(x)))

        ASB_to_BSA = SetMorphism(Hom(ASB, BSA, category=module_category), exalted_untwist)
        ASB_to_BSA.register_as_coercion()

        def exalted_twist(x):
            return ASB(twist(BAmod(x)))

        BSA_to_ASB = SetMorphism(Hom(BSA, ASB, category=module_category), exalted_twist)
        BSA_to_ASB.register_as_coercion()

        return ASB

    def __init__(self, A, B, twist):
        """
        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: B = A
            sage: BA = tensor([B,A],category=ModulesWithBasis(ZZ))
            sage: AB = BA
            sage: def t_map(x):
            ...       return AB.monomial((x[1],x[0]))
            sage: twist = BA.module_morphism(on_basis=t_map, codomain=AB, category=ModulesWithBasis(ZZ))
            sage: C = SmashProductAlgebra(A,B,twist); C # indirect doctest
            Smash product of A and A

        """
        # all type checking is done through :class:`SmashProductAlgebra`.
        self._factors = (A,B)
        # Note the difference between self._factors and the eventual self._sets coming from
        # the superclass CombinatorialFreeModule_Tensor: the latter
        # give the factors of the flattened tensor product if A or B happen to also be tensor products.
        self._twist = twist
        self._has_opposite = False # default: no "opposite" smash product is available
        self._id_keys = (A.one_basis(), B.one_basis())
        R = A.base_ring()
        module_category=ModulesWithBasis(R)
        tensor_category=module_category.TensorProducts()
        self._factor_is_tensor = tuple([x in tensor_category for x in self._factors])
        self._n_factors = tuple([len(x) if isinstance(x, tuple) else 1 for x in self._id_keys])
        if not all(x != 0 for x in self._n_factors):
            raise TypeError, "Tensor factors should be nontrivial"
        self._key_flattener = CartesianProductWithFlattening(self._factor_is_tensor)
        self._key_unflattener = CartesianProductWithUnflattening(self._factor_is_tensor, self._n_factors)
        
        # Define the product morphism using twist
        ItwistI = tensor([A._identity_map(), twist, B._identity_map()], category=module_category)
        mAmB = tensor([A._product_morphism(),B._product_morphism()], category=module_category)
        self._product_morphism_map = SetMorphism(Hom(ItwistI.domain(), mAmB.codomain(), category=module_category), mAmB * ItwistI)

        # the following imparts the structure of having a product and bilinear tensor product
        # TODO: for Hopf algebras, define the additional structure maps
        self._category = Category.join((AlgebrasWithBasis(R), tensor_category))

        CombinatorialFreeModule_Tensor.__init__(self, twist.codomain()._sets, category=self._category)

        # wanted to register self.factor_embedding(i) as coercions but this produces
        # incorrect results if the tensor factors are the same

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ)
            sage: B = A
            sage: BA = tensor([B,A],category=ModulesWithBasis(ZZ))
            sage: AB = BA
            sage: def t_map(x):
            ...       return AB.monomial((x[1],x[0]))
            sage: twist = BA.module_morphism(on_basis=t_map, codomain=AB, category=ModulesWithBasis(ZZ))
            sage: C = SmashProductAlgebra(A,B,twist); C
            Smash product of A and A

            Question: Should this name include twist?

        """
        return "Smash product of %s and %s"%(self.factors()[0], self.factors()[1])

    @cached_method
    def one_basis(self):
        r"""
        The tuple that indexes the unit element of the smash product.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ)
            sage: cat = ModulesWithBasis(ZZ)
            sage: AA = tensor([A,A], category=cat)
            sage: twist = AA.module_morphism(on_basis=lambda x: AA.monomial((x[1],x[0])), category=cat, codomain=AA)
            sage: A2 = SmashProductAlgebra(A,A,twist)
            sage: A2.one_basis()
            (1, 1)
            sage: AAA = tensor([A,A,A], category=cat)
            sage: twist2= AAA.module_morphism(on_basis=lambda x: AAA.monomial((x[2],x[0],x[1])), category=cat, codomain=AAA)
            sage: AA2 = SmashProductAlgebra(A,A2,twist2)
            sage: AA2.one_basis()
            (1, 1, 1)

        """
        return self._key_flattener(*self._id_keys)

    def factors(self):
        r"""
        The list of factors in the smash product.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: AA = tensor([A,A],category=ModulesWithBasis(ZZ))
            sage: def t_map(x):
            ...       return AA.monomial((x[1],x[0]))
            sage: twist = AA.module_morphism(on_basis=t_map, codomain=AA, category=ModulesWithBasis(ZZ))
            sage: C = SmashProductAlgebra(A,A,twist)
            sage: C.factors()
            (A, A)

        """
        return self._factors

    def twist(self):
        r"""
        The map that defines the smash product structure.

        EXAMPLES::

            sage: W = WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: AA = tensor([A,A],category=ModulesWithBasis(ZZ))
            sage: def t_map(x):
            ...       return AA.monomial((x[1],x[0]))
            sage: twist = AA.module_morphism(on_basis=t_map, codomain=AA, category=ModulesWithBasis(ZZ))
            sage: C = SmashProductAlgebra(A,A,twist)
            sage: C.twist()
            Generic endomorphism of A # A

        """
        return self._twist

    def from_direct_product(self, a, b):
        r"""
        Map the pair of elements into ``self``.

        EXAMPLES::

            sage: W=WeylGroup("A2",prefix="s")
            sage: A = W.algebra(ZZ)
            sage: AA = tensor([A,A],category=ModulesWithBasis(ZZ))
            sage: def t_map(x):
            ...       return AA.monomial((x[1],x[0]))
            sage: twist = AA.module_morphism(on_basis=t_map,category=ModulesWithBasis(ZZ),codomain=AA)
            sage: A2 = SmashProductAlgebra(A,A,twist)
            sage: a = A.an_element(); a
            B[s1*s2*s1] + 3*B[s1*s2] + 3*B[s2*s1]
            sage: b = A.monomial(W.from_reduced_word([1]))
            sage: ab = A2.from_direct_product(a,b); ab
            B[s1*s2*s1] # B[s1] + 3*B[s1*s2] # B[s1] + 3*B[s2*s1] # B[s1]
            sage: ab.parent()
            Smash product of A and A

        """
        return self(tensor([self.factors()[0](a),self.factors()[1](b)],category=ModulesWithBasis(self.base_ring())))

    def an_element(self):
        r"""
        Make an element of ``self``.

        EXAMPLES::

            sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: AA = tensor([A,A], category=ModulesWithBasis(ZZ))
            sage: twist = AA.module_morphism(on_basis=lambda x: AA.monomial((x[0]*x[1]*x[0]**(-1),x[0])),codomain=AA)
            sage: C = SmashProductAlgebra(A, A, twist)
            sage: C.an_element()
            B[s1*s2*s1] # B[s1*s2*s1] + 3*B[s1*s2*s1] # B[s1*s2] + 3*B[s1*s2*s1] # B[s2*s1] + 3*B[s1*s2] # B[s1*s2*s1] + 9*B[s1*s2] # B[s1*s2] + 9*B[s1*s2] # B[s2*s1] + 3*B[s2*s1] # B[s1*s2*s1] + 9*B[s2*s1] # B[s1*s2] + 9*B[s2*s1] # B[s2*s1]

        """
        return self.from_direct_product(self.factors()[0].an_element(),self.factors()[1].an_element())

    @cached_method
    def factor_embedding(self, i):
        r"""
        The algebra embedding from factor ``i`` to ``self``.

        INPUT:

        - ``i`` -- If 0, the left factor, and if 1, the right factor.

        EXAMPLES::

        sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
        sage: r = W.from_reduced_word
        sage: A = W.algebra(ZZ); A.rename("A")
        sage: AA = tensor([A,A], category=ModulesWithBasis(ZZ))
        sage: twist = AA.module_morphism(on_basis=lambda x: AA.monomial((x[0]*x[1]*x[0]**(-1),x[0])),codomain=AA)
        sage: C = SmashProductAlgebra(A, A, twist)
        sage: c = C.factor_embedding(0)(A.monomial(r([1,2]))); c
        B[s1*s2] # B[1]
        sage: d = C.factor_embedding(1)(A.monomial(r([2]))); d
        B[1] # B[s2]
        sage: c*d
        B[s1*s2] # B[s2]
        sage: d*c
        B[s2*s1] # B[s2]

        """
        if i == 0:
            def on_basis_left(x):
                return self.monomial(self._key_flattener(x,self._id_keys[1]))
            on_basis = on_basis_left
        elif i == 1:
            def on_basis_right(x):
                return self.monomial(self._key_flattener(self._id_keys[0],x))
            on_basis = on_basis_right
        else:
            raise ValueError, "Embedding Factor %s must be 0 or 1"%(i)
        return self.factors()[i].module_morphism(on_basis = on_basis, codomain = self, category=ModulesWithBasis(self.base_ring()))

    def _product_morphism(self):
        r"""
        The multiplication map on the smash product.

        This is a module homomorphism from the twofold tensor product of ``self``, to ``self``.
        """
        return self._product_morphism_map

    def product_on_basis(self, p1, p2):
        r"""
        The product of basis elements indexed by `p1` and `p2`.

        EXAMPLES::

            sage: W = WeylGroup(CartanType(['A',2]),prefix="s")
            sage: r = W.from_reduced_word
            sage: A = W.algebra(ZZ); A.rename("A")
            sage: AA = tensor([A,A], category=ModulesWithBasis(ZZ))
            sage: twist = AA.module_morphism(on_basis=lambda x: AA.monomial((x[0]*x[1]*x[0]**(-1),x[0])),codomain=AA)
            sage: C = SmashProductAlgebra(A, A, twist)
            sage: p1 = (r([1]),r([1,2]))
            sage: p2 = (r([1,2,1]),r([2]))
            sage: C.product_on_basis(p1,p2)
            B[1] # B[s1]

        """
        mult = self._product_morphism()
        return self(mult(mult.domain().monomial(p1+p2)))

    def opposite(self):
        r"""
        The smash product based on the other ordering of tensor factors.

        For a doctest, see :meth:`SmashProductAlgebraElement.to_opposite`
        """
        if not self._has_opposite:
            raise NotImplementedError, "Opposite smash product is not defined"
        return self._opposite_product

SmashProductAlgebra.Element = SmashProductAlgebraElement

