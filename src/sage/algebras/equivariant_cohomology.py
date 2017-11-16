r"""
Equivariant cohomology for Kac-Moody flag manifolds
"""
from sage.combinat.root_system.cartan_type import CartanType
from sage.sets.family import Family
from sage.algebras.symmetric_algebra import SymmetricAlgebra
from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule
from sage.rings.rational_field import QQ
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup

class EquivariantCohomologyPoint(SymmetricAlgebra):
    r"""
    The torus-equivariant cohomology of a point.

    INPUT:

    - ``cartan_type`` -- A cartan type
    - ``lattice`` -- one of 'weight', 'root', 'ambient' (default: 'weight')
    - ``prefix`` -- A string (default: None)

    This class implements a polynomial ring with distinguished generators
    coming from a distinguished basis of a lattice specified by a Cartan type::

        sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
        sage: EquivariantCohomologyPoint(['A',2]).algebra_generators()
        Finite family {1: o1, 2: o2}
        sage: EquivariantCohomologyPoint(['A',2], lattice='root').algebra_generators()
        Finite family {1: a1, 2: a2}
        sage: S = EquivariantCohomologyPoint(['A',2], lattice='ambient')
        sage: x = S.algebra_generators(); x
        Finite family {0: x0, 1: x1, 2: x2}

    The automorphisms induced by Weyl group elements are built in::

        sage: s = S.weyl_group().simple_reflections()
        sage: f = S.weyl_automorphism(s[1]*s[2])
        sage: [(t, f(t**2)) for t in x]
        [(x0, x1^2), (x1, x2^2), (x2, x0^2)]

    Divided difference operators are implemented::

        sage: S.divided_difference(x[0]**2, 1)
        x0 + x1

    Localization values for opposite Schubert classes are available::

        sage: S.schubert_loc(s[1]*s[2], s[1]*s[2]*s[1]).factor(proof=False)
        (-x0 + x2) * (-x0 + x1)
    """
    @staticmethod
    def __classcall__(cls, cartan_type, lattice='weight', prefix=None):
        cartan_type = CartanType(cartan_type)
        if lattice == 'ambient':
            if not cartan_type.is_finite() and not (cartan_type.is_affine() and cartan_type.is_untwisted_affine()):
                raise TypeError("Ambient space is only used for finite and untwisted affine Cartan types")
        elif lattice not in ['weight', 'root']:
            raise ValueError("lattice specification is not valid")
        if prefix is None:
            prefix = ""
        return super(EquivariantCohomologyPoint, cls).__classcall__(cls,cartan_type,lattice,prefix=prefix)

    def __init__(self, cartan_type, lattice, prefix):
        r"""                                                                    
        """
        self._cartan_type = cartan_type
        if lattice == 'ambient':
            lattice = cartan_type.root_system().ambient_space(base_ring=QQ)
            if prefix == "":
                prefix = 'x'
        elif lattice == 'root':
            lattice = cartan_type.root_system().root_space(base_ring=QQ)
            if prefix == "":
                prefix = 'a'
        else:
            if cartan_type.is_affine() and cartan_type.is_untwisted_affine():
                lattice = cartan_type.root_system().weight_space(base_ring=QQ, extended=True)
                if prefix == "":
                    prefix = 'L'
            else:
                lattice = cartan_type.root_system().weight_space(base_ring=QQ)
                if prefix == "":
                    prefix = 'o'
        self._lattice = lattice
        self._weyl_group = lattice.weyl_group(prefix="s")
        SymmetricAlgebra.__init__(self, self._lattice, prefix=prefix)
        self._forgetful_morphism = self.induced_algebra_endomorphism([lattice.zero() for x in lattice.basis()])
        self._nilCoxeter_action = HeckeAlgebraRepresentation(self, self.divided_difference, cartan_type, 0, 0)
        self._group_algebra = self._lattice.algebra(QQ)

    def _repr_(self):
        r"""                                                                    
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1])
            Equivariant cohomology of a point for the maximal torus of type ['A', 2, 1]
        """
        return "Equivariant cohomology of a point for the maximal torus of type %s"%self.cartan_type()

    def cartan_type(self):
        r"""                                                                    
        The Cartan type of ``self``.

        EXAMPLES::                                                              
                                                                                
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1]).cartan_type()
            ['A', 2, 1]                                                         
        """
        return self._cartan_type

    def lattice(self):
        r"""
        The lattice of which ``self`` is the symmetric algebra.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1]).lattice()
            Extended weight space over the Rational Field of the Root system of type ['A', 2, 1]
            sage: EquivariantCohomologyPoint(['A',2]).lattice()
            Weight space over the Rational Field of the Root system of type ['A', 2]
            sage: EquivariantCohomologyPoint(['C',2], lattice='ambient').lattice()
            Ambient space of the Root system of type ['C', 2]
            sage: EquivariantCohomologyPoint(['C',2], lattice='root').lattice()
            Root space over the Rational Field of the Root system of type ['C', 2]
        """
        return self._lattice

    def group_algebra(self):
        r"""
        The group algebra of the lattice of ``self``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1]).group_algebra()
            Group algebra of the Extended weight space over the Rational Field of the Root system of type ['A', 2, 1] over Rational Field
        """
        return self._group_algebra

    @cached_method
    def simple_roots(self):
        r"""
        The simple roots embedded into the polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1]).simple_roots()
            Finite family {0: 2*L0 - L1 - L2 + Ldelta, 1: -L0 + 2*L1 - L2, 2: -L0 - L1 + 2*L2}
            sage: EquivariantCohomologyPoint(['A',2]).simple_roots()
            Finite family {1: 2*o1 - o2, 2: -o1 + 2*o2}
            sage: EquivariantCohomologyPoint(['C',2], lattice='ambient').simple_roots()
            Finite family {1: x0 - x1, 2: 2*x1}
            sage: EquivariantCohomologyPoint(['C',2], lattice='root').simple_roots()
            Finite family {1: a1, 2: a2}
        """
        return Family(self.cartan_type().index_set(), lambda i: self.from_module_map()(self.base_module().simple_root(i)))

    @cached_method
    def simple_root(self, i):
        r"""
        The `i`-th simple root embedded into the polynomial ring.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1]).simple_root(0)
            2*L0 - L1 - L2 + Ldelta
            sage: EquivariantCohomologyPoint(['A',2]).simple_root(2)
            -o1 + 2*o2
            sage: EquivariantCohomologyPoint(['C',2], lattice='ambient').simple_root(2)
            2*x1
            sage: EquivariantCohomologyPoint(['C',2], lattice='root').simple_root(1)
            a1
        """
        return self.simple_roots()[i]

    def weyl_group(self):
        r"""
        The Weyl group of ``self``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: EquivariantCohomologyPoint(['A',2,1]).weyl_group()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the extended weight space)
        """
        return self._weyl_group

    @cached_method
    def weyl_automorphism(self, w):
        r"""
        The automorphism of ``self`` induced by the Weyl group element ``w``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',1],lattice='root')
            sage: s = H.weyl_group().gen(0)
            sage: S = H.weyl_automorphism(s)
            sage: S(H.simple_root(1))
            -a1

            sage: H = EquivariantCohomologyPoint(['A',2],lattice='ambient')
            sage: S = Family(H.cartan_type().index_set(), lambda i: H.weyl_automorphism(H.weyl_group().simple_reflection(i)))
            sage: [[v, i, S[i](v)] for v in H.base_module().basis() for i in H.cartan_type().index_set()]
            [[(1, 0, 0), 1, x1], [(1, 0, 0), 2, x0], [(0, 1, 0), 1, x0], [(0, 1, 0), 2, x2], [(0, 0, 1), 1, x2], [(0, 0, 1), 2, x1]]
        """
        return self.induced_algebra_endomorphism([w.action(x) for x in self.lattice().basis()])

    @cached_method
    def simple_reflection_automorphism(self, i):
        r"""
        The automorphism of the `i`-th simple reflection acting on ``self``.

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',2],lattice='root')
            sage: s = H.weyl_group().simple_reflections()
            sage: S = H.simple_reflection_automorphism(2)
            sage: S(H.simple_root(1))
            a1 + a2
            sage: S(H.simple_root(2))
            -a2
        """
        return self.weyl_automorphism(self.weyl_group().simple_reflection(i))

    def divided_difference(self, f, i):
        r"""
        The `i`-th divided difference of `f`.

        MATH::

            (f - s_i(f)) / alpha_i

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',2],lattice='root')
            sage: a = H.simple_roots()
            sage: H.divided_difference(a[1]*(a[1]+a[2]), 1)
            a1 + 2*a2
            sage: H.divided_difference(a[1]*(a[1]+a[2]), 2)
            0
        """
        return self.divide_elements(f - self.weyl_automorphism(self.weyl_group().simple_reflection(i))(f), self.simple_root(i))

    def forgetful_morphism(self):
        r"""
        The morphism to the base ring that forgets equivariance.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: S = EquivariantCohomologyPoint(['A',2,1])
            sage: f = S.forgetful_morphism()
            sage: p = S.an_element(); p
            L0
            sage: f((p+3)**2)
            9
        """
        return self._forgetful_morphism

    @cached_method
    def schubert_loc(self, v, w):
        r"""
        Localize the `v`-th Schubert class at the torus-fixed point `w`.

        The identity indexes the big cell.

        INPUT::

            - v, w -- Weyl group elements

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',2,1])
            sage: W = H.weyl_group()
            sage: r = W.from_reduced_word
            sage: H.schubert_loc(r([]), r([]))
            1
            sage: H.schubert_loc(r([]), r([1]))
            1
            sage: H.schubert_loc(r([1]), r([]))
            0
            sage: H.schubert_loc(r([0]), r([0,1]))
            2*L0 - L1 - L2 + Ldelta
            sage: H.schubert_loc(r([0]), r([0,1,0]))
            L0 + L1 - 2*L2 + Ldelta
            sage: H.schubert_loc(r([1]), r([0,1,0]))
            L0 + L1 - 2*L2 + Ldelta
            sage: H.schubert_loc(r([1,0]), r([0,1,0])).factor(proof=False)
            (-L0 + 2*L1 - L2) * (L0 + L1 - 2*L2 + Ldelta)
            sage: H = EquivariantCohomologyPoint(['C',2],lattice='ambient')
            sage: W = H.weyl_group()
            sage: r = W.from_reduced_word
            sage: H.algebra_generators()
            Finite family {0: x0, 1: x1}
            sage: H.schubert_loc(r([2]),r([1,2]))
            2*x0
            sage: H.schubert_loc(r([1,2]),r([1,2])).factor(proof=False)
            (-2) * (-x0 + x1) * x0
            sage: H.schubert_loc(r([1,2]),r([1,2,1,2])).factor(proof=False)
            (2) * x0 * (x0 + x1)

            sage: H = EquivariantCohomologyPoint(['C',2,1], lattice='root')
            sage: W = H.weyl_group()
            sage: r = W.from_reduced_word
            sage: H.schubert_loc(r([0,2]),r([0,2,1,2])).factor(proof=False)
            (2) * a0 * (a0 + a1 + a2)

        """
        if w == self.weyl_group().one():
            if v == w:
                return self.base_ring().one()
            else:
                return self.base_ring().zero()
        i = w.first_descent(side='right')
        wsi = w.apply_simple_reflection(i, side='right')
        ans = self.schubert_loc(v, wsi)
        if v.has_descent(i, side='right'):
            ans = ans + (self.from_module_map()(self.base_module().simple_root(i).weyl_action(wsi))) * self.schubert_loc(v.apply_simple_reflection(i,side='right'),wsi)
        return ans

    def exponent_to_weight(self, exp):
        r"""
        Change an exponent vector to a lattice element.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',2], lattice='ambient')
            sage: x = H.algebra_generators()
            sage: f = x[0]**2+3*x[1]*x[2]
            sage: f.dict()
            {(0, 1, 1): 3, (2, 0, 0): 1}
            sage: _.keys()
            [(0, 1, 1), (2, 0, 0)]
            sage: [H.exponent_to_weight(vector(ex)) for ex in _]
            [(0, 1, 1), (2, 0, 0)]
        """
        return self.lattice().sum_of_terms([(self.keys()[i], exp[i]) for i in range(len(exp))], distinct=True)

    def to_group_algebra(self, f):
        r"""
        Map an element to the group algebra.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',2,1])
            sage: x = H.an_element(); x
            L0
            sage: H.to_group_algebra(x*x+3*x)
            3*B[Lambda[0]] + B[2*Lambda[0]]
        """
        d = f.dict()
        return self.group_algebra().sum_of_terms([(self.exponent_to_weight(exp),d[exp]) for exp in d.keys()
], distinct=True)

    def from_group_algebra_on_basis(self, wt):
        r"""
        The monomial of ``self`` associated with a weight.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint("A2",lattice='ambient')
            sage: X = H.lattice()
            sage: H.from_group_algebra_on_basis(X((2,2,3)))
            x0^2*x1^2*x2^3
        """
        return self.prod([self.algebra_generators()[i] ** wt[i] for i in self.keys()])

    @cached_method
    def from_group_algebra_map(self):
        return self.group_algebra().module_morphism(on_basis=self.from_group_algebra_on_basis, codomain=self)

    def from_group_algebra(self, f):
        r"""
        Map an element from the group algebra to ``self``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyPoint
            sage: H = EquivariantCohomologyPoint(['A',2,1])
            sage: h = H.group_algebra().an_element(); h
            B[2*Lambda[0] + 2*Lambda[1] + 3*Lambda[2]]
            sage: H.from_group_algebra((h+1)**2)
            L0^4*L1^4*L2^6 + 2*L0^2*L1^2*L2^3 + 1
        """
        return self.from_group_algebra_map()(f)

class EquivariantCohomologyFlags(CombinatorialFreeModule):
    r"""
    The maximal-torus equivariant cohomology ring of the generalized flag manifold.

    INPUT:

    - ``cartan_type`` -- A cartan type

    EXAMPLES::

        sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
        sage: H = EquivariantCohomologyFlags(['A',2])
        sage: S = H.base_ring()
        sage: W = S.weyl_group()
        sage: s = W.simple_reflections()

    This is a free `S`-algebra where `S` is the polynomial ring given by the equivariant cohomology
    of a point. The distinguished basis is given by the equivariant classes of opposite
    Schubert varieties in the (thick) Kac-Moody flag manifolds

    The unit is the class indexed by the identity::

        sage: H.one()
        X[1]

    Here the square of the class of a simple reflection is computed.

        sage: H.monomial(s[1])*H.monomial(s[1])
        X[s2*s1] + (x0-x1)*X[s1]
    """
    @staticmethod
    def __classcall__(cls, cartan_type):
        cartan_type = CartanType(cartan_type)
        return super(EquivariantCohomologyFlags, cls).__classcall__(cls,cartan_type)

    def __init__(self, cartan_type):
        r"""
        """
        if cartan_type.is_finite() or cartan_type.is_affine():
            lattice = 'ambient'
        else:
            lattice = 'weight'
        S = EquivariantCohomologyPoint(cartan_type, lattice)
        CombinatorialFreeModule.__init__(self, S, S.weyl_group(), prefix="X", category=AlgebrasWithBasis(S))

    @cached_method
    def cartan_type(self):
        r"""                                                                    
        The Cartan type of ``self``.

        EXAMPLES::                                                              
                                                                                
            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
            sage: EquivariantCohomologyFlags(['A',2,1]).cartan_type()
            ['A', 2, 1]                                                         
        """
        return self.base_ring().cartan_type()

    def _repr_(self):
        r"""
        A string representing ``self``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
            sage: EquivariantCohomologyFlags(['A',2,1])
            Equivariant cohomology of the flag variety of type ['A', 2, 1]
        """
        return "Equivariant cohomology of the flag variety of type %s"%self.cartan_type()

    def localization_function(self, elt):
        r"""
        The localization function associated with a 
        linear combination of Weyl group elements.

        Each Weyl group element represents a Schubert class.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
            sage: H = EquivariantCohomologyFlags(['A',2])
            sage: S = H.base_ring()
            sage: W = S.weyl_group()
            sage: s = W.simple_reflections()
            sage: elt = S.simple_root(1)*H.monomial(s[1]) + H.monomial(s[2]*s[1])
            sage: f = H.localization_function(elt)
            sage: [(w, f(w)) for w in W]
            [(1, 0), (s1*s2*s1, x0^2 - 2*x0*x2 + x2^2), (s1*s2, x0^2 - 2*x0*x1 + x1^2), (s1, x0^2 - 2*x0*x1 + x1^2), (s2*s1, x0^2 - 2*x0*x2 + x2^2), (s2, 0)]
        """
        return lambda w: self.base_ring().sum([c*self.base_ring().schubert_loc(v,w) for (v,c) in elt])

    def schubert_expansion(self, f, support_lower_bounds=None, length_upper_bound=None, index_set=None):
        r"""
        Expand an equivariant class into Schubert classes.

        INPUT:

        - `f` -- A function from the Weyl group to ``self``, interpreted
        as an equivariant class in the flag manifold under
        localization at `T`-fixed points.

        - ``support_lower_bounds`` -- (default: None) a set of Weyl group elements `x` such that for every element `y` in the support, `x \le y`.

        - ``length_upper_bound`` -- (default: None) An upper bound for the
        length of elements in the support.

        - ``index_set`` -- (default: None, meaning []) Asserts that the support is contained in the set
        of minimum length coset representatives modulo the subgroup generated by the simple reflections
        in ``index_set``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
            sage: H = EquivariantCohomologyFlags(['A',2])
            sage: S = H.base_ring()
            sage: W = S.weyl_group()
            sage: s = W.simple_reflections()
            sage: elt = S.simple_root(1)*H.monomial(s[1]) + H.monomial(s[2]*s[1])
            sage: f = H.localization_function(elt)
            sage: H.schubert_expansion(f)
            X[s2*s1] + (x0-x1)*X[s1]
        """
        S = self.base_ring()
        W = S.weyl_group()
        if not support_lower_bounds:
            support_lower_bounds = [W.one()]
        min_length = min([x.length() for x in support_lower_bounds])
        if W.cartan_type().is_finite():
            if not length_upper_bound:
                length_upper_bound = W.long_element().length()
            else:
                length_upper_bound = min(length_upper_bound, W.long_element().length())
        elif not length_upper_bound:
            raise ValueError("A length upper bound must be specified for infinite type")
        answer = self.zero()
        for i in range(min_length, length_upper_bound+1):
            answer_loc_func = self.localization_function(answer)
            wlist = [w for w in W.elements_of_length(i) if all([x.bruhat_le(w) for x in support_lower_bounds])]
            if index_set:
                wlist = [w for w in wlist if w == w.coset_representative(index_set)]
            answer = answer + self.sum([S.divide_elements(f(w)-answer_loc_func(w),S.schubert_loc(w,w)) * self.monomial(w) for w in wlist])
        return answer

    @cached_method
    def one_basis(self):
        r"""
        The index element for the identity of ``self``.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
            sage: H = EquivariantCohomologyFlags(['A',2]).one_basis()
        """
        return self.base_ring().weyl_group().one()

    @cached_method
    def product_on_basis(self, u, v):
        r"""
        Product of two Schubert basis elements.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyFlags
            sage: H = EquivariantCohomologyFlags(['A',2])
            sage: S = H.base_ring()
            sage: W = S.weyl_group()
            sage: s = W.simple_reflections()
            sage: H.product_on_basis(s[1], s[1])
            X[s2*s1] + (x0-x1)*X[s1]
        """
        S = self.base_ring()
        f = lambda w: S.schubert_loc(u,w)*S.schubert_loc(v,w)
        du = u.descents()
        dv = v.descents()
        index_set = self.base_ring().cartan_type().index_set()
        index_set = [i for i in index_set if i not in du and i not in dv]
        return self.schubert_expansion(f, support_lower_bounds=[u,v], length_upper_bound = u.length()+v.length(), index_set=index_set)

class EquivariantCohomologyAffineFlags(EquivariantCohomologyFlags):
    r"""
    The torus equivariant cohomology algebra for the affine flag manifolds.

    This class has support for localizations at translation elements.
    """
    @staticmethod
    def __classcall__(cls, cartan_type):
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_untwisted_affine():
            raise ValueError("not an untwisted affine Cartan type")
        return super(EquivariantCohomologyAffineFlags, cls).__classcall__(cls,cartan_type)

    def __init__(self, cartan_type):
        r"""
        """
        EquivariantCohomologyFlags.__init__(self, cartan_type)
        self._WE = ExtendedAffineWeylGroup(cartan_type)
        self._PvW0 = self._WE.PvW0()
        self._X = self.base_ring().lattice()
        self._Xbasis = self._X.basis()
        self._X0 = self._X.classical()  # classical ambient space
        self._X0basis = self._X0.basis()
        self._Ia = self._X0basis.keys() # keys to basis of classical ambient space

    def extended_affine_weyl_group(self):
        return self._PvW0

    def loc_at_func(self, w):
        r"""
        Return the function that localizes finite Weyl-invariant polynomials
        at a given fixed point in the affine Weyl group.


        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyAffineFlags
            sage: H = EquivariantCohomologyAffineFlags(['A',2,1])
            sage: S = H.base_ring()
            sage: W = S.weyl_group()
            sage: w = W.from_reduced_word([1,2,1,0])
            sage: lf = H.loc_at_func(w)
            sage: lf(S.sum([S.simple_root(i)**2 for i in range(3)]))
            2*x0^2 - 2*x0*x1 - 2*x0*x2 + 4*x0*xdelta + 2*x1^2 - 2*x1*x2 + 2*x2^2 - 4*x2*xdelta + 3*xdelta^2
            sage: lf(S.gen(0)*(S.gen(1)+S.gen(2))+S.gen(1)*S.gen(2))
            x0*x1 + x0*x2 + x0*xdelta + x1*x2 + 2*x1*xdelta + 3*x2*xdelta + 2*xdelta^2
            sage: lf(S.gen(0))
            x0 + 2*xdelta
            sage: lf(S.gen(1))
            x1 + xdelta
            sage: lf(S.gen(2))
            x2
        """
        we = self.extended_affine_weyl_group().from_reduced_word(w.reduced_word())
        trans = we.to_dual_translation_left().to_ambient()
        return self.loc_at_func_trans(trans)

    @cached_method
    def loc_at_func_trans(self, trans):
        r"""
        Above function for a given translation element.

        EXAMPLES::

            sage: from sage.algebras.equivariant_cohomology import EquivariantCohomologyAffineFlags
            sage: H = EquivariantCohomologyAffineFlags(['A',2,1])
            sage: S = H.base_ring()
            sage: S.gens()
            (x0, x1, x2, xdelta, xdeltacheck)
            sage: X0 = S.lattice().classical()
            sage: lf = H.loc_at_func_trans(X0((-1,0,1)))
            sage: lf(S.sum([S.simple_root(i)**2 for i in range(3)]))
            2*x0^2 - 2*x0*x1 - 2*x0*x2 + 4*x0*xdelta + 2*x1^2 - 2*x1*x2 + 2*x2^2 - 4*x2*xdelta + 3*xdelta^2
            sage: lf(S.gen(0)*(S.gen(1)+S.gen(2))+S.gen(1)*S.gen(2))
            x0*x1 + x0*x2 - x0*xdelta + x1*x2 + x2*xdelta - xdelta^2
            sage: lf(S.gen(0))
            x0 + xdelta
            sage: lf(S.gen(1))
            x1
            sage: lf(S.gen(2))
            x2 - xdelta
        """
        S = self.base_ring()
        return S.hom([S(self._Xbasis[i] - trans.scalar(self._X0basis[i]) * self._Xbasis['delta']) for i in self._Ia] + [S(self._Xbasis['delta']), S(self._Xbasis['deltacheck'])], codomain=S)

