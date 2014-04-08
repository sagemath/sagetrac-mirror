"""
Affine Hecke algebra with several realizations

AUTHOR:

- Mark Shimozono (2014)

"""

#*****************************************************************************
#       Copyright (C) 2014 Mark Shimozono <mshimo at math dot vt dot edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.realizations import Category_realization_of_parent, Realizations
from sage.categories.morphism import SetMorphism
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.category import Category
from sage.categories.homset import Hom
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.tensor import tensor
from sage.sets.family import Family, FiniteFamily
from sage.misc.bindable_class import BindableClass
from sage.misc.abstract_method import abstract_method
from sage.algebras.multiparameter_hecke_algebra import MultiParameterHeckeAlgebra, ParameterFamilies
from sage.misc.cachefunc import cached_method
from sage.algebras.smash_product_algebra import SmashProductAlgebra
from sage.rings.integer_ring import ZZ


class AffineHeckeAlgebra(UniqueRepresentation, Parent):
    r"""
    The affine Hecke algebra with several realizations.

    ..warning: Only implemented properly for untwisted or dual untwisted affine type.

    INPUT:

        - ``cartan_type`` -- An affine or finite Cartan type (a finite Cartan type is an abbreviation for its untwisted affinization)
        - ``q1, q2`` -- parameters (default: None)

    The parameters `q1` and `q2` represent the pair of eigenvalues `q1[i]` and `q2[i]` of the algebra generator `T_i`
    for `i` in the affine Dynkin node set `I^X`. See :func:`sage.algebras.multiparameter_hecke_algebra.ParameterFamilies`.

    The supported bases are:

    - "TX" -- `T_w` for `w \in W_a(\tilde{X})`
    - "TY_Y" -- `T_w Y^\mu` for `w \in W(Y)` and `\mu \in Q^Y`
    - "Y_TY" -- `Y^\mu T_w` for `w \in W(Y)` and `\mu \in Q^Y`

    EXAMPLES::

        sage: H = AffineHeckeAlgebra("A2")
        sage: TX = H.TX()
        sage: TX.an_element()
        TX[0,1,2] + 3*TX[0,1] + 2*TX[0] + 1
        sage: TY_Y = H.TY_Y()
        sage: TY_Y.an_element()
        Ty[1,2,1] # Y[(2, 2, 3)] + 3*Ty[1,2] # Y[(2, 2, 3)] + 3*Ty[2,1] # Y[(2, 2, 3)]
        sage: Y_TY = H.Y_TY()
        sage: Y_TY.an_element()
        Y[(2, 2, 3)] # Ty[1,2,1] + 3*Y[(2, 2, 3)] # Ty[1,2] + 3*Y[(2, 2, 3)] # Ty[2,1]

    There are built-in coercions between the bases::

        sage: m = TX.monomial(TX.basis().keys().an_element()); m
        TX[0,1,2]
        sage: Y_TY(m)
        Y[(1, 0, -1)] # Ty[1] + ((-v^2+1)/v)*Y[(1, 0, -1)] # 1
        sage: TY_Y(m)
        Ty[1] # Y[(0, 1, -1)]
        sage: TX(TY_Y(m))
        TX[0,1,2]

    Here is an example with a nonreduced root system and unequal parameters.

        sage: K = QQ['v,vl,v0'].fraction_field()
        sage: v,vl,v0=K.gens()
        sage: H = AffineHeckeAlgebra(['D',3,2], q1=Family(dict([[0,v0],[1,vl],[2,v]])), doubled_parameters=Family(dict([[2,v0-1/v0]])))
        sage: TX = H.TX()
        sage: TY_Y = H.TY_Y()
        sage: Y_TY = H.Y_TY()
        sage: TX.an_element()
        TX[0,1,2] + 3*TX[0,1] + 2*TX[0] + 1
        sage: TY_Y(TX.an_element())
        2*Ty[1,2,1] # Y[(-1, 0)] + 3*Ty[1,2] # Y[(0, -1)] + ((3*v0^2-3)/v0)*Ty[1] # Y[(0, 0)] + Ty[1] # Y[(0, 1)] + ((2*v0^2+v0-2)/v0)*1 # Y[(0, 0)]        
        sage: Y_TY(TX.an_element())
        Y[(0, 0)] # 1 + 2*Y[(1, 0)] # Ty[1,2,1] + ((-2*vl^2+3*vl+2)/vl)*Y[(1, 0)] # Ty[1,2] + ((-2*vl^2+2)/vl)*Y[(1, 0)] # Ty[2,1] + ((2*v^2*vl^2-3*v^2*vl-2*v^2+v*vl-2*vl^2+3*vl+2)/(v*vl))*Y[(1, 0)] # Ty[1] + ((2*vl^4-3*vl^3-4*vl^2+3*vl+2)/vl^2)*Y[(1, 0)] # Ty[2] + ((-2*v^2*vl^4+3*v^2*vl^3+2*v^2*vl^2-v*vl^3+2*vl^4-3*v^2*vl-3*vl^3-2*v^2+v*vl-2*vl^2+3*vl+2)/(v*vl^2))*Y[(1, 0)] # 1
        sage: TY_Y(TX.an_element()) == TY_Y(Y_TY(TX.an_element()))
        True
        sage: TX[1,0,1,0] == TX[0,1,0,1]
        True
        sage: TY_Y(TX[1,0,1,0]) == TY_Y(TX[0,1,0,1])
        True
        sage: Y_TY(TX[1,2,1,2]) == Y_TY(TX[2,1,2,1])
        True

    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, q1=None, q2=None, doubled_parameters=None):
        from sage.combinat.root_system.cartan_type import CartanType
        cartan_type = CartanType(cartan_type)
        if isinstance(q1, dict):
            q1 = Family(q1)
        if isinstance(q2, dict):
            q2 = Family(q2)
        return super(AffineHeckeAlgebra, cls).__classcall__(cls, cartan_type, q1, q2, doubled_parameters)

    def __init__(self, cartan_type, q1, q2, doubled_parameters):
        if cartan_type.is_reducible():
            raise ValueError, "Cartan type should be irreducible"
        if cartan_type.is_finite(): # a finite Cartan type is an abbreviation for its untwisted affinization
            cartan_type = cartan_type.affine()
        elif not cartan_type.is_affine():
            raise ValueError, "Cartan type must be finite or affine"
        self._cartan_type = cartan_type

        from sage.combinat.root_system.extended_affine_weyl_group import ExtendedAffineWeylGroup
        self._PvW0 = ExtendedAffineWeylGroup(cartan_type, style="PvW0")
        self._E = self._PvW0.realization_of()
        self._W0Pv = self._E.W0Pv()
        self._FW = self._E.FW()
        self._WX = self._E.affine_weyl()
        if cartan_type.is_untwisted_affine():
            self._cartan_type_Y = cartan_type.classical().dual()
        else:
            self._cartan_type_Y = cartan_type.classical()
        self._Y = self._cartan_type_Y.root_system().ambient_space()
        self._W0 = self._Y.weyl_group()

        I = cartan_type.index_set()
        self._base_ring, self._q1, self._q2 = ParameterFamilies(I, q1, q2)
        self._doubled_parameters = doubled_parameters

        Parent.__init__(self, category = AlgebrasWithBasis(self._base_ring).WithRealizations())

        # create the realizations (they are cached)
        TX = self.TX()
        TY_Y = self.TY_Y()
        Y_TY = self.Y_TY()
        # register coercion between TY_Y and Y_TY
        TY_Y.register_opposite(Y_TY)

        # coercion of TY into the affine Hecke algebra
        TY = TY_Y.factors()[0]
        TY_to_TX = TY.module_morphism(on_basis=TX.monomial * self.classical_weyl_to_affine_morphism, category=ModulesWithBasis(self._base_ring), codomain=TX)
        TY_to_TX.register_as_coercion()

        # coercion of group algebra of Y into affine Hecke algebra
        KY = TY_Y.factors()[1]

        def T_signs(mu):
            pi, word, signs = self._FW(mu.to_weight_space(ZZ)).alcove_walk_signs()
            if pi != pi.parent().one():
                raise ValueError, "%s should be in the root lattice"%mu
            return TX.product_by_signed_generator_sequence(TX.one(), word, tuple([-x for x in signs]))

        KY_to_TX = KY.module_morphism(on_basis=T_signs, category=ModulesWithBasis(self._base_ring), codomain=TX)
        KY_to_TX.register_as_coercion()

        TY_Y_to_TX = TY_Y.module_morphism(on_basis = lambda (w,mu): TY_to_TX(TY.monomial(w))*KY_to_TX(KY.monomial(mu)), category=ModulesWithBasis(self._base_ring), codomain=TX)
        TY_Y_to_TX.register_as_coercion()

        Y_TY_to_TX = Y_TY.module_morphism(on_basis = lambda (mu,w): KY_to_TX(KY.monomial(mu))*TY_to_TX(TY.monomial(w)), category=ModulesWithBasis(self._base_ring), codomain=TX)
        Y_TY_to_TX.register_as_coercion()

        def TX_to_TY_Y_func(w):
            i = w.first_descent(side="left")
            if i is None:
                return TY_Y.one()
            return TY_Y.algebra_generators()[i] * TX_to_TY_Y_func(w.apply_simple_reflection(i, side="left"))

        TX_to_TY_Y = TX.module_morphism(on_basis=TX_to_TY_Y_func, category=ModulesWithBasis(self._base_ring),codomain=TY_Y)
        TX_to_TY_Y.register_as_coercion()

        def TX_to_Y_TY_func(w):
            i = w.first_descent(side="right")
            if i is None:
                return Y_TY.one()
            return TX_to_Y_TY_func(w.apply_simple_reflection(i, side="right")) * Y_TY.algebra_generators()[i]

        TX_to_Y_TY = TX.module_morphism(on_basis=TX_to_Y_TY_func, category=ModulesWithBasis(self._base_ring),codomain=Y_TY)
        TX_to_Y_TY.register_as_coercion()

        KY_to_Y_TY = SetMorphism(Hom(KY,Y_TY,ModulesWithBasis(self._base_ring)),lambda y: Y_TY.factor_embedding(0)(y))
        KY_to_Y_TY.register_as_coercion()
        KY_to_TY_Y = SetMorphism(Hom(KY,TY_Y,ModulesWithBasis(self._base_ring)),lambda y: TY_Y.factor_embedding(1)(y))
        KY_to_TY_Y.register_as_coercion()

    @cached_method
    def TX(self):
        r"""
        Realizes ``self`` in "TX"-style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").TX()
            TX basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.AffineHeckeAlgebraTX()

    @cached_method
    def TY_Y(self):
        r"""
        Realizes ``self`` in "TY_Y" style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").TY_Y()
            TY_Y basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        return self.AffineHeckeAlgebraTY_Y()

    @cached_method
    def Y_TY(self):
        r"""
        Realizes ``self`` in "Y_TY" style.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").Y_TY()
            Y_TY basis of The affine Hecke algebra of type ['A', 2, 1]            
        """
        return self.AffineHeckeAlgebraY_TY()

    def a_realization(self):
        r"""
        Returns the default realization.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").a_realization()
            TX basis of The affine Hecke algebra of type ['A', 2, 1]            

        """
        return self.TX()

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: AffineHeckeAlgebra("A2")
            The affine Hecke algebra of type ['A', 2, 1]

        """
        return "The affine Hecke algebra of type %s"%self.cartan_type()

    def base_ring(self):
        return self._base_ring

    def cartan_type(self):
        return self._cartan_type

    @cached_method
    def index_set(self):
        r"""
        Returns the index set of the affine Dynkin diagram.

        EXAMPLES::

            sage: AffineHeckeAlgebra("B2").index_set()
            (0, 1, 2)

        """
        return self.cartan_type().index_set()

    def lattice(self):
        r"""
        Returns the `Y` lattice.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).lattice()
            Ambient space of the Root system of type ['A', 2]
            sage: AffineHeckeAlgebra(['A',5,2]).lattice()
            Ambient space of the Root system of type ['C', 3]
            sage: AffineHeckeAlgebra(['C',2,1]).lattice()
            Ambient space of the Root system of type ['B', 2]

        """
        return self._Y

    @cached_method
    def Y_algebra(self):
        r"""
        The group algebra of the `Y` lattice.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").Y_algebra()
            Group algebra of the Ambient space of the Root system of type ['A', 2] over Fraction Field of Univariate Polynomial Ring in v over Rational Field
        """
        return self.lattice().algebra(self.base_ring(), prefix="Y")

    def classical_weyl(self):
        r"""
        Returns the classical Weyl group of ``self``.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).classical_weyl()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
            sage: AffineHeckeAlgebra(['A',5,2]).classical_weyl()
            Weyl Group of type ['C', 3] (as a matrix group acting on the ambient space)
            sage: AffineHeckeAlgebra(['C',2,1]).classical_weyl()
            Weyl Group of type ['B', 2] (as a matrix group acting on the ambient space)

        """

        return self._W0

    @cached_method
    def TY(self):
        r"""
        The finite Hecke algebra of type `Y`.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").TY()
            Hecke algebra of type ['A', 2]

        """
        I0 = self._cartan_type_Y.index_set()
        return MultiParameterHeckeAlgebra(self._W0, self._q1, self._q2, prefix="Ty")

    def affine_weyl(self):
        r"""
        Returns the affine Weyl group of ``self``.

        EXAMPLES::

            sage: AffineHeckeAlgebra(['A',2,1]).affine_weyl()
            Weyl Group of type ['A', 2, 1] (as a matrix group acting on the root lattice)
            sage: AffineHeckeAlgebra(['A',5,2]).affine_weyl()
            Weyl Group of type ['B', 3, 1]^* (as a matrix group acting on the root lattice)
            sage: AffineHeckeAlgebra(['C',2,1]).affine_weyl()
            Weyl Group of type ['C', 2, 1] (as a matrix group acting on the root lattice)

        """
        return self._WX

    def classical_weyl_to_affine_morphism(self, w):
        r"""
        The image of `w` under the homomorphism from the classical Weyl group into the affine Weyl group.

        EXAMPLES::

            sage: H = AffineHeckeAlgebra("A2")
            sage: H.classical_weyl_to_affine_morphism(H.classical_weyl().an_element())
            S1*S2

        """
        return self.affine_weyl().from_reduced_word(w.reduced_word())

    def q1(self):
        r"""
        The Family of first eigenvalues for the `T_i` operators.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").q1()
            Finite family {0: v, 1: v, 2: v}

        """
        return self._q1

    def q2(self):
        r"""
        The Family of first eigenvalues for the `T_i` operators.

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").q2()
            Finite family {0: -1/v, 1: -1/v, 2: -1/v}


        """
        return self._q2

    class _BasesCategory(Category_realization_of_parent):
        r"""
        The category of realizations of an affine Hecke algebra
        """
        def __init__(self, basis, prefix=None):
            r"""
            Initialize a basis.

            INPUT:

            - ``basis`` -- a basis
            - ``prefix`` -- a prefix

            TESTS::

                sage: H = AffineHeckeAlgebra("A2")
                sage: bases = H._BasesCategory()
                sage: H.TX() in bases
                True
            """
            Category_realization_of_parent.__init__(self, basis)
            basis._prefix = prefix

        def super_categories(self):
            r"""
            EXAMPLES::

                sage: AffineHeckeAlgebra("A2")._BasesCategory().super_categories()
                [Category of realizations of The affine Hecke algebra of type ['A', 2, 1], Category of algebras with basis over Fraction Field of Univariate Polynomial Ring in v over Rational Field, Category of monoids with realizations, Category of additive unital additive magmas with realizations]

            """
            return [Realizations(self.base())]+self.base().category().super_categories()

        def _repr_(self):
            r"""
            Return the representation of ``self``.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2")._BasesCategory()
                Category of bases of The affine Hecke algebra of type ['A', 2, 1]
            """
            return "Category of bases of %s" % self.base()

        class ParentMethods:

            @abstract_method(optional=False)
            def algebra_generators(self):
                r"""
                The family of generators `T_i` in the given realization.

                It should be a Family with key set the Dynkin node set.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: H.TX().algebra_generators()
                    Finite family {0: TX[0], 1: TX[1], 2: TX[2]}
                    sage: H.TY_Y().algebra_generators()
                    Finite family {0: Ty[1,2,1] # Y[(-1, 0, 1)] + ((v^2-1)/v)*1 # Y[(0, 0, 0)], 1: Ty[1] # Y[(0, 0, 0)], 2: Ty[2] # Y[(0, 0, 0)]}
                    sage: K = QQ['v,vl'].fraction_field(); v,vl=K.gens()
                    sage: H = AffineHeckeAlgebra(['C',3,1],q1=Family(dict([[0,vl],[1,v],[2,v],[3,vl]])))
                    sage: H.TY_Y().algebra_generators()
                    Finite family {0: Ty[1,2,3,2,1] # Y[(-1, 0, 0)] + ((vl^2-1)/vl)*1 # Y[(0, 0, 0)], 1: Ty[1] # Y[(0, 0, 0)], 2: Ty[2] # Y[(0, 0, 0)], 3: Ty[3] # Y[(0, 0, 0)]}

                """
                pass

            @abstract_method(optional=False)
            def Y_morphism(self, y):
                r"""
                The image of `y` under the morphism from the group algebra of `Y` into ``self``.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: TX = H.TX()
                    sage: KY = H.Y_algebra()
                    sage: y = KY.monomial(KY.basis().keys().simple_root(1)); y
                    Y[(1, -1, 0)]
                    sage: z = TX.Y_morphism(y); z
                    TX[0,2,0,1] + ((-v^2+1)/v)*TX[0,2,1]                   
                    sage: TY_Y = H.TY_Y()
                    sage: TY_Y(z)
                    1 # Y[(1, -1, 0)]
                    sage: TY_Y(z) == TY_Y.Y_morphism(y)
                    True
                """
                pass

            @abstract_method(optional=False)
            def classical_hecke_morphism(self, a):
                r"""
                Returns the image of `a` from the finite Hecke algebra into ``self``.

                EXAMPLES::

                    sage: H = AffineHeckeAlgebra("A2")
                    sage: H0 = H.TY()
                    sage: h = H0.an_element(); h
                    Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                    sage: H.TX().classical_hecke_morphism(h)
                    TX[1,2,1] + 3*TX[2,1] + 3*TX[1,2]                    
                """
                pass

            def from_reduced_word(self, word):
                r"""
                Converts an affine or finite reduced word into a group element.

                .. warning::

                    Must be implemented in style "TX".

                EXAMPLES::

                    sage: AffineHeckeAlgebra("A2").TY_Y().from_reduced_word([0,2,1])
                    Ty[2] # Y[(1, -1, 0)]

                """
                return self(self.realization_of().TX().from_reduced_word(word))

    class _Bases(UniqueRepresentation, BindableClass):
        r"""
        The class of realizations of an affine Hecke algebra.
        """

        def _repr_(self):
            r"""
            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").TX() # indirect doctest
                TX basis of The affine Hecke algebra of type ['A', 2, 1]

            """
            return "%s basis of the %s"%(self._prefix, self.realization_of())

    class AffineHeckeAlgebraTX(MultiParameterHeckeAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "TX" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").TX()
            TX basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            MultiParameterHeckeAlgebra.__init__(self, E.affine_weyl(), E.q1(), E.q2(), prefix="TX", category=Category.join((E._BasesCategory(),AlgebrasWithBasis(E.base_ring()))))
            self._style = "TX"

        def _repr_(self):
            E = self.realization_of()
            return "%s basis of %s"%(self._style, E._repr_())

        def from_reduced_word(self, word):
            r"""
            Converts an affine or finite reduced word into a group element.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").TX().from_reduced_word([0,2,1])
                TX[0,2,1]
            """
            return self.monomial(self.realization_of().affine_weyl().from_reduced_word(word))

        def classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: H0 = H.TY()
                sage: h = H0.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.TX().classical_hecke_morphism(h)
                TX[1,2,1] + 3*TX[2,1] + 3*TX[1,2]

            """
            return self(a)

        def Y_morphism(self, y):
            r"""
            The image of the `Y`-lattice group algebra element `y` into ``self``.

            EXAMPLES::

                sage: H=AffineHeckeAlgebra("A2")
                sage: KY = H.Y_algebra()
                sage: z = H.TX().Y_morphism(KY.monomial(KY.basis().keys().simple_root(2))); z
                TX[0,1,0,2] + ((-v^2+1)/v)*TX[0,1,2]
                sage: H.Y_TY()(z)
                Y[(0, 1, -1)] # 1
                sage: H.TY_Y()(z)
                1 # Y[(0, 1, -1)]
                sage: H.Y_TY()(H.TY_Y()(z))
                Y[(0, 1, -1)] # 1

            """
            return self(self.realization_of().TY_Y()(y))

    class AffineHeckeAlgebraTY_Y(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "TY_Y" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").TY_Y()
            TY_Y basis of The affine Hecke algebra of type ['A', 2, 1]
        """

        def __init__(self, E):
            KY = E.Y_algebra()
            H = E.TY()
            module_category = ModulesWithBasis(E.base_ring())
            self._YoTY = tensor([KY,H],category=module_category)
            self._TYoY = tensor([H,KY],category=module_category)
            self._YM = KY.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention="antidominant", doubled_parameters=E._doubled_parameters, side="right")
            def right_action_on_TY_Y((w,mu), i):
                smu = mu.simple_reflection(i)
                return tensor([H.monomial(w), self._YM[i](KY.monomial(mu)) - E.q1()[i]*KY.monomial(smu)]) + tensor([H.product_by_generator_on_basis(w,i), KY.monomial(smu)])

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._TYoYM = HeckeAlgebraRepresentation(self._TYoY, right_action_on_TY_Y, E.cartan_type(), E.q1(), E.q2())

            def twist_func((mu,w)):
                return self._TYoYM.Tw(w)(self._TYoY.monomial((H.one_basis(), mu)))

            twist = self._YoTY.module_morphism(on_basis=twist_func,codomain=self._TYoY,category=module_category)
            SmashProductAlgebra.__init__(self, H, KY, twist, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))
            self._style = "TY_Y"

        def _repr_(self):
            E = self.realization_of()
            return "%s basis of %s"%(self._style, E._repr_())

        def Y_morphism(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of `Y` into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: KY = H.Y_algebra()
                sage: y = KY.an_element(); y
                Y[(2, 2, 3)]
                sage: H.TY_Y().Y_morphism(y)
                1 # Y[(2, 2, 3)]

            """
            return self.factor_embedding(1)(y)

        def classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: H0 = H.TY()
                sage: h = H0.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.TY_Y().classical_hecke_morphism(h)
                Ty[1,2,1] # Y[(0, 0, 0)] + 3*Ty[1,2] # Y[(0, 0, 0)] + 3*Ty[2,1] # Y[(0, 0, 0)]
            """
            return self.factor_embedding(0)(a)

        def T0(self):
            r"""
            The operator `T_0`.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").TY_Y().T0()
                Ty[1,2,1] # Y[(-1, 0, 1)] + ((v^2-1)/v)*1 # Y[(0, 0, 0)]
            """
            return self.realization_of().Y_TY().T0().to_opposite()

        @cached_method
        def algebra_generators(self):
            r"""
            The family of generators `T_i` in the given realization.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").TY_Y().algebra_generators()
                Finite family {0: Ty[1,2,1] # Y[(-1, 0, 1)] + ((v^2-1)/v)*1 # Y[(0, 0, 0)], 1: Ty[1] # Y[(0, 0, 0)], 2: Ty[2] # Y[(0, 0, 0)]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(0)(self.factors()[0].algebra_generators()[i])] for i in self.realization_of().index_set()]))


    class AffineHeckeAlgebraY_TY(SmashProductAlgebra, _Bases):
        r"""
        Affine Hecke algebra in "Y_TY" style.

        INPUT:

        - `E` -- Affine Hecke algebra parent

        EXAMPLES::

            sage: AffineHeckeAlgebra("A2").Y_TY()
            Y_TY basis of The affine Hecke algebra of type ['A', 2, 1]

        """
        def __init__(self, E):
            KY = E.Y_algebra()
            H = E.TY()
            module_category = ModulesWithBasis(E.base_ring())
            self._YoTY = tensor([KY,H],category=module_category)
            self._TYoY = tensor([H,KY],category=module_category)
            self._YM = KY.nonreduced_demazure_lusztig_operators(E.q1(), E.q2(), convention="antidominant", doubled_parameters=E._doubled_parameters, side="left")

            def left_action_on_Y_TY((mu,w), i):
                smu = mu.simple_reflection(i)
                return tensor([self._YM[i](KY.monomial(mu)) - E.q1()[i]*KY.monomial(smu), H.monomial(w)]) + tensor([KY.monomial(smu),H.product_by_generator_on_basis(w,i,side="left")])

            from sage.combinat.root_system.hecke_algebra_representation import HeckeAlgebraRepresentation
            self._YoTYM = HeckeAlgebraRepresentation(self._YoTY, left_action_on_Y_TY, E.cartan_type(), E.q1(), E.q2(), side="left")

            def untwist_func((w,mu)):
                return self._YoTYM.Tw(w)(self._YoTY.monomial((mu, H.one_basis())))

            untwist = self._TYoY.module_morphism(on_basis=untwist_func,codomain=self._YoTY,category=module_category)
            SmashProductAlgebra.__init__(self, KY, H, untwist, category=Category.join((E._BasesCategory(), AlgebrasWithBasis(E.base_ring()).TensorProducts())))

            self._style = "Y_TY"

        def _repr_(self):
            E = self.realization_of()
            return "%s basis of %s"%(self._style, E._repr_())

        def Y_morphism(self, y):
            r"""
            The image of `y` under the morphism from the group algebra of `Y` into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: KY = H.Y_algebra()
                sage: y = KY.an_element(); y
                Y[(2, 2, 3)]
                sage: H.Y_TY().Y_morphism(y)
                Y[(2, 2, 3)] # 1
            """
            return self.factor_embedding(0)(y)

        def classical_hecke_morphism(self, a):
            r"""
            Returns the image of `a` from the finite Hecke algebra into ``self``.

            EXAMPLES::

                sage: H = AffineHeckeAlgebra("A2")
                sage: H0 = H.TY()
                sage: h = H0.an_element(); h
                Ty[1,2,1] + 3*Ty[1,2] + 3*Ty[2,1]
                sage: H.Y_TY().classical_hecke_morphism(h)
                Y[(0, 0, 0)] # Ty[1,2,1] + 3*Y[(0, 0, 0)] # Ty[1,2] + 3*Y[(0, 0, 0)] # Ty[2,1]
            """
            return self.factor_embedding(1)(a)

        def T0(self):
            r"""
            The `T_0` operator.

            Let `\varphi` be the short dominant root of type `Y`. Then

            ..MATH::

                T_0^X = Y^\varphi T^{-1}_{s_\varphi}.

            """
            E = self.realization_of()
            phi = E._cartan_type_Y.root_system().coroot_lattice().highest_root().associated_coroot()
            s_phi = phi.associated_reflection()
            TY = self.factors()[1]
            return self.from_direct_product(self.factors()[0].monomial(phi.to_ambient()), TY.product_by_signed_generator_sequence(TY.one(), s_phi, [-1 for i in range(len(s_phi))]))

        @cached_method
        def algebra_generators(self):
            r"""
            The algebra generators `T_i`.

            EXAMPLES::

                sage: AffineHeckeAlgebra("A2").Y_TY().algebra_generators()
                Finite family {0: Y[(1, 0, -1)] # Ty[1,2,1] + ((-v^2+1)/v)*Y[(1, 0, -1)] # Ty[1,2] + ((-v^2+1)/v)*Y[(1, 0, -1)] # Ty[2,1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] # Ty[1] + ((v^4-2*v^2+1)/v^2)*Y[(1, 0, -1)] # Ty[2] + ((-v^6+2*v^4-2*v^2+1)/v^3)*Y[(1, 0, -1)] # 1, 1: Y[(0, 0, 0)] # Ty[1], 2: Y[(0, 0, 0)] # Ty[2]}

            """
            return Family(dict([[i, self.T0() if i == 0 else self.factor_embedding(1)(self.factors()[1].algebra_generators()[i])] for i in self.realization_of().index_set()]))
