r"""
Lie Conformal Algebras

.. include:: ../../../algebras/lie_conformal_algebras/lie_conformal_algebra_desc.rst

.. SEEALSO::

    :mod:`sage.algebras.lie_conformal_algebras.lie_conformal_algebra`

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation
"""

#******************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.modules import Modules
from .category_types import Category_over_base_ring
from sage.categories.category_with_axiom import \
                                     CategoryWithAxiom_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.element import coerce_binop
from sage.categories.category_with_axiom import all_axioms as all_axioms
from sage.misc.misc_c import prod
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()

class LieConformalAlgebras(Category_over_base_ring):
    r"""
    The category of Lie conformal algebras.

    This is the base category for all Lie conformal algebras.
    Subcategories with axioms are ``FinitelyGenerated`` and
    ``WithBasis``. A *finitely generated* Lie conformal algebra is a
    Lie conformal algebra over `R` which is finitely generated as an
    `R[T]`-module. A Lie conformal algebra *with basis* is one with a
    preferred basis as an `R`-module.

    EXAMPLES:

    The base category::

        sage: C = LieConformalAlgebras(QQ); C
        Category of Lie conformal algebras over Rational Field
        sage: C.is_subcategory(VectorSpaces(QQ))
        True

    Some subcategories::

        sage: LieConformalAlgebras(QQbar).FinitelyGenerated().WithBasis()
        Category of finitely generated Lie conformal algebras with basis over Algebraic Field

    In addition we support functorial constructions ``Graded`` and
    ``Super``. These functors commute::

        sage: LieConformalAlgebras(AA).Graded().Super()
        Category of super H-graded Lie conformal algebras over Algebraic Real Field
        sage: LieConformalAlgebras(AA).Graded().Super() is LieConformalAlgebras(AA).Super().Graded()
        True

    That is, we only consider gradings on super Lie conformal algebras
    that are compatible with the `\mathbb{Z}/2\mathbb{Z}` grading.
    All four subcategories commute in this sense::

        sage: C = LieConformalAlgebras(AA).Graded().FinitelyGenerated().WithBasis().Super()
        sage: D = LieConformalAlgebras(AA).Super().WithBasis().Graded().FinitelyGenerated()
        sage: C is D
        True
        sage: C
        Category of finitely generated super H-graded Lie conformal algebras with basis over Algebraic Real Field

    Vertex algebras are Lie conformal algebras::

        sage: V = FreeBosonsVertexAlgebra(QQ)
        sage: V in LieConformalAlgebras(QQ)
        True

    So are Poisson vertex algebras::

        sage: P = V.classical_limit(); P
        The classical limit of The Free Bosons vertex algebra with generators (alpha_-1|0>,) over Rational Field
        sage: P in LieConformalAlgebras(QQ)
        True

    The base ring needs to be a commutative ring::

        sage: LieConformalAlgebras(QuaternionAlgebra(2))
        Traceback (most recent call last):
        ValueError: base must be a commutative ring got Quaternion Algebra (-1, -1) with base ring Rational Field
    """

    @staticmethod
    def __classcall_private__(cls,R,check=True):
        r"""
        INPUT:

        - `R` -- a commutative ring.
        - ``check`` -- a boolean (default: ``True``); whether to check
          that `R` is a commutative ring.

        EXAMPLES::

            sage: LieConformalAlgebras(QuaternionAlgebra(2))
            Traceback (most recent call last):
            ValueError: base must be a commutative ring got Quaternion Algebra (-1, -1) with base ring Rational Field
            sage: LieConformalAlgebras(ZZ)
            Category of Lie conformal algebras over Integer Ring
        """
        if check:
            if not (R in _CommutativeRings):
                    raise ValueError("base must be a commutative ring got {}"\
                                     .format(R))
        return super(LieConformalAlgebras,cls).__classcall__(cls,R)

    @cached_method
    def super_categories(self):
        """
        The list of super categories of this category.

        EXAMPLES::

            sage: C = LieConformalAlgebras(QQ)
            sage: C.super_categories()
            [Category of vector spaces over Rational Field]
            sage: C = LieConformalAlgebras(QQ).FinitelyGenerated(); C
            Category of finitely generated Lie conformal algebras over Rational Field
            sage: C.super_categories()
            [Category of Lie conformal algebras over Rational Field]
            sage: C.all_super_categories()
            [Category of finitely generated Lie conformal algebras over Rational Field,
             Category of Lie conformal algebras over Rational Field,
             Category of vector spaces over Rational Field,
             Category of modules over Rational Field,
             Category of bimodules over Rational Field on the left and Rational Field on the right,
             Category of right modules over Rational Field,
             Category of left modules over Rational Field,
             Category of commutative additive groups,
             Category of additive groups,
             Category of additive inverse additive unital additive magmas,
             Category of commutative additive monoids,
             Category of additive monoids,
             Category of additive unital additive magmas,
             Category of commutative additive semigroups,
             Category of additive commutative additive magmas,
             Category of additive semigroups,
             Category of additive magmas,
             Category of sets,
             Category of sets with partial maps,
             Category of objects]
        """
        return [Modules(self.base_ring())]

    def _repr_object_names(self):
        """
        The name of the objects of this category.

        EXAMPLES::

            sage: LieConformalAlgebras(QQ)
            Category of Lie conformal algebras over Rational Field
            sage: LieConformalAlgebras(QQ).Graded().FinitelyGenerated()
            Category of finitely generated H-graded Lie conformal algebras over Rational Field
        """
        return "Lie conformal algebras over {}".format(self.base_ring())

    class ParentMethods:

        def set_lift(self, liftmorphism):
            r"""
            Register ``liftmorphsm`` as a coercion between this Lie
            conformal algebra and its universal enveloping vertex
            algebra.

            INPUT:

            - ``liftmorphism`` a function; a function from this Lie
              conformal algebra to its universal enveloping vertex
              algebra.

            .. NOTE::

                This method should not be used directly, instead the
                user should use :meth:`register_lift<sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra.UniversalEnvelopingVertexAlgebra.register_lift>`

            EXAMPLES::

                sage: L = VirasoroLieConformalAlgebra(QQ);
                sage: V = L.universal_enveloping_algebra()
                sage: L.lift
                Generic morphism:
                  From: The Virasoro Lie conformal algebra over Rational Field
                  To:   The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
                sage: W = VirasoroVertexAlgebra(QQ,1/2)
                sage: L.lift
                Generic morphism:
                  From: The Virasoro Lie conformal algebra over Rational Field
                  To:   The Virasoro vertex algebra of central charge 1/2 over Rational Field
                sage: V.register_lift()
                sage: L.lift
                Generic morphism:
                  From: The Virasoro Lie conformal algebra over Rational Field
                  To:   The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

            TESTS::

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: R = NeveuSchwarzLieConformalAlgebra(QQ);
                sage: R.set_lift(V._lca.lift)
                Traceback (most recent call last):
                ...
                ValueError: liftmorphism must be a morphism between The Neveu-Schwarz super Lie conformal algebra over Rational Field and its universal enveloping algebra
            """
            from sage.categories.morphism import Morphism
            if not isinstance(liftmorphism, Morphism):
                raise TypeError("liftmorphism must be a morphism.")

            V = liftmorphism.codomain()
            from sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra\
                import UniversalEnvelopingVertexAlgebra as UEA
            if not isinstance(V, UEA):
                raise TypeError("codomain must be a universal enveloping "\
                                "vertex algebra.")

            if liftmorphism.domain() is not self or V._lca is not self:
                raise ValueError("liftmorphism must be a morphism between {}"\
                           " and its universal enveloping algebra".format(self))
            f = None
            try:
                f = self.lift
            except AttributeError:
                pass

            if f is not None and f.codomain() is V:
                return 

            self.lift = liftmorphism
            try: 
                self.lift.register_as_coercion()
            except AssertionError:
                #We have already constructed previously a coercion
                pass

        @abstract_method
        def universal_enveloping_algebra(self,
                                        central_parameters=None,
                                        names=None):
            r"""
            The universal enveloping vertex algebra of this Lie
            conformal algebra.

            INPUT:

            - ``central_parameters`` -- A family, a ``dict`` or
              ``None`` (default: ``None``); If not ``None``, this
              parameter must be a family of constants in the
              base ring of this Lie conformal algebra parametrized by
              the central elements. This family encodes a central
              character of its universal enveloping vertex algebra. If
              ``None``, the *zero* central character is used.

            - ``names`` -- The names of the generators of the universal
              enveloping vertex algebra.

            OUTPUT:

            If `L` is a Lie conformal algebra over `R` with some
            central elements `C_i \in L` indexed by a set `I`,
            ``central_parameters`` is a family of elements `c_i \in R`
            indexed by the same set, then this method returns the
            central quotient of the universal enveloping vertex algebra
            of `L` by the ideal generated by `C_i - c_i |0\rangle`.

            EXAMPLES:

            We construct the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra of central charge `0` over
            the real algebraic numbers::

                sage: Vir = VirasoroLieConformalAlgebra(AA)
                sage: V = Vir.universal_enveloping_algebra(); V
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Algebraic Real Field
                sage: V.0.bracket(V.0)
                {0: L_-3|0>, 1: 2*L_-2|0>}

            We construct the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra of central charge 2 over the
            rationals::

                sage: Vir = VirasoroLieConformalAlgebra(QQ);
                sage: Vir.gens()
                (L, C)
                sage: cp = Family({Vir.1:2})
                sage: V = Vir.universal_enveloping_algebra(cp)
                sage: V.0.nproduct(V.0,3)
                |0>

            The Virasoro Algebra is not defined over `\ZZ`::

                sage: Vir = VirasoroLieConformalAlgebra(ZZ)
                Traceback (most recent call last):
                ...
                ArithmeticError: inverse does not exist

            The universal enveloping vertex algebra of the Neveu-Schwarz
            Lie conformal algebra::

                sage: NS = NeveuSchwarzLieConformalAlgebra(QQbar)
                sage: NS.inject_variables()
                Defining L, G, C
                sage: V = NS.universal_enveloping_algebra({C:1})
                sage: G*G
                L_-3|0>

            The list of central parameters needs to be indexed by the
            central elements of the algebra::

                sage: Vir = VirasoroLieConformalAlgebra(QQ);
                sage: cp = Family({Vir.0 : 3})
                sage: V = Vir.universal_enveloping_algebra(cp)
                Traceback (most recent call last):
                ...
                ValueError: central_parameters must be parametrized by central elements
            """
            raise NotImplementedError("the universal enveloping vertex " +
                           "algebra of {} is not implemented yet".format(self))

        @abstract_method
        def ideal(self, *gens, **kwds):
            """
            The ideal of this conformal algebra generated by ``gens``.

            .. TODO::

                Ideals of Lie Conformal Algebras are not implemented
                yet.
            """
            raise NotImplementedError("Ideals of Lie Conformal algebras are "\
                                      "not implemented yet")

        def is_super(self):
            """
            Wether this Lie conformal algebra is a super Lie
            conformal algebra.

            EXAMPLES::

                sage: V = VirasoroLieConformalAlgebra(QQ)
                sage: V.is_super()
                False
                sage: NeveuSchwarzLieConformalAlgebra(QQbar).is_super()
                True

            Notice that we can force to have a *purely even* super Lie
            conformal algebra::

                sage: bosondict = {('a','a'):{1:{('K',0):1}}}
                sage: R = LieConformalAlgebra(QQ,bosondict,names=('a',),central_elements=('K',),super=True)
                sage: R.is_super()
                True
                sage: [g.is_even_odd() for g in R.gens()]
                [0, 0]
            """
            return self in LieConformalAlgebras(self.base_ring()).Super()

        def is_graded(self):
            """
            Wether this Lie conformal algebra is graded or not

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: Vir
                The Virasoro Lie conformal algebra over Rational Field
                sage: Vir.is_graded()
                True
                sage: WeylLieConformalAlgebra(QQ).is_graded()
                False
            """
            return self in LieConformalAlgebras(self.base_ring()).Graded()

    class ElementMethods:
        @coerce_binop
        def bracket(self,rhs):
            r"""
            The `\lambda`-bracket of these two elements.

            EXAMPLES:

            The brackets of the Virasoro Lie conformal Algebra::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L.bracket(L)
                {0: TL, 1: 2*L, 3: 1/2*C}
                sage: L.bracket(L.T())
                {0: 2*T^(2)L, 1: 3*TL, 2: 4*L, 4: 2*C}

            Now with a current algebra::

                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: V.gens()
                (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                sage: E = V.0; H = V.1; F = V.2;
                sage: H.bracket(H)
                {1: 2*B['K']}
                sage: E.bracket(F)
                {0: B[alphacheck[1]], 1: B['K']}

            With a super Lie Conformal Algebra::

                sage: R = FermionicGhostsLieConformalAlgebra(QQbar); R
                The Fermionic ghosts Lie conformal algebra with generators (b, c, K) over Algebraic Field
                sage: R.inject_variables()
                Defining b, c, K
                sage: b.bracket(c)
                {0: K}
                sage: c.bracket(b)
                {0: K}

            This methods coerces elements to a common parent::

                sage: R = NeveuSchwarzLieConformalAlgebra(QQbar)
                sage: R.inject_variables()
                Defining L, G, C
                sage: V = R.universal_enveloping_algebra({C:1})
                sage: L.bracket(L*G)
                {0: L_-3G_-3/2|0> + L_-2G_-5/2|0>,
                 1: 7/2*L_-2G_-3/2|0>,
                 2: 3*G_-5/2|0>,
                 3: 13/2*G_-3/2|0>}

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user
                needs to implement :meth:`_bracket_`
            """
            return self._bracket_(rhs)

        @abstract_method
        def _bracket_(self,rhs):
            r"""
            Returns the `\lambda`-bracket of these two elements.

            EXAMPLES:

            The brackets of the Virasoro Lie conformal Algebra::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L._bracket_(L)
                {0: TL, 1: 2*L, 3: 1/2*C}
                sage: L._bracket_(L.T())
                {0: 2*T^(2)L, 1: 3*TL, 2: 4*L, 4: 2*C}

            Now with a current algebra::

                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: V.gens()
                (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                sage: E = V.0; H = V.1; F = V.2;
                sage: H._bracket_(H)
                {1: 2*B['K']}
                sage: E._bracket_(F)
                {0: B[alphacheck[1]], 1: B['K']}

            .. NOTE::

                It is guaranteed that both are elements of the same
                parent.
            """
            raise NotImplementedError("Not implemented")

        @coerce_binop
        def nproduct(self,rhs,n):
            r"""
            Returns the n-th product of these two elements.

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L.nproduct(L,3)
                1/2*C
                sage: L.nproduct(L.T(),0)
                2*T^(2)L
                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: E = V.0; H = V.1; F = V.2;
                sage: E.nproduct(H,0) == - 2*E
                True
                sage: E.nproduct(F,1)
                B['K']

            when `n<0`  we obtain an element of the universal enveloping
            vertex algebra, however, the universal enveloping vertex
            algebra needs to be constructed first::

                sage: R = AffineLieConformalAlgebra(AA,'A1', names=('e','h','f'))
                sage: R.inject_variables()
                Defining e, h, f, K
                sage: e.nproduct(f,1)
                K
                sage: e.nproduct(f,-1)
                Traceback (most recent call last):
                ...
                NotImplementedError: In order to lift an element first need to construct the universal enveloping vertex algebra
                sage: V = R.universal_enveloping_algebra()
                sage: e.nproduct(f,-1)
                e_-1f_-1|0>

            The product of the two fermions in the `N=2` vertex algebra
            has conformal weight `3`::

                sage: R = N2LieConformalAlgebra(QQ)
                sage: R.inject_variables()
                Defining L, J, G1, G2, C
                sage: V = R.universal_enveloping_algebra({C:1})
                sage: L.nproduct(G1*G2,1)
                3*G1_-3/2G2_-3/2|0>

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user
                needs to implement :meth:`_nproduct_`
            """
            return self._nproduct_(rhs,n)

        def _nproduct_(self,rhs,n):
            r"""
            Returns the n-th product of this two elements.

            If `n\geq 0` it returns the element of this Lie conformal
            algebra. If `n < 0` then it first lifts this element to the
            universal  enveloping vertex algebra and returns the
            corresponding element

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L._nproduct_(L,3)
                1/2*C
                sage: L._nproduct_(L.T(),0)
                2*T^(2)L
                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: E = V.0; H = V.1; F = V.2;
                sage: E._nproduct_(H,0) == - 2*E
                True
                sage: E._nproduct_(F,1)
                B['K']

            when `n<0`  we obtain an element of the universal enveloping
            vertex algebra, however, the universal enveloping vertex
            algebra needs to be constructed first::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L._nproduct_(L,-3)
                L_-4L_-2|0>

            .. NOTE::

                It is guaranteed that both are elements of the same
                parent.
            """
            if n >= 0:
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                return self.lift().nproduct(rhs,n)

        def _mul_(self,right):
            """
            The normally ordered product of these two elements in the
            universal enveloping vertex algebra.

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: Vir.inject_variables()
                Defining L, C
                sage: V = Vir.universal_enveloping_algebra({C:1})
                sage: L*L.T()
                L_-3L_-2|0> + L_-5|0>

            Note that the universal enveloping algebra needs to be
            constructed first::

                sage: W = WeylLieConformalAlgebra(QQ,4)
                sage: W.inject_variables()
                Defining alpha0, alpha1, alpha2, alpha3, K
                sage: alpha0*alpha1
                Traceback (most recent call last):
                NotImplementedError: In order to lift an element first need to construct the universal enveloping vertex algebra
                sage: V = W.universal_enveloping_algebra({K:1})
                sage: alpha0*alpha1
                alpha0_(-1)alpha1_(-1)|0>
            """
            return self.lift()*right.lift()

        @abstract_method(optional=True)
        def lift(self):
            r"""
            The image of this element under the canonical lift
            to the universal enveloping vertex algebra.

            .. WARNING::

                The universal enveloping algebra needs to be constructed
                first for this morphism to be defined.

                This morphism is registered as a coercion between this
                Lie conformal algebra and its universal enveloping
                vertex algebra upon creation. Since we consider central
                quotients of the universal enveloping vertex algebras
                by fixed central parameters, each time a different
                universal enveloping vertex algebra is constructed,
                this lift morphism is changed. See the examples below
                and also :meth:`register_lift\
                <sage.algebras.vertex_algebras.universal_enveloping_\
                vertex_algebra.UniversalEnvelopingVertexAlgebra.\
                register_lift>`.


            EXAMPLES:

            We lift to the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra with central charge `0`::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L.lift()
                L_-2|0>
                sage: L.lift().__class__
                <class 'sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

            Notice that the target of the ``lift`` morphism changes when
            we construct another universal enveloping vertex algebra::

                sage: V = VirasoroVertexAlgebra(QQ,1);
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 1 over Rational Field
                sage: V = VirasoroVertexAlgebra(QQ,3)
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 3 over Rational Field

            Notice that recreation may not re-establish the right coercion
            depending on the method of construction::

                sage: cp = Family({Vir.1:1/2}); V = Vir.universal_enveloping_algebra(cp)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 1/2 over Rational Field
            """
            raise NotImplementedError("Not implemented")

        @abstract_method(optional=False)
        def T(self, n=1):
            r"""
            The n-th derivative of this element.

            INPUT:

            - ``n`` -- integer (default:``1``); how many times
              to apply `T` to this element.

            OUTPUT:

            `T^n a` where `a` is this element. Notice that We use the
            *divided powers* notation `T^{(j)} = \frac{T^j}{j!}`

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ)
                sage: Vir.inject_variables()
                Defining L, C
                sage: L.T()
                TL
                sage: L.T(3)
                6*T^(3)L
                sage: C.T()
                0
            """
            raise NotImplementedError("Not implemented")

    class SubcategoryMethods:

        def FinitelyGeneratedAsLieConformalAlgebra(self):
            """
            The subcategory of finitely generated Lie conformal
            algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated()
                Category of finitely generated Lie conformal algebras over Rational Field
            """
            return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")


        def FinitelyGenerated(self):
            """
            The subcategory of finitely generated Lie conformal algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated()
                Category of finitely generated Lie conformal algebras over Rational Field
            """
            return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

        def WithBasis(self):
            """
            The subcategory of Lie conformal algebras with a preferred
            basis.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).WithBasis()
                Category of Lie conformal algebras over Rational Field with basis
            """
            return self._with_axiom("WithBasis")

        def Graded(self, base_ring=None):
            """
            The subcategory of H-Graded Lie conformal algebras.

            EXAMPLES::

                sage: C = LieConformalAlgebras(ZZ).WithBasis().Graded(); C
                Category of H-graded Lie conformal algebras with basis over Integer Ring
                sage: D = LieConformalAlgebras(ZZ).Graded().WithBasis()
                sage: D is C
                True
            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self, CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsLieConformalAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Graded()._with_axioms(axioms)
            return GradedModulesCategory.category_of(self)

        def Super(self, base_ring=None):
            """
            The subcategory of super Lie conformal algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(AA).Super().WithBasis()
                Category of super Lie conformal algebras with basis over Algebraic Real Field
            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self,CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsLieConformalAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Super()._with_axioms(axioms)
            return SuperModulesCategory.category_of(self)

    class Super(SuperModulesCategory):
        """
        The subcategory of super Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).Super().WithBasis()
            Category of super Lie conformal algebras with basis over Algebraic Real Field
        """
        #Need to do all this to make Super commute with Graded.
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                [Category of finitely generated super Lie conformal algebras over Rational Field,
                 Category of super H-graded Lie conformal algebras over Rational Field,
                 Category of finitely generated H-graded Lie conformal algebras over Rational Field]
            """
            return [self.base_category(),]

        class SubcategoryMethods:

            def Graded(self, base_ring=None):
                """
                The subcategory of H-graded super Lie conformal
                algebras.

                EXAMPLES::

                    sage: LieConformalAlgebras(QQ).Super().Graded()
                    Category of super H-graded Lie conformal algebras over Rational Field
                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                        "FinitelyGeneratedAsLieConformalAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Graded()._with_axioms(axioms)

                return GradedModulesCategory.category_of(
                                         self.base_category()).Super()


        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of Super Lie conformal algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Super().WithBasis()
                Category of super Lie conformal algebras with basis over Integer Ring
                sage: LieConformalAlgebras(ZZ).Super().WithBasis() is LieConformalAlgebras(ZZ).WithBasis().Super()
                True
            """

            class FinitelyGeneratedAsLieConformalAlgebra(
                                                CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated super Lie
                conformal algebras with basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Super().FinitelyGenerated().WithBasis()
                    Category of finitely generated super Lie conformal algebras with basis over Integer Ring
                """
                pass

        class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated super Lie
            conformal algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Super().FinitelyGenerated()
                Category of finitely generated super Lie conformal algebras over Integer Ring
            """
            pass


    class Graded(GradedModulesCategory):
        """
        The subcategory of H-graded Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).Graded()
            Category of H-graded Lie conformal algebras over Algebraic Field
        """

        def _repr_object_names(self):
            """
            The names of the objects of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).Super().Graded()
                Category of super H-graded Lie conformal algebras over Algebraic Field
                sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
                Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())


        class SubcategoryMethods:

            def Super(self, base_ring=None):
                """
                The subcategory of H-graded super Lie conformal algebras

                EXAMPLES::

                    sage: LieConformalAlgebras(QQ).Super().Graded()
                    Category of super H-graded Lie conformal algebras over Rational Field
                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                      "FinitelyGeneratedAsLieConformalAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Super()._with_axioms(axioms)
                return SuperModulesCategory.category_of(self)

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of H-graded Lie conformal algebras with
            basis.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Graded().WithBasis()
                Category of H-graded Lie conformal algebras with basis over Integer Ring
            """

            class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated H-graded
                Lie conformal algebras with basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
                    Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
                """
                pass

        class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated H-graded Lie
            conformal algebras.

            EXAMPLES::

                sage: C = LieConformalAlgebras(ZZ).Graded().FinitelyGenerated(); C
                Category of finitely generated H-graded Lie conformal algebras over Integer Ring
                sage: C is LieConformalAlgebras(ZZ).FinitelyGenerated().Graded()
                True
            """
            pass

        class Super(SuperModulesCategory):
            """
            The subcategory of H-graded super Lie conformal algebras.

            EXAMPLES::

                sage: C = LieConformalAlgebras(QQbar).Graded().Super(); C
                Category of super H-graded Lie conformal algebras over Algebraic Field
                sage: C is LieConformalAlgebras(QQbar).Super().Graded()
                True
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: LieConformalAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                    [Category of finitely generated super Lie conformal algebras over Rational Field,
                     Category of super H-graded Lie conformal algebras over Rational Field,
                     Category of finitely generated H-graded Lie conformal algebras over Rational Field]
                """
                return [self.base_category(),]


        class ElementMethods:
            @abstract_method
            def degree(self):
                """
                The degree of this element.

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                    sage: L.degree()
                    2
                    sage: L.T(3).degree()
                    5

                    sage: L = N2LieConformalAlgebra(QQbar); L
                    The N=2 super Lie conformal algebra over Algebraic Field
                    sage: L.inject_variables()
                    Defining L, J, G1, G2, C
                    sage: G1.degree()
                    3/2
                """
                raise NotImplementedError("Not implemented")

    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The subcategory of Lie conformal algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).WithBasis()
            Category of Lie conformal algebras over Algebraic Field with basis
        """
        def _repr_object_names(self):
            """
            The names of objects of this category.
            """
            return "{} with basis".format(self.base_category().\
                                          _repr_object_names())

        class ElementMethods:

            def index(self):
                """
                The index of this basis element.

                EXAMPLES::

                    sage: V = NeveuSchwarzLieConformalAlgebra(QQ)
                    sage: V.inject_variables()
                    Defining L, G, C
                    sage: G.T(3).index()
                    ('G', 3)
                    sage: v = V.an_element(); v
                    L + G + C
                    sage: v.index()
                    Traceback (most recent call last):
                    ...
                    ValueError: index can only be computed for monomials, got L + G + C
                """
                if self.is_zero():
                    return None
                if not self.is_monomial():
                    raise ValueError ("index can only be computed for "\
                                      "monomials, got {}".format(self))

                return next(iter(self.monomial_coefficients()))

        class SubcategoryMethods:
            def FinitelyGenerated(self):
                """
                The subcategory of finitely generated Lie conformal
                algebras with basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                    Category of finitely generated Lie conformal algebras with basis over Algebraic Field
                """
                return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

        class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated Lie conformal
            algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                Category of finitely generated Lie conformal algebras with basis over Algebraic Field
            """
            pass


    class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
        """
        The subcategory of finitely generated Lie conformal
        algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
            Category of finitely generated Lie conformal algebras with basis over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of objects of this category.
            """
            return "finitely generated {}".format(self.base_category().\
                            _repr_object_names())

        class ParentMethods:
            def ngens(self):
                r"""
                The number of generators of this Lie conformal algebra.

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ)
                    sage: Vir.ngens()
                    2

                    sage: V = AffineLieConformalAlgebra(QQ, 'A2')
                    sage: V.ngens()
                    9

                    sage: L = N2LieConformalAlgebra(QQbar); L
                    The N=2 super Lie conformal algebra over Algebraic Field
                    sage: L.gens()
                    (L, J, G1, G2, C)
                    sage: L.ngens()
                    5
                """
                return len(self.gens())

            def gen(self,i):
                r"""
                The i-th generator of this Lie conformal algebra.

                EXAMPLES::

                    sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                    sage: V.gens()
                    (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                    sage: V.gen(0)
                    B[alpha[1]]
                    sage: V.1
                    B[alphacheck[1]]
                """
                return self.gens()[i]
