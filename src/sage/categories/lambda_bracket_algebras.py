r"""
Lambda Bracket Algebras

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation.
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

from .category_types import Category_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.categories.modules import Modules
from sage.structure.element import coerce_binop
from sage.misc.cachefunc import cached_method
from sage.categories.commutative_rings import CommutativeRings
from sage.misc.lazy_import import LazyImport
_CommutativeRings = CommutativeRings()


class LambdaBracketAlgebras(Category_over_base_ring):
    r"""
    The category of Lambda bracket algebras.

    This is an abstract base category for Lie conformal algebras and
    super Lie conformal algebras.

    """
    @staticmethod
    def __classcall_private__(cls, R, check=True):
        r"""
        INPUT:

        - `R` -- a commutative ring
        - ``check`` -- a boolean (default: ``True``); whether to check
          that `R` is a commutative ring

        EXAMPLES::

            sage: LieConformalAlgebras(QuaternionAlgebra(2))
            Traceback (most recent call last):
            ValueError: base must be a commutative ring got Quaternion Algebra (-1, -1) with base ring Rational Field
            sage: LieConformalAlgebras(ZZ)
            Category of Lie conformal algebras over Integer Ring
        """
        if check:
            if not (R in _CommutativeRings):
                    raise ValueError("base must be a commutative ring got {}".format(R))
        return super(LambdaBracketAlgebras, cls).__classcall__(cls, R)

    @cached_method
    def super_categories(self):
        """
        The list of super categories of this category.

        EXAMPLES::

            sage: from sage.categories.lambda_bracket_algebras import LambdaBracketAlgebras
            sage: LambdaBracketAlgebras(QQ).super_categories()
            [Category of vector spaces over Rational Field]
        """
        return [Modules(self.base_ring())]

    def _repr_object_names(self):
        """
        The name of the objects of this category.

        EXAMPLES::

            sage: from sage.categories.lambda_bracket_algebras import LambdaBracketAlgebras
            sage: LambdaBracketAlgebras(QQ)
            Category of Lambda bracket algebras over Rational Field
        """
        return "Lambda bracket algebras over {}".format(self.base_ring())

    class SubcategoryMethods:

        def FinitelyGeneratedAsLambdaBracketAlgebra(self):
            """
            The category of finitely generated Lambda bracket algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated()
                Category of finitely generated lie conformal algebras over Rational Field
            """
            return self._with_axiom("FinitelyGeneratedAsLambdaBracketAlgebra")

        def FinitelyGenerated(self):
            """
            The category of finitely generated Lambda bracket algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated()
                Category of finitely generated lie conformal algebras over Rational Field
            """
            return self._with_axiom("FinitelyGeneratedAsLambdaBracketAlgebra")

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

                sage: L = lie_conformal_algebras.Virasoro(QQ);
                sage: V = L.universal_enveloping_algebra()
                sage: L.lift
                Generic morphism:
                  From: The Virasoro Lie conformal algebra over Rational Field
                  To:   The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
                sage: W = vertex_algebras.Virasoro(QQ,1/2)
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

                sage: V = vertex_algebras.Virasoro(QQ,1/2)
                sage: R = lie_conformal_algebras.NeveuSchwarz(QQ);
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

        @abstract_method(optional=True)
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

                sage: Vir = lie_conformal_algebras.Virasoro(AA)
                sage: V = Vir.universal_enveloping_algebra(); V
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Algebraic Real Field
                sage: V.0.bracket(V.0)
                {0: L_-3|0>, 1: 2*L_-2|0>}

            We construct the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra of central charge 2 over the
            rationals::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ);
                sage: Vir.gens()
                (L, C)
                sage: cp = Family({Vir.1:2})
                sage: V = Vir.universal_enveloping_algebra(cp)
                sage: V.0.nproduct(V.0,3)
                |0>

            The Virasoro Algebra is not defined over `\ZZ`::

                sage: Vir = lie_conformal_algebras.Virasoro(ZZ)
                Traceback (most recent call last):
                ...
                ArithmeticError: inverse does not exist

            The universal enveloping vertex algebra of the Neveu-Schwarz
            Lie conformal algebra::

                sage: NS = lie_conformal_algebras.NeveuSchwarz(QQbar)
                sage: NS.inject_variables()
                Defining L, G, C
                sage: V = NS.universal_enveloping_algebra({C:1})
                sage: G*G
                L_-3|0>

            The list of central parameters needs to be indexed by the
            central elements of the algebra::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ);
                sage: cp = Family({Vir.0 : 3})
                sage: V = Vir.universal_enveloping_algebra(cp)
                Traceback (most recent call last):
                ...
                ValueError: central_parameters must be parametrized by central elements
            """
            raise NotImplementedError("the universal enveloping vertex " +
                           "algebra of {} is not implemented yet".format(self))

        def ideal(self, *gens, **kwds):
            r"""
            The ideal of this Lambda bracket algebra generated by ``gens``.

            .. TODO::

                Ideals of Lie Conformal Algebras are not implemented yet.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)
                sage: Vir.ideal()
                Traceback (most recent call last):
                ...
                NotImplementedError: ideals of Lie Conformal algebras are not implemented yet
            """
            raise NotImplementedError("ideals of Lie Conformal algebras are "
                                      "not implemented yet")
    class ElementMethods:

        @coerce_binop
        def bracket(self,rhs):
            r"""
            The `\lambda`-bracket of these two elements.

            EXAMPLES:

            The brackets of the Virasoro Lie conformal Algebra::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                sage: L.bracket(L)
                {0: TL, 1: 2*L, 3: 1/2*C}
                sage: L.bracket(L.T())
                {0: 2*T^(2)L, 1: 3*TL, 2: 4*L, 4: 2*C}

            Now with a current algebra::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
                sage: V.gens()
                (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                sage: E = V.0; H = V.1; F = V.2;
                sage: H.bracket(H)
                {1: 2*B['K']}
                sage: E.bracket(F)
                {0: B[alphacheck[1]], 1: B['K']}

            With a super Lie Conformal Algebra::

                sage: R = lie_conformal_algebras.FermionicGhosts(QQbar); R
                The Fermionic ghosts Lie conformal algebra with generators (b, c, K) over Algebraic Field
                sage: R.inject_variables()
                Defining b, c, K
                sage: b.bracket(c)
                {0: K}
                sage: c.bracket(b)
                {0: K}

            This methods coerces elements to a common parent::

                sage: R = lie_conformal_algebras.NeveuSchwarz(QQbar)
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
        def _bracket_(self, rhs):
            r"""
            The `\lambda`-bracket of these two elements.

            .. NOTE::

                It is guaranteed that both are elements of the same
                parent.

            EXAMPLES:

            The brackets of the Virasoro Lie conformal Algebra::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                sage: L._bracket_(L)
                {0: TL, 1: 2*L, 3: 1/2*C}
                sage: L._bracket_(L.T())
                {0: 2*T^(2)L, 1: 3*TL, 2: 4*L, 4: 2*C}

            Now with a current algebra::

                sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
                sage: V.gens()
                (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
                sage: E = V.0; H = V.1; F = V.2;
                sage: H._bracket_(H)
                {1: 2*B['K']}
                sage: E._bracket_(F)
                {0: B[alphacheck[1]], 1: B['K']}
            """

        @coerce_binop
        def nproduct(self,rhs,n):
            r"""
            The ``n``-th product of these two elements.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                sage: L.nproduct(L,3)
                1/2*C
                sage: L.nproduct(L.T(),0)
                2*T^(2)L
                sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
                sage: E = V.0; H = V.1; F = V.2;
                sage: E.nproduct(H,0) == - 2*E
                True
                sage: E.nproduct(F,1)
                B['K']

            when `n<0`  we obtain an element of the universal enveloping
            vertex algebra, however, the universal enveloping vertex
            algebra needs to be constructed first::

                sage: R = lie_conformal_algebras.Affine(AA,'A1', names=('e','h','f'))
                sage: R.inject_variables()
                Defining e, h, f, K
                sage: e.nproduct(f,1)
                K
                sage: e.nproduct(f,-1)
                Traceback (most recent call last):
                ...
                NotImplementedError: in order to lift an element first need to construct the universal enveloping vertex algebra
                sage: V = R.universal_enveloping_algebra()
                sage: e.nproduct(f,-1)
                e_-1f_-1|0>

            The product of the two fermions in the `N=2` vertex algebra
            has conformal weight `3`::

                sage: R = lie_conformal_algebras.N2(QQ)
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
            The ``n``-th product of these two elements.

            If `n\geq 0` it returns the element of this Lie conformal
            algebra. If `n < 0` then it first lifts this element to the
            universal  enveloping vertex algebra and returns the
            corresponding element

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                sage: L._nproduct_(L,3)
                1/2*C
                sage: L._nproduct_(L.T(),0)
                2*T^(2)L
                sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
                sage: E = V.0; H = V.1; F = V.2;
                sage: E._nproduct_(H,0) == - 2*E
                True
                sage: E._nproduct_(F,1)
                B['K']

            when `n<0`  we obtain an element of the universal enveloping
            vertex algebra, however, the universal enveloping vertex
            algebra needs to be constructed first::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
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

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)
                sage: Vir.inject_variables()
                Defining L, C
                sage: V = Vir.universal_enveloping_algebra({C:1})
                sage: L*L.T()
                L_-3L_-2|0> + L_-5|0>

            Note that the universal enveloping algebra needs to be
            constructed first::

                sage: W = lie_conformal_algebras.Weyl(QQ,4)
                sage: W.inject_variables()
                Defining alpha0, alpha1, alpha2, alpha3, K
                sage: alpha0*alpha1
                Traceback (most recent call last):
                NotImplementedError: in order to lift an element first need to construct the universal enveloping vertex algebra
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

                sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L.lift()
                L_-2|0>
                sage: L.lift().__class__
                <class 'sage.algebras.vertex_algebras.universal_enveloping_vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

            Notice that the target of the ``lift`` morphism changes when
            we construct another universal enveloping vertex algebra::

                sage: V = vertex_algebras.Virasoro(QQ,1);
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 1 over Rational Field
                sage: V = vertex_algebras.Virasoro(QQ,3)
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 3 over Rational Field

            Notice that recreation may not re-establish the right coercion
            depending on the method of construction::

                sage: cp = Family({Vir.1:1/2}); V = Vir.universal_enveloping_algebra(cp)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

                sage: V = vertex_algebras.Virasoro(QQ,1/2)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field

                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra of central charge 1/2 over Rational Field
            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def T(self, n=1):
            r"""
            The ``n``-th derivative of ``self``.

            INPUT:

            - ``n`` -- integer (default:``1``); how many times
              to apply `T` to this element

            OUTPUT:

            `T^n a` where `a` is this element. Notice that we use the
            *divided powers* notation `T^{(j)} = \frac{T^j}{j!}`.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)
                sage: Vir.inject_variables()
                Defining L, C
                sage: L.T()
                TL
                sage: L.T(3)
                6*T^(3)L
                sage: C.T()
                0
            """

    WithBasis = LazyImport("sage.categories.lambda_bracket_algebras_with_basis",
                           "LambdaBracketAlgebrasWithBasis", "WithBasis")

    FinitelyGeneratedAsLambdaBracketAlgebra = LazyImport(
        'sage.categories.finitely_generated_lambda_bracket_algebras',
        'FinitelyGeneratedLambdaBracketAlgebras')
