r"""
Lie Conformal Algebras

Let `R` be a commutative ring, a *super Lie conformal algebra*
[Kac1997]_ over `R`
(also known as a *vertex Lie algebra*) is an `R[T]` super module `L`
together with a `\mathbb{Z}/2\mathbb{Z}`-graded `R`-bilinear
operation (called the `\lambda`-bracket)
`L\otimes L \rightarrow L[\lambda]`
(polynomials in `\lambda` with
coefficients in `L`), `a \otimes b \mapsto [a_\lambda b]` satisfying

1. Sesquilinearity:

   .. MATH::

        [Ta_\lambda b] = - \lambda [a_\lambda b], \qquad [a_\lambda Tb] =
        (\lambda+ T) [a_\lambda b].

2. Skew-Symmetry:

   .. MATH::

        [a_\lambda b] = - (-1)^{p(a)p(b)} [b_{-\lambda - T} a],

   where `p(a)` is `0` if `a` is *even* and `1` if `a` is *odd*. The
   bracket in the RHS is computed as follows. First we evaluate
   `[b_\mu a]` with the formal
   parameter `\mu` to the *left*, then
   replace each appearance of the formal variable `\mu` by `-\lambda - T`.
   Finally apply `T` to the coefficients in `L`.

3. Jacobi identity:

   .. MATH::

       [a_\lambda [b_\mu c]] = [ [a_{\lambda + \mu} b]_\mu c] +
       (-1)^{p(a)p(b)} [b_\mu [a_\lambda c ]],

   which is understood as an equality in `L[\lambda,\mu]`.

   `T` is usually called the *translation operation* or the *derivative*.
   For an element `a \in L` we will say that `Ta` is the *derivative of*
   `a`. We define the *n-th products* `a_{(n)} b` for `a,b \in L` by

   .. MATH::

        [a_\lambda b] = \sum_{n \geq 0} \frac{\lambda^n}{n!} a_{(n)} b.

   A Lie conformal algebra is called *H-Graded* [DSK2006]_ if there exists
   a decomposition `L = \oplus L_n` such that the
   `\lambda`-bracket becomes graded of degree `-1`, that is:

   .. MATH::

        a_{(n)} b \in L_{p + q -n -1} \qquad
        a \in L_p, \: b \in L_q, \: n \geq 0.

   In particular this implies that the action of `T` increases
   degree by `1`.

.. NOTE::

    In the literature arbitrary gradings are allowed. In this
    implementation we only support non-negative rational gradings.


EXAMPLES:

1. The **Virasoro** Lie conformal algebra `Vir` over a ring `R`
   where `12` is invertible has two generators `L, C` as an `R[T]`-module.
   It is the direct sum of a free module of rank `1` generated by `L`, and
   a free rank one `R` module generated by `C` satisfying `TC = 0`.  `C`
   is central (the `\lambda`-bracket of `C` with any other vector
   vanishes). The remaining `\lambda`-bracket is given by

   .. MATH::

        [L_\lambda L] = T L + 2 \lambda L + \frac{\lambda^3}{12} C.

2. The **affine** or current Lie conformal algebra `L(\mathfrak{g})`
   associated to a finite dimensional Lie algebra `\mathfrak{g}` with
   non-degenerate, invariant `R`-bilinear form `(,)` is given as a central
   extension of the free
   `R[T]` module generated by `\mathfrak{g}` by a central element `K`. The
   `\lambda`-bracket of generators is given by

   .. MATH::

        [a_\lambda b] = [a,b] + \lambda (a,b) K, \qquad a,b \in \mathfrak{g}

3. The **Weyl** Lie conformal algebra, or `\beta-\gamma` system is
   given as the central extension of a free `R[T]` module with two
   generators `\beta` and `\gamma`, by a central element `K`.
   The only non-trivial brackets among generators are

   .. MATH::

        [\beta_\lambda \gamma] = - [\gamma_\lambda \beta] = K

4. The **Neveu-Schwarz** super Lie conformal algebra is a super Lie
   conformal algebra which is an extension of the Virasoro Lie conformal
   algebra. It consists of a Virasoro generator `L` as in example 1 above
   and an *odd* generator `G`. The remaining brackets are given by:

   .. MATH::

        [L_\lambda G] = \left( T + \frac{3}{2} \lambda \right) G \qquad
        [G_\lambda G] = 2 L + \frac{\lambda^2}{3} C

.. SEEALSO::

    - :mod:`sage.algebras.lie_conformal_algebras.lie_conformal_algebra`
    - :mod:`sage.algebras.lie_conformal_algebras.examples`

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

from sage.categories.modules import Modules
from .category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.element import coerce_binop
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.commutative_rings import CommutativeRings
from sage.misc.lazy_import import LazyImport

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

    def example(self):
        """
        An example of parent in this category.

        EXAMPLES::

            sage: LieConformalAlgebras(QQ).example()
            The Virasoro Lie conformal algebra over Rational Field
        """
        from sage.algebras.lie_conformal_algebras.virasoro_lie_conformal_algebra\
                                            import VirasoroLieConformalAlgebra
        return VirasoroLieConformalAlgebra(self.base_ring())

    def _repr_object_names(self):
        """
        The name of the objects of this category.

        EXAMPLES::

            sage: LieConformalAlgebras(QQ)
            Category of Lie conformal algebras over Rational Field
        """
        return "Lie conformal algebras over {}".format(self.base_ring())

    class ParentMethods:

        def ideal(self, *gens, **kwds):
            """
            The ideal of this conformal algebra generated by ``gens``.

            .. TODO::

                Ideals of Lie Conformal Algebras are not implemented
                yet.

            EXAMPLES::

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)
                sage: Vir.ideal()
                Traceback (most recent call last):
                ...
                NotImplementedError: ideals of Lie Conformal algebras are not implemented yet
            """
            raise NotImplementedError("ideals of Lie Conformal algebras are "\
                                      "not implemented yet")

        def is_super(self):
            """
            Wether this Lie conformal algebra is a super Lie
            conformal algebra.

            EXAMPLES::

                sage: V = lie_conformal_algebras.Virasoro(QQ)
                sage: V.is_super()
                False
                sage: lie_conformal_algebras.NeveuSchwarz(QQbar).is_super()
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

                sage: Vir = lie_conformal_algebras.Virasoro(QQ)
                sage: Vir
                The Virasoro Lie conformal algebra over Rational Field
                sage: Vir.is_graded()
                True
            """
            return self in LieConformalAlgebras(self.base_ring()).Graded()

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

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user
                needs to implement :meth:`_bracket_`
            """
            return self._bracket_(rhs)

        @abstract_method
        def _bracket_(self,rhs):
            r"""
            The `\lambda`-bracket of these two elements.

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

            .. NOTE::

                It is guaranteed that both are elements of the same
                parent.
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

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user
                needs to implement :meth:`_nproduct_`
            """
            return self._nproduct_(rhs,n)

        def _nproduct_(self,rhs,n):
            r"""
            The ``n``-th product of these two elements.

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

            .. NOTE::

                It is guaranteed that both are elements of the same
                parent.
            """
            if n >= 0:
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                raise NotImplementedError("vertex algebras are not implemented")

        @abstract_method
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

        def is_even_odd(self):
            """
            Return ``0`` if this element is *even* and ``1`` if it is
            *odd*.

            .. NOTE::

                This method returns ``0`` by default since every Lie
                conformal algebra can be thought as a purely even Lie
                conformal algebra. In order to
                implement a super Lie conformal algebra, the user
                needs to implement this method.

            EXAMPLES::

                sage: R = lie_conformal_algebras.NeveuSchwarz(QQ);
                sage: R.inject_variables()
                Defining L, G, C
                sage: G.is_even_odd()
                1
            """
            return 0

    class SubcategoryMethods:

        def FinitelyGeneratedAsLieConformalAlgebra(self):
            """
            The subcategory of finitely generated Lie conformal algebras.

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

    class Super(SuperModulesCategory):
        """
        The subcategory of super Lie conformal algebras.

        EXAMPLES::

            sage: LieConformalAlgebras(AA).Super()
            Category of super Lie conformal algebras over Algebraic Real Field
        """
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).Super().super_categories()
                [Category of super modules over Rational Field,
                 Category of Lie conformal algebras over Rational Field]
            """
            return [self.base_category(),]

        def example(self):
            """
            An example parent in this category.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).Super().example()
                The Neveu-Schwarz super Lie conformal algebra over Rational Field
            """
            from sage.algebras.lie_conformal_algebras.neveu_schwarz_lie_conformal_algebra\
                                          import NeveuSchwarzLieConformalAlgebra
            return NeveuSchwarzLieConformalAlgebra(self.base_ring())

        class SubcategoryMethods:
            def Graded(self, base_ring=None):
                """
                The subcategory of super H-graded Lie conformal algebras.

                EXAMPLES::

                    sage: LieConformalAlgebras(QQ).Super().Graded()
                    Category of super H-graded Lie conformal algebras over Rational Field
                """
                assert base_ring is None or base_ring is self.base_ring()
                return GradedModulesCategory.category_of(
                                         self.base_category()).Super()

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

                sage: LieConformalAlgebras(QQbar).Graded()
                Category of H-graded Lie conformal algebras over Algebraic Field
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class Super(SuperModulesCategory):
            """
            The subcategory of super H-graded Lie conformal algebras.

            EXAMPLES::

                sage: C = LieConformalAlgebras(QQbar)
                sage: C.Graded().Super()
                Category of super H-graded Lie conformal algebras over Algebraic Field
                sage: C.Graded().Super() is C.Super().Graded()
                True
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: C = LieConformalAlgebras(QQ).Graded().Super()
                    sage: C.super_categories()
                    [Category of super Lie conformal algebras over Rational Field,
                     Category of H-graded Lie conformal algebras over Rational Field]
                """
                return [self.base_category(),]

        class ElementMethods:
            @abstract_method
            def degree(self):
                """
                The degree of this element.

                EXAMPLES::

                    sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
                    sage: L.degree()
                    2
                    sage: L.T(3).degree()
                    5
                """

    WithBasis = LazyImport('sage.categories.lie_conformal_algebras_with_basis',
                           'LieConformalAlgebrasWithBasis')

    FinitelyGeneratedAsLieConformalAlgebra = LazyImport(
        'sage.categories.finitely_generated_lie_conformal_algebras',
        'FinitelyGeneratedAsLieConformalAlgebra')
