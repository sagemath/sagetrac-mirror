r"""
Vertex Algebras

AUTHORS:

- Reimundo Heluani (10-09-2019): Initial implementation.

.. include:: ../../../algebras/vertex_algebras/vertex_algebra_desc.inc

.. SEEALSO::

    :mod:`sage.algebras.vertex_algebras.vertex_algebra`
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
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.categories.quotients import QuotientsCategory
from sage.misc.abstract_method import abstract_method
from sage.structure.element import coerce_binop
from sage.functions.other import factorial
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.misc.cachefunc import cached_method
from sage.rings.all import QQ,ZZ
from sage.combinat.free_module import CombinatorialFreeModule

class VertexAlgebras(Category_over_base_ring):
    """
    The category of vertex algebras.

    EXAMPLES::

        sage: VertexAlgebras(QQ)
        Category of vertex algebras over Rational Field
        sage: VertexAlgebras(QQ).is_subcategory(LieConformalAlgebras(QQ))
        True

    TESTS::

        sage: VertexAlgebras(ZZ).super_categories()
        [Category of Lie conformal algebras over Integer Ring]
        sage: VertexAlgebras(ZZ).Super().WithBasis()
        Category of super vertex algebras with basis over Integer Ring
        sage: type(VertexAlgebras(ZZ).Super().WithBasis())
        <class 'sage.categories.vertex_algebras.VertexAlgebras.Super.WithBasis_with_category'>
        sage: type(VertexAlgebras(QQ).Super().FinitelyGenerated())
        <class 'sage.categories.vertex_algebras.VertexAlgebras.Super.FinitelyGeneratedAsVertexAlgebra_with_category'>
        sage: type(VertexAlgebras(QQ).Super().FinitelyGenerated().WithBasis())
        <class 'sage.categories.vertex_algebras.VertexAlgebras.Super.FinitelyGeneratedAsVertexAlgebra.WithBasis_with_category'>
        sage: type(VertexAlgebras(QQ).Graded().FinitelyGenerated().WithBasis())
        <class 'sage.categories.vertex_algebras.VertexAlgebras.Graded.FinitelyGeneratedAsVertexAlgebra.WithBasis_with_category'>
        sage: type(VertexAlgebras(QQ).Graded().WithBasis())
        <class 'sage.categories.vertex_algebras.VertexAlgebras.Graded.WithBasis_with_category'>
        sage: type(VertexAlgebras(QQ).Graded().WithBasis().Super())
        <class 'sage.categories.category.JoinCategory_with_category'>
        sage: VertexAlgebras(QQbar).Super().WithBasis() is VertexAlgebras(QQbar).WithBasis().Super()
        True
        sage: VertexAlgebras(QQbar).Super().WithBasis().Graded().FinitelyGenerated() is VertexAlgebras(QQbar).FinitelyGenerated().WithBasis().Graded().Super()
        True
    """
    @cached_method
    def super_categories(self):
        """
        The super categories of this category.

        EXAMPLES::

            sage: C = VertexAlgebras(QQ)
            sage: C.super_categories()
            [Category of Lie conformal algebras over Rational Field]

            sage: VertexAlgebras(AA).Super().Graded().super_categories()
            [Category of H-graded vertex algebras over Algebraic Real Field,
             Category of super H-graded Lie conformal algebras over Algebraic Real Field,
             Category of super vertex algebras over Algebraic Real Field]
        """
        return [LieConformalAlgebras(self.base_ring()),]

    def _repr_object_names(self):
        """
        The name of the objects of this category.

        EXAMPLES::

            sage: VertexAlgebras(QQbar)
            Category of vertex algebras over Algebraic Field
        """
        return "vertex algebras over {}".format(self.base_ring())

    class Quotients(QuotientsCategory):
        """
        The category of quotients of vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).Quotients()
            Category of quotients of vertex algebras over Algebraic Field
        """
        class ParentMethods:

            @abstract_method
            def defining_ideal(self):
                raise NotImplementedError()

            @abstract_method
            def cover_algebra(self):
                raise NotImplementedError()

            def _repr_(self):
                return "Quotient of {} by {}".format(self.cover_algebra(),
                                                     self.defining_ideal())

            def an_element(self):
                return self.retract(self.cover_algebra().an_element())

            @cached_method
            def zero(self):
                return self.retract(self.cover_algebra().zero())

            @cached_method
            def vacuum(self):
                return self.retract(self.cover_algebra().vacuum())

        class ElementMethods:

            @abstract_method
            def lift(self):
                raise NotImplementedError()

            def is_even_odd(self):
                return self.lift().is_even_odd()

    def example(self):
        """
        An example of parent in this category.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).example()
            The Virasoro vertex algebra of central charge 1 over Algebraic Field
        """
        from sage.algebras.vertex_algebras.virasoro_vertex_algebra import \
                                                        VirasoroVertexAlgebra
        return VirasoroVertexAlgebra(self.base_ring(),self.base_ring().one())

    class ParentMethods:
        @abstract_method(optional=True)
        def singular_support(self):
            r"""
            The algebra of functions of the singular support of this
            vertex algebra.

            The graded Poisson vertex algebra freely generated
            as a differential algebra by the `C_2` quotient of this
            vertex algebra.

            .. TODO::

                We only support arc algebras of universal enveloping
                vertex algebras and their quotients.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2); Q=V.quotient(V.ideal(V.find_singular(6)[0]))
                sage: Q.arc_algebra()
                Quotient of The arc algebra over Rational Field generated by ('L',) by the differential ideal generated by (L_2^3,)
                sage: V.arc_space()
                The arc algebra over Rational Field generated by ('L',)
            """
            raise NotImplementedError("singular_support is not implemented "\
                                      "for {}".format(self))


        @abstract_method(optional=True)
        def ideal(self, *gens):
            """
            The ideal of this vertex algebra generated by ``gens``.

            INPUT:

            - ``gens`` -- a list or tuple of elements of this vertex
               algebra.

            EXAMPLES:

            We construct the ideal defining the *Virasoro Ising*
            module::

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: I = V.ideal(v)
                sage: I
                ideal of The Virasoro vertex algebra of central charge 1/2 over Rational Field generated by (L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>,)
            """
            raise NotImplementedError("Ideals of {} are not implemented yet"\
                                        .format(self))

        @abstract_method(optional=True)
        def quotient(self, I):
            """
            The quotient of this vertex algebra by the ideal ``I``.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: I = V.ideal(v)
                sage: Q = V.quotient(I); Q
                Quotient of The Virasoro vertex algebra of central charge 1/2 over Rational Field by ideal of The Virasoro vertex algebra of central charge 1/2 over Rational Field generated by (L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>,)
                sage: Q(L*(L*L))
                -93/64*L_-3L_-3|0> + 33/8*L_-4L_-2|0> + 27/16*L_-6|0>
            """
            raise NotImplementedError("Quotients of {} are not implemented yet"\
                                        .format(self))

        def classical_limit(self):
            """
            The Poisson vertex algebra classical limit of this vertex
            algebra.

            EXAMPLES:

            We construct the classical limit of the universal Virasoro
            vertex algebra of central charge `1/2`::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2)
                sage: P = V.classical_limit()
                sage: V.inject_variables()
                Defining L
                sage: (L*L)*L == L*(L*L)
                False
                sage: (P(L)*P(L))*P(L) == P(L)*(P(L)*P(L))
                True
                sage: V(L).bracket(V(L))
                {0: L_-3|0>, 1: 2*L_-2|0>, 3: 1/4*|0>}
                sage: P(L).bracket(P(L))
                {}

            We construct the classical limit of the *Ising* model::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: Q = V.quotient(V.ideal(v)); P = Q.classical_limit()
                sage: L*(L*L)
                L_-2L_-2L_-2|0>
                sage: Q(L)*(Q(L)*Q(L))
                33/8*L_-4L_-2|0>-93/64*L_-3L_-3|0>+27/16*L_-6|0>
                sage: P(L)*(P(L)*P(L)) == P.zero()
                True
            """
            from sage.algebras.poisson_vertex_algebras.poisson_vertex_algebra \
                 import PoissonVertexAlgebra
            return PoissonVertexAlgebra(self.base_ring(), self)

        def is_finitely_generated(self):
            """
            If this vertex algebra is finitely generated.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V.is_finitely_generated()
                True
                sage: W = AffineVertexAlgebra(QQ, 'A1', 1)
                sage: W.is_finitely_generated()
                True
            """
            return self in VertexAlgebras(self.base_ring()).FinitelyGenerated()

        def is_super(self):
            """
            If this vertex algebra is a super vertex algebra.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,3)
                sage: V.is_super()
                False
                sage: NeveuSchwarzVertexAlgebra(QQ,1/2).is_super()
                True

            Note however that there are purely even vertex algebras::

                sage: R = AbelianLieConformalAlgebra(QQ,2,parity=(0,0))
                sage: V = R.universal_enveloping_algebra(); V.is_super()
                True
            """
            return self in VertexAlgebras(self.base_ring()).Super()

        def is_graded(self):
            """
            If this vertex algebra is H-Graded.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V.is_graded()
                True
                sage: W = AffineVertexAlgebra(QQ, 'A1', 1)
                sage: W.is_graded()
                True
            """
            return self in VertexAlgebras(self.base_ring()).Graded()

        @abstract_method
        def vacuum(self):
            """
            The vacuum vector of this vertex algebra.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V.vacuum()
                |0>
            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def zero(self):
            """
            The zero vector in this vertex algebra.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2);
                sage: V.register_lift(); L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: Q = V.quotient(V.ideal(v))
                sage: Q(0)
                0
                sage: V(0)
                0
                sage: V(0) == V.zero()
                True
                sage: Q(0) == Q.zero()
                True
            """
            raise NotImplementedError("Not Implemented")

        def central_charge(self):
            """
            The central charge of this vertex algebra.
            """
            try:
                return self._c
            except AttributeError:
                raise NotImplementedError("The central charge of {} is not "\
                                          "implemented".format(self))

    class ElementMethods:
        @coerce_binop
        def nproduct(self,rhs,n):
            """
            The ``n``-th product of these two elements.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                sage: L.nproduct(L,3)
                1/4*|0>
                sage: L.nproduct(L,-3)
                L_-4L_-2|0>

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement an n-th product the user
                needs to implement :meth:`_bracket_`.
            """
            return self._nproduct_(rhs,n)

        def _nproduct_(self,rhs,n):
            """
            The `n`-th product of these two elements.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                sage: L.nproduct(L,3)
                1/4*|0>
                sage: L.nproduct(L,-3)
                L_-4L_-2|0>
            """
            if n >= 0 :
                return self._bracket_(rhs).get(n,self.parent().zero())
            else:
                return self.T(-n-1)._mul_(rhs)/factorial(-1-n)

        @abstract_method
        def is_singular(self):
            """
            Return whether this vector is a singular vector.

            If `a \in V` is a vector in a finitely generated H-Graded
            vertex algebra, then `a` is singular if for each homogeneous
            vector `v in V` we have `v_n a = 0` whenever `n > 0`.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: v.is_singular()
                True
                sage: V = AffineVertexAlgebra(QQ, 'A1', 2); E = V.0
                sage: (E*E*E).is_singular()
                True
            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def li_filtration_degree(self):
            """
            The smallest space `F^p` in the Li filtration of this
            vertex algebra containing this element.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                sage: L.li_filtration_degree()
                0
                sage: (L.T(2)*L.T()).li_filtration_degree()
                3
            """
            raise NotImplementedError("Not implemented")

    class SubcategoryMethods:

        def FinitelyGeneratedAsVertexAlgebra(self):
            """
            The subcategory of finitely and strongly generated vertex
            algebras.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).FinitelyGenerated()
                Category of finitely generated vertex algebras over Algebraic Field
            """
            return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        def FinitelyGenerated(self):
            """
            The subcategory of finitely and strongly generated vertex
            algebras.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).FinitelyGenerated()
                Category of finitely generated vertex algebras over Algebraic Field
            """
            return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        def WithBasis(self):
            """
            The subcategory of vertex algebras with a preferred basis

            EXAMPLES::

                sage: VertexAlgebras(ZZ).WithBasis()
                Category of vertex algebras with basis over Integer Ring
            """
            return self._with_axiom("WithBasis")

        def Graded(self, base_ring=None):
            """
            The subcategory of H-Graded vertex algebras.

            EXAMPLES::

                sage: C = VertexAlgebras(ZZ).WithBasis().Graded(); C
                Category of H-graded vertex algebras with basis over Integer Ring
                sage: D = VertexAlgebras(ZZ).Graded().WithBasis()
                sage: D is C
                True
            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self,CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsVertexAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Graded()._with_axioms(axioms)
            return GradedModulesCategory.category_of(self)

        def Super(self, base_ring=None):
            """
            The subcategory of super vertex algebras.
            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self,CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsVertexAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Super()._with_axioms(axioms)
            return SuperModulesCategory.category_of(self)

    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The subcategory of vertex algebras with a preferred basis

        EXAMPLES::

            sage: VertexAlgebras(ZZ).WithBasis()
            Category of vertex algebras with basis over Integer Ring
        """
        class SubcategoryMethods:

            def FinitelyGenerated(self):
                """
                The subcategory of finitely and strongly generated vertex
                algebras with basis.

                EXAMPLES::

                    sage: VertexAlgebras(QQbar).FinitelyGenerated()
                    Category of finitely generated vertex algebras over Algebraic Field
                """
                return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        class ElementMethods:

            def is_monomial(self):
                """
                Whether this element is a monomial.

                .. WARNING::

                    The product on a vertex algebra is not associative
                    nor commutative, so this method relies on an
                    explicit PBW basis which is available for
                    example in universal enveloping vertex algebras.

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1); V.inject_variables()
                    Defining L
                    sage: L.T(2)*L
                    2*L_-4L_-2|0>
                    sage: (L.T(2)*L).is_monomial()
                    True
                    sage: L*L.T(2)
                    2*L_-4L_-2|0> + 4*L_-6|0>
                    sage: (L*L.T(2)).is_monomial()
                    False
                """
                return (len(self.monomial_coefficients()) == 1 or\
                        self.is_zero())

            def index(self):
                """
                The basis index parametrizing this monomial element.

                INPUT: 

                - ``self`` -- a monomial in this vertex algebra.

                EXAMPLES::

                    sage: V = AffineVertexAlgebra(QQ, 'A1', 1, names = ('e','h', 'f'));
                    sage: V.inject_variables()
                    Defining e, h, f
                    sage: v = e.T(3)*(e.T(2)*(e*(h.T()*(h*f)))); v
                    12*e_-4e_-3e_-1h_-2h_-1f_-1|0>
                    sage: v.index()
                    ([4, 3, 1], [2, 1], [1])

                    sage: v = f*e; v.index()
                    Traceback (most recent call last):
                    ...
                    ValueError: index can only be computed for monomials, got e_-1f_-1|0> - h_-2|0>
                """
                if self.is_zero():
                    return None
                if not self.is_monomial():
                    raise ValueError ("index can only be computed for monomials,"\
                                      " got {}".format(self))

                return next(iter(self.monomial_coefficients()))

    class Super(SuperModulesCategory):
        """
        The subcategory of super vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).Super()
            Category of super vertex algebras over Algebraic Field
            sage: VertexAlgebras(QQbar).Super().all_super_categories()
            [Category of super vertex algebras over Algebraic Field,
             Category of super Lie conformal algebras over Algebraic Field,
             Category of super modules over Algebraic Field,
             Category of graded modules over Algebraic Field,
             Category of vertex algebras over Algebraic Field,
             Category of filtered modules over Algebraic Field,
             Category of Lie conformal algebras over Algebraic Field,
             Category of vector spaces over Algebraic Field,
             Category of modules over Algebraic Field,
             Category of bimodules over Algebraic Field on the left and Algebraic Field on the right,
             Category of right modules over Algebraic Field,
             Category of left modules over Algebraic Field,
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
        #Need to do all this to make Super commute with Graded.
        def extra_super_categories(self):
            """
            The extra super categories of this category.

            EXAMPLES::

                sage: VertexAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                [Category of super H-graded vertex algebras over Rational Field,
                 Category of finitely generated H-graded vertex algebras over Rational Field,
                 Category of finitely generated super vertex algebras over Rational Field]
            """
            return [self.base_category(),]

        def example(self):
            """
            An example of a Super Vertex Algebra

            EXAMPLES::

                sage: VertexAlgebras(QQ).Super().example()
                The Neveu-Schwarz super vertex algebra of central charge 1 over Rational Field
            """
            from sage.algebras.vertex_algebras.neveu_schwarz_vertex_algebra\
                    import NeveuSchwarzVertexAlgebra
            return NeveuSchwarzVertexAlgebra(self.base_ring(),
                                             self.base_ring().one())

        class SubcategoryMethods:

            def Graded(self, base_ring=None):
                """
                The subcategory of H-graded super vertex algebras.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).Super().Graded()
                    Category of super H-graded vertex algebras over Rational Field
                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                        "FinitelyGeneratedAsVertexAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Graded()._with_axioms(axioms)

                return GradedModulesCategory.category_of(
                                         self.base_category()).Super()

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of super vertex algebras with basis.

            EXAMPLES::

                sage: VertexAlgebras(QQ).WithBasis().Super()
                Category of super vertex algebras with basis over Rational Field
                sage: _ is VertexAlgebras(QQ).Super().WithBasis()
                True
            """
            pass

        class FinitelyGeneratedAsVertexAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated super vertex algebras.

            EXAMPLES::

                sage: VertexAlgebras(GF(5)).FinitelyGenerated().Super()
                Category of finitely generated super vertex algebras over Finite Field of size 5

            TESTS::

                sage: VertexAlgebras(GF(5)).FinitelyGenerated().Super() is VertexAlgebras(GF(5)).Super().FinitelyGenerated()
                True
            """

            class WithBasis(CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated super vertex
                algebras with basis.

                EXAMPLES::

                    sage: VertexAlgebras(ZZ).WithBasis().FinitelyGeneratedAsVertexAlgebra()
                    Category of finitely generated vertex algebras with basis over Integer Ring
                """
                pass

    class Graded(GradedModulesCategory):
        """
        The subcategory of H-graded vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQ).Graded()
            Category of H-graded vertex algebras over Rational Field
        """
        def _repr_object_names(self):
            """
            The names of objects in this category

            EXAMPLES::

                sage: VertexAlgebras(QQ).Graded().FinitelyGenerated()
                Category of finitely generated H-graded vertex algebras over Rational Field
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class SubcategoryMethods:

            def Super(self, base_ring=None):
                """
                The subcategory of H-graded super Lie conformal algebras

                EXAMPLES::

                    sage: VertexAlgebras(CC).Super().Graded()
                    Category of super H-graded vertex algebras over Complex Field with 53 bits of precision
                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                        "FinitelyGeneratedAsVertexAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Super()._with_axioms(axioms)
                return SuperModulesCategory.category_of(self)

        class Super(SuperModulesCategory):
            """
            The subcategory of H-graded super Lie vertex algebras.

            EXAMPLES::

                sage: VertexAlgebras(AA).Graded().Super()
                Category of super H-graded vertex algebras over Algebraic Real Field
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                    [Category of super H-graded vertex algebras over Rational Field,
                     Category of finitely generated H-graded vertex algebras over Rational Field,
                     Category of finitely generated super vertex algebras over Rational Field]
                """
                return [self.base_category(),]

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of graded vertex algebras with basis.

            EXAMPLES::

                sage: VertexAlgebras(QQ).Graded().WithBasis()
                Category of H-graded vertex algebras with basis over Rational Field
            """
            pass

        class FinitelyGeneratedAsVertexAlgebra(
                                        CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated H-graded vertex
            algebras.

            EXAMPLES::

                sage: VertexAlgebras(QQ).Graded().FinitelyGenerated()
                Category of finitely generated H-graded vertex algebras over Rational Field
            """
            class ElementMethods:

                def _action_from_partition_tuple(self,pt,negative=True):
                    """
                    Helper function to apply elements of a vertex
                    algebra constructed from partitions.

                    INPUT:

                    - ``pt`` -- an index of the parent vertex algebra
                      as returned by 
                      :meth:`sage.structure.indexed_generators.indices`.

                    - ``negative`` -- boolean (default: `True`);

                    OUTPUT: the result of repeatedly applying
                    modes determined by ``pt`` of the generators of
                    this vertex algebra to the vector ``self``. By
                    default negative modes are applied. Thus if
                    `pt = [[1]]` and `L` is the unique generator of `V`,
                    the mode `L_{-1}` will be applied. If ``negative``
                    is `False`, non-negative modes are applied instead.
                    Thus in the example above `L_0` will be applied.

                    EXAMPLES::

                        sage: V = VirasoroVertexAlgebra(QQ,1/2); vac = V.vacuum()
                        sage: vac._action_from_partition_tuple([[3,2]])
                        L_-3L_-2|0>
                        sage: V = NeveuSchwarzVertexAlgebra(QQ,1); vac = V.vacuum()
                        sage: vac._action_from_partition_tuple([[2,2],[5/2,3/2]])
                        L_-2L_-2G_-5/2G_-3/2|0>
                        sage: V = AffineVertexAlgebra(QQ,'A1',2, names = ('e','h','f'));
                        sage: pt = V.indices().an_element(); pt.to_list()
                        [[1, 1, 1, 1], [2, 1, 1], [3, 1]]
                        sage: pt
                        ([1, 1, 1, 1], [2, 1, 1], [3, 1])
                        sage: pt = V.indices().an_element(); pt = pt.to_list()
                        sage: V.vacuum()._action_from_partition_tuple(pt)
                        e_-1e_-1e_-1e_-1h_-2h_-1h_-1f_-3f_-1|0>
                    """
                    p = self.parent()
                    if not isinstance(pt,(tuple,list)) or len(pt) != p.ngens():
                        raise ValueError("pt needs to be a list of length {}".\
                                         format(p.ngens()))

                    ret = self
                    pt = list(pt)
                    pt.reverse()
                    for j,mu in enumerate(pt):
                        mu.reverse()
                        g = p.gen(-j-1)
                        for n in mu:
                            if negative:
                                ret = g.nmodeproduct(ret,-n)
                            else:
                                ret = g.nmodeproduct(ret,n-1)
                    return ret

                def is_singular(self):
                    r"""
                    Return whether this vector is a singular vector.

                    If `a \in V` is a vector in a finitely generated
                    H-Graded vertex algebra, then `a` is singular if
                    for each homogeneous vector `v in V` we have
                    `v_n a = 0` whenever `n > 0`.

                    EXAMPLES::

                        sage: V = VirasoroVertexAlgebra(QQ,1/2); V.register_lift(); L = V.0
                        sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                        sage: v.is_singular()
                        True
                        sage: V = AffineVertexAlgebra(QQ, 'A1', 2); E = V.0;
                        sage: (E*E*E).is_singular()
                        True

                    """
                    p = self.parent()

                    if not self.is_homogeneous():
                        raise ValueError("{} is not homogeneous".format(self))

                    weight = self.weight()
                    for g in p.gens():
                        gw = g.weight()
                        br = g._bracket_(self)
                        if any(not br.get(n+gw-1, p.zero()).is_zero() for
                                   n in range(1,weight+2)):
                            return False
                    return True



            class WithBasis(CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated H-graded vertex
                algebras with basis.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).Graded().FinitelyGenerated().WithBasis()
                    Category of finitely generated H-graded vertex algebras with basis over Rational Field
                """
                class ElementMethods: #Graded.FinitelyGenerated.WithBasis

                    @abstract_method
                    def _li_filtration_monomial_degree(self):
                        raise NotImplementedError()

                    def li_filtration_lt(self):
                        """
                        The leading terms of this element with respect
                        to the Li filtration.

                        EXAMPLES::

                            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.inject_variables(); V.register_lift()
                            Defining L
                            sage: v = L*(L*L) - 33/16*L.T(2)*L + 93/64*L.T()*L.T() - 27/384*L.T(4); v
                            L_-2L_-2L_-2|0> - 33/8*L_-4L_-2|0> + 93/64*L_-3L_-3|0> - 27/16*L_-6|0>
                            sage: v.li_filtration_lt()
                            L_-2L_-2L_-2|0>
                            sage: v.li_filtration_degree()
                            0
                        """
                        if self.is_zero():
                            return self
                        lt = [(m._li_filtration_monomial_degree(),m) for m in self.terms()]
                        lideg = min(k for k,v in lt)
                        return sum(v for k,v in lt if k==lideg)

                    def li_filtration_degree(self):
                        r"""
                        The minimum `p` such that this element belongs
                        to the `p`-th filtered part of this vertex
                        algebra.

                        EXAMPLES::

                            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.inject_variables()
                            Defining L
                            sage: v = L*(L*L) - 33/16*L.T(2)*L + 93/64*L.T()*L.T() - 27/384*L.T(4); v
                            L_-2L_-2L_-2|0> - 33/8*L_-4L_-2|0> + 93/64*L_-3L_-3|0> - 27/16*L_-6|0>
                            sage: v.li_filtration_lt()
                            L_-2L_-2L_-2|0>
                            sage: v.li_filtration_degree()
                            0
                        """
                        if self.is_zero():
                            return Infinity
                        return min(m._li_filtration_monomial_degree() for m in self.monomials())

                class ParentMethods: #Graded.FinitelyGenerated.WithBasis

                    def get_weight_less_than(self,n):
                        """
                        The sub-vector space of this vertex algebra of
                        vectors with
                        conformal weight less than or equal to ``n``.

                        .. NOTE::

                            To implement a finitely generated H-graded
                            vertex algebra with basis. The basis indices
                            should implement a method `subset` which
                            admits at least the keyword parameter
                            ``energy`` such that
                            ``_indices.subset(energy=n)`` enumerates
                            the indices of the basis of monomials with
                            conformal weight ``n``. See for example
                            :meth:`EnergyPartitionTuples_all.subset`

                        INPUT:

                        - ``n`` -- a non-negative rational number;

                        OUTPUT: a submodule of this vertex algebra.

                        EXAMPLES::

                            sage: V = FreeFermionsVertexAlgebra(QQ); M = V.get_weight_less_than(5/2); M
                            Free module generated by {0, 1, 2, 3} over Rational Field
                            sage: [v.lift() for v in M.basis()]
                            [|0>, psi_-1/2|0>, psi_-3/2|0>, psi_-3/2psi_-1/2|0>]

                        TESTS::

                            sage: V = FreeFermionsVertexAlgebra(QQ); M = V.get_weight_less_than(5/2);
                            sage: M.reduce(V.vacuum())
                            0
                        """
                        if any(g.weight() not in QQ or g.weight == 0 for g in
                               self.gens()):

                            raise NotImplementedError("get_weight_less_than is"\
                                         " not implemented for {}".format(self))

                        if n not in QQ or n < 0:
                            raise ValueError("n needs to be a non-negative"\
                                             " rational number")

                        basis = []
                        for i in self._indices:
                            v = self(i)
                            if v.weight() < n:
                                basis.append(v)
                            else:
                                break
                        return self.submodule(basis)

                    def get_weight(self,n):
                        r"""
                        The sub-vector space of this vertex algebra of
                        vectors with conformal weight equal to ``n``.

                        INPUT:

                        - ``n`` -- a non-negative rational number;

                        OUTPUT: 

                        A submodule of this vertex algebra.

                        EXAMPLES::

                            sage: V = NeveuSchwarzVertexAlgebra(QQ, 1/2); M = V.get_weight(7/2)
                            sage: [v.lift() for v in M.basis()]
                            [G_-7/2|0>, L_-2G_-3/2|0>]
                        """
                        if any(g.weight() not in QQ or g.weight == 0 for g in\
                               self.gens()):
                            raise NotImplementedError("get_weight is not "\
                                            "implemented for {}".format(self))

                        if n not in QQ or n < 0:
                            raise ValueError("n needs to be a non-negative "\
                                             "rational number")

                        return self.submodule([self(v) for v in \
                                                self._indices.subset(energy=n)])

                    def dimension_at_weight(self,n):
                        """
                        The dimension of the space of conformal weight
                        ``n`` of this vertex algebra.

                        INPUT:

                        - ``n`` -- a non-negative rational number.

                        EXAMPLES::

                            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.dimension_at_weight(8)
                            7
                            sage: Q = V.quotient(V.ideal(V.find_singular(6))); Q.dimension_at_weight(8)
                            5
                            sage: NeveuSchwarzVertexAlgebra(QQ,3/2).dimension_at_weight(11/2)
                            5
                        """
                        if n not in QQ or n < 0:
                            raise ValueError("n must be a non-negative "\
                                             "rational number")
                        return self._indices.subset(energy=n).cardinality()

                    def find_singular(self,n):
                        """
                        Return a basis of the vector space of singular
                        vectors of weight `n`.

                        EXAMPLES:

                        We find the singular vector of the Virasoro
                        vertex algebra of central charge 1/2::

                            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.find_singular(6)
                            (L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>,)

                        Finding singular vectors work for quotients of
                        vertex algebras as well::

                            sage: V = AffineVertexAlgebra(QQ,'A1',1,names = ('e','h','f')); V.register_lift()
                            sage: sing = V.find_singular(2); sing
                            (f_-1f_-1|0>,
                             f_-2|0> + h_-1f_-1|0>,
                             h_-1h_-1|0> + h_-2|0> - 2*e_-1f_-1|0>,
                             e_-1h_-1|0> + e_-2|0>,
                             e_-1e_-1|0>)
                            sage: Q = V.quotient(V.ideal([sing[0]],check=False))
                            sage: Q.find_singular(2)
                            (f_-2|0> + h_-1f_-1|0>,
                             h_-1h_-1|0> + h_-2|0> - 2*e_-1f_-1|0>,
                             e_-1h_-1|0> + e_-2|0>,
                             e_-1e_-1|0>)

                        And for super vertex algebras::

                            sage: V = NeveuSchwarzVertexAlgebra(QQ, 7/10); V
                            The Neveu-Schwarz super vertex algebra of central charge 7/10 over Rational Field
                            sage: V.find_singular(4)
                            (G_-5/2G_-3/2|0> - 1/3*L_-2L_-2|0> - 1/10*L_-4|0>,)

                        TESTS::

                            sage: V = AffineVertexAlgebra(QQ,'A1',1); V.register_lift(); sing = V.find_singular(2)
                            sage: Q = V.quotient(V.ideal(sing))
                            sage: Q.find_singular(2)
                            ()
                            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.find_singular(0)
                            (|0>,)
                        """
                        if any(g.weight() not in QQ or g.weight == 0 for g in \
                               self.gens()):
                            raise NotImplementedError("find_singular is not "\
                                             "implemented for {}".format(self))

                        if n not in QQ or n < 0:
                            raise ValueError("n needs to be a non-negative "
                                             "rational number")
                        M = self.get_weight(n)
                        B = M.basis()
                        N = self.get_weight_less_than(n)
                        W = CombinatorialFreeModule(self.base_ring(),
                                                    self.gens())
                        from sage.categories.tensor import tensor
                        Z = tensor([N,W])
                        br = {(g,v): g.bracket(v.lift()) for v in B for g in \
                              self.gens()}
                        on_basis = lambda v: sum(tensor([N.retract(sum(\
                                br[(g,B[v])][k] for k in br[(g,B[v])] if \
                                k > g.weight()-1)),W(g)]) for g in self.gens())
                        f = M.module_morphism(on_basis, codomain=Z)
                        return tuple([v.lift() for v in f.kernel_basis()])


                class Quotients(QuotientsCategory):
                    """
                    The category of quotients of finitely generated
                    H-graded vertex algebras with basis.

                    EXAMPLES::

                        sage: VertexAlgebras(QQ).Graded().FinitelyGenerated().WithBasis().Quotients()
                        Category of quotients of finitely generated H-graded vertex algebras with basis over Rational Field
                    """
                    class ElementMethods:

                        def _li_filtration_monomial_degree(self):
                            return self.lift()._li_filtration_monomial_degree()

                        def weight(self):
                            return self.lift().weight()

                        def T(self, n=1):
                            return self.parent().retract(self.lift().T(n))

                        def _mul_(self, other):
                            return self.parent().retract(
                                                self.lift()._mul_(other.lift()))

                        def _bracket_(self, other):
                            p = self.parent()
                            sl = self.lift()
                            ol = other.lift()
                            return {k:p.retract(v) for k,v in sl._bracket_(ol).\
                                    items() if not p.retract(v).is_zero()}


        class ElementMethods:   #VertexAlgebras.Graded
            """
            Base class for elements of an H-graded vertex algebra
            """
            def degree(self):
                """
                The conformal weight of this element.

                .. NOTE::

                    This method is an alias of :meth:`weight`. Parents
                    of this category need to implement ``weight``.

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                    sage: L.degree()
                    2
                    sage: W = AffineVertexAlgebra(QQ, 'A1', 2); E = W.0
                    sage: E.degree()
                    1
                    sage: L.T().degree()
                    3
                    sage: (L + L.T()).degree()
                    Traceback (most recent call last):
                    ...
                    ValueError: L_-2|0> + L_-3|0> is not homogeneous!
                """
                return self.weight()

            @abstract_method
            def weight(self):
                """
                The conformal weight of this element.

                This method is an alias of :meth:`degree`
                """
                raise NotImplementedError("Not implemented")

            def homogeneous_terms(self):
                """
                Return a tuple with the homogeneous terms of this
                element.

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2)
                    sage: V.inject_variables()
                    Defining L
                    sage: (L+L.T(3) + L*L.T() + L.T(7)/factorial(7)).homogeneous_terms()
                    (L_-2|0>, 7*L_-5|0> + L_-3L_-2|0>, L_-9|0>)

                TESTS::

                    sage: V = VirasoroVertexAlgebra(QQ,1); V(0).homogeneous_terms()
                    (0,)
                """
                if self.is_zero():
                    return tuple([self])
                S = {}
                p = self.parent()
                for m in self.terms():
                    w = m.weight()
                    S[w] = S.get(w,p.zero()) + m
                return tuple(S.values())

            def is_homogeneous(self):
                """
                Whether this element is homogeneous with respect to conformal
                weight

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0;
                    sage: L.is_homogeneous()
                    True
                    sage: (L + L.T()).is_homogeneous()
                    False
                """
                return self.is_zero() or len(self.homogeneous_terms()) == 1


            def nmodeproduct(self, other, n):
                r"""
                The shifted product of these two elements.

                INPUT:

                - ``other`` an element of this vertex algebra
                - ``n`` an integer number.

                OUTPUT: 

                The shifted `n`-product of ``self``
                with ``other``, which is defined as follows.

                For an element `a` of degree `p`, that is `a \in V_p`,
                then `a_n b` is defined as `a_{(n+p-1)}b`.

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ, 1/2)
                    sage: V.inject_variables()
                    Defining L
                    sage: L.nmodeproduct(L.T(),0)
                    3*L_-3|0>

                    sage: (L + V.vacuum()).nmodeproduct(L,0)
                    Traceback (most recent call last):
                    ...
                    ValueError: Couldn't compute weight of L_-2|0> + |0>, it's not homogeneous?

                """
                try:
                    weight = self.weight()
                except ValueError:
                    raise ValueError("Couldn't compute weight of {}, "\
                                    "is it not homogeneous?".format(self))
                return self._nproduct_(other, n+weight-1)

    class FinitelyGeneratedAsVertexAlgebra(CategoryWithAxiom_over_base_ring):
        """
        The subcategory of finitely and strongly generated vertex
        algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQ).FinitelyGenerated()
            Category of finitely generated vertex algebras over Rational Field
        """

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The category of finitely generated vertex algebras with
            basis.

            EXAMPLES::

                sage: VertexAlgebras(QQ).FinitelyGenerated().WithBasis()
                Category of finitely generated vertex algebras with basis over Rational Field
            """
            pass

        class ParentMethods:
            @abstract_method
            def gens(self):
                """
                The list of generators of this vertex algebra

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2)
                    sage: V.gens()
                    (L_-2|0>,)
                    sage: L = V.0; L in V
                    True
                    sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: Q.gens()
                    (L_-2|0>,)
                    sage: Q.gens()[0] in Q
                    True
                    sage: V = AffineVertexAlgebra(QQ, 'A1', 1); V
                    The universal affine vertex algebra of CartanType ['A', 1] at level 1 over Rational Field
                    sage: V.gens()
                    (alpha[1]_-1|0>, alphacheck[1]_-1|0>, -alpha[1]_-1|0>)
                """
                raise NotImplementedError("Not implemented")

            def ngens(self):
                """
                The number of generators of this vertex algebra

                EXAMPLES::

                    sage: VirasoroVertexAlgebra(QQ,1/2).ngens()
                    1
                    sage: AffineVertexAlgebra(QQ,'A2',1).ngens()
                    8
                """
                return len(self.gens())

            def gen(self,i):
                return self.gens()[i]

            def hilbert_series(self,ord):
                r"""
                The graded dimension of this algebra

                INPUT:

                - ``ord`` -- a positive rational number; the precision
                  order of the result.

                OUTPUT:

                The sum

                .. MATH::

                    \sum_{n = 0}^{ord} q^n \mathrm{dim} V_n

                where `n` runs over all rationals such that `V_n \neq 0`

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2)
                    sage: V.hilbert_series(8)
                    1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 4*q^6 + 4*q^7 + 7*q^8
                    sage: L = V.0; v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: Q.hilbert_series(9)
                    1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 3*q^6 + 3*q^7 + 5*q^8 + 5*q^9

                    sage: V = FreeFermionsVertexAlgebra(QQ)
                    sage: V.hilbert_series(11/2)
                    1 + q^(1/2) + q^(3/2) + q^2 + q^(5/2) + q^3 + q^(7/2) + 2*q^4 + 2*q^(9/2) + 2*q^5 + O(q^(11/2))
                """
                from sage.arith.functions import lcm
                from sage.functions.other import floor
                weights = [g.weight() for g in self.gens()]
                if any([w not in QQ or w < 0 for w in weights]):
                    raise NotImplementedError("hilbert_series is not "\
                                              "implemented for {}".format(self))
                if ord not in QQ or ord < 0:
                    raise ValueError("ord must be a positive rational number")
                l = lcm([g.weight().denominator() for g in self.gens()])
                if l==1:
                    from sage.rings.power_series_ring import PowerSeriesRing
                    q = PowerSeriesRing(ZZ,'q', default_prec = ord).gen()
                    return sum(self.dimension_at_weight(n)*q**n for\
                               n in range(floor(ord))).O(floor(ord))
                else:
                    from sage.rings.puiseux_series_ring import PuiseuxSeriesRing
                    q = PuiseuxSeriesRing(ZZ,'q').gen()
                    ord = floor(ord*l)
                    f = sum(self.dimension_at_weight(n/l)*q**(n/l) for n in \
                            range(ord))
                    return f.add_bigoh(ord/l)


        class Quotients(QuotientsCategory):
            class ParentMethods:

                def gens(self):
                    return tuple(self.retract(g) for g in\
                                 self.cover_algebra().gens()\
                                 if not self.retract(g).is_zero())
