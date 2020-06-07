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

class VertexAlgebras(Category_over_base_ring):
    """
    The category of vertex algebras.

    EXAMPLES::

        sage: VertexAlgebras(QQ)
        Category of Vertex algebras over Rational Field
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
        sage: VertexAlgebras(QQbar).Super().WithBasis().Graded().FinitelyGenerated() is VertexAlgebras(QQbar).FinitelyGenerate
        ....: d().WithBasis().Graded().Super()
        True
    """
    def super_categories(self):
        """
        The super categories of this category.

        EXAMPLES::

            sage: C = VertexAlgebras(QQ)
            sage: C.super_categories()
            [Category of Lie conformal algebras over Rational Field]

            sage: VertexAlgebras(AA).Super().Graded().super_categories()
            [Category of super vertex algebras over Algebraic Real Field,
             Category of super H-graded Lie conformal algebras over Algebraic Real Field]
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
        pass

    class ParentMethods:
        @abstract_method(optional=True)
        def arc_algebra(self, termorder='wdegrevlex'):
            r"""
            The algebra of functions of the arc space of the `C_2`
            quotient of this Vertex algebra.

            INPUT:

            - ``termorder`` a string (default: ``'wdegrevlex'``); the
              monomial ordering of the algebra.

            OUTPUT: The graded Poisson vertex algebra freely generated
            as a differential algebra by the `C_2` quotient of this
            vertex algebra.

            TODO: we only support arc algebras of universal enveloping
            vertex algebras and their quotients.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2); Q=V.quotient(V.ideal(V.find_singular(6)[0]))
                sage: Q.arc_algebra()
                Quotient of The arc algebra over Rational Field generated by ('L',) by the differential ideal generated by (L_2^3,)
                sage: V.arc_space()
                The arc algebra over Rational Field generated by ('L',)
            """
            raise NotImplementedError("The arc algebra of {} is not " +
                                      "implemented yet".format(self))

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
                ideal of The Virasoro vertex algebra at central charge 1/2 generated by (L_-2L_-2L_-2|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>-27/16*L_-6|0>,)

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
                Quotient of The Virasoro vertex algebra at central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>-27/16*L_-6|0>,)
                sage: Q(L*(L*L))
                33/8*L_-4L_-2|0>-93/64*L_-3L_-3|0>+27/16*L_-6|0>
            """
            raise NotImplementedError("Quotients of {} are not implemented yet"\
                                        .format(self))

        @abstract_method(optional=True)
        def classical_limit(self):
            """
            The Poisson vertex algebra classical limit of this vertex
            algebra.

            EXAMPLES:

            We construct the classical limit of the universal Virasoro
            vertex algebra at central charge `1/2`::

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
            raise NotImplementedError("Classical limit of {} is not " +
                                     "implemented yet".format(self))

        def is_finitely_generated(self):
            """
            If this vertex algebra is finitely generated.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V.is_strongly_generated()
                True
                sage: W = AffineVertexAlgebra(QQ, 'A1', 1)
                sage: W.is_strongly_generated()
                True
            """
            return self in VertexAlgebras(self.base_ring()).FinitelyGenerated()

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
        def module(self):
            """
            The underlying module of this vertex algebra.

            EXAMPLES:

            For universal enveloping vertex algebras we get a
            :class:`CombinatorialFreeModule<sage.combinat.free_module.CombinatorialFreeModule>`::

                sage: W = AffineVertexAlgebra(QQ, 'A1', 2)
                sage: W.module()
                Free module generated by Partition tuples of level 3 over Rational Field

            For quotient algebras we get the algebra itself::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2)
                sage: V.register_lift()
                sage: L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: Q = V.quotient(V.ideal(v))
                sage: Q.module()
                Quotient of The Virasoro vertex algebra at central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>-27/16*L_-6|0>,)
                sage: Q.module() is Q
                True
            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def zero(self):
            """
            The zero vector in this vertex algebra.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2); L = V.0
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


        @abstract_method(optional=True)
        def find_singular(self,n):
            """
            Return the vector space of singular vectors of weight `n`.
            """
            raise NotImplementedError("Not implemented")

        def central_charge(self):
            """
            The central charge of this vertex algebra.
            """
            try:
                return self._c
            except AttributeError:
                raise NotImplementedError("The central charge of {} is not "+
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
        def monomials(self):
            """
            The tuple of monomials in this element.

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0; w = (L*L)*L;
                sage: w.monomials()
                (2*L_-3L_-3|0>, 4*L_-4L_-2|0>, 1/2*L_-6|0>, L_-2L_-2L_-2|0>)
            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def filtered_degree(self):
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

                sage: C = LieConformalAlgebras(ZZ).WithBasis().Graded(); C
                Category of H-graded Lie conformal algebras with basis over Integer Ring
                sage: D = LieConformalAlgebras(ZZ).Graded().WithBasis()
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

            EXAMPLES::

                sage: LieConformalAlgebras(AA).Super().WithBasis()
                Category of super Lie conformal algebras with basis over Algebraic Real Field

            TESTS::

                sage: C = LieConformalAlgebras(ZZ).Graded().WithBasis(); C
                Category of H-graded Lie conformal algebras with basis over Integer Ring
                sage: D = LieConformalAlgebras(ZZ).WithBasis().Graded(); D
                Category of H-graded Lie conformal algebras with basis over Integer Ring
                sage: C is D
                True
                sage: C = LieConformalAlgebras(ZZ).WithBasis().Graded().FinitelyGenerated(); C
                Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
                sage: D = LieConformalAlgebras(ZZ).FinitelyGenerated().WithBasis().Graded(); D
                Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
                sage: C is D
                True
                sage: C = LieConformalAlgebras(AA).Super().WithBasis(); C
                Category of super Lie conformal algebras with basis over Algebraic Real Field
                sage: D = LieConformalAlgebras(AA).WithBasis().Super(); D
                Category of super Lie conformal algebras with basis over Algebraic Real Field
                sage: C is D
                True
                sage: C = LieConformalAlgebras(AA).Graded().FinitelyGenerated().Super().WithBasis(); C
                Category of finitely generated super H-graded Lie conformal algebras with basis over Algebraic Real Field
                sage: D = LieConformalAlgebras(AA).WithBasis().FinitelyGenerated().Super().Graded(); D
                Category of finitely generated super H-graded Lie conformal algebras with basis over Algebraic Real Field
                sage: C is D
                True
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
        pass

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
                [Category of finitely generated super vertex algebras over Rational Field,
                 Category of super H-graded vertex algebras over Rational Field]
            """
            return [self.base_category(),]
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

        class ParentMethods:

            def is_super(self):
                """
                Wether this vertex algebra is a super vertex algebra.
                """
                return True

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
            pass

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

            class WithBasis(CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated H-graded vertex
                algebras with basis.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).Graded().FinitelyGenerated().WithBasis()
                    Category of finitely generated H-graded vertex algebras with basis over Rational Field
                """
                pass


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
                    ValueError: L_-2|0>+L_-3|0> is not homogeneous!
                """
                return self.weight()

            @abstract_method
            def weight(self):
                """
                The conformal weight of this element.

                This method is an alias of :meth:`degree`
                """
                raise NotImplementedError("Not implemented")

            @abstract_method
            def is_homogeneous(self):
                """
                Whether this element is homogeneous with respect to
                conformal weight.

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                    sage: L.is_homogeneous()
                    True
                    sage: (L + L.T()).is_homogeneous()
                    False
                """
                raise NotImplementedError("Not implemented")

            def nmodeproduct(self, other, n):
                r"""
                The shifted product of these two elements.

                INPUT:

                - ``other`` an element of this vertex algebra
                - ``n`` an integer number.

                OUTPUT: returns the shifted `n`-product of ``self``
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
                    ValueError: Couldn't compute weight of L_-2|0>+|0>, it's not homogeneous?

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
                    [L_-2|0>]
                    sage: Q.gens()[0] in Q
                    True
                    sage: V = AffineVertexAlgebra(QQ, 'A1', 1); V
                    The universal affine vertex algebra of CartanType ['A', 1] at level 1
                    sage: V.gens()
                    (E(alpha[1])_-1|0>, E(alphacheck[1])_-1|0>, E(-alpha[1])_-1|0>)
                """
                raise NotImplementedError("Not implemented")

            @abstract_method
            def ngens(self):
                """
                The number of generators of this vertex algebra

                EXAMPLES::

                    sage: VirasoroVertexAlgebra(QQ,1/2).ngens()
                    1
                    sage: AffineVertexAlgebra(QQ,'A2',1).ngens()
                    8
                """
                raise NotImplementedError("Not implemented")

            def hilbert_series(self,ord):
                """
                The graded dimension of this algebra

                INPUT:

                - ``ord`` -- integer; the precision order of the result.

                OUTPUT: The sum

                .. MATH::

                    \sum_{n = 0}^{ord} q^n \mathrm{dim} V_n

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2)
                    sage: V.hilbert_series(8)
                    1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 4*q^6 + 4*q^7 + 7*q^8
                    sage: L = V.0; v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: Q.hilbert_series(9)
                    1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 3*q^6 + 3*q^7 + 5*q^8 + 5*q^9
                """
                from sage.rings.power_series_ring import PowerSeriesRing
                q = PowerSeriesRing(self.base_ring(),'q',
                                    default_prec = ord+1).gen()
                return sum(self.dimension(n)*q**n for n in range(ord+1 ))


