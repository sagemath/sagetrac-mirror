r"""
Poisson Vertex Algebras.

AUTHORS:

- Reimundo Heluani (10-09-2019): Initial implementation.

.. include:: ../../../algebras/vertex_algebras/poisson_vertex_algebra_desc.rst

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
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.graded_modules import GradedModulesCategory
from .lie_conformal_algebras import LieConformalAlgebras
from .super_modules import SuperModulesCategory
from sage.rings.all import QQ
from sage.categories.quotients import QuotientsCategory

class PoissonVertexAlgebras(Category_over_base_ring):

    @cached_method
    def super_categories(self):
        """
        The super categories of this category

        EXAMPLES::

            sage: C = PoissonVertexAlgebras(QQ); C
            Category of Poisson vertex algebras over Rational Field
            sage: C.super_categories()
            [Category of commutative algebras over Rational Field,
             Category of Lie conformal algebras over Rational Field]
            sage: C.Graded()
            Category of H-graded Poisson vertex algebras over Rational Field
            sage: C.Graded().super_categories()
            [Category of graded algebras over Rational Field,
             Category of Poisson vertex algebras over Rational Field,
             Category of H-graded Lie conformal algebras over Rational Field]

        """
        return [LieConformalAlgebras(self.base_ring()),
                CommutativeAlgebras(self.base_ring())]


    def _repr_object_names(self):
        """
        The names of objects of this category
        """
        return "Poisson vertex algebras over {}".format(self.base_ring())

    class ElementMethods:

        @abstract_method
        def lift(self):
            raise NotImplementedError()

        def is_even_odd(self):
            return self.lift().is_even_odd()

    class ParentMethods:
        @abstract_method
        def ideal(self, gens):
            raise NotImplementedError()

        def is_super(self):
            return self in PoissonVertexAlgebras(self.base_ring()).Super()

    class Quotients(QuotientsCategory):
        """
        The category of quotients of Poisson vertex algebras.

        EXAMPLES::
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
            def one(self):
                return self.retract(self.cover_algebra().one())

        class ElementMethods:

            @abstract_method
            def lift(self):
                raise NotImplementedError()

            def is_even_odd(self):
                return self.lift().is_even_odd()

    class SubcategoryMethods:
        def FinitelyGenerated(self):
            """
            The subcategory of finitely generated Poisson vertex algebras

            This returns ``self.FinitelyGeneratedAsPoissonVertexAlgebra()`` implicitly using
            only the algebra structure over `R[T]`.
            """
            return self._with_axiom('FinitelyGeneratedAsPoissonVertexAlgebra')

        def WithBasis(self):
            """
            The subcategory of Poisson vertex algebras with a preferred basis
            """
            return self._with_axiom('WithBasis')

        def Graded(self, base_ring=None):
            """
            The subcategory of H-Graded Poisson vertex algebras.

            EXAMPLES::
            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self,CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsPoissonVertexAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Graded()._with_axioms(axioms)
            return GradedModulesCategory.category_of(self)

        def Super(self, base_ring=None):
            """
            The subcategory of super Poisson vertex algebras.

            """
            assert base_ring is None or base_ring is self.base_ring()
            if isinstance(self,CategoryWithAxiom_over_base_ring):
                axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsPoissonVertexAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Super()._with_axioms(axioms)
            return SuperModulesCategory.category_of(self)


    class WithBasis(CategoryWithAxiom_over_base_ring):
        class SubcategoryMethods:
            def FinitelyGenerated(self):
                return self._with_axiom(
                                     "FinitelyGeneratedAsPoissonVertexAlgebra")

        class ParentMethods:
            @abstract_method
            def basis(self):
                """
                Return a basis of this Poisson vertex algebra
                """
                raise NotImplementedError("Not implemented")

    class FinitelyGeneratedAsPoissonVertexAlgebra(
                                            CategoryWithAxiom_over_base_ring):

        def _repr_object_names(self):
            return "finitely generated {}".format(self.base_category().\
                                                  _repr_object_names())

        class WithBasis(CategoryWithAxiom_over_base_ring):
            pass

        class ParentMethods:
            @abstract_method
            def gens(self):
                raise NotImplementedError("Not implemented")

            def ngens(self):
                return len(self.gens())

            def gen(self,i):
                return self.gens()[i]

        class Quotients(QuotientsCategory):
            class ParentMethods:
                def gens(self):
                    return tuple(self.retract(g) for g in\
                                 self.cover_algebra().gens()\
                                 if not self.retract(g).is_zero())

    class Super(SuperModulesCategory):
        def extra_super_categories(self):
            return [self.base_category(),]

        class SubcategoryMethods:

            def Graded(self, base_ring=None):
                """
                The subcategory of H-graded super Poisson vertex
                algebras.
                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsPoissonVertexAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Graded()._with_axioms(axioms)

                return GradedModulesCategory.category_of(
                                         self.base_category()).Super()

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of super Poisson vertex algebras with basis.
            """
            pass

        class FinitelyGeneratedAsPoissonVertexAlgebra(
                                            CategoryWithAxiom_over_base_ring):

            class WithBasis(CategoryWithAxiom_over_base_ring):
                pass

    class Graded(GradedModulesCategory):
        def _repr_object_names(self):
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class ElementMethods:

            @abstract_method
            def weight(self):
                raise NotImplementedError("Not Implemented")

        class SubcategoryMethods:

            def Super(self, base_ring=None):
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                    "FinitelyGeneratedAsPoissonVertexAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Super()._with_axioms(axioms)
                return SuperModulesCategory.category_of(self)

            def FinitelyGenerated(self):
                return self._with_axiom(
                                      'FinitelyGeneratedAsPoissonVertexAlgebra')


        class Super(SuperModulesCategory):
            def extra_super_categories(self):
                return [self.base_category(),]

        class WithBasis(CategoryWithAxiom_over_base_ring):

            class FinitelyGeneratedAsPoissonVertexAlgebra(
                                            CategoryWithAxiom_over_base_ring):

                class Quotients(QuotientsCategory):
                    """
                    The category of quotients of finitely generated
                    H-graded Poisson vertex algebras with basis.
                    """
                    class ElementMethods:

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

                    class ParentMethods:
                        @cached_method
                        def one(self):
                            #I need this here instead than PVA.Quotients()
                            #because otherwise mro picks the one from
                            #monoids I do not know why
                            return self.retract(self.cover_algebra().one())

                class ParentMethods: #Graded.WithBasis.FinitelyGenerated

                    def get_weight(self,n):
                        r"""
                        The sub-vector space of this poisson vertex
                        algebra of vectors with weight equal to ``n``.

                        INPUT:

                        - ``n`` -- a non-negative rational number;

                        OUTPUT: a submodule of this Poisson vertex
                        algebra.

                        EXAMPLES::
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

                    def hilbert_series(self,ord):
                        r"""
                        The graded dimension of this algebra

                        INPUT:

                        - ``ord`` -- a positive rational number; the precision
                        order of the result.

                        OUTPUT: The sum

                        .. MATH::

                            \sum_{n = 0}^{ord} q^n \mathrm{dim} P_n

                        where `n` runs over all rationals such that `P_n \neq 0`

                        EXAMPLES::

                            sage: V = NeveuSchwarzVertexAlgebra(QQ,1); P = V.classical_limit()
                            sage: P.hilbert_series(11/2)
                            1 + q^(3/2) + q^2 + q^(5/2) + q^3 + 2*q^(7/2) + 3*q^4 + 3*q^(9/2) + 3*q^5 + O(q^(11/2))
                            sage: V = VirasoroVertexAlgebra(QQ,1); P = V.singular_support()
                            sage: P.hilbert_series(10)
                            1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 4*q^6 + 4*q^7 + 7*q^8 + 8*q^9 + O(q^10)
                        """
                        from sage.arith.functions import lcm
                        from sage.functions.other import floor
                        from sage.rings.all import QQ, ZZ
                        weights = [g.weight() for g in self.gens()]
                        if any([w not in QQ or w < 0 for w in weights]):
                            raise NotImplementedError("hilbert_series is not "\
                                              "implemented for {}".format(self))
                        if ord not in QQ or ord < 0:
                            raise ValueError("ord must be a positive rational "\
                                             "number")
                        l = lcm([g.weight().denominator() for g in self.gens()])
                        if l==1:
                            from sage.rings.power_series_ring import\
                                                                PowerSeriesRing
                            q = PowerSeriesRing(ZZ,'q', default_prec=ord).gen()
                            return sum(self.dimension_at_weight(n)*q**n for\
                                       n in range(floor(ord))).O(floor(ord))
                        else:
                            from sage.rings.puiseux_series_ring import\
                                                              PuiseuxSeriesRing
                            q = PuiseuxSeriesRing(ZZ,'q').gen()
                            ord = floor(ord*l)
                            f = sum(self.dimension_at_weight(n/l)*q**(n/l) for\
                                    n in range(ord))
                            return f.add_bigoh(ord/l)


        class FinitelyGeneratedAsPoissonVertexAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            pass
