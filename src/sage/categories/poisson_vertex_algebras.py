r"""
Poisson vertex algebras
AUTHORS

- Reimundo Heluani (10-09-2019): Initial implementation

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
from sage.categories.category_with_axiom import all_axioms as all_axioms
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.categories.graded_modules import GradedModulesCategory

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
        return [CommutativeAlgebras(self.base_ring()), 
                LieConformalAlgebras(self.base_ring())]

    def _repr_object_names(self):
        """
        The names of objects of this category
        """
        return "Poisson vertex algebras over {}".format(self.base_ring())

    class SubcategoryMethods:
        def FinitelyGenerated(self):
            """
            The subcategory of finitely generated Poisson vertex algebras

            This returns ``self.FinitelyGeneratedAsMagma()`` implicitly using 
            only the algebra structure over `R[T]`.
            """
            return self.FinitelyGeneratedAsMagma()

        def WithBasis(self):
            """
            The subcategory of Poisson vertex algebras with a preferred basis
            """
            return self._with_axiom('WithBasis')

    class WithBasis(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            """
            The names of objects of this category
            """
            return "{} with basis".format(self.base_category().\
                                          _repr_object_names())

        class ParentMethods:
            @abstract_method
            def basis(self):
                """
                Return a basis of this Poisson vertex algebra
                """
                raise NotImplementedError("Not implemented")

        class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                """
                The names of objects of this category
                """
                return "finitely generated {}".format(self.base_category().\
                                                      _repr_object_names())


    class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            """
            The names of objects of this category
            """
            return "finitely generated {}".format(self.base_category().\
                                                  _repr_object_names())

        class ParentMethods:
            @abstract_method
            def gens(self):
                """
                The list of generators of this Poisson vertex algebra

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2);
                    sage: P = PoissonVertexAlgebra(QQ,V)
                    sage: P.gens()
                    (L_-2|0>,)
                    sage: P.0 in P
                    True
                    sage: P.inject_variables()
                    Defining L
                    sage: L*L*L
                    L_-2L_-2L_-2|0>
                    sage: L*L*L == L*(L*L)
                    True
                    sage: L in V
                    False

                    sage: L = V.0
                    sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: P = PoissonVertexAlgebra(QQ,Q)
                    sage: P.gens()
                    ({2: {0: B[0]}},)
                    sage: P.inject_variables()
                    Defining L
                    sage: L in P
                    True
                    sage: (L*L*L).is_zero()
                    True

                """
                raise NotImplementedError("Not implemented")

            def ngens(self):
                """
                The number of generators of this Poisson vertex algebra
                
                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2);
                    sage: P = PoissonVertexAlgebra(QQ,V)
                    sage: P.ngens()
                    1

                 """
                return len(self.gens())
        
        class ElementMethods:
            def monomials(self):
                """
                The monomials in this element

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0
                    sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: P = PoissonVertexAlgebra(QQ,Q)
                    sage: Q.inject_variables()
                    Defining L
                    sage: mon = P(L*L*L).monomials(); mon
                    [{6: {2: 65/8*B[0]}}, {6: {2: 35/64*B[1]}}, {6: {4: 35/16*B[0]}}]
                    sage: [m.lift() for m in mon]
                    [65/8*L_-4L_-2|0>, 35/64*L_-3L_-3|0>, 35/16*L_-6|0>]

                """
                return tuple(v[1]*v[0] for v in 
                             self.monomial_coefficients().items())

    class Graded(GradedModulesCategory):
        def _repr_object_names(self):
            """
            The names of objects of this category

            TESTS::
                sage: PoissonVertexAlgebras(QQ).WithBasis().Graded().__class__
                <class 'sage.categories.poisson_vertex_algebras.PoissonVertexAlgebras.Graded.WithBasis_with_category'>
                sage: PoissonVertexAlgebras(QQ).WithBasis().Graded() is PoissonVertexAlgebras(QQ).Graded().WithBasis()
                True

            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                """
                The names of objects of this category
                """
                return "finitely generated {}".format(self.base_category().\
                                                      _repr_object_names())

            class WithBasis(CategoryWithAxiom_over_base_ring):
                def _repr_object_names(self):
                    """
                    The names of objects of this category
                    """
                    return "{} with basis".format(self.base_category().\
                                                  _repr_object_names())


        class WithBasis(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                """
                The names of objects of this category
                """
                return "{} with basis".format(self.base_category().\
                                              _repr_object_names())


        class ElementMethods:
            @abstract_method
            def weight(self):
                """
                The conformal weight of this element

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2)
                    sage: P = PoissonVertexAlgebra(QQ,V)
                    sage: P.0.weight()
                    2
                    sage: L = V.0; v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: P = PoissonVertexAlgebra(QQ,Q)
                    sage: P.0.weight()
                    2

                """
                raise NotImplementedError("Not implemented")
