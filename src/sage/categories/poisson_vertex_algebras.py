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
            return self._with_axiom('FinitelyGeneratedAsMagma')

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
                                    "FinitelyGeneratedAsMagma"])
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
                                    "FinitelyGeneratedAsVertexAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Super()._with_axioms(axioms)
            return SuperModulesCategory.category_of(self)


    class WithBasis(CategoryWithAxiom_over_base_ring):
        class SubcategoryMethods:
            def FinitelyGenerated(self):
                return self._with_axiom("FinitelyGeneratedAsMagma")

        class ParentMethods:
            @abstract_method
            def basis(self):
                """
                Return a basis of this Poisson vertex algebra
                """
                raise NotImplementedError("Not implemented")

        class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
            pass

    class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
        
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
                                        "FinitelyGeneratedAsMagma"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Graded()._with_axioms(axioms)

                return GradedModulesCategory.category_of(
                                         self.base_category()).Super()

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of super Poisson vertex algebras with basis.
            """
            pass

        class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
            class WithBasis(CategoryWithAxiom_over_base_ring):
                pass

    class Graded(GradedModulesCategory):
        def _repr_object_names(self):
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class SubcategoryMethods:

            def Super(self, base_ring=None):
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",
                                        "FinitelyGeneratedAsMagma"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Super()._with_axioms(axioms)
                return SuperModulesCategory.category_of(self)

        class Super(SuperModulesCategory):
            def extra_super_categories(self):
                return [self.base_category(),]

        class WithBasis(CategoryWithAxiom_over_base_ring):
            pass

        class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
            class WithBasis(CategoryWithAxiom_over_base_ring):
                pass
 
