r"""
Poisson vertex algebras
AUTHORS

- Reimundo Heluani (08-09-2019): Initial implementation
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
from lie_conformal_algebras import LieConformalAlgebras
from sage.categories.commutative_algebras import CommutativeAlgebras
all_axioms += ("HGraded",)


class PoissonVertexAlgebras(Category_over_base_ring):

    @cached_method
    def super_categories(self):
        return [CommutativeAlgebras(self.base_ring()), 
            LieConformalAlgebras(self.base_ring())]

    def _repr_object_names(self):
        return "Poisson vertex algebras over {}".format(self.base_ring())

    class SubcategoryMethods:
        def Graded(self):
            return self.HGraded()

        def HGraded(self):
            return self._with_axiom('HGraded')

        def FinitelyGenerated(self):
            return self.FinitelyGeneratedAsMagma()

        def WithBasis(self):
            return self._with_axiom('WithBasis')

    class WithBasis(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            return "poisson vertex algebra with basis over {}".format(
                                                            self.base_ring())

        class ParentMethods:
            @abstract_method
            def basis(self):
                raise NotImplementedError("Not implemented")

    class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            return "finitely generated poisson vertex algebras over {}"\
                        .format(self.base_ring())

        class ParentMethods:
            @abstract_method
            def gens(self):
                return

            def ngens(self):
                return len(self.gens())
        
        class ElementMethods:
            def monomials(self):
                return tuple( v[1]*v[0] for v in self.monomial_coefficients().\
                            items() )

    class HGraded(CategoryWithAxiom_over_base_ring):
        def extra_super_categories(self):
            return [CommutativeAlgebras(self.base_ring()).Graded(),]

        def _repr_object_names(self):
            return "H-graded poisson vertex algebra over {}"\
                    .format(self.base_ring())

        class FinitelyGeneratedAsMagma(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                return "Finitely generated H-graded poisson vertex algebra "\
                    "over {}".format(self.base_ring())

        class ElementMethods:
            @abstract_method
            def weight(self):
                """
                weight
                """
 
        
    


