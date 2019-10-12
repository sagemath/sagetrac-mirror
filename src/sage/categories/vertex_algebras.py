r"""
Vertex algebras
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

from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import all_axioms as all_axioms
from sage.categories.quotients import QuotientsCategory
from lie_conformal_algebras import LieConformalAlgebras
from sage.algebras.vertex_algebras.vertex_algebra_quotient import VertexAlgebraQuotient
from sage.functions.other import factorial

all_axioms += ("FinitelyGeneratedAsVertexAlgebra","HGraded")
class VertexAlgebras(Category_over_base_ring):
    """
    The category of vertex algebras.
    """

    @cached_method
    def super_categories(self):
        return [LieConformalAlgebras(self.base_ring()),]

    def _repr_object_names(self):
        return "Vertex algebras over {}".format(self.base_ring())

    class Quotients(QuotientsCategory):
        pass

    class ParentMethods:
        def ideal(self, *gens, **kwds):
            from sage.algebras.vertex_algebras.vertex_algebra_ideal import VertexAlgebraIdeal
            return VertexAlgebraIdeal(self,gens)

        def quotient(self, I):
            return VertexAlgebraQuotient(I)

        def classical_limit(self):
            raise NotImplementedError("General classical limit is only"+
                " implemented for H-graded vertex agebras")

        def _element_constructor_(self,x):
            if x in self.base_ring():
                if x != 0 :
                    raise ValueError("can only convert the scalar 0 into a vertex algebra element")
                return self.zero()
            return self.element_class(self,x)

        def is_strongly_generated(self):
            return self in VertexAlgebras(self.base_ring()).FinitelyGenerated()

        def is_graded(self):
            return self in VertexAlgebras(self.base_ring()).HGraded()

        @abstract_method
        def vacuum(self):
            return

        @abstract_method
        def module(self):
            return

        @abstract_method
        def zero(self):
            raise NotImplementedError("Not Implemented")

    class ElementMethods:
        def _nproduct_(self,rhs,n):
            if n >= 0 :
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                return factorial(-1 -n)**(-1 )*self.T(-n-1 )._mul_(rhs)


    class SubcategoryMethods:
        def Graded(self):
            return self.HGraded()

        def HGraded(self):
            return self._with_axiom('HGraded')

        def FinitelyGeneratedAsVertexAlgebra(self):
            return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        def FinitelyGenerated(self):
            return self.FinitelyGeneratedAsVertexAlgebra()

        def WithBasis(self):
            return self._with_axiom("WithBasis")

    class WithBasis(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            return "vertex algebras with basis over {}".format(self.base_ring())

        class ElementMethods:
            def monomials(self):
                return tuple( v[1 ]*v[0 ] for v in
                                            self.monomial_coefficients().items() )
 

    class HGraded(CategoryWithAxiom_over_base_ring):

        def _repr_object_names(self):
            return "H-graded vertex algebras over {}".format(self.base_ring())

        class SubcategoryMethods:
            def FinitelyGenerated(self):
                return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        class WithBasis(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                return "H-graded vertex algebras with basis "\
                            "over {}".format(self.base_ring())


        class ParentMethods:
            def classical_limit(self):
                """
                Returns the Poisson vertex algebra which is the classical
                limit of ``self``.
                """
                from sage.algebras.vertex_algebras.poisson_vertex_algebra import PoissonVertexAlgebra
                return PoissonVertexAlgebra(self.base_ring(), self)

        class ElementMethods:
            #for compatibility with LCA as `self` is in LCA
            def degree(self):
                return self.weight()

            def filtered_degree(self):
                return max(m.weight() for m in self.monomial_coefficients())

            @abstract_method
            def weight(self):
                """Conformal weight"""

            def is_homogeneous(self):
                try:
                    self.weight()
                except ValueError:
                    return False
                return True

            def nmodeproduct(self,other, n):
                try:
                    weight = self.weight()
                except ValueError:
                    raise ValueError( "Couldn't compute weight of {}, "\
                        "it's not homogeneous?".format(self) )
                return self.nproduct(other, n + weight - 1  )

    class FinitelyGeneratedAsVertexAlgebra(CategoryWithAxiom_over_base_ring):

        def _repr_object_names(self):
            return "finitely and strongly generated" \
                        " vertex algebras over {}".format(self.base_ring())
        
        class ParentMethods:
            @abstract_method
            def gens(self):
                return

            def ngens(self):
                return len(self.gens())

            @abstract_method
            def central_parameters(self):
                return
        
            def hilbert_series(self,ord):
                from sage.rings.power_series_ring import PowerSeriesRing
                q = PowerSeriesRing(self.base_ring(), 'q', 
                                                default_prec = ord).gen()
                return sum(self.dimension(n)*q**n for n in range(ord+1 ))     

       
        class SubcategoryMethods:
            def HGraded(self):
                return self._with_axiom('HGraded')

        class HGraded(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                return "H-graded finitely and strongly generated vertex"\
                    " algebras over {}".format(self.base_ring())

            class SubcategoryMethods:
                def WithBasis(self):
                    return self._with_axiom("WithBasis")

            class WithBasis(CategoryWithAxiom_over_base_ring):
                def _repr_object_names(self):
                    return "H-graded finitely and strongly generated vertex"\
                        " algebras with basis over {}".format(self.base_ring())

            class ElementMethods:
                def is_singular(self):
                    p = self.parent()
                    try:
                        weight = self.weight()
                    except ValueError:
                        raise ValueError( "Couldn't compute weight of {}, "\
                            "it's not homogeneous?".format(self) )

                    return all ( p(g).nmodeproduct(self,n).is_zero() for 
                        n in range(1 ,weight+2 ) for g in p.gens() )

                def _action_from_partition_tuple(self,p,negative=True):
                    """
                    helper function. From a partition tuple `p` applies the
                    corresponding basis element from `V` to self.
                    """
                    ngens = self.parent().ngens()
                    if len(p) != ngens:
                        raise ValueError("p has to be a partition tuple of "
                            "level {0}, got {1}".format(ngens,p))

                    ret = self
                    p = p.to_list()
                    p.reverse()
                    for j in range(ngens):
                        p[j].reverse()
                        g = self.parent()(self.parent().gen(ngens-j-1 ))
                        for n in p[j]:
                            if negative:
                                ret = g.nmodeproduct(ret,-n)
                            else:
                                ret = g.nmodeproduct(ret, n-1 )
                    return ret


