r"""
Lie conformal algebras
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

from sage.categories.modules import Modules
from sage.categories.category_types import Category_over_base_ring
from sage.categories.category_with_axiom import \
                                     CategoryWithAxiom_over_base_ring 
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.element import coerce_binop
from sage.categories.morphism import Morphism
from sage.categories.category_with_axiom import all_axioms as all_axioms
from sage.misc.misc_c import prod
from sage.combinat.partition import Partition

all_axioms += ("FinitelyGeneratedAsLieConformalAlgebra", "HGraded")

class LieConformalAlgebras(Category_over_base_ring):
    """
    The category of Lie conformal algebras.

    EXAMPLES::

        sage: LieConformalAlgebras(QQ)
        Category of Lie conformal algebras over Rational Field
        sage: LieConformalAlgebras(QQ).super_categories()
        [Category of vector spaces over Rational Field]
        sage: R = PolynomialRing(QQ, 'x')
        sage: LieConformalAlgebras(R).super_categories()
        [Category of modules over Univariate Polynomial Ring in x over Rational Field]
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::
            sage: LieConformalAlgebras(QQ)
            Category of Lie conformal algebras over Rational Field
        """
        return [Modules(self.base_ring())]

    def _repr_object_names(self):
        return "Lie conformal algebras over {}".format(self.base_ring())

    class ParentMethods:
        def universal_enveloping_algebra(self, central_parameters=None):
            if hasattr(self, 'lift'):
                if self.lift.codomain().central_parameters() ==\
                        central_parameters:
                    return self.lift.codomain()
            self._construct_UEA(central_parameters)
            return self.lift.codomain()

        def _construct_UEA(self, central_parameters=None):
            """ 
            Returns the universal enveloping vertex algebra of ``self``. 

            Unlike :meth:`universal_enveloping_algebra`, this method does not 
            (usually) construct the canonical lift homomorphism from ``self`` to 
            the universal enveloping vertex algebra. 
            """
            from sage.algebras.vertex_algebras.vertex_algebra import \
                                                        VertexAlgebra
            return VertexAlgebra(self.base_ring(), self, 
                central_parameters=central_parameters)

        def ideal(self, *gens, **kwds):
            raise NotImplementedError("ideals are not implemented yet")

 
    class ElementMethods:
        @coerce_binop
        def bracket(self,rhs):
            """
            Returns the lambda bracket ``[self_\lambda rhs]``.
            """
            return self._bracket_(rhs)

        @abstract_method
        def _bracket_(self,rhs):
            """
            Returns the lambda bracket ```[self_\lambda rhs]``, where 
            rhs is an element of the same Lie conformal algebra of ``self``.
            """

        @coerce_binop
        def nproduct(self,rhs,n):
            return self._nproduct_(rhs,n)

        def _nproduct_(self,rhs,n):
            if n >= 0 :
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                return self.lift().nproduct(rhs,n)


        @abstract_method(optional=True)
        def lift(self):
            """
            Returns the image of ``self`` under the canonical lift to the
            universal enveloping vertex algebra of ``self.parent()``.
            """

    class SubcategoryMethods:
        def Graded(self):
            return self.HGraded()

        def HGraded(self):
            return self._with_axiom('HGraded')

        def FinitelyGeneratedAsLieConformalAlgebra(self):
            return self._with_axiom('FinitelyGeneratedAsLieConformalAlgebra')

        def FinitelyGenerated(self):
            return self.FinitelyGeneratedAsLieConformalAlgebra()

        def WithBasis(self):
            return self._with_axiom("WithBasis")

    class HGraded(CategoryWithAxiom_over_base_ring):
        r"""
        Category of H-graded Lie conformal algebras

        TESTS::

            sage: C = LieConformalAlgebras(QQ).Graded()
            sage: TestSuite(C).run()
        """

        def extra_super_categories(self):
            from sage.categories.graded_modules import GradedModules
            return (GradedModules(self.base_ring()),)

        def _repr_object_names(self):
            return "H-graded Lie conformal algebra over {}".format(
                self.base_ring())

        class SubcategoryMethods:
            def WithBasis(self):
                return self._with_axiom("WithBasis")

        class ParentMethods:
            def is_graded(self):
                return True

        class ElementMethods:
            @abstract_method
            def degree(self):
                pass

        class WithBasis(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                return "H-graded Lie conformal algebra with basis "\
                    "over {}".format(self.base_ring())

            class SubcategoryMethods:
                def FinitelyGenerated(self):
                    return self._with_axiom(
                                    "FinitelyGeneratedAsLieConformalAlgebra")

            class FinitelyGeneratedAsLieConformalAlgebra( 
                                        CategoryWithAxiom_over_base_ring):
                def _repr_object_names(self):
                    return "Finitely generated H-graded Lie conformal algebra"\
                                " with basis over {}".format(self.base_ring())

    class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            return "finitely generated Lie conformal algebra over {}"\
                    .format(self.base_ring())

        class ParentMethods:
            def ngens(self):
                return len(self.gens())

            def gen(self,i):
                return self.gens()[i]

    class WithBasis(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            return "Lie conformal algebras with basis over {}".format(
                    self.base_ring())

        class SubcategoryMethods:
            def FinitelyGenerated(self):
                return self.FinitelyGeneratedAsLieConformalAlgebra()

            def FinitelyGeneratedAsLieConformalAlgebra(self):
                return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

        class FinitelyGeneratedAsLieConformalAlgebra(CategoryWithAxiom_over_base_ring):
            def _repr_object_names(self):
                return "Finitely Generated Lie conformal algebras with basis over {}".format(
                    self.base_ring())

            class ParentMethods:
                @abstract_method
                def gens(self):
                    return

                @abstract_method(optional=True)
                def central_elements(self):
                    return tuple()

            class ElementMethods:
                #Do we want this in a category "With Generators"?
                def T(self,n=1 ):
                    if n==0 :
                        return self
                    coef = self.monomial_coefficients()
                    p = self.parent()
                    ret = p.zero()
                    for k in coef.keys():
                        if (k[0 ],k[1 ]+1 ) in p._indices:
                            ret += prod(range(k[1 ]+1 ,k[1 ]+n+1 ))*coef[k]*p.monomial((k[0 ],k[1 ]+n))
                    return ret

                def lift(self):
                    p = self.parent()
                    V = p.lift.codomain()
                    ret = V.zero()
                    for c in self.monomial_coefficients().items():
                        if p.monomial(c[0 ]) in p.central_elements():
                            ret += c[1 ]*V.central_parameters()[
                                        p.monomial(c[0 ])]*V.vacuum()
                        else:
                            l = [Partition([])]*V.ngens()
                            l[p._index_to_pos[c[0 ][0 ]]] = Partition(
                                                            [c[0 ][1 ]+1 ])
                            ret += c[1 ]*V(l)
                    return ret
                        
        class ElementMethods:

            def is_monomial(self):
                return ( len(self.monomial_coefficients()) == 1  or self.is_zero() )

            def index(self):
                if self.is_zero():
                    return tuple()
                if not self.is_monomial():
                    raise ValueError ("index can only be computed for monomials")
                return self.monomial_coefficients().keys()[0 ]


class LiftMorphism(Morphism):
    def _call_(self,x):
        return x.lift()    

