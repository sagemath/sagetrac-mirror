r"""
Lie conformal algebras

AUTHORS:

- Reimundo Heluani (2019-10-05): Initial implementation

.. include:: ../../../algebras/vertex_algebras/lie_conformal_algebra_desc.rst

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
from sage.categories.category_with_axiom import \
                                     CategoryWithAxiom_over_base_ring 
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.element import coerce_binop
from sage.categories.morphism import Morphism
from sage.categories.category_with_axiom import all_axioms as all_axioms
from sage.misc.misc_c import prod
from sage.combinat.partition import Partition
from sage.categories.graded_modules import GradedModulesCategory

all_axioms += ("FinitelyGeneratedAsLieConformalAlgebra",)

class LieConformalAlgebras(Category_over_base_ring):
    r"""
    The category of Lie conformal algebras. 

    EXAMPLES::

        sage: C = LieConformalAlgebras(QQ); C
        Category of Lie conformal algebras over Rational Field
        sage: C.is_subcategory(VectorSpaces(QQ))
        True
    """ 

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

    def _repr_object_names(self):
        """
        The name of the objects of this category

        EXAMPLES::

            sage: LieConformalAlgebras(QQ)
            Category of Lie conformal algebras over Rational Field
            sage: LieConformalAlgebras(QQ).Graded().FinitelyGenerated()
            Category of finitely generated as lie conformal algebra H-graded Lie
            conformal algebras over Rational Field
        """
        return "Lie conformal algebras over {}".format(self.base_ring())

    class ParentMethods:
        def universal_enveloping_algebra(self, 
                                        central_parameters=None, 
                                        names=None):
            r"""
            The universal enveloping vertex algebra of this Lie conformal
            algebra.

            INPUT:

            - ``central_parameters`` -- A family of constants in the base ring
              of this Lie conformal algebra parametrized by the central
              elements.

            - ``names`` -- The names of the generators of the universal
              enveloping vertex algebra.

            OUTPUT: 

            If `L` is a Lie conformal algebra over `R` with some 
            central elements
            `C_i \in L` indexed by a set `I`, ``central_parameters`` is a
            family of elements `c_i \in R` indexed by the same set, then this
            method returns the central quotient of the universal enveloping
            vertex algebra of `L` by the ideal generated by `C_i - c_i
            |0\rangle`. 

        
            EXAMPLES:

            We construct the universal enveloping vertex algebra of the Virasoro
            Lie conformal algebra at central charge `0` over the complex
            numbers::

                sage: Vir = VirasoroLieConformalAlgebra(CC)
                sage: V = Vir.universal_enveloping_algebra(); V
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Complex Field with 53 bits of precision.
                sage: V.0.bracket(V.0)
                {0: L_-3|0>, 1: 2.00000000000000*L_-2|0>}

            We construct the universal enveloping vertex algebra of the Virasoro
            Lie conformal algebra at central charge 2 over the rationals::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); 
                sage: Vir.gens()
                (L, C)
                sage: cp = Family({Vir.1:2}); V = Vir.universal_enveloping_algebra(cp)
                sage: V.0.nproduct(V.0,3)
                |0>

            The Virasoro Algebra is not defined over `\ZZ`::

                sage: Vir = VirasoroLieConformalAlgebra(ZZ)
                Traceback (most recent call last):
                ...
                ArithmeticError: inverse does not exist

            The list of central parameters needs to be indexed by the central
            elements of the algebra::

                sage: Vir = VirasoroLieConformalAlgebra(QQ);
                sage: cp = Family({Vir.0 : 3})
                sage: V = Vir.universal_enveloping_algebra(cp)
                Traceback (most recent call last):
                ...
                ValueError: central_parameters must be parametrized by central elements

            """  
            if central_parameters == None:
                from sage.sets.family import Family
                central_parameters = Family({ce:0 for ce in
                                            self.central_elements()})
            if hasattr(self, 'lift'):
                cp = self.lift.codomain().central_parameters()
                if cp == central_parameters:
                    return self.lift.codomain()
            V = self._construct_UEA(central_parameters, 
                                    names=names)

            from sage.categories.homset import Hom
            self.lift = LiftMorphism(Hom(self, V, category = 
                            LieConformalAlgebras(self.base_ring())))
            try: 
                self.lift.register_as_coercion()
            except AssertionError:
                #we already constructed this morphisms and its fine
                pass

            return V

        def _construct_UEA(self, central_parameters=None, names=None):
            """ 
            Returns the universal enveloping vertex algebra of ``self``. 

            see :meth:`universal_enveloping_algebra`
            """
            from sage.algebras.vertex_algebras.vertex_algebra import VertexAlgebra
            return VertexAlgebra(self.base_ring(), self, 
                                central_parameters=central_parameters, 
                                names=names)

        @abstract_method
        def ideal(self, *gens, **kwds):
            """ 
            The ideal of this conformal algebra generated by `gens`

            """
            raise NotImplementedError("Ideals are not implemented yet")

 
    class ElementMethods:
        @coerce_binop
        def bracket(self,rhs):
            r"""
            Returns the `\lambda`-bracket of these two elements.

            EXAMPLES:

            The brackets of the Virasoro Lie conformal Algebra::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L.bracket(L)
                {0: TL, 1: 2*L, 3: 1/2*C}
                sage: L.bracket(L.T())
                {0: 2*T^(2)L, 1: 3*TL, 2: 4*L, 4: 2*C}

            Now with a current algebra::

                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: V.gens()
                (alpha[1], alphacheck[1], -alpha[1], K)
                sage: E = V.0; H = V.1; F = V.2;
                sage: H.bracket(H)
                {1: 2*K}
                sage: E.bracket(F)
                {0: alphacheck[1], 1: K}


            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user needs
                to implement :meth:`_bracket_`
            """
            return self._bracket_(rhs)

        @abstract_method
        def _bracket_(self,rhs):
            r"""
            Returns the `\lambda`-bracket of these two elements.

            EXAMPLES:

            The brackets of the Virasoro Lie conformal Algebra::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L._bracket_(L)
                {0: TL, 1: 2*L, 3: 1/2*C}
                sage: L._bracket_(L.T())
                {0: 2*T^(2)L, 1: 3*TL, 2: 4*L, 4: 2*C}

            Now with a current algebra::

                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: V.gens()
                (alpha[1], alphacheck[1], -alpha[1], K)
                sage: E = V.0; H = V.1; F = V.2;
                sage: H._bracket_(H)
                {1: 2*K}
                sage: E._bracket_(F)
                {0: alphacheck[1], 1: K}

            .. NOTE::

                It is guaranteed that both are elements of the same parent.
        
            """

        @coerce_binop
        def nproduct(self,rhs,n):
            r"""
            Returns the n-th product of this two elements.

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L.nproduct(L,3)
                1/2*C
                sage: L.nproduct(L.T(),0)
                2*T^(2)L
                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: E = V.0; H = V.1; F = V.2;
                sage: E.nproduct(H,0) == - 2*E
                True
                sage: E.nproduct(F,1)
                K     
            
            when `n<0`  we obtain an element of the universal enveloping
            vertex algebra, however, the universal enveloping vertex algebra
            needs to be constructed first::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L.nproduct(L,-3)
                L_-4L_-2|0>

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user needs
                to implement :meth:`_nproduct_`

            """
            return self._nproduct_(rhs,n)

        def _nproduct_(self,rhs,n):
            r"""
            Returns the n-th product of this two elements.

            If `n\geq 0` it returns the element of this Lie conformal algebra.
            If `n < 0` then it first lifts this element to the universal 
            enveloping vertex algebra and returns the corresponding element

            EXAMPLES::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: L._nproduct_(L,3)
                1/2*C
                sage: L._nproduct_(L.T(),0)
                2*T^(2)L
                sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                sage: E = V.0; H = V.1; F = V.2;
                sage: E._nproduct_(H,0) == - 2*E
                True
                sage: E._nproduct_(F,1)
                K     

            when `n<0`  we obtain an element of the universal enveloping
            vertex algebra, however, the universal enveloping vertex algebra
            needs to be constructed first::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L._nproduct_(L,-3)
                L_-4L_-2|0>

            .. NOTE::

                It is guaranteed that both are elements of the same parent.
            """ 
            if n >= 0 :
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                return self.lift().nproduct(rhs,n)

        @abstract_method(optional=True)
        def lift(self):
            r"""
            Returns the image of this element under the canonical lift to the
            universal enveloping vertex algebra. 

            .. WARNING::

                The universal enveloping algebra needs to be constructed first for
                this morphism to be defined. 

                This morphism is registered as a coercion between this Lie conformal
                algebra and its universal enveloping vertex algebra upon creation.
                Since we consider central quotients of the universal enveloping
                vertex algebras by fixed central parameters, each time a different
                universal enveloping vertex algebra is constructed, this lift
                morphism is changed. See the examples below and also 
                :meth:`register_lift()<sage.algebras.vertex_algebras.vertex_algebra.UniversalEnvelopingVertexAlgebra.register_lift>`.


            EXAMPLES:

            We lift to the universal enveloping vertex algebra of the Virasoro
            Lie conformal algebra with central charge `0`::
               
                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L.lift()
                L_-2|0>
                sage: L.lift().__class__
                <class 'sage.algebras.vertex_algebras.vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                sage: L.lift().parent()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.

            Notice that the target of the ``lift`` morphism changes when we
            construct another universal enveloping vertex algebra::

                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra at central charge 1/2
                sage: V = VirasoroVertexAlgebra(QQ,3)
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra at central charge 3

            Notice that recreation may not re-establish the right coercion
            depending on the method of construction::

                sage: cp = Family({Vir.1:1/2}); V = Vir.universal_enveloping_algebra(cp)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra at central charge 1/2                

            """
            raise NotImplementedError("Not implemented")


    class SubcategoryMethods:

        def FinitelyGeneratedAsLieConformalAlgebra(self):
            """
            The subcategory of finitely generated Lie conformal algebras

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated()
                Category of finitely generated Lie conformal algebras over
                Rational Field

            """
            return self._with_axiom('FinitelyGeneratedAsLieConformalAlgebra')

        def FinitelyGenerated(self):
            """
            The subcategory of finitely generated Lie conformal algebras

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated()
                Category of finitely generated Lie conformal algebras over Rational Field

            """
            return self.FinitelyGeneratedAsLieConformalAlgebra()

        def WithBasis(self):
            """
            The subcategory of Lie conformal algebras with a preferred
            basis

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).WithBasis()
                Category of Lie conformal algebras over Rational Field with basis

            """
            return self._with_axiom("WithBasis")

    class Graded(GradedModulesCategory):

        def _repr_object_names(self):
            """
            The names of the objects of this category
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())

        class ParentMethods:
            def is_graded(self):
                """
                Wether this Lie conformal algebra is graded or not

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ)
                    sage: Vir
                    Lie conformal algebra on 2 generators (L, C) over Rational Field.
                    sage: Vir.is_graded()
                    True
                """
                return True

        class ElementMethods:
            @abstract_method
            def degree(self):
                """
                The degree of this element

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                    sage: L.degree()
                    2
                    sage: L.T(3).degree()
                    5

                """
                raise NotImplementedError("Not implemented")

        class WithBasis(CategoryWithAxiom_over_base_ring):
            def _repr_object_bames(self):
                return "pingo"


    class WithBasis(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            """
            The names of objects of this category
            """
            return "{} with basis".format(self.base_category().\
                                          _repr_object_names())

        class SubcategoryMethods:

            def Graded(self):
                axioms = self.axioms()
                return self.base_category().Graded()._with_axioms(axioms)
            
            def FinitelyGenerated(self):
                """
                The subcategory of finitely generated Lie conformal algebras
                with a preferred basis
                """
                return self.FinitelyGeneratedAsLieConformalAlgebra()

            def FinitelyGeneratedAsLieConformalAlgebra(self):
                """
                The subcategory of finitely generated Lie conformal algebras
                with a preferred basis
                """
                return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")


        class FinitelyGeneratedAsLieConformalAlgebra(CategoryWithAxiom_over_base_ring):

            def _repr_object_names(self):
                """
                The names of objects of this category
                """
                return "Finitely Generated {}".format(self.base_category().\
                                                      _repr_object_names())
            class ParentMethods:
                @abstract_method(optional=True)
                def central_elements(self):
                    r"""
                    The central generators of this Lie conformal algebra

                    EXAMPLES::

                        sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                        sage: V.gens()
                        (alpha[1], alphacheck[1], -alpha[1], K)
                        sage: V.central_elements()
                        (K,)
                        sage: Vir = VirasoroLieConformalAlgebra(QQ); Vir.central_elements()
                        (C,)

                    """
                    raise NotImplementedError("Not implemented")

            class ElementMethods:
                def T(self,n=1):
                    r"""
                    The n-th derivative of this element

                    INPUT:

                    - ``n`` -- integer (default:`1`); How many times to apply
                      `T` to this element. 

                    We use the notation `T^{(j)} = \frac{T^j}{j!}` 

                    EXAMPLES::

                        sage: Vir = VirasoroLieConformalAlgebra(QQ)
                        sage: L = Vir.0; C = Vir.1
                        sage: L.T()
                        TL
                        sage: L.T(3)
                        6*T^(3)L
                        sage: C.T()
                        0

                    """
                    if n == 0 :
                        return self
                    coef = self.value.monomial_coefficients()
                    p = self.parent()
                    ret = p.zero()
                    for k in coef.keys():
                        if (k[0],k[1]+1) in p._indices:
                            ret += prod(range(k[1]+1,k[1 ]+n+1))*coef[k]\
                                        *p.monomial((k[0],k[1]+n))
                    return ret

                def lift(self):
                    r"""
                    Returns the image of this element under the canonical lift to the
                    universal enveloping vertex algebra. 

                    .. WARNING::

                        The universal enveloping algebra needs to be constructed first for
                        this morphism to be defined. 

                        This morphism is registered as a coercion between this Lie conformal
                        algebra and its universal enveloping vertex algebra upon creation.
                        Since we consider central quotients of the universal enveloping
                        vertex algebras by fixed central parameters, each time a different
                        universal enveloping vertex algebra is constructed, this lift
                        morphism is changed. See the examples below and also 
                        :meth:`register_lift()<sage.algebras.vertex_algebras.vertex_algebra.UniversalEnvelopingVertexAlgebra.register_lift>`.


                    EXAMPLES:

                    We lift to the universal enveloping vertex algebra of the Virasoro
                    Lie conformal algebra with central charge `0`::
                       
                        sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                        sage: V = Vir.universal_enveloping_algebra()
                        sage: L.lift()
                        L_-2|0>
                        sage: L.lift().__class__
                        <class 'sage.algebras.vertex_algebras.vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                        sage: L.lift().parent()
                        The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.

                    Notice that the target of the ``lift`` morphism changes when we
                    construct another universal enveloping vertex algebra::

                        sage: Vir.lift.codomain()
                        The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                        sage: V = VirasoroVertexAlgebra(QQ,1/2);
                        sage: V.register_lift()
                        sage: Vir.lift.codomain()
                        The Virasoro vertex algebra at central charge 1/2

                    Notice that recreation may not re-establish the right coercion
                    depending on the method of construction::

                        sage: Vir = VirasoroLieConformalAlgebra(QQ)
                        sage: cp = Family({Vir.1:1/3}); V = Vir.universal_enveloping_algebra(cp)
                        sage: Vir.lift.codomain()
                        The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                        sage: V = VirasoroVertexAlgebra(QQ,1/2)
                        sage: Vir.lift.codomain()
                        The universal enveloping vertex algebra of Lie conformal algebra on 2 generators (L, C) over Rational Field.
                        sage: V.register_lift()
                        sage: Vir.lift.codomain()
                        The Virasoro vertex algebra at central charge 1/2                

                    """
                    p = self.parent()
                    if not hasattr(p, 'lift'):
                        raise NotImplementedError(
                            "In order to lift an element first need to "\
                            "construct the universal enveloping vertex "\
                            "algebra")
                    V = p.lift.codomain()
                    ret = V.zero()
                    for c in self.value.monomial_coefficients().items():
                        if p.monomial(c[0]) in p.central_elements():
                            ret += c[1]*V.central_parameters()[
                                        p.monomial(c[0])]*V.vacuum()
                        else:
                            l = [Partition([])]*V.ngens()
                            l[p._index_to_pos[c[0][0]]] = Partition(
                                                            [c[0][1]+1])
                            ret += c[1]*V(l)
                    return ret
                        
        class ElementMethods:

            def is_monomial(self):
                r"""
                Whether this element is a monomial

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                    sage: (L + L.T()).is_monomial()
                    False
                    sage: L.T().is_monomial()
                    True

                """
                return (len(self.value.monomial_coefficients()) == 1  or self.is_zero())

            def index(self):
                r"""
                The index parametrizing this monomial element

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                    sage: L.index()
                    ('L', 0)
                    sage: L.T(4).index()
                    ('L', 4)
                    sage: (L + L.T()).index()
                    Traceback (most recent call last):
                    ...
                    ValueError: index can only be computed for monomials

                """
                if self.is_zero():
                    return tuple()
                if not self.is_monomial():
                    raise ValueError ("index can only be computed for monomials")
                return self.value.monomial_coefficients().keys()[0]

    class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            """
            The names of objects of this category
            """
            return "finitely generated {}".format(self.base_category().\
                            _repr_object_names())

        class ParentMethods:
            def ngens(self):
                r"""
                The number of generators of this Lie conformal algebra
                
                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ)
                    sage: Vir.ngens()
                    2
                    sage: V = AffineLieConformalAlgebra(QQ, 'A2')
                    sage: V.ngens()
                    9

                """
                return len(self.gens())

            def gen(self,i):
                r"""
                The generator of this Lie conformal algebra indexed by ``i``

                EXAMPLES::

                    sage: V = AffineLieConformalAlgebra(QQ, 'A1')
                    sage: V.gens()
                    (alpha[1], alphacheck[1], -alpha[1], K)
                    sage: V.gen(0)
                    alpha[1]
                    sage: V.1
                    alphacheck[1]

                """
                return self.gens()[i]


class LiftMorphism(Morphism):
    """
    The lift morphism from this Lie conformal algebra to its universal
    enveloping vertex algebra
    """
    def _call_(self,x):
        return x.lift()    

