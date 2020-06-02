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
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.commutative_rings import CommutativeRings
_CommutativeRings = CommutativeRings()

class LieConformalAlgebras(Category_over_base_ring):
    r"""
    The category of Lie conformal algebras. 

    EXAMPLES::

        sage: C = LieConformalAlgebras(QQ); C
        Category of Lie conformal algebras over Rational Field
        sage: C.is_subcategory(VectorSpaces(QQ))
        True

    The base ring needs to be a commutative ring::

        sage: LieConformalAlgebras(QuaternionAlgebra(2))
        Traceback (most recent call last):
        ValueError: base must be a commutative ring got Quaternion Algebra (-1, -1) with base ring Rational Field

    """ 
        
    @staticmethod
    def __classcall_private__(cls,R,check=True):
        """
        INPUT:

        - `R` -- a commutative ring.
        - ``check`` -- a boolean (default: True) wether to check that `R` is 
          a commutative ring

        EXAMPLES::

            sage: LieConformalAlgebras(QuaternionAlgebra(2))
            Traceback (most recent call last):
            ValueError: base must be a commutative ring got Quaternion Algebra (-1, -1) with base ring Rational Field
            sage: LieConformalAlgebras(ZZ)
            Category of Lie conformal algebras over Integer Ring

        """
        if check:
            if not (R in _CommutativeRings): 
                    raise ValueError("base must be a commutative ring" + 
                                     " got {}".format(R))
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

    def _repr_object_names(self):
        """
        The name of the objects of this category.

        EXAMPLES::

            sage: LieConformalAlgebras(QQ)
            Category of Lie conformal algebras over Rational Field
            sage: LieConformalAlgebras(QQ).Graded().FinitelyGenerated()
            Category of finitely generated H-graded Lie conformal algebras over Rational Field
        """
        return "Lie conformal algebras over {}".format(self.base_ring())

    class ParentMethods:
        def universal_enveloping_algebra(self, 
                                        central_parameters=None, 
                                        names=None):
            r"""
            The universal enveloping vertex algebra of this Lie 
            conformal algebra.

            INPUT:

            - ``central_parameters`` -- A family of constants in the 
              base ring of this Lie conformal algebra parametrized by 
              the central elements.

            - ``names`` -- The names of the generators of the universal
              enveloping vertex algebra.

            OUTPUT: 

            If `L` is a Lie conformal algebra over `R` with some 
            central elements `C_i \in L` indexed by a set `I`, 
            ``central_parameters`` is a family of elements `c_i \in R`
            indexed by the same set, then this method returns the 
            central quotient of the universal enveloping vertex algebra
            of `L` by the ideal generated by `C_i - c_i |0\rangle`. 
        
            EXAMPLES:

            We construct the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra at central charge `0` over 
            the complex numbers::

                sage: Vir = VirasoroLieConformalAlgebra(CC)
                sage: V = Vir.universal_enveloping_algebra(); V
                The universal enveloping vertex algebra of The Virasoro Lie conformal algebra over Complex Field with 53 bits of precision
                sage: V.0.bracket(V.0)
                {0: L_-3|0>, 1: 2.00000000000000*L_-2|0>}

            We construct the universal enveloping vertex algebra of the
            Virasoro Lie conformal algebra at central charge 2 over the
            rationals::

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

            The list of central parameters needs to be indexed by the
            central elements of the algebra::

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
            self.lift = _LiftMorphism(Hom(self, V, category = 
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
            from sage.algebras.vertex_algebras.vertex_algebra\
                                                        import VertexAlgebra
            return VertexAlgebra(self.base_ring(), self, 
                                central_parameters=central_parameters, 
                                names=names)

        @abstract_method
        def ideal(self, *gens, **kwds):
            """ 
            The ideal of this conformal algebra generated by `gens`.

            .. TODO::

                Not implemented
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
                in order to implement a Lie conformal algebra the user 
                needs to implement :meth:`_bracket_`
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

                It is guaranteed that both are elements of the same 
                parent.
            """
            raise NotImplementedError("Not implemented")

        @coerce_binop
        def nproduct(self,rhs,n):
            r"""
            Returns the n-th product of these two elements.

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
            vertex algebra, however, the universal enveloping vertex
            algebra needs to be constructed first::

                sage: R = AffineLieConformalAlgebra(AA,'A1'); e = R.0; f = R.2;
                sage: e.nproduct(f,1)
                K
                sage: e.nproduct(f,-1)
                Traceback (most recent call last):
                ...
                NotImplementedError: In order to lift an element first need to construct the universal enveloping vertex algebra
                sage: V = R.universal_enveloping_algebra()
                sage: e.nproduct(f,-1)
                alpha[1]_-1-alpha[1]_-1|0>

            .. NOTE::

                This method coerces both elements to the same parent
                in order to implement a Lie conformal algebra the user
                needs to implement :meth:`_nproduct_`
            """
            return self._nproduct_(rhs,n)

        def _nproduct_(self,rhs,n):
            r"""
            Returns the n-th product of this two elements.

            If `n\geq 0` it returns the element of this Lie conformal 
            algebra. If `n < 0` then it first lifts this element to the
            universal  enveloping vertex algebra and returns the 
            corresponding element

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
            vertex algebra, however, the universal enveloping vertex
            algebra needs to be constructed first::

                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L._nproduct_(L,-3)
                L_-4L_-2|0>

            .. NOTE::

                It is guaranteed that both are elements of the same
                parent.
            """ 
            if n >= 0:
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                return self.lift().nproduct(rhs,n)

        @abstract_method(optional=True)
        def lift(self):
            r"""
            Returns the image of this element under the canonical lift
            to the universal enveloping vertex algebra. 

            .. WARNING::

                The universal enveloping algebra needs to be constructed 
                first for this morphism to be defined. 

                This morphism is registered as a coercion between this 
                Lie conformal algebra and its universal enveloping 
                vertex algebra upon creation. Since we consider central
                quotients of the universal enveloping vertex algebras 
                by fixed central parameters, each time a different 
                universal enveloping vertex algebra is constructed, 
                this lift morphism is changed. See the examples below
                and also :meth:`register_lift()\
                <sage.algebras.vertex_algebras.vertex_algebra.\
                UniversalEnvelopingVertexAlgebra.register_lift>`.


            EXAMPLES:

            We lift to the universal enveloping vertex algebra of the 
            Virasoro Lie conformal algebra with central charge `0`::
               
                sage: Vir = VirasoroLieConformalAlgebra(QQ); L = Vir.0
                sage: V = Vir.universal_enveloping_algebra()
                sage: L.lift()
                L_-2|0>
                sage: L.lift().__class__
                <class 'sage.algebras.vertex_algebras.vertex_algebra.UniversalEnvelopingVertexAlgebra_with_category.element_class'>
                sage: L.lift().parent()
                The universal enveloping vertex algebra of The Virasoro Lie conformal algebra over Rational Field

            Notice that the target of the ``lift`` morphism changes when
            we construct another universal enveloping vertex algebra::

                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of The Virasoro Lie conformal algebra over Rational Field

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
                The universal enveloping vertex algebra of The Virasoro Lie conformal algebra over Rational Field

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: Vir.lift.codomain()
                The universal enveloping vertex algebra of The Virasoro Lie conformal algebra over Rational Field

                sage: V.register_lift()
                sage: Vir.lift.codomain()
                The Virasoro vertex algebra at central charge 1/2                
            """
            raise NotImplementedError("Not implemented")
        
        @abstract_method(optional=False)
        def T(self,n=1):
            r"""
            The n-th derivative of this element.

            INPUT:

            - ``n`` -- integer (default:`1`); How many times
              to apply `T` to this element. 

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
            raise NotImplementedError("Not implemented")

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

        def WithBasis(self):
            """
            The subcategory of Lie conformal algebras with a preferred
            basis.

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).WithBasis()
                Category of Lie conformal algebras over Rational Field with basis
            """
            return self._with_axiom("WithBasis")

        def Graded(self, base_ring=None):
            """
            The subcategory of H-Graded Lie conformal algebras.

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
                                    "FinitelyGeneratedAsLieConformalAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Graded()._with_axioms(axioms)
            return GradedModulesCategory.category_of(self)

        def Super(self, base_ring=None):
            """
            The subcategory of super Lie conformal algebras.

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
                                    "FinitelyGeneratedAsLieConformalAlgebra"])
                axioms = axioms_whitelist.intersection(self.axioms())
                return self.base_category().Super()._with_axioms(axioms)
            return SuperModulesCategory.category_of(self)

    class Super(SuperModulesCategory):
        """
        The subcategory of super Lie conformal algebras.

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
        #Need to do all this to make Super commute with Graded. 
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: LieConformalAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                [Category of finitely generated super Lie conformal algebras over Rational Field,
                 Category of super H-graded Lie conformal algebras over Rational Field]

            """
            return [self.base_category(),]

        class SubcategoryMethods:
            
            def Graded(self, base_ring=None):
                """
                The subcategory of H-graded super Lie conformal algebras

                EXAMPLES::
                
                    sage: LieConformalAlgebras(QQ).Super().Graded()
                    Category of super H-graded Lie conformal algebras over Rational Field

                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",  
                                        "FinitelyGeneratedAsLieConformalAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Graded()._with_axioms(axioms)

                return GradedModulesCategory.category_of(
                                         self.base_category()).Super()

            def WithBasis(self):
                """
                The subcategory of Super Lie conformal algebras
                with basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Super().WithBasis()
                    Category of super Lie conformal algebras with basis over Integer Ring
                    sage: LieConformalAlgebras(ZZ).Super().WithBasis() is LieConformalAlgebras(ZZ).WithBasis().Super()
                    True
                """
                return self._with_axiom("WithBasis")

            def FinitelyGeneratedAsLieConformalAlgebra(self):
                """
                The subcategory of finitely generated super Lie 
                conformal algebras.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Super().FinitelyGenerated()
                    Category of finitely generated super Lie conformal algebras over Integer Ring
                """
                return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of Super Lie conformal algebras with basis.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Super().WithBasis()
                Category of super Lie conformal algebras with basis over Integer Ring
                sage: LieConformalAlgebras(ZZ).Super().WithBasis() is LieConformalAlgebras(ZZ).WithBasis().Super()
                True
            """
            
            class SubcategoryMethods:

                def FinitelyGeneratedAsLieConformalAlgebra(self):
                    """
                    The subcategory of finitely generated super Lie
                    conformal algebras with basis.

                    EXAMPLES::

                        sage: LieConformalAlgebras(ZZ).Super().FinitelyGenerated().WithBasis()
                        Category of finitely generated super Lie conformal algebras with basis over Integer Ring
                    """
                    return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

            class FinitelyGeneratedAsLieConformalAlgebra(
                                                CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated super Lie
                conformal algebras with basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Super().FinitelyGenerated().WithBasis()
                    Category of finitely generated super Lie conformal algebras with basis over Integer Ring
                """
                pass

        class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated super Lie 
            conformal algebras.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Super().FinitelyGenerated()
                Category of finitely generated super Lie conformal algebras over Integer Ring
            """
            pass

        class ParentMethods:

            def is_super(self):
                """ 
                Wether this vertex algebra is a super vertex algebra.

                EXAMPLES::

                    sage: V = VirasoroLieConformalAlgebra(QQ)
                    sage: V.is_super()
                    Traceback (most recent call last):
                    AttributeError: 'VirasoroLieConformalAlgebra_with_category' object has no attribute 'is_super'
                """
                return True

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

                sage: LieConformalAlgebras(QQbar).Super().Graded()
                Category of super H-graded Lie conformal algebras over Algebraic Field
                sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
                Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
            """
            return "H-graded {}".format(self.base_category().\
                                        _repr_object_names())


        class SubcategoryMethods:

            def Super(self, base_ring=None):
                """
                The subcategory of H-graded super Lie conformal algebras

                EXAMPLES::
                
                    sage: LieConformalAlgebras(QQ).Super().Graded()
                    Category of super H-graded Lie conformal algebras over Rational Field

                """
                assert base_ring is None or base_ring is self.base_ring()
                if isinstance(self,CategoryWithAxiom_over_base_ring):
                    axioms_whitelist = frozenset(["WithBasis",  
                                        "FinitelyGeneratedAsLieConformalAlgebra"])
                    axioms = axioms_whitelist.intersection(self.axioms())
                    return self.base_category().Super()._with_axioms(axioms)
                return SuperModulesCategory.category_of(self)

            def WithBasis(self):
                """
                The subcategory of H-graded Lie conformal algebras with 
                basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Graded().WithBasis()
                    Category of H-graded Lie conformal algebras with basis over Integer Ring
                """
                return self._with_axiom("WithBasis")

            def FinitelyGeneratedAsLieConformalAlgebra(self):
                """
                The subcategory of finitely generated H-graded Lie 
                conformal algebras

                EXAMPLES::

                    sage: C = LieConformalAlgebras(ZZ).Graded().FinitelyGenerated(); C
                    Category of finitely generated H-graded Lie conformal algebras over Integer Ring
                    sage: C is LieConformalAlgebras(ZZ).FinitelyGenerated().Graded()
                    True
                """
                return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of H-graded Lie conformal algebras with 
            basis.

            EXAMPLES::

                sage: LieConformalAlgebras(ZZ).Graded().WithBasis()
                Category of H-graded Lie conformal algebras with basis over Integer Ring
            """
            
            class SubcategoryMethods:

                def FinitelyGeneratedAsLieConformalAlgebra(self):
                    """
                    The subcategory of finitely generated H-graded
                    Lie conformal algebras with basis.

                    EXAMPLES::

                        sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
                        Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
                    """
                    return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")

            class FinitelyGeneratedAsLieConformalAlgebra(
                                                CategoryWithAxiom_over_base_ring):
                """
                The subcategory of finitely generated H-graded
                Lie conformal algebras with basis.

                EXAMPLES::

                    sage: LieConformalAlgebras(ZZ).Graded().FinitelyGenerated().WithBasis()
                    Category of finitely generated H-graded Lie conformal algebras with basis over Integer Ring
                    """
                pass

        class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated H-graded Lie 
            conformal algebras

            EXAMPLES::

                sage: C = LieConformalAlgebras(ZZ).Graded().FinitelyGenerated(); C
                Category of finitely generated H-graded Lie conformal algebras over Integer Ring
                sage: C is LieConformalAlgebras(ZZ).FinitelyGenerated().Graded()
                True
            """
            pass

        class Super(SuperModulesCategory):
            """
            The subcategory of H-graded super Lie conformal algebras.

            EXAMPLES::

                sage: C = LieConformalAlgebras(QQbar).Graded().Super(); C
                Category of super H-graded Lie conformal algebras over Algebraic Field
                sage: C is LieConformalAlgebras(QQbar).Super().Graded()
                True
            """
            pass

        class ParentMethods:
            def is_graded(self):
                """
                Wether this Lie conformal algebra is graded or not

                EXAMPLES::

                    sage: Vir = VirasoroLieConformalAlgebra(QQ)
                    sage: Vir
                    The Virasoro Lie conformal algebra over Rational Field
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
        """
        The subcategory of Lie conformal algebras with basis.

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).WithBasis()
            Category of Lie conformal algebras over Algebraic Field with basis
        """
        def _repr_object_names(self):
            """
            The names of objects of this category
            """
            return "{} with basis".format(self.base_category().\
                                          _repr_object_names())

        class SubcategoryMethods:
            def FinitelyGenerated(self):
                """
                The subcategory of finitely generated Lie conformal 
                algebras with basis. 

                EXAMPLES::

                    sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                    Category of finitely generated Lie conformal algebras with basis over Algebraic Field
                """
                return self._with_axiom("FinitelyGeneratedAsLieConformalAlgebra")
            
        class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
            """
            The subcategory of finitely generated Lie conformal 
            algebras with basis. 

            EXAMPLES::

                sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
                Category of finitely generated Lie conformal algebras with basis over Algebraic Field
            """
            pass


    class FinitelyGeneratedAsLieConformalAlgebra(
                                            CategoryWithAxiom_over_base_ring):
        """
        The subcategory of finitely generated Lie conformal 
        algebras with basis. 

        EXAMPLES::

            sage: LieConformalAlgebras(QQbar).WithBasis().FinitelyGenerated()
            Category of finitely generated Lie conformal algebras with basis over Algebraic Field
        """
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


class _LiftMorphism(Morphism):
    """
    The lift morphism from this Lie conformal algebra to its universal
    enveloping vertex algebra

    EXAMPLES::
    
        sage: Vir = VirasoroLieConformalAlgebra(AA)
        sage: cp = Family({Vir.1:2})
        sage: V = Vir.universal_enveloping_algebra(cp)
        sage: Vir.lift
        Generic morphism:
          From: The Virasoro Lie conformal algebra over Algebraic Real Field
          To:   The universal enveloping vertex algebra of The Virasoro Lie conformal algebra over Algebraic Real Field
        sage: W = VirasoroVertexAlgebra(AA,1)
        sage: Vir.lift
        Generic morphism:
          From: The Virasoro Lie conformal algebra over Algebraic Real Field
          To:   The Virasoro vertex algebra at central charge 1
    """
    def _call_(self,x):
        return x.lift()    

