"""
Finitely Generated Vertex Algebras

AUTHORS:

- Reimundo Heluani (2019-10-09): Initial implementation.
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

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.rings.all import QQ,ZZ
from sage.categories.graded_modules import GradedModulesCategory
from sage.categories.super_modules import SuperModulesCategory
from sage.categories.quotients import QuotientsCategory
from sage.misc.abstract_method import abstract_method
from sage.combinat.free_module import CombinatorialFreeModule

class FinitelyGeneratedAsVertexAlgebra(CategoryWithAxiom_over_base_ring):
    """
    The subcategory of finitely generated vertex algebras.

    EXAMPLES::

        sage: VertexAlgebras(QQ).FinitelyGenerated()
        Category of finitely generated vertex algebras over Rational Field
    """
    class ParentMethods:
        @abstract_method
        def gens(self):
            """
            The list of generators of this vertex algebra.

            EXAMPLES::

                sage: V = vertex_algebras.Virasoro(QQ, 1/2);
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
                sage: V = vertex_algebras.Affine(QQ, 'A1', 1); V
                The universal affine vertex algebra of CartanType ['A', 1] at level 1 over Rational Field
                sage: V.gens()
                (E(alpha[1])_-1|0>, E(alphacheck[1])_-1|0>, E(-alpha[1])_-1|0>)
            """

        def ngens(self):
            """
            The number of generators of this vertex algebra.

            EXAMPLES::

                sage: vertex_algebras.Virasoro(QQ, 1/2).ngens()
                1
                sage: vertex_algebras.Affine(QQ, 'A2', 1).ngens()
                8
            """
            return len(self.gens())

        def gen(self,i):
            """
            The ``i``-th generator of this vertex algebra.

            EXAMPLES::

                sage: V = vertex_algebras.N2(QQbar, 1); V
                The N=2 super vertex algebra of central charge 1 over Algebraic Field
                sage: V.gen(1)
                J_-1|0>
            """
            return self.gens()[i]
    
        def some_elements(self):
            """
            A list of elements in this vertex algebra.

            EXAMPLES::

                sage: V = vertex_algebras.FreeFermions(QQ, 4)
                sage: V.some_elements()
                (psi_0_-1/2|0>, psi_1_-1/2|0>, psi_2_-1/2|0>, psi_3_-1/2|0>)
            """
            #Arithmetics are so slow that if we include any cuadratic term
            #here the testsuite will not succeed on time. 
            return self.gens()

    class Quotients(QuotientsCategory):
        """
        Base class of quotients of finitely generated vertex algebras.
        """
        class ParentMethods:

            def gens(self):
                """
                A list of generators for this quotient.

                EXAMPLES::

                    sage: V = vertex_algebras.NeveuSchwarz(QQ,7/10)
                    sage: Q = V.quotient(V.ideal(V.find_singular(4)))
                    sage: Q.gens()
                    (L_-2|0>, G_-3/2|0>)
                """
                return tuple(self.retract(g) for g in\
                             self.cover_algebra().gens()\
                             if not self.retract(g).is_zero())

    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The category of finitely generated vertex algebras with basis.

        EXAMPLES::

            sage: VertexAlgebras(QQ).FinitelyGenerated().WithBasis()
            Category of finitely generated vertex algebras with basis over Rational Field
        """
        class ElementMethods: #FinitelyGenerated.WithBasis

            @abstract_method
            def _li_filtration_monomial_degree(self):
                r"""
                The maximum `p` such that this monomial belongs
                to the `p`-th filtered part of this vertex
                algebra.

                EXAMPLES::

                    sage: V = vertex_algebras.Weyl(QQ,4)
                    sage: V.inject_variables()
                    Defining alpha0, alpha1, alpha2, alpha3
                    sage: v = alpha0.T(2)*(alpha1.T()*alpha2); v
                    2*alpha0_(-3)alpha1_(-2)alpha2_(-1)|0>
                    sage: v._li_filtration_monomial_degree()
                    3
                """

            def li_filtration_lt(self):
                """
                The leading terms of this element with respect
                to the Li filtration.

                EXAMPLES:

                A vertex algebra does not need to be H-graded
                to work with its Li filtration::

                    sage: V = vertex_algebras.Weyl(QQ,4)
                    sage: v = V.an_element()
                    sage: V = vertex_algebras.Weyl(QQ,4)
                    sage: v = V.an_element(); v
                    alpha1_(-1)alpha2_(-2)alpha3_(-3)|0> + |0> + 2*alpha0_(-1)|0> + 3*alpha1_(-1)|0>
                    sage: v.li_filtration_degree()
                    0
                    sage: w = v.li_filtration_lt(); w
                    |0> + 2*alpha0_(-1)|0> + 3*alpha1_(-1)|0>
                    sage: [m._li_filtration_monomial_degree() for m in w.monomials()]
                    [0, 0, 0]

                For graded algebras it works the same way::

                    sage: V = vertex_algebras.Virasoro(QQ,1)
                    sage: v = V.an_element(); v
                    |0> + 2*L_-2|0> + 3*L_-3|0> + L_-2L_-2L_-2L_-2|0>
                    sage: v.li_filtration_lt()
                    |0> + 2*L_-2|0> + L_-2L_-2L_-2L_-2|0>
                """
                if self.is_zero():
                    return self
                lt = [(m._li_filtration_monomial_degree(),m) for m in self.terms()]
                lideg = min(k for k,v in lt)
                return sum(v for k,v in lt if k==lideg)

            def li_filtration_degree(self):
                r"""
                The maximum `p` such that this element belongs
                to the `p`-th filtered part of this vertex
                algebra.

                EXAMPLES::

                    sage: V = vertex_algebras.Virasoro(QQ,1/2); V.register_lift()
                    sage: v = V.find_singular(6)[0]; v
                    L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>
                    sage: v.li_filtration_degree()
                    0
                    sage: v.li_filtration_lt()
                    L_-2L_-2L_-2|0>
                """
                if self.is_zero():
                    from sage.rings.infinity import Infinity
                    return Infinity
                return min(m._li_filtration_monomial_degree() for m in self.monomials())

        class Super(SuperModulesCategory):
            """
            The subcategory of super finitely generated vertex algebras
            with basis.

            EXAMPLES::

                sage: VertexAlgebras(AA).FinitelyGenerated().WithBasis().Super()
                Category of super finitely generated vertex algebras with basis over Algebraic Real Field
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category

                EXAMPLES::

                    sage: VertexAlgebras(QQ).FinitelyGenerated().WithBasis().Super().super_categories()
                    [Category of super vertex algebras with basis over Rational Field,
                    Category of finitely generated vertex algebras with basis over Rational Field,
                    Category of super finitely generated vertex algebras over Rational Field]
                """
                return [self.base_category(),]

        class Graded(GradedModulesCategory):
            """
            The subcategory of H-graded finitely generated vertex
            algebras with basis.

            .. NOTE::

                To implement a finitely generated H-graded
                vertex algebra with basis. The basis indices
                should implement a method `subset` which
                admits at least the keyword parameter
                ``energy`` such that
                ``_indices.subset(energy=n)`` enumerates
                the indices of the basis of monomials with
                conformal weight ``n``. See for example
                :mod:`EnergyPartitionTuples<sage.algebras.vertex_algebras.energy_partition_tuples>`.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).FinitelyGenerated().WithBasis().Graded()
                Category of H-graded finitely generated vertex algebras with basis over Algebraic Field
            """
            def _repr_object_names(self):
                """
                The names of the objects of this category.

                EXAMPLES::

                    sage: VertexAlgebras(QQbar).FinitelyGenerated().Graded()
                    Category of H-graded finitely generated vertex algebras over Algebraic Field
                """
                return "H-graded {}".format(self.base_category()._repr_object_names())

            class ParentMethods:

                def get_weight_less_than(self,n):
                    """
                    The sub-vector space of this vertex algebra of
                    vectors with conformal weight less than ``n``.

                    INPUT:

                    - ``n`` -- a non-negative rational number;

                    OUTPUT:

                    A submodule of this vertex algebra.

                    EXAMPLES::

                        sage: V = vertex_algebras.FreeFermions(QQ); M = V.get_weight_less_than(5/2); M
                        Free module generated by {0, 1, 2, 3} over Rational Field
                        sage: [v.lift() for v in M.basis()]
                        [|0>, psi_-1/2|0>, psi_-3/2|0>, psi_-3/2psi_-1/2|0>]

                    TESTS::

                        sage: V = vertex_algebras.FreeFermions(QQ); M = V.get_weight_less_than(5/2);
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
                    The sub-module of this vertex algebra of
                    vectors with conformal weight equal to ``n``.

                    INPUT:

                    - ``n`` -- a non-negative rational number;

                    OUTPUT:

                    A submodule of this vertex algebra.

                    EXAMPLES::

                        sage: V = vertex_algebras.NeveuSchwarz(QQ, 1/2); M = V.get_weight(7/2)
                        sage: [v.lift() for v in M.basis()]
                        [L_-2G_-3/2|0>, G_-7/2|0>]
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

                        sage: V = vertex_algebras.Virasoro(QQ,1/2); V.register_lift(); V.dimension_at_weight(8)
                        7
                        sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                        sage: Q.dimension_at_weight(8)      # long time (2 seconds)
                        5
                        sage: vertex_algebras.NeveuSchwarz(QQ,3/2).dimension_at_weight(11/2)
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

                    INPUT:

                    - ``n`` -- a positive rational number

                    OUTPUT:

                    A list of vectors in this vertex algebra that
                    form a basis for the subspace of conformal
                    weight ``n``.

                    .. SEEALSO::

                        :meth:`get_weight`

                    EXAMPLES:

                    We find the singular vector of the Virasoro
                    vertex algebra of central charge 1/2::

                        sage: V = vertex_algebras.Virasoro(QQ,1/2); V.find_singular(6)
                        (L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>,)

                    Finding singular vectors work for quotients of
                    vertex algebras as well::

                        sage: V = vertex_algebras.Affine(QQ,'A1',1,names = ('e','h','f'));
                        sage: V.register_lift()
                        sage: sing = V.find_singular(2); sing
                        (e_-1e_-1|0>,
                         e_-1h_-1|0> + e_-2|0>,
                         h_-1h_-1|0> - 2*e_-1f_-1|0> + h_-2|0>,
                         h_-1f_-1|0> + f_-2|0>,
                         f_-1f_-1|0>)
                        sage: Q = V.quotient(V.ideal([sing[0]],check=False))
                        sage: Q.find_singular(2)        # long time (1 second)
                        (e_-1h_-1|0> + e_-2|0>,
                         h_-1h_-1|0> - 2*e_-1f_-1|0> + h_-2|0>,
                         h_-1f_-1|0> + f_-2|0>,
                         f_-1f_-1|0>)

                    And for super vertex algebras::

                        sage: V = vertex_algebras.NeveuSchwarz(QQ, 7/10); V.register_lift();  V
                        The Neveu-Schwarz super vertex algebra of central charge 7/10 over Rational Field
                        sage: V.find_singular(4)
                        (L_-2L_-2|0> - 3/2*G_-5/2G_-3/2|0> + 3/10*L_-4|0>,)

                    TESTS::

                        sage: V = vertex_algebras.Affine(QQ,'A1',1); V.register_lift(); sing = V.find_singular(2)
                        sage: Q = V.quotient(V.ideal(sing))
                        sage: Q.find_singular(2)
                        ()
                        sage: V = vertex_algebras.Virasoro(QQ,1/2); V.find_singular(0)
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

            class Super(SuperModulesCategory):
                """
                The subcategory of super H-graded finitely generated
                vertex algebras with basis.

                EXAMPLES::

                    sage: C = VertexAlgebras(QQbar).FinitelyGenerated().WithBasis()
                    sage: C.Graded().Super()
                    Category of super H-graded finitely generated vertex algebras with basis over Algebraic Field
                    sage: C.Graded().Super() is C.Super().Graded()
                    True
                """
                def extra_super_categories(self):
                    """
                    The extra super categories of this category.

                    EXAMPLES::

                        sage: VertexAlgebras(QQ).FinitelyGenerated().WithBasis().Graded().Super().super_categories()
                        [Category of super finitely generated vertex algebras with basis over Rational Field,
                        Category of super H-graded vertex algebras with basis over Rational Field,
                        Category of H-graded finitely generated vertex algebras with basis over Rational Field,
                        Category of super H-graded finitely generated vertex algebras over Rational Field]
                    """
                    return [self.base_category(),]

            class Quotients(QuotientsCategory):
                """
                The category of quotients of finitely generated
                H-graded vertex algebras with basis.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).FinitelyGenerated().WithBasis().Graded().Quotients()
                    Category of quotients of H-graded finitely generated vertex algebras with basis over Rational Field
                """
                class ElementMethods:

                    def _li_filtration_monomial_degree(self):
                        r"""
                        The maximum `p` such that this monomial belongs
                        to the `p`-th filtered part of this vertex
                        algebra.

                        EXAMPLES::

                            sage: V = vertex_algebras.Virasoro(QQ,1/2)
                            sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                            sage: v = Q([[4,2,1]]); v._li_filtration_monomial_degree()  # long time (2 seconds)
                            4
                        """
                        return self.lift()._li_filtration_monomial_degree()

                    def weight(self):
                        """
                        The conformal weight of this element.

                        EXAMPLES::

                            sage: V = vertex_algebras.Virasoro(QQ,1/2)
                            sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                            sage: v = Q([[4,2,1]]); v.weight()  # long time (2 seconds)
                            10
                        """
                        return self.lift().weight()

                    def T(self, n=1):
                        """
                        The ``n``-th derivative of this element.

                        INPUT:

                        - ``n`` -- a non-negative integer (default:
                          ``1``); how many times to apply `T`.

                        EXAMPLES::

                            sage: V = vertex_algebras.Virasoro(QQ,1/2)
                            sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                            sage: v = Q([[4,2,1]]); v.T()   # long time (3 seconds)
                            -18*L_-5L_-4L_-2|0> - 4*L_-6L_-3L_-2|0> + 8*L_-7L_-2L_-2|0> - 6*L_-6L_-5|0> - 9/2*L_-7L_-4|0> + 45/2*L_-8L_-3|0> + 81*L_-9L_-2|0> + 16*L_-11|0>
                        """
                        return self.parent().retract(self.lift().T(n))

                    def _mul_(self, other):
                        """
                        The normally ordered product of these two
                        elements.

                        EXAMPLES::

                            sage: V = vertex_algebras.Virasoro(QQ,1/2)
                            sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                            sage: L = Q.0
                            sage: L*(L*L)
                            -93/64*L_-3L_-3|0> + 33/8*L_-4L_-2|0> + 27/16*L_-6|0>
                        """
                        return self.parent().retract(
                                            self.lift()._mul_(other.lift()))

                    def _bracket_(self, other):
                        r"""
                        The `\lambda`-bracket of these two elements.

                        EXAMPLES::

                            sage: V = vertex_algebras.Virasoro(QQ,1/2)
                            sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                            sage: L = Q.0;
                            sage: L._bracket_(L*L)
                            {0: 2*L_-3L_-2|0> + L_-5|0>,
                             1: 4*L_-2L_-2|0>,
                             2: 3*L_-3|0>,
                             3: 17/2*L_-2|0>,
                             5: 3/2*|0>}
                        """
                        p = self.parent()
                        sl = self.lift()
                        ol = other.lift()
                        return {k:p.retract(v) for k,v in sl._bracket_(ol).\
                                items() if not p.retract(v).is_zero()}

    class Super(SuperModulesCategory):
        """
        The subcategory of super finitely generated vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(AA).FinitelyGenerated().Super()
            Category of super finitely generated vertex algebras over Algebraic Real Field
        """
        def extra_super_categories(self):
            """
            The extra super categories of this category

            EXAMPLES::

                sage: VertexAlgebras(QQ).FinitelyGenerated().Super().super_categories()
                [Category of super vertex algebras over Rational Field,
                Category of finitely generated vertex algebras over Rational Field]
            """
            return [self.base_category(),]

    class Graded(GradedModulesCategory):
        """
        The subcategory of H-graded finitely generated vertex algebras.

        EXAMPLES::

            sage: VertexAlgebras(QQbar).FinitelyGenerated().Graded()
            Category of H-graded finitely generated vertex algebras over Algebraic Field
        """
        def _repr_object_names(self):
            """
            The names of the objects of this category.

            EXAMPLES::

                sage: VertexAlgebras(QQbar).FinitelyGenerated().Graded()
                Category of H-graded finitely generated vertex algebras over Algebraic Field
            """
            return "H-graded {}".format(self.base_category()._repr_object_names())

        class ParentMethods:

            @abstract_method
            def get_weight_less_than(self,n):
                """
                The sub-vector space of this vertex algebra of
                vectors with conformal weight less than ``n``.

                INPUT:

                - ``n`` -- a non-negative rational number;

                OUTPUT:

                A submodule of this vertex algebra.

                EXAMPLES::

                    sage: V = vertex_algebras.FreeFermions(QQ); M = V.get_weight_less_than(5/2); M
                    Free module generated by {0, 1, 2, 3} over Rational Field
                    sage: [v.lift() for v in M.basis()]
                    [|0>, psi_-1/2|0>, psi_-3/2|0>, psi_-3/2psi_-1/2|0>]

                TESTS::

                    sage: V = vertex_algebras.FreeFermions(QQ); M = V.get_weight_less_than(5/2);
                    sage: M.reduce(V.vacuum())
                    0
                """

            @abstract_method
            def get_weight(self,n):
                r"""
                The sub-module of this vertex algebra of
                vectors with conformal weight equal to ``n``.

                INPUT:

                - ``n`` -- a non-negative rational number;

                OUTPUT:

                A submodule of this vertex algebra.

                EXAMPLES::

                    sage: V = vertex_algebras.NeveuSchwarz(QQ, 1/2); M = V.get_weight(7/2)
                    sage: [v.lift() for v in M.basis()]
                    [L_-2G_-3/2|0>, G_-7/2|0>]
                """

            @abstract_method
            def dimension_at_weight(self,n):
                """
                The dimension of the space of conformal weight
                ``n`` of this vertex algebra.

                INPUT:

                - ``n`` -- a non-negative rational number.

                EXAMPLES::

                    sage: V = vertex_algebras.Virasoro(QQ,1/2); V.register_lift(); V.dimension_at_weight(8)
                    7
                    sage: Q = V.quotient(V.ideal(V.find_singular(6)))
                    sage: Q.dimension_at_weight(8)  # long time (2 seconds)
                    5
                    sage: vertex_algebras.NeveuSchwarz(QQ,3/2).dimension_at_weight(11/2)
                    5
                """

            @abstract_method
            def find_singular(self,n):
                """
                Return a basis of the vector space of singular
                vectors of weight `n`.

                INPUT:

                - ``n`` -- a positive rational number

                OUTPUT:

                A list of vectors in this vertex algebra that
                form a basis for the subspace of conformal
                weight ``n``.

                .. SEEALSO::

                    :meth:`get_weight`

                EXAMPLES:

                We find the singular vector of the Virasoro
                vertex algebra of central charge 1/2::

                    sage: V = vertex_algebras.Virasoro(QQ,1/2); V.find_singular(6)
                    (L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>,)

                Finding singular vectors work for quotients of
                vertex algebras as well::

                    sage: V = vertex_algebras.Affine(QQ,'A1',1,names = ('e','h','f'));
                    sage: V.register_lift()
                    sage: sing = V.find_singular(2); sing
                    (e_-1e_-1|0>,
                     e_-1h_-1|0> + e_-2|0>,
                     h_-1h_-1|0> - 2*e_-1f_-1|0> + h_-2|0>,
                     h_-1f_-1|0> + f_-2|0>,
                     f_-1f_-1|0>)
                    sage: Q = V.quotient(V.ideal([sing[0]],check=False))
                    sage: Q.find_singular(2)    # long time (1 second)
                    (e_-1h_-1|0> + e_-2|0>,
                     h_-1h_-1|0> - 2*e_-1f_-1|0> + h_-2|0>,
                     h_-1f_-1|0> + f_-2|0>,
                     f_-1f_-1|0>)

                And for super vertex algebras::

                    sage: V = vertex_algebras.NeveuSchwarz(QQ, 7/10); V.register_lift();  V
                    The Neveu-Schwarz super vertex algebra of central charge 7/10 over Rational Field
                    sage: V.find_singular(4)
                    (L_-2L_-2|0> - 3/2*G_-5/2G_-3/2|0> + 3/10*L_-4|0>,)

                TESTS::

                    sage: V = vertex_algebras.Affine(QQ,'A1',1); V.register_lift(); sing = V.find_singular(2)
                    sage: Q = V.quotient(V.ideal(sing))
                    sage: Q.find_singular(2)
                    ()
                    sage: V = vertex_algebras.Virasoro(QQ,1/2); V.find_singular(0)
                    (|0>,)
                """

            def hilbert_series(self,ord):
                r"""
                The graded dimension of this algebra.

                INPUT:

                - ``ord`` -- a positive rational number; the precision
                  order of the result.

                OUTPUT:

                The sum

                .. MATH::

                    \sum_{n = 0}^{ord} q^n \mathrm{dim} V_n

                where `n` runs over all rationals such that `V_n \neq 0`

                EXAMPLES::

                    sage: V = vertex_algebras.Virasoro(QQ,1/2); V.register_lift()
                    sage: V.hilbert_series(8)
                    1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 4*q^6 + 4*q^7 + O(q^8)
                    sage: L = V.0; v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: Q = V.quotient(V.ideal(v))
                    sage: Q.hilbert_series(9)   # long time (2 seconds)
                    1 + q^2 + q^3 + 2*q^4 + 2*q^5 + 3*q^6 + 3*q^7 + 5*q^8 + O(q^9)

                    sage: V = vertex_algebras.FreeFermions(QQ)
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

        class ElementMethods:

            def _action_from_partition_tuple(self, pt, negative=True):
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

                    sage: V = vertex_algebras.Virasoro(QQ,1/2); vac = V.vacuum()
                    sage: vac._action_from_partition_tuple([[3,2]])
                    L_-3L_-2|0>
                    sage: V = vertex_algebras.NeveuSchwarz(QQ,1); vac = V.vacuum()
                    sage: vac._action_from_partition_tuple([[2,2],[5/2,3/2]])
                    L_-2L_-2G_-5/2G_-3/2|0>
                    sage: V = vertex_algebras.Affine(QQ,'A1',2, names = ('e','h','f'));
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
                for each homogeneous vector `v \in V` we have
                `v_n a = 0` whenever `n > 0`.

                EXAMPLES::

                    sage: V = vertex_algebras.Virasoro(QQ,1/2); V.register_lift(); L = V.0
                    sage: v =  L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                    sage: v.is_singular()
                    True
                    sage: V = vertex_algebras.Affine(QQ, 'A1', 2); E = V.0;
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

        class Super(SuperModulesCategory):
            """
            The subcategory of super H-graded finitely generated vertex
            algebras.

            EXAMPLES::

                sage: C = VertexAlgebras(QQbar).FinitelyGenerated()
                sage: C.Graded().Super()
                Category of super H-graded finitely generated vertex algebras over Algebraic Field
                sage: C.Graded().Super() is C.Super().Graded()
                True
            """
            def extra_super_categories(self):
                """
                The extra super categories of this category.

                EXAMPLES::

                    sage: VertexAlgebras(QQ).FinitelyGenerated().Graded().Super().super_categories()
                    [Category of super finitely generated vertex algebras over Rational Field,
                    Category of super H-graded vertex algebras over Rational Field,
                    Category of H-graded finitely generated vertex algebras over Rational Field]
                """
                return [self.base_category(),]
