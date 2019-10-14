r"""
Vertex algebras
AUTHORS

- Reimundo Heluani (10-09-2019): Initial implementation

Let `R` be a ring. A *Vertex algebra* [Kac1997]_ over `R` is the datum of:

- an `R` module `V` called *the space of states*
- a vector `|0\rangle \in V` called the *vacuum* vector.
- an endomorphism `T \in End(V)` called the *translation operator*
- a bilinear multiplication map 

.. MATH:: 
    V \otimes V \rightarrow V( (z)), \qquad a \otimes b \mapsto Y(a,z)b =:
    \sum_{n \in \ZZ} a_{(n)}b z^{-1-n}.
    :label: statefield

Subject to the following set of axioms:

- The vacuum axioms:

.. MATH::
    Y(a,z)|0\rangle = a + Ta \cdot z + O(z^2), \qquad Y(|0\rangle,z)a = a, 
    \qquad a \in V.

- Translation invariance:

.. MATH::
    [T,Y(a,z)] = \frac{d}{dz} Y(a,z)

- Locality:

.. MATH::

    (z-w)^n Y(a,z)Y(b,w) = (z-w)^n Y(b,w)Y(a,z), \qquad a,b \in V, \: n \gg 0

The `\ZZ`-many products `a_{(n)}b` defined in :eq:`statefield` are called the
*nth*-products. A vertex algebra toghether with its non-negative nth-products
and its translation operator `T` is a 
:mod:`Lie Conformal Algebra<sage.categories.lie_conformal_algebras>`. 
The generating function for these non-negative products

.. MATH::
    [a_\lambda b] = \sum_{n \geq 0} \frac{\lambda^n}{n!} a_{(n)} b

is called the *OPE* or the `\lambda`-bracket of `a` and `b`. This forgetful 
functor admits a
left adjoint: to each Lie conformal algebra `L` we attach a universal enveloping
vertex algebra `U(L)`. This vertex algebra admits an embedding (when `R` is a
field of characteristic zero) `L \hookrightarrow U(L)` of Lie conformal
algebras. 

A vertex algebra is called H-Graded [DSK2006]_ if there exists a diagonalizable
operator `H \in End(V)` such that 

.. MATH::
    [H,Y(a,z)] = Y(Ha,z) + z \frac{d}{dz} Y(a,z).

Equivalently, there exists a decomposition `V = \oplus_{\Delta \in \CC}
V_\Delta` such that the nth product becomes graded of degree `-1-n`, that is

.. MATH::
    a_{(n)}b \in V_{p + q - 1 - n} \qquad a \in V_p, b \in V_q, n \in \ZZ.

.. NOTE:: 

    Although arbitrary gradings are allowed in the literature, we implement here
    only non-negative integer gradings. 

In this situation, for `a \in V_p` we define the *shifted nth-product* `a_n b =
a_{(n+p-1)}b`. With this convention, the *shifted nth-product* map `a \otimes b
\mapsto a_n b`  is graded of degree `-n`. 

A vertex algebra is called *strongly generated* by a collection of vectors `a^i
\in V` indexed by an ordered set `I`, 
if every vector of `V` can be written as a linear combination of vectors
of the form 

.. MATH::
    a^{i_1}_{(-j_{1,1})} \cdots a^{i_1}_{(-j_{1,n_1})} a^{i_2}_{(-j_{2,1})}
    \cdots a^{i_k}_{(-j_{k,n_j})} |0\rangle, \qquad i_1 < \ldots < i_k,
    \: j_{i,l} > 0. 

A vertex algebra is called *finitely strongly generated* if there is such a set
with a finite `I`.

EXAMPLES:

Typical examples of finitely and strongly generated H-Graded vertex algebra
arise as twisted universal enveloping vertex algebras of a finitely generated
H-Graded :mod:`Lie conformal algebra<sage.categories.lie_conformal_algebras>`.
They are described explicitly in terms of the OPE (`\lambda`-brackets) of 
their generators:

    - The **Virasoro** vertex algebra of central charge `c` is generated by one
      vector `L` with `\lambda`-bracket

        .. MATH::
            [L_\lambda L] = TL + 2 \lambda L + \frac{\lambda^3}{12} c |0\rangle

    - The **universal affine** vertex algebra `V^\kappa(\mathfrak{g})` with level
      `\kappa` associated to a finite dimensional 
      Lie algebra `\mathfrak{g}` with non-degenerate,
      invariant `R`-bilinear form `\kappa(,)` is generated by `\mathfrak{g}` with
      `\lambda`-bracket of generators given by

        .. MATH::
            [a_\lambda b] = [a,b] + \lambda \kappa(a,b) |0\rangle, 
            \qquad a,b \in \mathfrak{g}


    - The **Weyl** vertex algebra, or `\beta-\gamma` system has two generators 
      `\beta` and
      `\gamma` The only non-trivial brackets among
      generators are

        .. MATH::
            [\beta_\lambda \gamma] = - [\gamma_\lambda \beta] = |0\rangle

.. SEEALSO::
    :mod:`sage.algebras.vertex_algebras.vertex_algebra`

Every vertex algebra carries a decreasing filtration 
called the *Li standard filtration*
[Li2005]_. It is defined as follows, we define `F^pV` to be the subspace spanned
by vectors of the form

.. MATH::
    a^1_{(-n_1-1)} \cdots a^r_{(-n_r -1)} b, \qquad a^i,b \in V, n_i \in
    \ZZ_{\geq 0}, \: n_1 + \cdots n_r \geq p.

The associated graded `\mathrm{gr}_FV` is a
:mod:`Poisson vertex algebra<sage.categories.poisson_vertex_algebras>` known as
the *quasi-classical limit of* `V`. 
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
    
    EXAMPLES::

        sage: VertexAlgebras(QQ)
        Category of Vertex algebras over Rational Field
        sage: VertexAlgebras(QQ).is_subcategory(LieConformalAlgebras(QQ))
        True

    """

    @cached_method
    def super_categories(self):
        """
        The super categories of this category

        EXAMPLES::

            sage: C = VertexAlgebras(QQ)
            sage: C.super_categories()
            [Category of Lie conformal algebras over Rational Field]

        """
        return [LieConformalAlgebras(self.base_ring()),]

    def _repr_object_names(self):
        """
        The name of the objects of this category
        """
        return "Vertex algebras over {}".format(self.base_ring())

    class Quotients(QuotientsCategory):
        """
        The category of quotients of vertex algebras
        """
        pass

    class ParentMethods:
        def ideal(self, *gens):
            """
            The ideal of this vertex algebra generated by ``gens``

            INPUT: 

            - ``gens`` -- a list or tuple of elements of this vertex algebra.
              They must be singular vectors. 

            EXAMPLES:

            We construct the ideal defining the *Virasoro Ising* module::

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: I = V.ideal(v)
                sage: I
                ideal of The Virasoro vertex algebra at central charge 1/2 generated by (L_-2L_-2L_-2|0>-27/16*L_-6|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>,)

            If we instead use a non-singular vector::

                sage: V.ideal(L*L)
                Traceback (most recent call last):
                ...
                ValueError: Generators must be singular vectors of The Virasoro vertex algebra at central charge 1/2

            """
            from sage.algebras.vertex_algebras.vertex_algebra_ideal import VertexAlgebraIdeal
            return VertexAlgebraIdeal(self,gens)

        def quotient(self, I):
            """ 
            The quotient of this vertex algebra by the ideal ``I``

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ,1/2)
                sage: L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: I = V.ideal(v)
                sage: Q = V.quotient(I); Q
                Quotient of The Virasoro vertex algebra at central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>-27/16*L_-6|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>,)
                sage: Q(L*(L*L))
                33/8*L_-4L_-2|0>+27/16*L_-6|0>-93/64*L_-3L_-3|0>
            
            """
            return VertexAlgebraQuotient(I)

        def classical_limit(self):
            """
            The Poisson vertex algebra classical limit of this vertex algebra

            EXAMPLES:

            We construct the classical limit of the universal Virasoro vertex
            algebra at central charge `1/2`::

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
                33/8*L_-4L_-2|0>+27/16*L_-6|0>-93/64*L_-3L_-3|0>
                sage: P(L)*(P(L)*P(L)) == P.zero()
                True

            """
            raise NotImplementedError("General classical limit is only"+
                " implemented for H-graded vertex agebras")

        def _element_constructor_(self,x):
            """ 
            Constructs elements of this vertex algebra

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V([[3]])
                L_-4|0>
                sage: V.0
                L_-2|0>
                sage: V.zero()
                0
                sage: V.zero().__class__
                <class 'sage.algebras.vertex_algebras.vertex_algebra.VirasoroVertexAlgebra_with_category.element_class'>
                sage: W = AffineVertexAlgebra(QQ, 'A1', 1)
                sage: W([[2,1],[],[3]])
                E(alpha[1])_-2E(alpha[1])_-1E(-alpha[1])_-3|0>
                sage: 3*W.0
                3*E(alpha[1])_-1|0>

            """
            if x in self.base_ring():
                if x != 0 :
                    raise ValueError("can only convert the scalar 0 "\
                                     "into a vertex algebra element")
                return self.zero()
            return self.element_class(self,x)

        def is_strongly_generated(self):
            """
            If this vertex algebra is strongly generated

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
            If this vertex algebra is H-Graded

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V.is_graded()
                True
                sage: W = AffineVertexAlgebra(QQ, 'A1', 1)
                sage: W.is_graded()
                True

            """
            return self in VertexAlgebras(self.base_ring()).HGraded()

        @abstract_method
        def vacuum(self):
            """
            The vacuum vector of this vertex algebra

            EXAMPLES::

                sage: V = VirasoroVertexAlgebra(QQ, 3)
                sage: V.vacuum()
                |0>

            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def module(self):
            """
            The underlying module of this vertex algebra

            EXAMPLES:

            For universal enveloping vertex algebras we get a
            :class:`CombinatorialFreeModule<sage.combinat.free_module.CombinatorialFreeModule>`::

                sage: W = AffineVertexAlgebra(QQ, 'A1', 2)
                sage: W.module()
                Free module generated by Partition tuples of level 3 over Rational Field

            For quotient algebras we get the algebra itself::

                sage: V = VirasoroVertexAlgebra(QQ, 1/2)
                sage: L = V.0
                sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
                sage: Q = V.quotient(V.ideal(v))
                sage: Q.module()
                Quotient of The Virasoro vertex algebra at central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>-27/16*L_-6|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>,)
                sage: Q.module() is Q
                True

            """
            raise NotImplementedError("Not implemented")

        @abstract_method
        def zero(self):
            """
            The zero vector in this vertex algebra

            EXAMPLES::

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

    class ElementMethods:
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
                return self.bracket(rhs).get(n,self.parent().zero())
            else:
                return factorial(-1-n)**(-1)*self.T(-n-1)._mul_(rhs)


    class SubcategoryMethods:
        def Graded(self):
            """
            The subcategory of H-Graded vertex algebras
            """
            return self.HGraded()

        def HGraded(self):
            """
            The subcategory of H-Graded vertex algebras
            """
            return self._with_axiom('HGraded')

        def FinitelyGeneratedAsVertexAlgebra(self):
            """
            The subcategory of finitely and strongly generated vertex algebras
            """
            return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        def FinitelyGenerated(self):
            """
            The subcategory of finitely and strongly generated vertex algebras
            """
            return self.FinitelyGeneratedAsVertexAlgebra()

        def WithBasis(self):
            """
            The subcategory of vertex algebras with a preferred basis
            """
            return self._with_axiom("WithBasis")

    class WithBasis(CategoryWithAxiom_over_base_ring):
        def _repr_object_names(self):
            """
            The names of objects in this category
            """
            return "vertex algebras with basis over {}".format(self.base_ring())

        class ElementMethods:
            def monomials(self):
                """
                The tuple of monomials in this element

                EXAMPLES::

                    sage: V = VirasoroVertexAlgebra(QQ,1/2); L = V.0; w = (L*L)*L; 
                    sage: w.monomials()
                    (4*L_-4L_-2|0>, 2*L_-3L_-3|0>, 1/2*L_-6|0>, L_-2L_-2L_-2|0>)

                """
                return tuple(v[1]*v[0] for v in 
                            self.monomial_coefficients().items())
 

    class HGraded(CategoryWithAxiom_over_base_ring):

        def _repr_object_names(self):
            """
            The names of objects in this category
            """
            return "H-graded vertex algebras over {}".format(self.base_ring())

        class SubcategoryMethods:
            def FinitelyGenerated(self):
                """
                The subcategory of finitely and strongly generated H-Graded
                vertex algebras
                """
                return self._with_axiom("FinitelyGeneratedAsVertexAlgebra")

        class WithBasis(CategoryWithAxiom_over_base_ring):
            """
            The subcategory of H-Graded vertex algebras with a preferred
            basis
            """
            def _repr_object_names(self):
                """
                The names of objects in this category
                """
                return "H-graded vertex algebras with basis "\
                            "over {}".format(self.base_ring())


        class ParentMethods:
            def classical_limit(self):
                """
                The Poisson vertex algebra classical limit of this vertex algebra

                EXAMPLES:

                We construct the classical limit of the universal Virasoro vertex
                algebra at central charge `1/2`::

                    sage: V = VirasoroVertexAlgebra(QQ, 1/2)
                    sage: P = V.classical_limit()
                    sage: V.inject_variables()
                    Defining L
                    sage: (L*L)*L == L*(L*L)
                    False
                    sage: (P(L)*P(L))*P(L) == P(L)*(P(L)*P(L))
                    True
                    sage: L.bracket(L)
                    {0: TL, 1: 2*L, 3: 1/2*C}
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
                    33/8*L_-4L_-2|0>+27/16*L_-6|0>-93/64*L_-3L_-3|0>
                    sage: P(L)*(P(L)*P(L)) == P.zero()
                    True

                """
                from sage.algebras.vertex_algebras.poisson_vertex_algebra import PoissonVertexAlgebra
                return PoissonVertexAlgebra(self.base_ring(), self)

        class ElementMethods:
            #for compatibility with LCA as `self` is in LCA
            def degree(self):
                """
                The conformal weight of this element

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

            def filtered_degree(self):
                """
                The smallest space `F^p` in the Li filtration of this vertex
                algebra containing this element

                EXAMPLES::

                    sage: L.li_filtration_degree()
                    0
                    sage: (L.T(2)*L.T()).li_filtration_degree()
                    3

                """
                return max(m.weight() for m in self.monomial_coefficients())

            @abstract_method
            def weight(self):
                """
                The conformal weight of this element

                This method is an alias of :meth:`degree`
                """

            def is_homogeneous(self):
                """
                Whether this element is homogeneous with respect to conformal
                weight

                EXAMPLES::

                    sage: L.is_homogeneous()
                    True
                    sage: (L + L.T()).is_homogeneous()
                    False

                """
                try:
                    self.weight()
                except ValueError:
                    return False
                return True

            def nmodeproduct(self, other, n):
                try:
                    weight = self.weight()
                except ValueError:
                    raise ValueError("Couldn't compute weight of {}, "\
                                    "it's not homogeneous?".format(self))
                return self.nproduct(other, n+weight-1)

    class FinitelyGeneratedAsVertexAlgebra(CategoryWithAxiom_over_base_ring):

        def _repr_object_names(self):
            return "finitely and strongly generated" \
                        " vertex algebras over {}".format(self.base_ring())
        
        class ParentMethods:
            @abstract_method
            def gens(self):
                return

            @abstract_method
            def ngens(self):
                return

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

                    return all (p(g).nmodeproduct(self,n).is_zero() for 
                                n in range(1,weight+2) for g in p.gens())

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


