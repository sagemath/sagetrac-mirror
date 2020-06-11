"""
Universal Enveloping Vertex Algebra

AUTHORS:

- Reimundo Heluani (08-09-2019): Initial implementation.
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


from .vertex_algebra import VertexAlgebra
from sage.categories.lie_conformal_algebras import LieConformalAlgebras
from sage.categories.vertex_algebras import VertexAlgebras
from sage.sets.family import Family
from sage.combinat.free_module import CombinatorialFreeModule
from sage.misc.cachefunc import cached_method
from .vertex_algebra_element import UniversalEnvelopingVertexAlgebraElement
from sage.rings.rational_field import QQ
from .energy_partition_tuples import EnergyPartitionTuples
from sage.functions.other import floor


class UniversalEnvelopingVertexAlgebra(VertexAlgebra,CombinatorialFreeModule):

    def __init__(self, R, L, category=None,
                 central_parameters=None,
                 names=None, latex_names=None):
        """
        The (central quotient of the) universal enveloping vertex
        algebra of the Lie conformal algebra `L` over the ring `R`.

        INPUT:

        - ``R`` -- a commutative ring; the base ring of this vertex
          algebra. Undefined
          behaviour if this ring is not a field of characteristic zero.

        - ``L`` -- a Lie conformal algebra.

        - ``central_parameters`` -- a finite family parametrized by
          central elements of this vertex algebra (default: ``0`` for
          each central element of ``L``);
          a family describing the action of the central
          elements in this vertex algebra.
        """
        if L not in LieConformalAlgebras(R).WithBasis().FinitelyGenerated():
            raise ValueError ( "L needs to be a finitely generated " \
                "Lie conformal algebra with basis, got {}".format(L) )

        category = VertexAlgebras(R).FinitelyGenerated().WithBasis().\
           or_subcategory(category)

        if L in LieConformalAlgebras(R).Graded():
            category = category.Graded()

        if L in LieConformalAlgebras(R).Super():
            category = category.Super()

        VertexAlgebra.__init__(self, R, category=category, names=names,
                                latex_names=latex_names)

        self._lca = L
        if not central_parameters:
            central_parameters = Family({i:0  for i in L.central_elements()})
        if set(central_parameters.keys()) != set(L.central_elements()):
            raise ValueError ("central_parameters must be parametrized by "\
                              "central elements")

        self._central_parameters = central_parameters
        self._ngens = L.ngens() - central_parameters.cardinality()

        #need to call directly this because of 1 generator.
        #Also:self._module is needed for self.ngens()
        regular = tuple([2*g.is_even_odd() for g in L.gens()\
                    if g not in L.central_elements()])
        if self.is_graded():
            weights = tuple([g.degree() for g in L.gens() if g not in\
                             L.central_elements()])
            basis = EnergyPartitionTuples(weights,self._ngens,regular=regular)
        else:
            from sage.combinat.partition_tuple import PartitionTuples
            basis = PartitionTuples(level=self._ngens,regular=regular)
        CombinatorialFreeModule.__init__(self, R, basis_keys=basis,
                        category=category,
                        element_class=UniversalEnvelopingVertexAlgebraElement) 
        self.register_lift()


    def _repr_(self):
        """
        The name of this vertex algebra.

        EXAMPLES::

            sage: L = NeveuSchwarzLieConformalAlgebra(QQ)
            sage: L.universal_enveloping_algebra()
            The universal enveloping vertex algebra of the Neveu-Schwarz super Lie conformal algebra over Rational Field
        """
        lcafirst,lcaleft = format(self._lca).split(' ',1)
        if lcafirst == "The":
            return "The universal enveloping vertex algebra of the " + lcaleft
        return "The universal enveloping vertex algebra of " + lcafirst +\
                " " + lcaleft

    def register_lift(self):
        """
        Register a new coercion from its Lie conformal algebra.

        If this vertex algebra is the universal enveloping vertex
        algebra of the Lie conformal algebra `L`, register a new
        coercion from `L`. 

        .. SEEALSO::

            :meth:`sage.algebras.lie_conformal_algebras.lie_conformal_algebra_element.lift`

        EXAMPLES::

            sage: L = VirasoroLieConformalAlgebra(QQ); 
            sage: V = L.universal_enveloping_algebra()
            sage: L.lift
            Generic morphism:
              From: The Virasoro Lie conformal algebra over Rational Field
              To:   The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
            sage: W = VirasoroVertexAlgebra(QQ,1/2)
            sage: L.lift
            Generic morphism:
              From: The Virasoro Lie conformal algebra over Rational Field
              To:   The Virasoro vertex algebra at central charge 1/2
            sage: V.register_lift()
            sage: L.lift
            Generic morphism:
              From: The Virasoro Lie conformal algebra over Rational Field
              To:   The universal enveloping vertex algebra of the Virasoro Lie conformal algebra over Rational Field
        """
        from sage.categories.homset import Hom
        from sage.algebras.lie_conformal_algebras.\
            lie_conformal_algebra_with_structure_coefs import _LiftMorphism
        newlift = _LiftMorphism(Hom(self._lca, self,
                        category = LieConformalAlgebras(self._lca.base_ring())))
        self._lca.set_lift(newlift)

    def gens(self):
        """
        The generators of this vertex algebra.

        EXAMPLES::

            sage: V = NeveuSchwarzVertexAlgebra(QQ,1); V.gens()
            (L_-2|0>, G_-3/2|0>)
            sage: V = AffineVertexAlgebra(QQ,'A1', 1, names =('e','h', 'f')); V.gens()
            (e_-1|0>, h_-1|0>, f_-1|0>)
            sage: V = AffineVertexAlgebra(QQ,'A1', 1); V.gens()
            (alpha[1]_-1|0>, alphacheck[1]_-1|0>, -alpha[1]_-1|0>)
        """
        return tuple([self.gen(i) for i in range(self._ngens)])

    def ngens(self):
        """
        The number of generators of this vertex algebra.

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'B3', 1); V.ngens()
            21
            sage: V = N2VertexAlgebra(QQ,1); V.ngens()
            4
        """
        return self._ngens

    def gen(self,i):
        r"""
        The `i`-th generator of this vertex algebra.

        EXAMPLES::

            sage: V = NeveuSchwarzVertexAlgebra(QQ); V.gen(1)
            G_-3/2|0>
        """
        l = [[]]*self._ngens
        l[i] = [1]
        return self(l)

    def central_parameters(self):
        """
        The central character used to construct this universal
        enveloping vertex algebra. 

        EXAMPLES::

            sage: V = FreeFermionsVertexAlgebra(QQ); V.central_parameters()
            Finite family {K: 1}
        """
        return self._central_parameters
    
    def get_weight_less_than(self,n):
        """
        The sub-vector space of this vertex algebra of vectors with
        conformal weight less than or equal to ``n``.

        INPUT:

        - ``n`` -- a non-negative rational number; 

        OUTPUT: a submodule of this vertex algebra. 

        EXAMPLES::

            sage: V = FreeFermionsVertexAlgebra(QQ); M = V.get_weight_less_than(5/2); M
            Free module generated by {0, 1, 2, 3} over Rational Field
            sage: [v.lift() for v in M.basis()]
            [|0>, psi_-1/2|0>, psi_-3/2|0>, psi_-3/2psi_-1/2|0>]

        TESTS::
           
            sage: V = FreeFermionsVertexAlgebra(QQ); M = V.get_weight_less_than(5/2);
            sage: M.reduce(V.vacuum())
            0  
        """
        if not self.is_graded() or any(g.weight() not in QQ or g.weight == 0\
                                                        for g in self.gens()):
            raise NotImplementedError("get_weight_less_than is not implemented"\
                                      " for {}".format(self))

        if n not in QQ or n < 0:
            raise ValueError("n needs to be a non-negative rational number")
         
        weights = tuple([g.degree() for g in self._lca.gens() if g not in\
                         self._lca.central_elements()])
        regular = tuple([2*g.is_even_odd() for g in self._lca.gens()\
                    if g not in self._lca.central_elements()])
        basis = []
        for i in EnergyPartitionTuples(weights,self._ngens,regular=regular):
            if i.energy() < n:
                basis.append(self(i))
            else:
                break
        return self.submodule(basis)

    def get_weight(self,n):
        r"""
        The sub-vector space of this vertex algebra of vectors with
        conformal weight equal to ``n``.

        INPUT:

        - ``n`` -- a non-negative rational number; 

        OUTPUT: a submodule of this vertex algebra. 

        EXAMPLES::

            sage: V = NeveuSchwarzVertexAlgebra(QQ, 1/2); M = V.get_weight(7/2)
            sage: [v.lift() for v in M.basis()]
            [G_-7/2|0>, L_-2G_-3/2|0>]
        """
        if not self.is_graded() or any(g.weight() not in QQ or g.weight == 0\
                                                        for g in self.gens()):
            raise NotImplementedError("get_weight is not implemented for {}"\
                                        .format(self))

        if n not in QQ or n < 0:
            raise ValueError("n needs to be a non-negative rational number")
        
        weights = tuple([g.degree() for g in self._lca.gens() if g not in\
                         self._lca.central_elements()])
        regular = tuple([2*g.is_even_odd() for g in self._lca.gens()\
                    if g not in self._lca.central_elements()])
        basis = EnergyPartitionTuples(weights,self._ngens,n,regular=regular)
        return self.submodule([self(v) for v in basis])

    def dimension(self,n):
        """
        The dimension of the degree `n` part of this vertex algebra.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.dimension(4)
            2
            sage: V = AffineVertexAlgebra(QQ,'A1', 1); V.dimension(5)
            108
            sage: V = FreeFermionsVertexAlgebra(QQ); V.dimension(13/2)
            6
        """
        if not self.is_graded() or any(g.weight() not in QQ or g.weight == 0\
                                                        for g in self.gens()):
            raise NotImplementedError("dimension is not implemented for {}"\
                                        .format(self))

        if n not in QQ or n < 0:
            raise ValueError("n must be a non-negative rational number")

        weights = tuple([g.degree() for g in self._lca.gens() if g not in\
                         self._lca.central_elements()])
        regular = tuple([2*g.is_even_odd() for g in self._lca.gens()\
                    if g not in self._lca.central_elements()])
        return EnergyPartitionTuples(weights,self._ngens,n,regular=regular).\
                                                                cardinality()

    @cached_method
    def vacuum(self):
        """
        The vacuum vector of this vertex algebra

        EXAMPLES::

            sage: VirasoroVertexAlgebra(QQ,1).vacuum()
            |0>
        """
        vac = [[],]*self.ngens()
        return self(vac)

    def find_singular(self,n):
        """
        Return a basis of the vector space of singular vectors of
        weight `n`.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.find_singular(6)
            (L_-2L_-2L_-2|0> + 93/64*L_-3L_-3|0> - 33/8*L_-4L_-2|0> - 27/16*L_-6|0>,)
            sage: V = AffineVertexAlgebra(QQ,'A1',1,names=('e', 'h', 'f')); V
            The universal affine vertex algebra of CartanType ['A', 1] at level 1
            sage: V.find_singular(2)
            (f_-1f_-1|0>,
             f_-2|0> + h_-1f_-1|0>,
             h_-1h_-1|0> + h_-2|0> - 2*e_-1f_-1|0>,
             e_-1h_-1|0> + e_-2|0>,
             e_-1e_-1|0>)

        TESTS::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.find_singular(0)
            (|0>,)
        """
        if not self.is_graded() or any(g.weight() not in QQ or g.weight == 0\
                                                        for g in self.gens()):
            raise NotImplementedError("find_singular is not implemented"\
                                      " for {}".format(self))

        if n not in QQ or n < 0:
            raise ValueError("n needs to be a non-negative rational number")
        M = self.get_weight(n)
        B = M.basis()
        N = self.get_weight_less_than(n)
        W = CombinatorialFreeModule(self.base_ring(),self.gens())
        from sage.categories.tensor import tensor
        Z = tensor([N,W])
        br = {(g,v): g.bracket(v.lift()) for v in B for g in self.gens()}
        on_basis = lambda v: sum(tensor([N.retract(sum(br[(g,B[v])][k] for k in\
                br[(g,B[v])] if k > g.weight()-1)),W(g)]) for g in self.gens())
        f = M.module_morphism(on_basis, codomain=Z)
        return tuple([v.lift() for v in f.kernel_basis()])

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
            33/8*L_-4L_-2|0>-93/64*L_-3L_-3|0>+27/16*L_-6|0>
            sage: P(L)*(P(L)*P(L)) == P.zero()
            True

        """
        from .poisson_vertex_algebra import PoissonVertexAlgebra
        return PoissonVertexAlgebra(self.base_ring(), self)

    def li_filtration(self,n,k=None):
        r"""Let `V` be this vertex algebra and `V_n` its conformal weight `n`
        part. This method returns the filtered vector space `F_\bullet V_n` with respect
        to the Li filtration. If ``k`` is specified it returns `F_k V_n`.

        EXAMPLES::

            sage: V = AffineVertexAlgebra(QQ, 'A1', 1);
            sage: V.li_filtration(2)
            {0: Free module generated by {0, 1, 2, 3, 4, 5, 6, 7, 8} over Rational Field,
             1: Free module generated by {0, 1, 2} over Rational Field,
             2: Free module generated by {} over Rational Field,
             3: Free module generated by {} over Rational Field}
            sage: V.li_filtration(2,1)
            Free module generated by {0, 1, 2} over Rational Field
            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.dimension(12)
            21
            sage: V.li_filtration(12,7)
            Free module generated by {0, 1, 2, 3, 4, 5} over Rational Field

        """
        A = self.get_graded_part(n)
        ret = {}
        for m in range(n+2) if k==None else range(k,k+1):
            basis = [b for b in A.basis() if self._from_dict(
                b.monomial_coefficients()).li_filtration_degree()>=m]
            ret[m] = A.submodule(basis)
        return ret if k==None else ret[k]

    def ideal(self, gens, check=True):
        """
        The ideal of this vertex algebra generated by ``gens``.

        INPUT:

        - ``gens`` -- a list or tuple of elements of this vertex
        algebra.

        - ``check`` -- a boolean (default: ``True``); whether to
        check that the generators are singular vectors.

        EXAMPLES:

        We construct the ideal defining the *Virasoro Ising* module::

            sage: V = VirasoroVertexAlgebra(QQ,1/2)
            sage: L = V.0
            sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
            sage: I = V.ideal(v)
            sage: I
            ideal of The Virasoro vertex algebra at central charge 1/2 generated by (L_-2L_-2L_-2|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>-27/16*L_-6|0>,)

        If we instead use a non-singular vector::

            sage: V.ideal(L*L)
            Traceback (most recent call last):
            ...
            ValueError: Generators must be singular vectors of The Virasoro vertex algebra at central charge 1/2

        NOTE::

        We only implement ideals of universal enveloping vertex
        algebras and their quotients, generated by singular
        vectors.
        """
        from .vertex_algebra_ideal import VertexAlgebraIdeal
        return VertexAlgebraIdeal(self, gens, check=check)

    def quotient(self, I):
        """
        The quotient of this vertex algebra by the ideal ``I``.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2)
            sage: L = V.0
            sage: v = L*(L*L) + 93/64*L.T()*L.T() - 33/16*L.T(2)*L - 9/128*L.T(4)
            sage: I = V.ideal(v)
            sage: Q = V.quotient(I); Q
            Quotient of The Virasoro vertex algebra at central charge 1/2 by the ideal generated by (L_-2L_-2L_-2|0>+93/64*L_-3L_-3|0>-33/8*L_-4L_-2|0>-27/16*L_-6|0>,)
            sage: Q(L*(L*L))
            33/8*L_-4L_-2|0>-93/64*L_-3L_-3|0>+27/16*L_-6|0>
        """
        from .vertex_algebra_quotient import VertexAlgebraQuotient
        return VertexAlgebraQuotient(I)

    def arc_algebra(self, termorder='wdegrevlex'):
        r"""
        The algebra of functions of the arc space of the `C_2`
        quotient of this Vertex algebra.


        INPUT:

        - ``termorder`` a string (default: ``'wdegrevlex'``); the
          monomial ordering of the algebra.

        OUTPUT: The graded Poisson vertex algebra freely generated
        as a differential algebra by the `C_2` quotient of this
        vertex algebra.

        TODO: we only support arc algebras of universal enveloping
        vertex algebras and their quotients.

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ, 1/2); Q=V.quotient(V.ideal(V.find_singular(6)[0]))
            sage: Q.arc_algebra()
            Quotient of The arc algebra over Rational Field generated by ('L',) by the differential ideal generated by (L_2^3,)
            sage: V.arc_space()
            The arc algebra over Rational Field generated by ('L',)
        """
        from sage.algebras.vertex_algebras.poisson_vertex_algebra \
                import VertexAlgebraArcAlgebra
        return VertexAlgebraArcAlgebra(self, termorder)



