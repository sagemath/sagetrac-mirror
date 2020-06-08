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
from sage.combinat.partition_tuple import PartitionTuples, PartitionTuples_level
from sage.combinat.partition import Partition
from .vertex_algebra_element import UniversalEnvelopingVertexAlgebraElement

class UniversalEnvelopingVertexAlgebra(VertexAlgebra,CombinatorialFreeModule):

    def __init__(self, R, L, category=None,
                 central_parameters=None,
                 names=None, latex_names=None):
        """
        The (central quotient of the) universal enveloping vertex
        algebra of the Lie conformal algebra `L` over the ring `R`.

        INPUT:

        - ``R`` a commutative ring; the base ring of this vertex
          algebra. Undefined
          behaviour if this ring is not a field of characteristic zero.

        - ``L`` a Lie conformal algebra.

        - ``central_parameters`` a finite family parametrized by
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
        if central_parameters:
            cp = central_parameters
        else:
            cp = Family({i:0  for i in L.central_elements()})
        if set(cp.keys()) != set(L.central_elements()):
            raise ValueError ("central_parameters must be parametrized by "\
                              "central elements")

        self._central_parameters = cp
        self._ngens = L.ngens() - cp.cardinality()

        #need to call directly this because of 1 generator.
        #Also:self._module is needed for self.ngens()
        regular = tuple([2*g.is_even_odd() for g in L.gens()\
                    if g not in L.central_elements()])
        basis = PartitionTuples(level=self._ngens,regular=regular)
        CombinatorialFreeModule.__init__(self, R, basis_keys=basis,
                        category=category,
                        element_class=UniversalEnvelopingVertexAlgebraElement) 
        self.register_lift()


    def _repr_(self):
        lcafirst,lcaleft = format(self._lca).split(' ',1)
        if lcafirst == "The":
            return "The universal enveloping vertex algebra of the " + lcaleft
        return "The universal enveloping vertex algebra of " + lcafirst +\
                " " + lcaleft

    def register_lift(self):
        from sage.categories.homset import Hom
        from sage.algebras.lie_conformal_algebras.\
            lie_conformal_algebra_with_structure_coefs import _LiftMorphism
        newlift = _LiftMorphism(Hom(self._lca, self,
                        category = LieConformalAlgebras(self._lca.base_ring())))
        self._lca.set_lift(newlift)

    def gens(self):
        """The generators of this vertex algebra"""
        return tuple([self.gen(i) for i in range(self._ngens)])

    def ngens(self):
        """
        The number of generators of this vertex algebra
        """
        return self._ngens

    def gen(self,i):
        r"""The `i`-th generator of this vertex algebra"""
        l = [[]]*self._ngens
        l[i] = [1]
        return self(l)

    def central_elements(self):
        r"""If this vertex algebra is the universal enveloping vertex algebra of
        the Lie conformal algebra `L`. This method returns the family of central
        elements of the `L`."""
        return tuple([i.lift() for i in self._lca.central_elements()])

    def central_parameters(self):
        """Return the central character used to construct this universal
        enveloping vertex algebra"""
        return self._central_parameters

    def get_weight(self,n):
        r"""
        return the degree n filtered part of `self`. This is a
        `CombinatorialFreeModule` with the same dictionary keys as
        `self.module()`.
        """
        if not self.is_graded() or any(g.weight() == 0 for g in self.gens()):
            raise NotImplementedError("get_weight is not implemented for {}"\
                                        .format(self))

        if n == 0:
            self.submodule(self.vacuum())
        else:
            basis = [PartitionTuples_level(self.ngens())(p) for m in range(1 ,n+1)
            for p in PartitionTuples(self.ngens(),m) if
            self(PartitionTuples_level(self.ngens())(p)).weight() == n ]
        return CombinatorialFreeModule(self.base_ring(), basis)

    def dimension(self,n):
        r"""The dimension of the degree `n` part of this vertex algebra

        EXAMPLES::

            sage: V = VirasoroVertexAlgebra(QQ,1/2); V.dimension(4)
            2
            sage: V.dimension(6)
            4
            sage: V.dimension(0)
            1
            sage: V.dimension(1)
            0
            sage: V = AffineVertexAlgebra(QQ, 'A1', 1); V.dimension(1)
            3
            sage: V.dimension(2)
            9

        """

        return self.get_graded_part(n).dimension()

    @cached_method
    def vacuum(self):
        """The vacuum vector of this vertex algebra"""
        vac = [Partition([]),]*self.ngens()
        return self(vac)

    def find_singular(self,n):
        """
        Return the vector space of singular vectors of weight `n`
        """
        M = self.get_graded_part(n)
        B = M.basis()
        from sage.matrix.constructor import Matrix
        ret = Matrix(self.base_ring(),B.cardinality(),0,0)
        for g in self.gens():
            br = {v:g.bracket(self._from_dict(
                  v.monomial_coefficients())) for v in B}
            w = g.weight()
            for j in range(1,n+1):
                Mj = self.get_graded_part(n-j)
                ret = ret.augment(Matrix([Mj._from_dict(
                    br[v].get(j+w-1,self.zero()).value\
                    .monomial_coefficients()).to_vector()
                    for v in B ]))
        myker = ret.kernel().basis()
        return [self._from_dict(M.from_vector(v)\
            .monomial_coefficients()) for v in myker]

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

        def ideal(self, *gens, check=True):
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
            return VertexAlgebraIdeal(self,gens, check=check)

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



