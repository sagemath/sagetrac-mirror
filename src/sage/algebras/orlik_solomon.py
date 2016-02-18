r"""
Orlik-Solomon Algebras
"""

#*****************************************************************************
#       Copyright (C) 2015 William Slofstra
#                          Travis Scrimshaw <tscrimsh at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.algebras import Algebras
from sage.sets.family import Family

class OrlikSolomonAlgebra(CombinatorialFreeModule):
    r"""
    An Orlik-Solomon algebra.

    Let `R` be a commutative ring. Let `M` be a matroid with ground set
    `X`. Let `C(M)` denote the set of circuits of `M`. Let `E` denote
    the exterior algera over `R` generated by `\{ e_x \mid x \in X \}`.
    The *Orlik-Solomon ideal* `J(M)` is the ideal of `E` generated by

    .. MATH::

        \partial e_S := \sum_{i=1}^t (-1)^{i-1} e_{j_1} \wedge e_{j_2}
        \wedge \cdots \wedge \widehat{e}_{j_i} \wedge \cdots \wedge e_{j_t}

    for all `S = \left\{ j_1 < j_2 < \cdots < j_t \right\} \in C(M)`,
    where `\widehat{e}_{j_i}` means that the term `e_{j_i}` is being
    omitted. The notation `\partial e_S` is not a coincidence, as
    `\partial e_S` is actually the image of
    `e_S := e_{j_1} \wedge e_{j_2} \wedge \cdots \wedge e_{j_t}` under the
    unique derivation `\partial` of `E` which sends all `e_x` to `1`.

    It is easy to see that `\partial e_S \in J(M)` not only for circuits
    `S`, but also for any dependent set `S` of `M`. Moreover, every
    dependent set `S` of `M` satisfies `e_S \in J(M)`.

    The *Orlik-Solomon algebra* `A(M)` is the quotient `E / J(M)`. This is
    a graded finite-dimensional skew-commutative `R`-algebra. Fix
    some ordering on `X`; then, the NBC sets of `M` (that is, the subsets
    of `X` containing no broken circuit of `M`) form a basis of `A(M)`.
    (Here, a *broken circuit* of `M` is defined to be the result of
    removing the smallest element from a circuit of `M`.)

    In the current implementation, the basis of `A(M)` is indexed by the
    NBC sets, which are implemented as frozensets.

    INPUT:

    - ``R`` -- the base ring
    - ``M`` -- the defining matroid
    - ``ordering`` -- (optional) an ordering of the ground set

    EXAMPLES:

    We create the Orlik-Solomon algebra of the uniform matroid `U(3, 4)`
    and do some basic computations::

        sage: M = matroids.Uniform(3, 4)
        sage: OS = M.orlik_solomon_algebra(QQ)
        sage: OS.dimension()
        14
        sage: G = OS.algebra_generators()
        sage: M.broken_circuits()
        frozenset({frozenset({1, 2, 3})})
        sage: G[1] * G[2] * G[3]
        OS{0, 1, 2} - OS{0, 1, 3} + OS{0, 2, 3}

    REFERENCES:

    .. [CE01] Raul Cordovil and Gwihen Etienne.
       *A note on the Orlik-Solomon algebra*.
       Europ. J. Combinatorics. **22** (2001). pp. 165-170.
       http://www.math.ist.utl.pt/~rcordov/Ce.pdf

    - :wikipedia:`Arrangement_of_hyperplanes#The_Orlik-Solomon_algebra`
    """
    @staticmethod
    def __classcall_private__(cls, R, M, ordering=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: from sage.algebras.orlik_solomon import OrlikSolomonAlgebra
            sage: OS1 = OrlikSolomonAlgebra(QQ, M)
            sage: OS2 = OrlikSolomonAlgebra(QQ, M, ordering=(0,1,2,3,4,5))
            sage: OS3 = OrlikSolomonAlgebra(QQ, M, ordering=[0,1,2,3,4,5])
            sage: OS1 is OS2 and OS2 is OS3
            True
        """
        if ordering is None:
            ordering = sorted(M.groundset())
        return super(OrlikSolomonAlgebra, cls).__classcall__(cls, R, M, tuple(ordering))

    def __init__(self, R, M, ordering=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: TestSuite(OS).run()

        We check on the matroid associated to the graph with 3 vertices and
        2 edges between each vertex::

            sage: G = Graph([[1,2],[1,2],[2,3],[2,3],[1,3],[1,3]], multiedges=True)
            sage: M = Matroid(G)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: elts = OS.some_elements() + list(OS.algebra_generators())
            sage: TestSuite(OS).run(elements=elts)
        """
        self._M = M
        self._sorting = {x:i for i,x in enumerate(ordering)}

        # set up the dictionary of broken circuits
        self._broken_circuits = dict()
        for c in self._M.circuits():
            L = sorted(c, key=lambda x: self._sorting[x])
            self._broken_circuits[frozenset(L[1:])] = L[0]

        cat = Algebras(R).FiniteDimensional().WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, R, M.no_broken_circuits_sets(ordering),
                                         prefix='OS', bracket='{',
                                         generator_cmp=self._cmp_term,
                                         category=cat)

    def _cmp_term(self, x, y):
        """
        Compare the terms indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS._cmp_term(frozenset({1, 2}), frozenset({1, 4}))
            -1
            sage: OS._cmp_term(frozenset({0, 1, 2}), frozenset({1, 4}))
            -1
            sage: OS._cmp_term(frozenset({}), frozenset({1, 4}))
            1
        """
        c = cmp(len(x), len(y))
        if c != 0:
            return -c
        return cmp(sorted(x), sorted(y))

    def _repr_term(self, m):
        """
        Return a string representation of the basis element indexed by `m`.

        EXAMPLES::

            sage: M = matroids.Uniform(3, 4)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS._repr_term(frozenset([0]))
            'OS{0}'
        """
        return "OS{{{}}}".format(str(list(m))[1:-1])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: M.orlik_solomon_algebra(QQ)
            Orlik-Solomon algebra of Wheel(3): Regular matroid of rank 3
             on 6 elements with 16 bases
        """
        return "Orlik-Solomon algebra of {}".format(self._M)

    @cached_method
    def one_basis(self):
        """
        Return the index of the basis element corresponding to `1`
        in ``self``.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.one_basis() == frozenset([])
            True
        """
        return frozenset({})

    @cached_method
    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.

        These form a family indexed by the ground set `X` of `M`. For
        each `x \in X`, the `x`-th element is `e_x`.

        EXAMPLES::

            sage: M = matroids.Uniform(2, 2)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.algebra_generators()
            Finite family {0: OS{0}, 1: OS{1}}

            sage: M = matroids.Uniform(1, 2)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.algebra_generators()
            Finite family {0: OS{0}, 1: OS{0}}

            sage: M = matroids.Uniform(1, 3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.algebra_generators()
            Finite family {0: OS{0}, 1: OS{0}, 2: OS{0}}
        """
        return Family(sorted(self._M.groundset()),
                      lambda i: self.subset_image(frozenset([i])))

    @cached_method
    def product_on_basis(self, a, b):
        """
        Return the product in ``self`` of the basis elements
        indexed by ``a`` and ``b``.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.product_on_basis(frozenset([2]), frozenset([3,4]))
            OS{0, 1, 2} - OS{0, 1, 4} + OS{0, 2, 3} + OS{0, 3, 4}

        ::

            sage: G = OS.algebra_generators()
            sage: prod(G)
            0
            sage: G[2] * G[4]
            -OS{1, 2} + OS{1, 4}
            sage: G[3] * G[4] * G[2]
            OS{0, 1, 2} - OS{0, 1, 4} + OS{0, 2, 3} + OS{0, 3, 4}
            sage: G[2] * G[3] * G[4]
            OS{0, 1, 2} - OS{0, 1, 4} + OS{0, 2, 3} + OS{0, 3, 4}
            sage: G[3] * G[2] * G[4]
            -OS{0, 1, 2} + OS{0, 1, 4} - OS{0, 2, 3} - OS{0, 3, 4}

        TESTS:

        Let us check that `e_{s_1} e_{s_2} \cdots e_{s_k} = e_S` for any
        subset `S = \{ s_1 < s_2 < \cdots < s_k \}` of the ground set::

            sage: G = Graph([[1,2],[1,2],[2,3],[3,4],[4,2]], multiedges=True)
            sage: M = Matroid(G)
            sage: E = M.groundset_list()
            sage: OS = M.orlik_solomon_algebra(ZZ)
            sage: G = OS.algebra_generators()
            sage: import itertools
            sage: def test_prod(F):
            ....:     LHS = OS.subset_image(frozenset(F))
            ....:     RHS = OS.prod([G[i] for i in sorted(F)])
            ....:     return LHS == RHS
            sage: all( test_prod(F) for k in range(len(E)+1)
            ....:                   for F in itertools.combinations(E, k) )
            True
        """
        if not a:
            return self.basis()[b]
        if not b:
            return self.basis()[a]

        if not a.isdisjoint(b):
            return self.zero()

        R = self.base_ring()
        # since a is disjoint from b, we can just multiply the generator
        if len(a) == 1:
            i = list(a)[0]
            # insert i into nbc, keeping track of sign in coeff
            ns = b.union({i})
            ns_sorted = sorted(ns, key=lambda x: self._sorting[x])
            coeff = (-1)**ns_sorted.index(i)

            return R(coeff) * self.subset_image(ns)

        # r is the accumalator
        # we reverse a in the product, so add a sign
        # note that l>=2 here
        if len(a) % 4 < 2:
            sign = R.one()
        else:
            sign = - R.one()
        r = self._from_dict({b: sign}, remove_zeros=False)

        # now do the multiplication generator by generator
        G = self.algebra_generators()
        for i in sorted(a, key=lambda x: self._sorting[x]):
            r = G[i] * r

        return r

    @cached_method
    def subset_image(self, S):
        """
        Return the element `e_S` of `A(M)` (``== self``) corresponding to
        a subset `S` of the ground set of `M`.

        INPUT:

        - ``S`` -- a frozenset which is a subset of the ground set of `M`

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: BC = sorted(M.broken_circuits(), key=sorted)
            sage: for bc in BC: (sorted(bc), OS.subset_image(bc))
            ([1, 3], -OS{0, 1} + OS{0, 3})
            ([1, 4, 5], OS{0, 1, 4} - OS{0, 1, 5} - OS{0, 3, 4} + OS{0, 3, 5})
            ([2, 3, 4], OS{0, 1, 2} - OS{0, 1, 4} + OS{0, 2, 3} + OS{0, 3, 4})
            ([2, 3, 5], OS{0, 2, 3} + OS{0, 3, 5})
            ([2, 4], -OS{1, 2} + OS{1, 4})
            ([2, 5], -OS{0, 2} + OS{0, 5})
            ([4, 5], -OS{3, 4} + OS{3, 5})

            sage: M4 = matroids.CompleteGraphic(4)
            sage: OS = M4.orlik_solomon_algebra(QQ)
            sage: OS.subset_image(frozenset({2,3,4}))
            OS{0, 2, 3} + OS{0, 3, 4}

        An example of a custom ordering::

            sage: G = Graph([[3, 4], [4, 1], [1, 2], [2, 3], [3, 5], [5, 6], [6, 3]])
            sage: M = Matroid(G)
            sage: s = [(5, 6), (1, 2), (3, 5), (2, 3), (1, 4), (3, 6), (3, 4)]
            sage: sorted([sorted(c) for c in M.circuits()])
            [[(1, 2), (1, 4), (2, 3), (3, 4)],
             [(3, 5), (3, 6), (5, 6)]]
            sage: OS = M.orlik_solomon_algebra(QQ, ordering=s)
            sage: OS.subset_image(frozenset([]))
            OS{}
            sage: OS.subset_image(frozenset([(1,2),(3,4),(1,4),(2,3)]))
            0
            sage: OS.subset_image(frozenset([(2,3),(1,2),(3,4)]))
            OS{(1, 2), (3, 4), (2, 3)}
            sage: OS.subset_image(frozenset([(1,4),(3,4),(2,3),(3,6),(5,6)]))
            -OS{(1, 2), (5, 6), (2, 3), (1, 4), (3, 6)}
             + OS{(1, 2), (5, 6), (3, 4), (1, 4), (3, 6)}
             - OS{(1, 2), (5, 6), (3, 4), (2, 3), (3, 6)}
            sage: OS.subset_image(frozenset([(1,4),(3,4),(2,3),(3,6),(3,5)]))
            OS{(1, 2), (5, 6), (2, 3), (1, 4), (3, 5)}
             - OS{(1, 2), (5, 6), (2, 3), (1, 4), (3, 6)}
             + OS{(1, 2), (5, 6), (3, 4), (1, 4), (3, 5)}
             + OS{(1, 2), (5, 6), (3, 4), (1, 4), (3, 6)}
             - OS{(1, 2), (5, 6), (3, 4), (2, 3), (3, 5)}
             - OS{(1, 2), (5, 6), (3, 4), (2, 3), (3, 6)}

        TESTS::

            sage: G = Graph([[1,2],[1,2],[2,3],[2,3],[1,3],[1,3]], multiedges=True)
            sage: M = Matroid(G)
            sage: sorted([sorted(c) for c in M.circuits()])
            [[0, 1], [0, 2, 4], [0, 2, 5], [0, 3, 4],
             [0, 3, 5], [1, 2, 4], [1, 2, 5], [1, 3, 4],
             [1, 3, 5], [2, 3], [4, 5]]
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.subset_image(frozenset([]))
            OS{}
            sage: OS.subset_image(frozenset([1, 2, 3]))
            0
            sage: OS.subset_image(frozenset([1, 3, 5]))
            0
            sage: OS.subset_image(frozenset([1, 2]))
            OS{0, 2}
            sage: OS.subset_image(frozenset([3, 4]))
            -OS{0, 2} + OS{0, 4}
            sage: OS.subset_image(frozenset([1, 5]))
            OS{0, 4}

            sage: G = Graph([[1,2],[1,2],[2,3],[3,4],[4,2]], multiedges=True)
            sage: M = Matroid(G)
            sage: sorted([sorted(c) for c in M.circuits()])
            [[0, 1], [2, 3, 4]]
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.subset_image(frozenset([]))
            OS{}
            sage: OS.subset_image(frozenset([1, 3, 4]))
            -OS{0, 2, 3} + OS{0, 2, 4}

        We check on a non-standard ordering::

            sage: M = matroids.Wheel(3)
            sage: o = [5,4,3,2,1,0]
            sage: OS = M.orlik_solomon_algebra(QQ, ordering=o)
            sage: BC = sorted(M.broken_circuits(ordering=o), key=sorted)
            sage: for bc in BC: (sorted(bc), OS.subset_image(bc))
            ([0, 1], OS{0, 3} - OS{1, 3})
            ([0, 1, 4], OS{0, 3, 5} - OS{0, 4, 5} - OS{1, 3, 5} + OS{1, 4, 5})
            ([0, 2], OS{0, 5} - OS{2, 5})
            ([0, 2, 3], -OS{0, 3, 5} + OS{2, 3, 5})
            ([1, 2], OS{1, 4} - OS{2, 4})
            ([1, 2, 3], -OS{1, 3, 5} + OS{1, 4, 5} + OS{2, 3, 5} - OS{2, 4, 5})
            ([3, 4], OS{3, 5} - OS{4, 5})
        """
        if not isinstance(S, frozenset):
            raise ValueError("S needs to be a frozenset")
        for bc in self._broken_circuits:
            if bc.issubset(S):
                i = self._broken_circuits[bc]
                if i in S:
                    # ``S`` contains not just a broken circuit, but an
                    # actual circuit; then `e_S = 0`.
                    return self.zero()
                coeff = self.base_ring().one()
                # Now, reduce ``S``, and build the result ``r``:
                r = self.zero()
                switch = False
                Si = S.union({i})
                Ss = sorted(Si, key=lambda x: self._sorting[x])
                for j in Ss:
                    if j in bc:
                        r += coeff * self.subset_image(Si.difference({j}))
                    if switch:
                        coeff *= -1
                    if j == i:
                        switch = True
                return r
        else: # So ``S`` is an NBC set.
            return self.monomial(S)

    def degree_on_basis(self, m):
        """
        Return the degree of the basis element indexed by ``m``.

        EXAMPLES::

            sage: M = matroids.Wheel(3)
            sage: OS = M.orlik_solomon_algebra(QQ)
            sage: OS.degree_on_basis(frozenset([1]))
            1
            sage: OS.degree_on_basis(frozenset([0, 2, 3]))
            3
        """
        return len(m)

