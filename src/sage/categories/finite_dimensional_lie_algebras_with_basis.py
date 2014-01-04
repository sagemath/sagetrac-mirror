r"""
Finite Dimensional Lie Algebras With Basis

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras
from sage.rings.all import ZZ
from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import LieBracket
from sage.sets.family import Family
from sage.matrix.constructor import matrix

class FiniteDimensionalLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of finite dimensional Lie algebras with a basis.
    """
    _base_category_class_and_axiom = [LieAlgebras.FiniteDimensional, "WithBasis"]

    class ParentMethods:
        @cached_method
        def _construct_UEA(self):
            """
            Construct the universal enveloping algebra of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
                sage: L._construct_UEA()
                Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
            """
            # Create the UEA relations
            # We need to get names for the basis elements, not just the generators
            names = self.variable_names()
            F = FreeAlgebra(self.base_ring(), len(names), names)
            gens = F.gens()
            d = F.gens_dict()
            rels = {}
            S = self.structure_coefficients()
            for k in S.keys():
                g0 = d[k._left._name]
                g1 = d[k._right._name]
                if g0 < g1:
                    rels[g1*g0] = g0*g1 - sum(val*d[g._name] for g, val in S[k])
                else:
                    rels[g0*g1] = g1*g0 + sum(val*d[g._name] for g, val in S[k])
            return F.g_algebra(rels)

        def ambient_module(self):
            """
            Return the ambient module of ``self``.
            """
            return FreeModule(self.base_ring(), self.dimension())

        def killing_matrix(self, x, y):
            r"""
            Return the Killing matrix of ``x`` and ``y``.

            The Killing form is defined as the matrix corresponding to the
            action of `\mathrm{ad}_x \circ \mathrm{ad}_y` in the basis
            of ``self``.

            EXAMPLES::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.killing_matrix(x, y)
                [ 0  0]
                [-1  0]
            """
            return x.adjoint_matrix() * y.adjoint_matrix()

        def killing_form(self, x, y):
            r"""
            Return the Killing form on ``x`` and ``y``.

            The Killing form is defined as

            .. MATH::

                \langle x \mid y \rangle = \mathrm{tr}\left( \mathrm{ad}_x
                \circ \mathrm{ad}_y \right).
            """
            return self.killing_matrix(x, y).trace()

        @cached_method
        def killing_form_matrix(self):
            """
            Return the matrix of the Killing form of ``self``.
            """
            B = self.basis()
            m = matrix([[self.killing_form(x, y) for x in B] for y in B])
            m.set_immutable()
            return m

        @cached_method
        def structure_coefficients(self):
            """
            Return the non-trivial structure coefficients of ``self``.
            In particular, if `[x, y] = 0`, then we don't include it in the
            output.
            """
            d = {}
            B = list(self.basis())
            one = self.base_ring().one()
            by_basis = lambda a: (self._from_dict({a[0]: one}, remove_zeros=False), a[1])
            for i,x in enumerate(B):
                for y in B[i+1:]:
                    if self._basis_cmp(x, y) > 0:
                        d[LieBracket(y, x)] = map(by_basis, self.bracket(y, x))
                    else:
                        d[LieBracket(x, y)] = map(by_basis, self.bracket(x, y))
            return Family(d)

        def centralizer(self, S):
            """
            Return the centralizer of ``S`` in ``self``.
            """
            from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
            if isinstance(S, LieSubalgebra) or S is self:
                K = S
            else:
                K = self.subalgebra(S)

            m = K.basis_matrix()
            sc = {k: v.to_vector() for k,v in self.structure_coefficients().items()}
            X = self.basis()
            d = self.dimension()
            c_mat = matrix([[sum(r[j]*sc[x,X[j]][k] for j in range(d)) for x in X]
                            for r in m for k in range(d)])
            C = c_mat.right_kernel()
            return self.subalgebra(C) # TODO: convert C back to elements of ``self``

        def center(self):
            """
            Return the center of ``self``.
            """
            return self.centralizer(self)

        def normalizer(self, V):
            """
            Return the normalizer of ``V`` in ``self``.
            """
            from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
            if not isinstance(V, LieSubalgebra) and V is not self:
                V = self.subalgebra(V)

            m = V.basis_matrix()
            sc = {k: v.to_vector() for k,v in self.structure_coefficients().items()}
            X = self.basis()
            d = self.dimension()
            t = m.nrows()
            n_mat = matrix([[sum(r[j]*sc[x,X[j]][k] for j in range(d)) for x in X]
                            + [0]*(t*l) + [m[j][k] for j in range(t)] + [0]*(t*(t-l-1))
                            for l,r in enumerate(m.rows()) for k in range(d)])
            N = n_mat.right_kernel()
            # TODO: convert N back to elements of ``self`` by taking the first ``n`` coefficients
            return self.subalgebra(N)

        def derived_subalgebra(self):
            """
            Return the derived subalgebra of ``self``.

            EXAMPLES::
            """
            return self.product_space(self)
            #basis = self.basis()
            # We only need to do [a, b] when a < b (on some total order of the basis elements)
            # We echelonize the matrix here
            #b_mat = matrix(self.base_ring(), [self.bracket(a, b).to_vector()
            #                 for i,a in enumerate(basis) for b in basis[i+1:]])
            #b_mat.echelonize()
            #r = b_mat.rank()
            #names = self.variable_names()
            #elt = lambda row: {names[i]: v for i, v in enumerate(row) if v != 0}
            #gens = [self.element_class(self, elt(row)) for row in b_mat.rows()[:r]]
            #return self.subalgebra(gens)

        @cached_method
        def derived_series(self):
            """
            Return the derived series `(\mathfrak{g}_i)_i` of ``self`` where
            the rightmost `\mathfrak{g}_k = \mathfrak{g}_{k+1} = \cdots`.

            EXAMPLES::
            """
            L = [self]
            while L[-1].dimension() > 0:
                p = L[-1].derived_subalgebra()
                if L[-1].dimension() == p.dimension():
                    break
                L.append(p)
            return tuple(L)

        @cached_method
        def lower_central_series(self):
            """
            Return the lower central series `(\mathfrak{g}_i)_i` of ``self``
            where the rightmost `\mathfrak{g}_k = \mathfrak{g}_{k+1} = \cdots`.

            EXAMPLES::
            """
            L = [self]
            while L[-1].dimension() > 0:
                s = L[-1].product_space(L[-1], self)
                if L[-1].dimension() == s.dimension():
                    break
                L.append(s)
            return tuple(L)

        @cached_method
        def upper_central_series(self):
            """
            Return the upper central series `(C_i)_i` of ``self``
            where the rightmost `C_k = C_{k+1} = \cdots`.

            EXAMPLES::
            """
            C = [self.center()]
            d = self.dimension()
            while C[-1].dimension() < d:
                s = self.centralizer(self.quotient(C[-1]))
                s = C[-1].product_space(s, self)
                if C[-1].dimension() == s.dimension():
                    break
                C.append(s)
            return tuple(C)

        def hypercenter(self):
            """
            Return the hypercenter of ``self``.
            """
            return upper_central_series()[-1]

        def is_abelian(self):
            """
            Return if ``self`` is an abelian Lie algebra.

            EXAMPLES::
            """
            return len(self.structure_coefficients()) == 0

        def is_solvable(self):
            r"""
            A Lie algebra `\mathfrak{g}` is solvable if the derived series

            .. MATH::

                \mathfrak{g} \subseteq [\mathfrak{g}, \mathfrak{g}] \subseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], [\mathfrak{g},
                \mathfrak{g}] \bigr] \subseteq
                \biggl[ \bigl[ [\mathfrak{g}, \mathfrak{g}], [\mathfrak{g},
                \mathfrak{g}] \bigr], \bigl[ [\mathfrak{g}, \mathfrak{g}],
                [\mathfrak{g}, \mathfrak{g}] \bigr] \biggr] \subseteq \cdots

            eventually becomes 0.
            """
            return self.derived_series()[-1].dimension() == 0

        def is_nilpotent(self):
            r"""
            Return if ``self`` is a nilpotent Lie algebra.

            A Lie algebra `\mathfrak{g}` is nilpotent if the lower central
            series

            .. MATH::

                \mathfrak{g} \subseteq [\mathfrak{g}, \mathfrak{g}] \subseteq
                \bigl[ [\mathfrak{g}, \mathfrak{g}], \mathfrak{g} \bigr]
                \subseteq\biggl[\bigl[ [\mathfrak{g}, \mathfrak{g}],
                \mathfrak{g} \bigr], \mathfrak{g}\biggr] \subseteq \cdots

            eventually becomes 0.
            """
            return g.lower_central_series()[-1].dimension() != 0

        def is_semisimple(self):
            """
            Return if ``self`` if a semisimple Lie algebra.

            A Lie algebra is semisimple if the solvable radical is zero. This
            is equivalent to saying the Killing form is non-degenerate
            (in characteristic 0).
            """
            return not self.killing_form_matrix().is_singular()

        def is_isomorphic(self, g):
            """
            Check to see if ``self`` is isomorphic to ``g`` (as Lie algebras).

            It is clear that `\mathfrak{g}` and `\mathfrak{g}^{\prime}` are
            isomorphic only if the structure graphs in are isomorphic.
            Since we strip all singletons for simplicity, we also need to check
            that the dimensions agree.

            .. SEEALSO::

                :meth:`structure_graph()` for a definition/construction of the
                structure graph of a Lie algebra.
            """
            # This is not quite sufficient...
            return self.dimension() == g.dimension() \
                and self.structure_graph().is_isomorphic(g.structure_graph(), edge_labels=True)

        def dimension(self):
            """
            Return the dimension of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
                sage: L.dimension()
                2
            """
            return ZZ(len(self.basis()))

        def _dense_free_module(self, base_ring=None):
            """
            Return ``self`` as a free module (i.e. forgetting the Lie
            bracket structure).

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
                sage: L._dense_free_module()
                Vector space of dimension 2 over Rational Field
            """
            if base_ring is None:
                base_ring = self.base_ring()
            from sage.modules.free_module import FreeModule
            return FreeModule(base_ring, self.dimension())

    class ElementMethods:
        def to_vector(self):
            """
            Return ``self`` as a dense vector.
            """
            M = self.parent()._dense_free_module()
            if not self:
                return M.zero()
            gd = self.parent().gens_dict()
            basis = self.parent().basis()
            ret = [0]*len(basis)
            # Since we always expand in the basis when bracketing,
            #   the monomials should be basis elements
            for k,c in self:
                ret[basis.index(gd[k._name])] = c
            return M(ret)

        def adjoint_matrix(self): # In #11111 (more or less) by using matrix of a mophism
            """
            Return the matrix of the adjoint action of ``self``.

            EXAMPLES::

                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: x.adjoint_matrix()
                [0 0]
                [1 0]
                sage: y.adjoint_matrix()
                [-1  0]
                [ 0  0]
            """
            P = self.parent()
            basis = P.basis()
            return matrix([P.bracket(self, b).to_vector() for b in basis])

