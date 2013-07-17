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
        def structure_coefficients(self):
            """
            Return the structure coefficients of ``self``.
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

        def structure_graph(self):
            r"""
            Return the edge labelled directed graph associated with ``self``.

            We construct the directed graph with vertex set being the basis
            elements `x \in B` and all pairs of basis elements `(x, y)` such
            that `x < y` for some total ordering on `B`. There is an egde
            between `(x, y)` and `z` labelled by the structure coefficient
            `c_{xy}^z` if and only if `c_{xy}^z \neq 0`.

            Here we return a modified form of the graph above by stripping all
            singletons, vertices which do not have either an incoming or
            outgoing edge.
            """
            from sage.graphs.digraph import DiGraph
            G = DiGraph()
            S = self.structure_coefficients()
            for k in S.keys():
                for g,c in S[k]:
                    G.add_edge(k, g, c)
            return G

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
            return len(self.basis())

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

        # I think we need some ambient module to do this computation in...
        @cached_method
        def basis(self):
            """
            Return the basis of ``self``.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
                sage: L.basis()
                (x, y)
            """
            basis = list(self.gens())

            # Setup the basic names for the generators
            basis_names = self.variable_names()

            # Only need to determine the structure coefficients from gens
            # In the process, we will also determine a basis as well
            coeffs = {}
            R = self.base_ring()
            i = 0

            # Modified form of to_vector to take our current basis
            from sage.modules.free_module import FreeModule
            def to_vector(v, cur_names):
                M = FreeModule(R, len(cur_names))
                if not v:
                    return M.zero()
                ret = [0]*len(cur_names)
                # Since we always expand in the basis when bracketing,
                #   the monomials should be basis elements
                for k,c in v:
                    ret[cur_names.index(k)] = c
                return M(ret)

            b_vecs = map(lambda x: to_vector(x, basis), basis_names)
            from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, LieBracket
            while i < len(basis):
                j = i + 1
                while j < len(basis):
                    b = self.bracket(basis[i], basis[j])
                    if not b: # If b == 0, nothing to do
                        continue

                    v = to_vector(b, basis_names)
                    dep = M.linear_dependence(b_vecs + [v]) # b_vecs are linearly indep
                    if len(dep) != 0:
                        dep = dep[0]
                        # If there is a dependency, rewrite as a linear combination of basis
                        s = -dep[-1]
                        d = {basis_names[k]:R(c/s) for k,c in enumerate(dep[:-1]) if c != 0}
                        coeffs[LieBracket(basis_names[i], basis_names[j])] = d
                    else:
                        # We need to add a new basis vector
                        if names is None:
                            basis_names.append(LieGenerator( 'x%s'%len(basis) ))
                        else:
                            comb_name = basis_names[i]._name + basis_names[j]._name
                            basis_names.append(LieGenerator(comb_name))
                        coeffs[LieBracket(basis_names[i], basis_names[j])] = {basis_names[-1]:1}
                        basis.append(b)
                        b_vecs.append(v)
                    j += 1
                i += 1
            return basis

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

