r"""
Lie Algebras With Basis

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
from sage.algebras.lie_algebras.lie_algebra_element import LieBracket

class LieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of Lie algebras with a basis.
    """
    _base_category_class_and_axiom = [LieAlgebras, "WithBasis"]

    class ParentMethods:
        def _basis_cmp(self, x, y):
            """
            Compare two basis element indices. The default is to call ``cmp``.
            """
            return cmp(x, y)

        @abstract_method(optional=True)
        def bracket_on_basis(self, x, y):
            """
            Return the bracket of basis elements indexed by ``x`` and ``y``
            where ``x < y``. If this is not implemented, then the method
            ``_bracket_()`` for the elements must be overwritten.
            """

        @cached_method
        def basis(self, names=None):
            """
            Return the basis of ``self``.

            INPUT:

            - ``gens`` -- a set of generators for the Lie subalgebra a linear
              combination of basis elements of ``self``
            - ``names`` -- (optional) a list of names for the generators

            EXAMPLES::

                sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
                sage: L.basis()
                (x, y)
                sage: L.<x,y> = LieAlgebra(QQ, {('x','y'):{'x':1}})
                sage: L.lie_subalgebra([x, y+3*x]).structure_coefficients()
                Finite family {[x0, x1]: ((x0, 1),)}
                sage: L.lie_subalgebra([2*x+3*y, y+3*x], ['a','b']).structure_coefficients()
                Finite family {[a, b]: ((a, 1), (b, -3))}

                sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
                sage: L.lie_subalgebra([x, y], ['a', 'b']).structure_coefficients()
                Finite family {[a, b]: ((ab, 1),), [a, ab]: ((b, -1),), [ab, b]: ((a, -1),)}
            """
            basis = list(self.gens())

            # Setup the basic names for the generators
            # This needs to be more generic
            if names is None:
                basis_names = map(LieGenerator, self.variable_names())
            else:
                basis_names = map(LieGenerator, names[:len(basis)])

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
                        if names is not None and len(names) < len(basis):
                            basis_names.append(LieGenerator( names[len(basis)] ))
                        else:
                            comb_name = basis_names[i]._name + basis_names[j]._name
                            basis_names.append(LieGenerator(comb_name))
                        coeffs[LieBracket(basis_names[i], basis_names[j])] = {basis_names[-1]:1}
                        basis.append(b)
                        b_vecs.append(v)
                    j += 1
                i += 1
            return basis

        @abstract_method(optional=True)
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

        def pbw_basis(self, basis_indices=None, basis_cmp=None, **kwds):
            """
            Return the Poincare-Birkhoff-Witt basis corresponding to ``self``.

            .. TODO::

                Currently only works for finite dimensional Lie algebras.
            """
            if basis_cmp is None:
                basis_cmp = self._basis_cmp
            from sage.algebras.lie_algebras.poincare_birkhoff_witt \
                import PoincareBirkhoffWittBasis
            return PoincareBirkhoffWittBasis(self, basis_indices, basis_cmp, **kwds)

        poincare_birkhoff_witt_basis = pbw_basis

        _construct_UEA = pbw_basis

    class ElementMethods:
        def _bracket_(self, y):
            """
            Return the Lie bracket of ``[self, y]``.
            """
            P = self.parent()
            def term(ml,mr):
                comp = P._basis_cmp(ml,mr)
                if comp == 0:
                    return P.zero()
                if comp < 0:
                    return P.bracket_on_basis(ml, mr)
                return -P.bracket_on_basis(mr, ml)
            return P.sum(cl*cr * term(ml,mr)
                         for ml,cl in self for mr,cr in y)

