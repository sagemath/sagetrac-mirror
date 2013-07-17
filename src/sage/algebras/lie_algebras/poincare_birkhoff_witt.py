"""
Poincare-Birkhoff-Witt Basis For the Universal Enveloping Algebras

AUTHORS:

- Travis Scrimshaw (2013-11-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.monoids.indexed_monoid import IndexedFreeAbelianMonoid
from sage.combinat.free_module import CombinatorialFreeModule
from sage.sets.family import Family

class PoincareBirkhoffWittBasis(CombinatorialFreeModule):
    r"""
    The Poincare-Birkhoff-Witt basis of the universal enveloping algebra of
    a Lie algebra.
    """
    @staticmethod
    def __classcall_private__(cls, g, basis_indices=None, basis_cmp=None,
                              prefix='PBW', **kwds):
        """
        Normalize input to ensure a unique representation.
        """
        if basis_indices is None:
            basis_indices = g.basis().keys()
        elif len(basis_indices) != len(g.basis()):
            raise ValueError("the number of basis indices does not"
                             "match the number of basis elements")
        return super(PoincareBirkhoffWittBasis, cls).__classcall__(cls,
                            g, tuple(basis_indices), basis_cmp, prefix, **kwds)

    def __init__(self, g, basis_indices, basis_cmp, prefix, **kwds):
        """
        Initialize ``self``.
        """
        if basis_cmp is None:
            lst = list(basis_indices)
            basis_cmp = lambda x,y: cmp(lst.index(x), lst.index(y))

        R = g.base_ring()
        self._g = g

        basis = IndexedFreeAbelianMonoid(basis_indices, prefix, monomial_cmp=basis_cmp, **kwds)
        self._basis_map = {x: basis.gen(basis_indices[i]) for i,x in enumerate(g.basis().keys())}
        self._basis_map_inv = {basis_indices[i]: x for i,x in enumerate(g.basis())}

        CombinatorialFreeModule.__init__(self, R, basis,
                                         prefix='', bracket=False, latex_bracket=False,
                                         category=AlgebrasWithBasis(R))

        # Default lifting coercion
        basis_function = lambda x: self.monomial(self._basis_map[x])
        g.module_morphism(basis_function, codomain=self,
                          triangular='upper', unitriangular=True).register_as_coercion()

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Poincare-Birkhoff-Witt basis corresponding to {}".format(self._g)

    def lie_algebra(self):
        """
        Return the underlying Lie algebra of ``self``.
        """
        return self._g

    def algebra_generators(self):
        """
        Return the algebra generators of ``self``.
        """
        G = self._indices.gens()
        return Family(self._indices._indices, lambda x: self.monomial(G[x]))

    @cached_method
    def one_basis(self):
        """
        Return the basis element indexing `1`.
        """
        return self._indices.one()

    def product_on_basis(self, lhs, rhs):
        """
        Return the product of the two basis elements ``lhs`` and ``rhs``.
        """
        # Some trivial base cases
        if lhs == self.one_basis():
            return self.monomial(rhs)
        if rhs == self.one_basis():
            return self.monomial(lhs)

        I = self._indices
        basis_cmp = I.print_options()['monomial_cmp']
        trail = lhs.trailing_support()
        lead = rhs.leading_support()
        if basis_cmp(trail, lead) <= 0:
            return self.monomial(lhs * rhs)

        # Create the commutator
        # We have xy - yx = [x, y] -> xy = yx + [x, y] and we have x > y
        terms = self._basis_map_inv[trail].bracket(self._basis_map_inv[lead])
        lead = I.gen(lead)
        trail = I.gen(trail)
        terms = self.sum_of_terms((self._basis_map[t], c) for t,c in terms)
        terms += self.monomial(lead * trail)
        return self.monomial(lhs.cancel(trail)) * terms * self.monomial(rhs.cancel(lead))

