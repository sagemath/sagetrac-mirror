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
        elif basis_indices.cardinality() != g.basis().cardinality():
            raise ValueError("the number of basis indices does not"
                             " match the number of basis elements")
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

        monoid_cmp = lambda x,y: basis_cmp(x[0],y[0])
        self._basis_cmp = basis_cmp
        monomials = IndexedFreeAbelianMonoid(basis_indices, prefix,
                                             monomial_cmp=monoid_cmp, **kwds)
        self._basis_map = {x: monomials.gen(basis_indices[i])
                           for i,x in enumerate(g.basis().keys())}
        self._basis_map_inv = {basis_indices[i]: x for i,x in enumerate(g.basis())}

        CombinatorialFreeModule.__init__(self, R, basis,
                                         prefix='', bracket=False, latex_bracket=False,
                                         category=AlgebrasWithBasis(R))

        # Default lifting coercion

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Poincare-Birkhoff-Witt basis corresponding to {}".format(self._g)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion map from ``R`` to ``self``.
        """
        if R == self._g:
            basis_function = lambda x: self.monomial(self._basis_map[x])
            # TODO: this diagonal, but with a smaller indexing set...
            return g.module_morphism(basis_function, codomain=self,
                                     triangular='upper', unitriangular=True)
        return super(PoincareBirkhoffWittBasis, self)._coerce_map_from_(R)

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
        trail = lhs.trailing_support()
        lead = rhs.leading_support()
        if self._basis_cmp(trail, lead) <= 0:
            return self.monomial(lhs * rhs)

        # Create the commutator
        # We have xy - yx = [x, y] -> xy = yx + [x, y] and we have x > y
        terms = self._basis_map_inv[trail].bracket(self._basis_map_inv[lead])
        lead = I.gen(lead)
        trail = I.gen(trail)
        terms = self.sum_of_terms((self._basis_map[t], c) for t,c in terms)
        terms += self.monomial(lead * trail)
        return self.monomial(lhs.cancel(trail)) * terms * self.monomial(rhs.cancel(lead))

class IdealPBW: # TODO: inherit from the correct classes
    """
    An ideal of a universal enveloping algebra in the PBW basis.
    """
    def __init__(self, UEA, gens):
        """
        Initialize ``self``.
        """
        self._UEA = UEA
        self._gens = map(lambda x: x / x.leading_coefficient(), gens) # Make everything monic
        # TODO: passing up __init__

    def gens(self):
        """
        Return the generators of ``self``.
        """
        return self._gens

    @cached_method
    def groebner_basis(self):
        """
        Compute a Groebner basis of ``self``.
        """
        G = self.gens()
        GU = self._UEA.gens()
        D = [(x,y) for x in G for y in G] + [(x,ZZ(i)) for x in G for i in range(len(GU))]
        while len(D) != 0:
            p = D.pop()
            if p[1].parent() is ZZ:
                h = p[0] * GU[p[1]] - GU[p[1]] * p[0]
            else:
                S = [p[0].support(), p[1].support()]
                h = self._UEA.monomial(prod([x for x in S[1] if x not in S[0]])) * p[0]
                h -= self._UEA.monomial(prod([x for x in S[0] if x not in S[1]])) * p[1]

            h = self._reduce(h, G)
            if h != 0:
                h /= h.leading_coefficient()
                for g in G:
                    D.append((g, h))
                for i in range(len(GU)):
                    D.append((h, i))
                G.append(h)
        return G

    def reduce(self, x, G=None):
        """
        Return ``x`` modulo ``self``.
        """
        ret = self._UEA.zero()
        R = self._UEA.base_ring()
        if G is None:
            G = self.groebner_basis()
        while x != 0:
            la,u = x.leading_item()
            found = False
            for g in G:
                mu, lm = g.leading_item()
                f = u.factorize_by(lm)
                if f is not None:
                    x -= la / mu * self._UEA.monomial(f[0]) * g * self._UEA.monomial(f[1])
                    found = True
                    break
            if not found:
                t = x.leading_term()
                x -= t
                ret += t
        return ret

class QuotientPBW: # TODO: inherit from the correct classes
    """
    A quotient algebra of a universal enveloping algebra in the PBW basis.

    INPUT:

    - ``UEA`` -- the universal enveloping algebra in the PBW basis
    - ``I`` -- the ideal
    """
    def __init__(self, UEA, I):
        """
        Initialize ``self``.
        """
        self._UEA = UEA
        self._I = I
        # TODO: passing up __init__

