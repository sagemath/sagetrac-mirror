"""
Finitely Presented Algebras

Finitely presented algebras realized as quotients of the free algebra.
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
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
from sage.algebras.algebra import Algebra
from sage.categories.algebras import Algebras
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.noncommutative_ideals import Ideal_nc

class TwoSidedAlgebraIdeal(Ideal_nc):
    """
    A two-sided ideal of an algebra.
    """
    def __init__(self, free_algebra, gens):
        """
        Initialize ``self``.
        """
        PBW = free_algebra.pbw_basis()
        self._gb = [1, map(lambda x: x / x.leading_coefficient(), map(PBW, gens))]
        Ideal_nc.__init__(self, free_algebra, gens, "twosided")

    def free_algebra(self):
        """
        Return the ambient free algebra of ``self``.
        """
        return self.ring()

    @cached_method
    def groebner_basis(self, d=infinity):
        """
        Return a Groebner basis of ``self``.
        """
        F = self.ring()
        if self._gb[0] is None:
            return tuple(map(F, self._gb[1]))
        if d <= self._gb[0]:
            return tuple(map(F, filter(lambda x: len(x.leading_support()) <= d, self._gb[1])))

        # Run Groebner basis algorithm on the PBW gens
        gens = self._gb[1] # The gens are in the PBW basis
        n = F.ngens()
        D = [(g, h) for g in gens for h in gens]
        D += [(g, ZZ(i)) for g in gens for i in range(n)]

        PBW = F.pbw_basis()
        FM = PBW.basis().keys()
        PBW_from_exp = lambda x: PBW.monomial(FM(enumerate(x)))

        while len(D) > 0:
            p = D.pop()
            if len(p[0].leading_support()) >= d:
                continue

            if p[1] in ZZ:
                h = p[0]*PBW.gen(p[1]) - PBW.gen(p[1])*p[0]
            else:
                if len(p[0].leading_support()) + len(p[1].leading_support()) > d:
                    continue

                f = map(lambda x: x.leading_support().to_word().evaluation(), p)
                lcm = list(f[0])
                for i,v in enumerate(f[1]):
                    if i > len(lcm):
                        lcm.append(v)
                    else:
                        lcm[i] = max(lcm[i], v)
                f = map(lambda x: [lcm[i] - v for i,v in enumerate(x)], f)
                h = PBW_from_exp(f[0]) * p[0] - PBW_from_exp(f[1]) * p[1]

            h = self._reduce_pbw(h, gens)
            if h != PBW.zero():
                h /= h.leading_coefficient()
                for g in gens:
                    D.append((g, h))
                for i in range(n):
                    D.append((g, ZZ(i)))
                gens.append(h)

        self._gb = [d, gens]
        # Convert back
        if d is None:
            return tuple(map(F, gens))
        # Remove all elements with larger degree
        return tuple(map(F, filter(lambda x: len(x.leading_support()) <= d, gens)))

    def _reduce_pbw(self, x, G):
        """
        Return ``x`` modulo ``G`` where ``x`` is in the PBW basis.
        """
        F = self.ring()
        PBW = F.pbw_basis()
        ret = PBW.zero()
        FM = PBW.basis().keys()

        while x != 0:
            print x
            u, la = x.leading_item()
            found = False
            for g in G:
                LM, mu = g.leading_item()
                #print LM, g

                # Factor u by LM if possible
                w = u.to_word()
                f = LM.to_word().first_pos_in(w)
                print w, LM.to_word(), f
                if f is not None:
                    print FM(w[:f]), g, FM(w[f+len(LM):])
                    print la / mu * PBW(FM(w[:f])) * g * PBW(FM(w[f+len(LM):]))
                    found = True
                    x -= la / mu * PBW(FM(w[:f])) * g * PBW(FM(w[f+len(LM):]))
                    break
            if not found:
                ret += la * F(u)
                x -= la * F(u)
            print "-------------------"
        return ret

    def reduce(self, x):
        """
        Return ``x`` modulo ``self``.
        """
        F = self.ring()
        PBW = F.pbw_basis()
        p = PBW(x)
        G = self.groebner_basis( len(p.leading_support()) )
        return F(self._reduce_pbw(p, map(PBW, G)))

class FinitelyPresentedAlgebra(Algebra, UniqueRepresentation):
    """
    A finitely presented algebra realized as a quotient of the free algebra.

    INPUT:

    - ``free_algebra`` -- the ambient free algebra
    - ``ideal`` -- the defining ideal
    - ``category`` -- (optional) the category
    """
    def __init__(self, free_algebra, ideal, names=None, category=None):
        """
        Initialize ``self``.
        """
        self._free_alg = free_algebra
        R = self._free_alg.base_ring()
        self._ideal = ideal
        category = Algebras(R).or_subcategory(category)
        if names is None:
            names = free_algebra.variable_names()
        #Parent.__init__(self, base=R, category=category)
        Algebra.__init__(self, R, names, normalize=True, category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Algebra generated by {} with relations {} over {}".format(
               self._free_alg.gens(),
               self._ideal.gens(), self.base_ring())

    def _latex_(self):
        """
        Return a latex representation of ``self``.
        """
        from sage.misc.latex import latex
        ret = latex(self.base_ring()) + "\\langle " + latex(self._free_alg.gens())
        ret += " \mid " + latex(self._ideal.gens()) + "\\rangle"
        return ret

    def relations(self):
        """
        Return the relations of ``self`` in the ambient free algebra.
        """
        return self._ideal.gens()

    def defining_ideal(self):
        """
        Return the defining ideal of ``self``.
        """
        return self._ideal

    class Element(ElementWrapper):
        """
        An element in a finitely presented algebra.
        """

