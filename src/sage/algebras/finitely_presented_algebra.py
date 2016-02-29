"""
Finitely Presented Algebras

Finitely presented algebras are realized as quotients of the free algebra.
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
from sage.algebras.algebra_element import AlgebraElement
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.noncommutative_ideals import Ideal_nc

class FreeAlgebraIdeal(Ideal_nc):
    """
    An ideal of a free algebra.
    """
    def __init__(self, F, gens, coerce=True, side="twosided"):
        """
        Initialize ``self``.
        """
        Ideal_nc.__init__(self, F, gens, coerce, side)

        gens = self.gens()
        if F.zero() in gens:
            self._gb = (F.zero(),)
            self._gb_todo = []
            return
        if F.one() in gens:
            self._gb = tuple(F.gens())
            self._gb_todo = []
            return

        gb = map(lambda x: x / x.leading_coefficient(FreeAlgebraIdeal._lead_cmp), gens)
        self._gb = gb
        self._gb_todo = [(g, h) for g in gb for h in gb]

    def free_algebra(self):
        """
        Return the ambient free algebra of ``self``.
        """
        return self.ring()

    @staticmethod
    def _lead_cmp(x, y):
        """
        Compare ``x`` and ``y`` by degree then reverse lex so
        :meth:`leading_item()` returns the largest degree and smallest lex.
        """
        x = x.to_word()
        y = y.to_word()
        c = cmp(len(x), len(y))
        if c != 0:
            return c
        return cmp(y, x)

    def groebner_basis(self, max_steps=infinity):
        """
        Return a Groebner basis of ``self``.

        INPUT:

        - ``max_steps`` -- (default: infinity) the maximum number of steps to
          do before terminating

        .. WARNING::

            This will run forever if the Groebner basis is infinite and
            ``max_steps`` is not specified.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ)
            sage: I = F.ideal(x^2 - x*y)
            sage: I.groebner_basis(5)
            WARNING: returning an incomplete Groebner basis
            (x^2 - x*y,
             x*y*x - x*y^2,
             x*y^2*x - x*y^3,
             x*y^3*x - x*y^4,
             x*y^4*x - x*y^5)
            sage: I = F.ideal(z^2 - z*y)
            sage: I.groebner_basis()
            (z*y - z^2,)
            sage: I = F.ideal(x^2 - x*y, x*y*x - y*x*y)
            sage: I.groebner_basis()
            (x^2 - x*y, x*y*x - y*x*y, x*y^2 - y*x*y)
            sage: I = F.ideal(z^2 - z*y, z*y*z - y*z*y)
            sage: I.groebner_basis()
            (z*y - z^2, y*z*y - z*y*z, y*z^4 - z^5)
        """
        side = self.side()
        if side != "twosided":
            raise NotImplementedError

        # Setup variables and functions
        l_cmp = FreeAlgebraIdeal._lead_cmp
        F = self.ring()

        n = len(F.gens())
        FM = F._basis_keys

        # Run Groebner basis algorithm

        zero = F.zero()
        from sage.combinat.words.word import Word
        while len(self._gb_todo) > 0 and len(self._gb) < max_steps:
            p = self._gb_todo.pop(0)

            # Compute the essential common multiples of the leading terms
            f = map(lambda x: list(x.leading_support(l_cmp).to_word()), p)
            ell = min(len(f[0]), len(f[1]))
            for k in range(1, ell):
                if f[0][:k] == f[1][-k:]:
                    w0 = Word(f[1][:-k])
                    w1 = Word(f[0][k:])

                    # Compute the S-polynomial
                    h = F(FM(w0)) * p[0] - p[1] * F(FM(w1))
                    h = self._normal_form(h, self._gb)

                    if h != zero:
                        h = h / h.leading_coefficient(l_cmp)
                        self._gb_todo.extend([(g, h) for g in self._gb])
                        self._gb.append(h)

        if len(self._gb_todo) != 0:
            print "WARNING: returning an incomplete Groebner basis"

        return tuple(self._gb)

    def _normal_form(self, x, G):
        """
        Return ``x`` modulo ``G`` (i.e. the normal form with respect
        to ``G``).
        """
        F = self.ring()
        ret = F.zero()
        FM = F._basis_keys
        l_cmp = FreeAlgebraIdeal._lead_cmp
        side = self.side()

        while x != 0:
            u, la = x.leading_item(l_cmp)
            found = False
            for g in G:
                LM, mu = g.leading_item(l_cmp)

                # Factor u by LM if possible
                w = u.to_word()
                lmw = LM.to_word()
                if side == 'twosided':
                    f = lmw.first_pos_in(w)
                    if f is not None:
                        found = True
                        x -= la / mu * F(FM(w[:f])) * g * F(FM(w[f+len(LM):]))
                        break
                elif side == 'left':
                    if lmw.is_proper_suffix(w):
                        x -= la / mu * F(FM(w[len(lmw):])) * g
                        break
                elif side == 'right':
                    if lmw.is_proper_prefix(w):
                        x -= la / mu * g * F(FM(w[:-len(lmw)]))
                        break
            if not found:
                ret += la * F(u)
                x -= la * F(u)
        return ret

    def reduce(self, x, max_steps=infinity):
        """
        Return ``x`` modulo ``self``.

        INPUT:

        - ``max_steps`` -- (default: infinity) the maximum number of steps to
          compute in the Groebner basis before terminating

        .. WARNING::

            If this terminates after ``max_steps``, the resulting reduction
            may not fully reduce the element.
        """
        if x == self.ring().zero():
            return x
        G = self.groebner_basis(max_steps)
        return self._normal_form(x, G)

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
               self.gens(), self._ideal.gens(), self.base_ring())

    def _latex_(self):
        """
        Return a latex representation of ``self``.
        """
        from sage.misc.latex import latex
        ret = latex(self.base_ring()) + "\\langle " + latex(self._free_alg.gens())
        ret += " \mid " + latex(self._ideal.gens()) + "\\rangle"
        return ret

    def ngens(self):
        """
        Return the number of generators of ``self``.
        """
        return self._free_alg.ngens()

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.
        """
        return tuple(map(lambda x: self.element_class(self, x), self._free_alg.gens()))

    def gen(self, i):
        return self.gens()[i]

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

    def quotient(self, I, names=None):
        """
        Return a quotient of ``self`` by the ideal ``I``.
        """
        if I.side() == 'twosided':
            Ip = self._free_alg.ideal(self._ideal.gens() + I.gens())
            return self._free_alg.quotient(Ip, names)
        return super(FinitelyPresentedAlgebra, self).quotient(I, names)

    @cached_method
    def basis(self):
        """
        Return a monomial basis of ``self`` if finite dimensional; otherwise
        this runs forever.

        EXAMPLES::

            sage: F.<x,y,z> = FreeAlgebra(QQ, 3)
            sage: PBW = F.pbw_basis()
            sage: xp,yp,zp = PBW.gens()
            sage: Ip = PBW.ideal(xp^2,yp^2,zp^4, xp*yp + yp*zp, zp*yp, yp*zp, xp*zp - zp*xp)
            sage: Qp = PBW.quotient(Ip)
            sage: Qp.basis() # long time
            (PBW[x],
             PBW[y],
             PBW[z],
             PBW[y]*PBW[x],
             PBW[z]*PBW[x],
             PBW[z]^2,
             PBW[z]^2*PBW[x],
             PBW[z]^3,
             PBW[z]^3*PBW[x])
        """
        mon = reduce(lambda x,y: x.union(set(y.monomials())), self.gens(), set([]))
        todo = set([(x,y) for x in mon for y in mon])
        while len(todo) != 0:
            x,y = todo.pop()
            n = x * y
            m = set(n.monomials())
            add = m.difference(mon)
            mon = mon.union(add)
            for x in mon:
                for y in add:
                    todo.add((x,y))
                    todo.add((y,x))
        mon.discard(self.zero())
        return tuple(sorted(mon, key=lambda x: x.value))

    class Element(AlgebraElement):
        """
        An element in a finitely presented algebra.
        """
        def __init__(self, parent, value, reduce=True):
            if reduce:
                value = parent._ideal.reduce(value)
            self.value = value
            AlgebraElement.__init__(self, parent)

        def _repr_(self):
            return repr(self.value)
        def _latex_(self):
            return latex(self.value)

        def __eq__(self,other):
            return isinstance(other, FinitelyPresentedAlgebra.Element) \
                and self.parent() == other.parent() and self.value == other.value

        def _add_(self, rhs):
            return self.__class__(self.parent(), self.value + rhs.value)

        def _sub_(self, rhs):
            return self.__class__(self.parent(), self.value - rhs.value)

        def _mul_(self, rhs):
            return self.__class__(self.parent(), self.value * rhs.value)

        def _rmul_(self, rhs):
            return self.__class__(self.parent(), rhs * self.value)

        def _lmul_(self, rhs):
            return self.__class__(self.parent(), self.value * rhs)

        def monomials(self):
            """
            Return the monomials of ``self``.
            """
            P = self.parent()
            return map(P, self.value.monomials())

