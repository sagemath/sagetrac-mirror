"""
Poincare-Birkhoff-Witt Algebras

Algebras in a PBW(-type) basis.
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
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.rings.noncommutative_ideals import Ideal_nc
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement

class PBWBasisOfFreeAlgebra(CombinatorialFreeModule):
    """
    The Poincare-Birkhoff-Witt basis of the free algebra.

    EXAMPLES::

        sage: F.<x,y> = FreeAlgebra(QQ, 2)
        sage: PBW = F.pbw_basis()
        sage: px, py = PBW.gens()
        sage: px * py
        PBW[x*y] + PBW[y]*PBW[x]
        sage: py * px
        PBW[y]*PBW[x]
        sage: px * py^3 * px - 2*px * py
        -2*PBW[x*y] - 2*PBW[y]*PBW[x] + PBW[x*y^3]*PBW[x] + PBW[y]*PBW[x*y^2]*PBW[x]
         + PBW[y]^2*PBW[x*y]*PBW[x] + PBW[y]^3*PBW[x]^2

    We can convert between the two bases::

        sage: p = PBW(x*y - y*x + 2); p
        2*PBW[1] + PBW[x*y]
        sage: F(p)
        2 + x*y - y*x
        sage: f = F.pbw_element(x*y*x + x^3*y + x + 3)
        sage: F(PBW(f)) == f
        True
        sage: p = px*py + py^4*px^2
        sage: F(p)
        x*y + y^4*x^2
        sage: PBW(F(p)) == p
        True

    Note that multiplication in the PBW basis agrees with multiplication
    as monomials::

        sage: F(px * py^3 * px - 2*px * py) == x*y^3*x - 2*x*y
        True

    TESTS:

    Check that going between the two bases is the identity::

        sage: F = FreeAlgebra(QQ, 2, 'x,y')
        sage: PBW = F.pbw_basis()
        sage: M = F.monoid()
        sage: L = [j.to_monoid_element() for i in range(6) for j in Words('xy', i)]
        sage: all(PBW(F(PBW(m))) == PBW(m) for m in L)
        True
        sage: all(F(PBW(F(m))) == F(m) for m in L)
        True
    """
    @staticmethod
    def __classcall_private__(cls, R, n=None, names=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.algebras.free_algebra import PBWBasisOfFreeAlgebra
            sage: PBW1 = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: PBW2.<x,y> = PBWBasisOfFreeAlgebra(QQ)
            sage: PBW3 = PBWBasisOfFreeAlgebra(QQ, 2, ['x','y'])
            sage: PBW1 is PBW2 and PBW2 is PBW3
            True
        """
        if n is None and names is None:
            from sage.algebras.free_algebra import FreeAlgebra_generic
            if not isinstance(R, FreeAlgebra_generic):
                raise ValueError("{} is not a free algebra".format(R))
            alg = R
        else:
            if n is None:
                n = len(names)
            from sage.algebras.free_algebra import FreeAlgebra
            alg = FreeAlgebra(R, n, names)
        return super(PBWBasisOfFreeAlgebra, cls).__classcall__(cls, alg)

    def __init__(self, alg):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: PBW = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: TestSuite(PBW).run()
        """
        R = alg.base_ring()
        self._alg = alg
        category = AlgebrasWithBasis(R)
        CombinatorialFreeModule.__init__(self, R, alg.monoid(), prefix='PBW',
                                         category=category)
        self._assign_names(alg.variable_names())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            The Poincare-Birkhoff-Witt basis of Free Algebra on 2 generators (x, y) over Rational Field
        """
        return "The Poincare-Birkhoff-Witt basis of %s"%(self._alg)

    def _repr_term(self, w):
        """
        Return a representation of term indexed by ``w``.

        EXAMPLES::

            sage: PBW = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: x,y = PBW.gens()
            sage: x*y # indirect doctest
            PBW[x*y] + PBW[y]*PBW[x]
            sage: y*x
            PBW[y]*PBW[x]
            sage: x^3
            PBW[x]^3
            sage: PBW.one()
            PBW[1]
            sage: 3*PBW.one()
            3*PBW[1]
        """
        if len(w) == 0:
            return super(PBWBasisOfFreeAlgebra, self)._repr_term(w)
        ret = ''
        p = 1
        cur = None
        for x in w.to_word().lyndon_factorization():
            if x == cur:
                p += 1
            else:
                if len(ret) != 0:
                    if p != 1:
                        ret += "^{}".format(p)
                    ret += "*"
                ret += super(PBWBasisOfFreeAlgebra, self)._repr_term(x.to_monoid_element())
                cur = x
                p = 1
        if p != 1:
            ret += "^{}".format(p)
        return ret

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: F.<x,y> = FreeAlgebra(QQ, 2)
            sage: R = F.pbw_basis()
            sage: R(3)
            3*PBW[1]
            sage: R(x*y)
            PBW[x*y] + PBW[y]*PBW[x]
        """
        from sage.algebras.free_algebra import FreeAlgebraElement
        if isinstance(x, FreeAlgebraElement):
            return self._alg.pbw_element(self._alg(x))
        return CombinatorialFreeModule._element_constructor_(self, x)

    def _coerce_map_from_(self, R):
        """
        Return ``True`` if there is a coercion from ``R`` into ``self`` and
        ``False`` otherwise.  The things that coerce into ``self`` are:

        - Anything that coerces into the associated free algebra of ``self``

        TESTS::

            sage: F = FreeAlgebra(ZZ, 3, 'x,y,z').pbw_basis()
            sage: G = FreeAlgebra(QQ, 3, 'x,y,z').pbw_basis()
            sage: H = FreeAlgebra(ZZ, 1, 'y').pbw_basis()
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            False
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
            sage: F._coerce_map_from_(G._alg.monoid())
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
            sage: F.has_coerce_map_from(FreeAlgebra(ZZ, 3, 'x,y,z'))
            True
        """
        return self._alg.has_coerce_map_from(R)

    def one_basis(self):
        """
        Return the index of the basis element for `1`.

        EXAMPLES::

            sage: PBW = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: PBW.one_basis()
            1
            sage: PBW.one_basis().parent()
            Free monoid on 2 generators (x, y)
        """
        return self._basis_keys.one()

    def algebra_generators(self):
        """
        Return the generators of ``self`` as an algebra.

        EXAMPLES::

            sage: PBW = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: gens = PBW.algebra_generators(); gens
            (PBW[x], PBW[y])
            sage: all(g.parent() is PBW for g in gens)
            True
        """
        return tuple(self.monomial(x) for x in self._basis_keys.gens())

    gens = algebra_generators

    def gen(self, i):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: PBW = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: PBW.gen(0)
            PBW[x]
            sage: PBW.gen(1)
            PBW[y]
        """
        return self.algebra_generators()[i]

    def ideal(self, *args, **kwds):
        """
        Return a side ``side`` ideal of ``self`` given by ``gens``.
        """
        if kwds.get('side', 'twosided') == 'twosided':
            return TwoSidedPBWIdeal(self, args)
        return super(PBWBasisOfFreeAlgebra, self).ideal(self, *args, **kwds)

    def quotient(self, I, names=None, category=None):
        """
        Return the quotient of ``self`` by ``I``.
        """
        from sage.algebras.finitely_presented_algebra import FinitelyPresentedAlgebra
        if not isinstance(I, TwoSidedPBWIdeal):
            raise TypeError("must be a two sided algebra ideal")
        return FinitelyPresentedAlgebra(self, I, names)

    def free_algebra(self):
        """
        Return the associated free algebra of ``self``.

        EXAMPLES::

            sage: PBW = FreeAlgebra(QQ, 2, 'x,y').pbw_basis()
            sage: PBW.free_algebra()
            Free Algebra on 2 generators (x, y) over Rational Field
        """
        return self._alg

    def product(self, u, v):
        """
        Return the product of two elements ``u`` and ``v``.

        EXAMPLES::

            sage: F = FreeAlgebra(QQ, 2, 'x,y')
            sage: PBW = F.pbw_basis()
            sage: x, y = PBW.gens()
            sage: PBW.product(x, y)
            PBW[x*y] + PBW[y]*PBW[x]
            sage: PBW.product(y, x)
            PBW[y]*PBW[x]
            sage: PBW.product(y^2*x, x*y*x)
            PBW[y]^2*PBW[x^2*y]*PBW[x] + PBW[y]^2*PBW[x*y]*PBW[x]^2 + PBW[y]^3*PBW[x]^3

        TESTS:

        Check that multiplication agrees with the multiplication in the
        free algebra::

            sage: F = FreeAlgebra(QQ, 2, 'x,y')
            sage: PBW = F.pbw_basis()
            sage: x, y = PBW.gens()
            sage: F(x*y)
            x*y
            sage: F(x*y*x)
            x*y*x
            sage: PBW(F(x)*F(y)*F(x)) == x*y*x
            True
        """
        return self(self.expansion(u) * self.expansion(v))

    def expansion(self, t):
        """
        Return the expansion of the element ``t`` of the Poincare-Birkhoff-Witt
        basis in the monomials of the free algebra.

        EXAMPLES::

            sage: F = FreeAlgebra(QQ, 2, 'x,y')
            sage: PBW = F.pbw_basis()
            sage: x,y = F.monoid().gens()
            sage: PBW.expansion(PBW(x*y))
            x*y - y*x
            sage: PBW.expansion(PBW.one())
            1
            sage: PBW.expansion(PBW(x*y*x) + 2*PBW(x) + 3)
            3 + 2*x + x*y*x - y*x^2

        TESTS:

        Check that we have the correct parent::

            sage: PBW.expansion(PBW(x*y)).parent() is F
            True
            sage: PBW.expansion(PBW.one()).parent() is F
            True
        """
        return sum([i[1] * self._alg.lie_polynomial(i[0]) for i in list(t)],
                   self._alg.zero())

    class Element(CombinatorialFreeModuleElement):
        def expand(self):
            """
            Expand ``self`` in the monomials of the free algebra.

            EXAMPLES::

                sage: F = FreeAlgebra(QQ, 2, 'x,y')
                sage: PBW = F.pbw_basis()
                sage: x,y = F.monoid().gens()
                sage: f = PBW(x^2*y) + PBW(x) + PBW(y^4*x)
                sage: f.expand()
                x + x^2*y - x*y*x + y^4*x
            """
            return self.parent().expansion(self)

class TwoSidedPBWIdeal(Ideal_nc):
    """
    A two-sided ideal of an algebra in a PBW(-type) basis.
    """
    def __init__(self, pbw_algebra, gens):
        """
        Initialize ``self``.
        """
        gb = map(lambda x: x / x.leading_coefficient(), gens)
        self._gb = [1, gb]
        n = len(pbw_algebra.algebra_generators())
        self._gb_todo = [(g, h) for g in gb for h in gb]
        self._gb_todo += [(g, ZZ(i)) for g in gb for i in range(n)]

        Ideal_nc.__init__(self, pbw_algebra, gens, "twosided")

    def pbw_basis(self):
        """
        Return the ambient PBW basis of ``self``.
        """
        return self.ring()

    @staticmethod
    def _lead_cmp(x, y):
        """
        Compare ``x`` and ``y`` by degree then reverse lex so
        :meth:`leading_item()` returns the largest degree and smallest lex.
        """
        c = cmp(len(x), len(y))
        if c != 0:
            return c
        return cmp(y, x)

    #@cached_method
    def groebner_basis(self, d=infinity):
        """
        Return a Groebner basis of ``self``.
        """
        l_cmp = TwoSidedPBWIdeal._lead_cmp
        PBW = self.ring()
        if self._gb[0] is None:
            return tuple(map(PBW, self._gb[1]))
        if d <= self._gb[0]:
            return tuple(map(PBW, filter(lambda x: len(x.leading_support(l_cmp)) <= d,
                                         self._gb[1])))

        # Run Groebner basis algorithm on the PBW gens
        gens = self._gb[1] # The gens are in the PBW basis
        n = len(PBW.algebra_generators())
        FM = PBW.basis().keys()
        PBW_from_exp = lambda x: PBW.monomial(FM( [a for a in enumerate(x) if a[1] > 0] ))

        cur = 0
        while cur < len(self._gb_todo):
            p = self._gb_todo[cur]
            if len(p[0].leading_support(l_cmp)) >= d:
                cur += 1
                continue

            if p[1] in ZZ:
                h = p[0]*PBW.gen(p[1]) - PBW.gen(p[1])*p[0]
            else:
                if len(p[0].leading_support(l_cmp)) + len(p[1].leading_support(l_cmp)) > d:
                    cur += 1
                    continue

                f = map(lambda x: x.leading_support(l_cmp).to_word().evaluation(), p)
                lcm = list(f[0])
                for i,v in enumerate(f[1]):
                    if i > len(lcm):
                        lcm.append(v)
                    else:
                        lcm[i] = max(lcm[i], v)
                f = map(lambda x: [lcm[i] - v for i,v in enumerate(x)], f)
                h = PBW_from_exp(f[0]) * p[0] - p[1] * PBW_from_exp(f[1])

            self._gb_todo.pop(cur) # We have taken care of the pairing
            h = self._normal_form(h, gens)

            if h != PBW.zero():
                h /= h.leading_coefficient(l_cmp)
                for g in gens:
                    self._gb_todo.append((g, h))
                for i in range(n):
                    self._gb_todo.append((g, ZZ(i)))
                gens.append(h)

        if len(self._gb_todo) == 0:
            d = infinity
        self._gb = [d, gens]
        # Convert back
        if d is None:
            return tuple(gens)
        # Remove all elements with larger degree
        return tuple(filter(lambda x: len(x.leading_support(l_cmp)) <= d, gens))

    def _normal_form(self, x, G):
        """
        Return ``x`` modulo ``G`` (i.e. the normal form) where ``x`` is in
        the PBW basis.
        """
        PBW = self.ring()
        ret = PBW.zero()
        FM = PBW.basis().keys()
        l_cmp = TwoSidedPBWIdeal._lead_cmp

        while x != 0:
            u, la = x.leading_item(l_cmp)
            found = False
            for g in G:
                LM, mu = g.leading_item(l_cmp)

                # Factor u by LM if possible
                w = u.to_word()
                f = LM.to_word().first_pos_in(w)
                if f is not None:
                    found = True
                    x -= la / mu * PBW(FM(w[:f])) * g * PBW(FM(w[f+len(LM):]))
                    break
            if not found:
                ret += la * PBW(u)
                x -= la * PBW(u)
        return ret

    def reduce(self, x):
        """
        Return ``x`` modulo ``self``.
        """
        if x == self.ring().zero():
            return x
        G = self.groebner_basis( len(x.leading_support()) )
        return self._normal_form(x, G)

