"""
Ideals of Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
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
from sage.structure.parent import Parent
from sage.structure.element import MonoidElement
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import ZZ
from sage.rings.infinity import infinity
from sage.categories.monoids import Monoids
from sage.categories.finite_dimensional_lie_algebras_with_basis import FiniteDimensionalLieAlgebrasWithBasis
from sage.combinat.permutation import Permutations
from sage.combinat.composition import Compositions
from sage.algebras.lie_algebras.free_lie_algebra import is_lyndon
from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator
from sage.algebras.lie_algebras.subalgebra import LieSubalgebra

class LieIdealMonoid(Parent, UniqueRepresentation):
    r"""
    The monoid of ideals in a Lie algebra.
    """
    def __init__(self, lie):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, 'x,y')
            sage: M = sage.algebras.lie_algebras.ideal.LieIdealMonoid(R)
            sage: TestSuite(M).run()
        """
        self.__lie = lie
        Parent.__init__(self, base=ZZ, category=Monoids())
        self._populate_coercion_lists_()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, 'x,y')
            sage: sage.algebras.lie_algebras.ideal.LieIdealMonoid(R)
            Monoid of ideals of Lie algebra
        """
        return "Monoid of ideals of %s"%self.__lie

    def lie_algebra(self):
        r"""
        Return the Lie algebra of which this is the ideal monoid.

        EXAMPLES::
        """
        return self.__lie

    def _element_constructor_(self, x):
        r"""
        Create an ideal in this monoid from ``x``.

        EXAMPLES::
        """
        y = self.__lie.ideal(x)
        y._set_parent(self)
        return y

    def _coerce_map_from_(self, x):
        r"""
        Used by coercion framework.

        EXAMPLES::
        """
        if isinstance(x, LieIdealMonoid):
            return self.lie_algebra().has_coerce_map_from(x.ring())
        else:
            return self.lie_algebra().has_coerce_map_from(x)

    def __cmp__(self, other):
        r"""
        Comparison function.

        EXAMPLES::
        """
        if not isinstance(other, LieIdealMonoid):
            return cmp(type(self), type(other))
        else:
            return cmp(self.lie_algebra(), other.lie_algebra())

class LieAlgebraIdeal(LieSubalgebra): #, MonoidElement): # FIXME: layout conflict
    r"""
    A generic ideal of a Lie algebra `\mathfrak{g}`.
    """
    def __init__(self, lie_algebra, gens, coerce=True):
        """
        Initialize this ideal.

        INPUT:

        - ``lie_algebra`` -- a Lie algebra

        - ``gens`` -- the generators for this ideal

        - ``coerce`` -- (default: ``True``) if ``gens`` needs to be coerced
          into ``lie_algebra``

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: x,y,z = Lyn.gens()
            sage: Lyn.ideal([Lyn([x, y]), z - x])
        """
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [lie_algebra(x) for x in gens]

        if len(gens) == 0:
            gens = (lie_algebra.zero_element(),)
        else:
            # Make sure the leading term has a coefficient of 1 (?)
            gens = tuple(map(lambda x: x/x.leading_coefficient(term_cmp), gens))
        LieSubalgebra.__init__(self, lie_algebra, gens)
        #MonoidElement.__init__(self, LieIdealMonoid(lie_algebra))

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Ideal {} of {}".format(self._gens, self._ambient)

    def _repr_short(self):
        """
        Short representation for the list of generators.

        EXAMPLES:

        If the string representation of a generator contains a line break,
        the generators are not represented from left to right but from
        top to bottom. This is the case, e.g., for matrices::
        """
        L = []
        has_return = False
        for x in self._gens:
            s = repr(x)
            if '\n' in s:
                has_return = True
                s = s.replace('\n','\n  ')
            L.append(s)
        if has_return:
            return '\n(\n  {}\n)\n',format(',\n\n  '.join(L))
        return '({})'.format(', '.join(L))

class LieAlgebraIdealFinitelyPresented(LieAlgebraIdeal):
    """
    An ideal of a finitely presented Lie algebra.
    """
    def __init__(self, lie_algebra, gens, coerce=True):
        """
        Initialize ``self``.
        """
        LieAlgebraIdeal.__init__(self, lie_algebra, gens, coerce=True)
        gb = self._inner_reduce(self._gens)
        self._gb = set(gb)
        gb = list(gb)
        self._gb_todo = set([(g,h) for i,g in enumerate(gb) for h in gb[i+1:]])

    def groebner_basis(self, max_size=infinity):
        """
        Return a Groebner-Shirshov basis for ``self``.
        """
        zero = self._ambient.zero()
        while len(self._gb_todo) != 0 and len(self._gb) < max_size:
            p = self._gb_todo.pop()
            for c in self._useful_compositions(p[0], p[1]):
                g = self._normal_from(self, self._gb)
                if g == zero: # Completely reduced, nothing to do
                    continue

                # Normalize and add more things to check
                g /= g.leading_coefficient(term_cmp)
                for f in M:
                    D.add((f, g))

                M, orig_to_red = self._inner_reduce_with_map(self._basis + [g])
                D = set()
                for k in orig_to_red:
                    # Filter out the entries (f, k)
                    for d in D:
                        if d[1] != k:
                            D.add(d)
                    # Add new entries as needed
                    val = orig_to_red[k]
                    if val is not None:
                        for f in M:
                            D.add((f, val))
                self._gb_todo = D
        if len(self._gb_todo) != 0:
            print("WARNING: only a partial Groebner basis has been computed,"
                  " the result may not be completely reduced")
        return tuple(self._gb)

    def is_groebner_basis(self, G):
        """
        Return if ``G`` is a Groebner basis.
        """
        # Currently we just test that we are a Groebner basis
        zero = self._ambient.zero()
        for g in G:
            for h in G:
                b = self._normal_form(g.bracket(h), G) # This should be a composition
                if b != zero: # It is not a Groebner basis
                    return False
        return True

    def _inner_reduce(self, G):
        """
        Return a self-reduced generating set ``M`` from a (finite)
        generating set ``G``.
        """
        reduced = False
        zero = self._ambient.zero()
        while not reduced:
            M = set()
            reduced = False
            for i,g in enumerate(G):
                red = self._normal_form(g, G[:i] + G[i+1:])
                red /= red.leading_coefficient(term_cmp)
                if red != g:
                    reduced = True
                if red != zero:
                    M.add(red)
            G = list(M)
        return G

    def _inner_reduce_with_map(self, G):
        """
        Return a self-reduced generating set from a (finite)
        generating set ``G`` along with a map from the original element to the
        final reduced element.
        """
        reduced = False
        zero = self._ambient.zero()
        MG = list(G)
        ret_map = {g:g for g in G}
        while not reduced:
            M = set()
            for i,g in enumerate(MG):
                red = self._normal_form(g, MG[:i] + MG[i+1:])
                red /= red.leading_coefficient(term_cmp)
                if red != g:
                    reduced = True
                    ret_map[red] = ret_map[g]
                    del ret_map[g]
                if red != zero:
                    M.add(red)
            MG = list(M)

        # Add in all the generators which were completely reduced
        ret_map = {v:k for k,v in ret_map.items()}
        for g in G:
            if g not in ret_map:
                ret_map[g] = None
        return G, ret_map

    def _useful_compositions(self, l, r):
        """
        Return all useful compositions of ``l`` and ``r`` with respect
        to the currently Groebner basis of ``self``.
        """
        lm,lc = l.leading_item(term_cmp)
        rm,rc = l.leading_item(term_cmp)
        lw = lm.to_word()
        rw = lm.to_word()

        L = self._ambient
        std = L._standard_bracket
        for i in range(min(len(lw), len(rw))):
            # Determine if there is an overlapping factor of ``lw`` and ``rw``.
            # In other words, can we write ``lw = ab`` and ``rw = bc`` for ``b`` non-empty?
            if lw[-i:] != rw[:i]:
                continue

            lhs = lw[:-i]
            #mid = lw[-i:]
            rhs = rw[i:]

            h = std(lw + rhs)
            app_l = self._appliance(lm, h)
            app_r = self._appliance(rm, h)

            lhs = l
            for key,side in app_l:
                if side: # Left
                    lhs = lhs.bracket(L.monomial(key))
                else:
                    lhs = L.monomial(key).bracket(lhs)
            lhs = r
            for key,side in app_r:
                if side: # Left
                    rhs = rhs.bracket(L.monomial(key))
                else:
                    rhs = L.monomial(key).bracket(rhs)

            yield lhs - rhs

    def reduce(self, x, max_size=infinity):
        """
        Return ``x`` modulo ``self``.
        """
        return self._normal_form(x, self.groebner_basis(max_size))

    def _normal_form(self, x, G):
        """
        Return the normal form of ``x`` modulo ``G``.
        """
        L = self._ambient
        zero = L.zero()
        cur = x
        ret = zero

        # TODO: this could be improved by sorting the generators by their
        #   leading monomials and using the descending chain condition
        while cur != zero:
            was_reduced = False
            for g in G:
                gm = g.leading_support(term_cmp)
                hm, hc = cur.leading_item(term_cmp)

                app = self._appliance(gm, hm)
                if app is None:
                    continue

                for key,side in app:
                    if side: # Left
                        g = g.bracket(L.monomial(key))
                    else:
                        g = L.monomial(key).bracket(g)
                g = g * hc / g.leading_coefficient(term_cmp)

                cur -= g
                was_reduced = True
                break

            if not was_reduced:
                ret += cur.leading_term(term_cmp)
                cur -= cur.leading_term(term_cmp)
        return ret

    def _appliance(self, gm, hm):
        """
        Return the output of an applicance, which when applied to ``gm`` has a
        leading term equal to the leading term of ``hm``. If no such
        applicance exists, then return ``None``.
        """
        # Trivial case of the leading monomials are equal
        if gm == hm:
            return []

        # Trivial case where hm is a single letter
        if isinstance(hm, LieGenerator) and gm != hm:
            return None

        # Check to see if gm exists as a subword of hm
        gw = gm.to_word()
        hw = hm.to_word()
        if not is_subword(gw, hw): # It does not exist as a subword
            return None

        #print "\nCurrent:", x
        #print "Reducing {} by {}".format(xm, g)

        L = self._ambient

        # app is the appliance as pairs (key, side)
        # the side is True for left, False for right
        app = []
        cur = hm
        while cur is not None:
            a = cur[0].to_word()
            b = cur[1].to_word()

            if gw == a:
                app.append((cur[1], False))
                break
            if gw == b:
                app.append((cur[0], True))
                break

            if is_subword(gw, a):
                cur = cur[0]
                app.append((cur[1], False))
                continue

            if is_subword(gw, b):
                cur = cur[1]
                app.append((cur[0], True))
                continue

            rf = cur.to_word()
            lyn_fact = lyndon_factorization(rf)
            brackets = reversed(map(L._standard_bracket, lyn_fact))
            app.extend([(x, False) for x in brackets])

            cur = None # This breaks us out of the loop

        return reversed(app)

LieIdealMonoid.Element = LieAlgebraIdeal

class FiniteDimensionalLieAlgebraIdeal(LieAlgebraIdeal):
    """
    The ideal of a finite dimensional Lie algebra `\mathfrak{g}`.
    """
    def groebner_basis(self):
        """
        Return a Groebner-Shirshov basis for ``self``.
        """
        return self.basis()

    def reduce(self, x):
        """
        Return ``x`` modulo ``self``.
        """
        M = self._dense_free_module()
        b = self.basis()
        comp = b.basis_matrix().right_kernel() # The complement basis TODO - this implementation is wrong
        # TODO: convert comp to vectors in self._ambient
        v = M.linear_dependence(comp + b + [x.to_vector()])
        return self._ambient.sum([c * comp[i] for i,c in enumerate(v)])

#####################################################################
## Some helper functions

def term_cmp(x, y):
    """
    Compare the Lyndon word terms ``x`` and ``y`` with grading.
    """
    wx = x.to_word()
    wy = y.to_word()
    if len(wx) > len(wy):
        return 1
    if len(wx) < len(wy):
        return -1
    return cmp(wy, wx)

def is_subword(u, v):
    """
    Helper function to see if ``v`` can be written as ``lur``.
    """
    lu = len(u)
    for i in range(len(v) - lu + 1):
        if all(v[i+j] == val for j,val in enumerate(u)):
            return True
    return False

def lyndon_factorization(w):
    """
    Return the Lyndon factorization of ``w``.
    """
    # We compute the indexes of the factorization.
    n = len(w)
    k = -1
    F = [0]
    while k < n-1:
        i = k+1
        j = k+2
        while j < n:
            if w[i] < w[j]:
                i = k+1
                j += 1
            elif w[i] == w[j]:
                i += 1
                j += 1
            else:
                break
        while k < i:
            F.append(k + j - i + 1)
            k = k + j - i
    return [w[F[i]:F[i+1]] for i in range(len(F)-1)]

def lyndon_factorization_old(w, c):
    """
    Return the Lyndon factorization of ``w`` given by the composition ``c``
    of ``len(w)`` or ``None`` if it does not give a Lyndon factorization.
    """
    prev = 0
    ret = []
    for i in c:
        next = prev + i
        sw = tuple(w[prev:next])
        if not is_lyndon(sw):
            return None
        ret.append(sw)
        prev = next
    return ret

