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
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import ZZ
from sage.categories.monoids import Monoids
from sage.monoids.free_monoid_element import MonoidElement
from sage.combinat.permutation import Permutations
from sage.combinat.composition import Compositions
from sage.algebras.lie_algebras.free_lie_algebra import is_lyndon

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

class LieAlgebraIdeal(MonoidElement):
    r"""
    The ideal of a Lie algebra `\mathfrak{g}`.
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
        self.__lie = lie_algebra
        if not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [lie_algebra(x) for x in gens]

        if len(gens) == 0:
            gens = (lie_algebra.zero_element(),)
        else:
            # Make sure the leading term has a coefficient of 1 (?)
            gens = tuple(map(lambda x: x/x.leading_coefficient(), gens))
        self.__gens = gens
        MonoidElement.__init__(self, LieIdealMonoid(lie_algebra))

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Ideal %s of %s"%(self.__gens, self.__lie)

    def _repr_short(self):
        """
        Short representation for the list of generators.

        EXAMPLES::

        If the string representation of a generator contains a line break,
        the generators are not represented from left to right but from
        top to bottom. This is the case, e.g., for matrices::
        """
        L = []
        has_return = False
        for x in self.gens():
            s = repr(x)
            if '\n' in s:
                has_return = True
                s = s.replace('\n','\n  ')
            L.append(s)
        if has_return:
            return '\n(\n  %s\n)\n'%(',\n\n  '.join(L))
        return '(%s)'%(', '.join(L))

    def gens(self):
        """
        Return the generators of ``self``.
        """
        return self.__gens

    @cached_method
    def groebner_basis(self):
        """
        Return a Groebner-Shirshov basis for ``self``.
        """
        # Currently we just test that we are a Groebner basis
        zero = self.__lie.zero()
        for x in self.__gens:
            for y in self.__gens:
                b = self.reduce(x.bracket(y))
                if b != zero: # It is not a Groebner basis
                    print x
                    print y
                    return False
        return True

    def reduce(self, x):
        """
        Return ``x`` modulo ``self``.
        """
        ret = x
        reducing = True
        zero = self.__lie.zero()

        while reducing:
            reducing = False
            for g in self.gens():
                if ret == zero:
                    return ret
                cur = self._appliance(g, ret)
                if cur is None:
                    continue
                reducing = True
                ret -= cur
                print "ret:", ret
        return ret

    def _appliance(self, g, x):
        """
        Return the applicance which applied to ``cur`` from the associated
        ``split`` has a leading term of ``xm``.
        """
        gm, gc = g.leading_item(term_cmp)
        xm, xc = x.leading_item(term_cmp)

        # Check to see if gm exists as a subword of xm
        l_factor = None
        xw = xm.to_word()
        gw = gm.to_word()
        l_gw = len(gw)
        for i in range(len(xw) - l_gw + 1):
            if all(xw[i+j] == gw[j] for j in range(l_gw)):
                l_factor,r_factor = xw[:i], xw[i+l_gw:]
                break
        if l_factor is None: # It does not exist as a subword
            return None

        #print "\nCurrent:", x
        #print "Reducing {} by {}".format(xm, g)

        L = self.__lie
        cur = gc**-1 * g
        one = L.base_ring().one()
        std = lambda x: L.element_class(L, {L._standard_bracket(x): one})

        # Go through all Lyndon factorizations of the left factors
        for l_comp in Compositions(len(l_factor)):
            l_fact = lyndon_factorization(l_factor, l_comp)
            if l_fact is None:
                continue
            l_fact = map(std, reversed(l_fact)) # Need to apply these in reversed order

            # Go through all Lyndon factorizations of the right factors
            for r_comp in Compositions(len(r_factor)):
                r_fact = lyndon_factorization(r_factor, r_comp)
                if r_fact is None:
                    continue
                r_fact = map(std, r_fact)

                # Loop through all left/right bracketing patterns
                for p in Permutations([0]*len(l_fact) + [1]*len(r_fact)):
                    il,ir = 0,0 # indexing counters
                    temp = cur
                    for i in p:
                        if i == 0: # left side
                            temp = L.bracket(l_fact[il], temp)
                            il += 1
                        else: # i == 1: right side
                            temp = L.bracket(temp, r_fact[ir])
                            ir += 1
                    #print "temp:", temp
                    #print "leading support:", temp.leading_support(term_cmp)

                    # If the leading term is xm, we have found an appliance
                    if temp.leading_support(term_cmp) == xm:
                        return xc / temp.leading_coefficient(term_cmp) * temp
        assert False, "could not find an appliance for {}".format(xm)

LieIdealMonoid.Element = LieAlgebraIdeal

## Some helper functions
########################

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

def lyndon_factorization(w, c):
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

