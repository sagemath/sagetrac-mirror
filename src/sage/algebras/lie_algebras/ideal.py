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


LieIdealMonoid.Element = LieAlgebraIdeal

