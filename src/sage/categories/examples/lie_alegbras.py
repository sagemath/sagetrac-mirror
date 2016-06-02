r"""
Examples of a Lie algebra
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import LieAlgebras
from sage.combinat.free_module import CombinatorialFreeModule

class AbelianLieAlgebra(CombinatorialFreeModule):
    r"""
    An example of a Lie algebra: the abelian Lie algebra

    This class illustrates a minimal implementation of a Lie algebra.
    """
    def __init__(self, R, gens = ("a", "b", "c")):
        """
        EXAMPLES::

            sage: A = LieAlgebras(QQ).example(); A
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
            sage: TestSuite(A).run()
        """
        CombinatorialFreeModule.__init__(self, R, gens, category=LieAlgebras(R).WithBasis())
        self._assign_names(gens)

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).example()
            sage: L._construct_UEA()
            Multivariate Polynomial Ring in a, b, c over Rational Field
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(self.base_ring(), self.variable_names())

    def _repr_(self):
        """
        EXAMPLES::

            sage: LieAlgebras(QQ).example()
            An example of a Lie algebra: the abelian Lie algebra on the generators (B['a'], B['b'], B['c']) over Rational Field
        """
        return "An example of a Lie algebra: the abelian Lie algebra on the" \
               " generators {} over {}".format(self.gens(), self.base_ring())

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).example()
            sage: L.gens()
            (B['a'], B['b'], B['c'])
        """
        return tuple(self.basis())

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket on basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L = LieAlgebras(QQ).example()
            sage: L.bracket_on_basis('a', 'c')
            0
        """
        return self.zero()

    class Element(CombinatorialFreeModule.Element):
        def lift(self):
            """
            Return the lift of ``self`` to the universal enveloping algebra.

            EXAMPLES::

                sage: L = LieAlgebras(QQ).example()
                sage: elt = L.an_element()
                sage: elt.lift()
                2*a + 2*b + 3*c
            """
            UEA = self.parent().universal_enveloping_algebra()
            gens_dict = UEA.gens_dict()
            return UEA.sum(c * gens_dict[t] for t, c in self)

Example = AbelianLieAlgebra

