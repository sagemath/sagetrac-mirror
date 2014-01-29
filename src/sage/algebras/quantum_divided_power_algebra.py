r"""
A minimal implementation of the quantum divided power algebra as a
graded algebra with basis.

AUTHOR:

- Bruce Westbury
"""
#*****************************************************************************
#  Copyright (C) 2011 Bruce W. Westbury <brucewestbury@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import GradedAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.q_analogues import q_binomial
from sage.categories.rings import Rings
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ


class QuantumDividedPowerAlgebra(CombinatorialFreeModule):
    r"""
    An example of a graded algebra with basis: the quantum divided power algebra.

    This class illustrates a minimal implementation of the quantum
    divided power algebra.
    """

    def __init__(self, R=None, q=None):
        if q is None:
            if R is None:
                self.q = ZZ['q'].gen()
            else:
                if not R in Rings():
                    raise ValueError('R is not a ring')
                self.q = R['q'].gen()
        else:
            if not R is None:
                if not R in Rings() or q not in R:
                    raise ValueError('q is not in the chosen ring')
                else:
                    self.q = q
            else:
                self.q = q

        S = self.q.parent()

        CombinatorialFreeModule.__init__(self, S, NonNegativeIntegers(),
                                         category=GradedAlgebrasWithBasis(S))

    def _repr_(self):
        return "The quantum divided power algebra over %s, %s" % (self.base_ring(),self.q)

    @cached_method
    def one(self):
        """
        Returns the unit of the algebra
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: from sage.algebras.all import QuantumDividedPowerAlgebra
            sage: A = QuantumDividedPowerAlgebra(ZZ)
            sage: A.one()
            B[0]
        """
        u = NonNegativeIntegers.from_integer(0)
        return self.monomial(u)

    def product_on_basis(self, left, right):
        r"""
        Product, on basis elements, as per :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

        INPUT:

        - ``left``, ``right`` - non-negative integers determining monomials (as the
          exponents of the generators) in this algebra

        OUTPUT: the product of the two corresponding monomials, as an
        element of ``self``.

        EXAMPLES::

            sage: from sage.algebras.all import QuantumDividedPowerAlgebra
            sage: A = QuantumDividedPowerAlgebra(ZZ)
            sage: B = A.basis()
            sage: B[2]*B[3]
            (q^6+q^5+2*q^4+2*q^3+2*q^2+q+1)*B[5]
         """
        return self.term(left + right, q_binomial(left + right, left, self.q))

    def degree_on_basis(self, t):
        """
        The degree of the element determined by the integer ``t`` in
        this graded module.

        INPUT:

        - ``t`` -- the index of an element of the basis of this module,
          i.e. a non-negative integer

        OUTPUT: an integer, the degree of the corresponding basis element

        EXAMPLES::

            sage: from sage.algebras.all import QuantumDividedPowerAlgebra
            sage: A = QuantumDividedPowerAlgebra(ZZ)
            sage: A.degree_on_basis(3)
            3
            sage: type(A.degree_on_basis(2))
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(t)

    @cached_method
    def algebra_generators(self):
        r"""
        The generators of this algebra, as per :meth:`Algebras.ParentMethods.algebra_generators`.

        """
        return Family(NonNegativeIntegers())
