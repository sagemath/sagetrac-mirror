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


class UnivariateQuantumDividedPowerAlgebra(CombinatorialFreeModule):
    r"""
    An example of a graded algebra with basis: the quantum divided power
    algebra in one variable.

    This class illustrates a minimal implementation of the quantum
    divided power algebra.

    Let `R` be a commutative ring, and let `q \in R`. The quantum
    divided power algebra `DPA(R, q)` in one variable `t` over `R`
    with quantum parameter `q` is the free `R`-module with basis
    `(t_0, t_1, t_2, t_3, \ldots)`. This free `R`-module `DPA(R, q)`
    is made into an `R`-algebra by setting

    .. MATH::

        t_i t_j = \binom{i+j}{i}_q t_{i+j} ,

    where `\binom{n}{k}_q` denotes the `q`-binomial coefficient
    (:func:`sage.combinat.q_analogues.q_binomial`). This `R`-algebra
    is commutative and has `t_0` as its unity; it is furthermore
    graded by giving `t_n` the degree `n`.

    When `1!_q`, `2!_q`, `3!_q`, ... are invertible in `R` (this
    holds, for example, when `R` is a field and `q` is not a root
    of unity), the quantum divided power algebra `DPA(R, q)` is
    isomorphic to the polynomial ring `R[t]` by the algebra
    isomorphism which sends every `t_i` to `t^i / i!_q`.

    .. NOTE::

        Due to this being a toy implementation, the univariate
        quantum divided power algebra is not available in the
        global namespace for immediate interactive use. Instead,
        it needs to be explicitly imported before using::

            sage: from sage.algebras.quantum_divided_power_algebra import UnivariateQuantumDividedPowerAlgebra
            sage: A = UnivariateQuantumDividedPowerAlgebra(Zmod(9), Zmod(9)(2)); A
            The quantum divided power algebra over Ring of integers modulo 9
             with quantum parameter 1

        If you want to implement an algebra which needs not be
        imported in order to be called, you need to add a
        ``lazy_import`` statement to ``src/sage/algebras/all.py``.
        In the case of the univariate quantum divided power
        algebra, it would look as follows::

            lazy_import('sage.algebras.quantum_divided_power_algebra', 'UnivariateQuantumDividedPowerAlgebra')

    INPUT:

    - ``R`` (default: `\ZZ`): base ring (a commutative ring).

    - ``q`` (default: the polynomial variable `q` of the
      polynomial ring `R[q]`): an element of `R`.

    OUTPUT:

    The univariate quantum divided power algebra `DPA(R, q)` over
    `R` with quantum parameter `q`, as a graded algebra with basis.

    EXAMPLES::

        sage: from sage.algebras.quantum_divided_power_algebra import UnivariateQuantumDividedPowerAlgebra
        sage: A = UnivariateQuantumDividedPowerAlgebra(ZZ, 1); A
        The quantum divided power algebra over Integer Ring
         with quantum parameter 1
        sage: TestSuite(A).run()
        sage: A_bas = A.basis()
        sage: A.one()
        B[0]
        sage: A_bas[2] * A_bas[5]
        21*B[7]

        sage: A = UnivariateQuantumDividedPowerAlgebra(ZZ, 2); A
        The quantum divided power algebra over Integer Ring
         with quantum parameter 2
        sage: TestSuite(A).run()
        sage: A_bas = A.basis()
        sage: A.one()
        B[0]
        sage: A_bas[2] * A_bas[5]
        2667*B[7]

        sage: A = UnivariateQuantumDividedPowerAlgebra(QQ, -1); A
        The quantum divided power algebra over Rational Field
         with quantum parameter -1
        sage: TestSuite(A).run()
        sage: A.base_ring()
        Rational Field
        sage: A_bas = A.basis()
        sage: A.one()
        B[0]
        sage: A_bas[2] * A_bas[5]
        3*B[7]

        sage: A = UnivariateQuantumDividedPowerAlgebra(QQ); A
        The quantum divided power algebra over Univariate
         Polynomial Ring in q over Rational Field with quantum
         parameter q
        sage: TestSuite(A).run()
        sage: A.base_ring()
        Univariate Polynomial Ring in q over Rational Field
        sage: A_bas = A.basis()
        sage: A.one()
        B[0]
        sage: A_bas[2] * A_bas[5]
        (q^10+q^9+2*q^8+2*q^7+3*q^6+3*q^5+3*q^4+2*q^3+2*q^2+q+1)*B[7]

        sage: R = LaurentPolynomialRing(QQ, ['x', 'y'])
        sage: x, y = R.gens()
        sage: A = UnivariateQuantumDividedPowerAlgebra(R, x-y); A
        The quantum divided power algebra over Multivariate
         Laurent Polynomial Ring in x, y over Rational Field with
         quantum parameter x - y
        sage: TestSuite(A).run()  # not tested -- seems to freeze
        sage: A.base_ring()
        Multivariate Laurent Polynomial Ring in x, y over Rational
         Field
        sage: A_bas = A.basis()
        sage: A.one()
        B[0]
        sage: A_bas[2] * A_bas[1]
        (x^2-2*x*y+y^2+x-y+1)*B[3]
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
                    self.q = R(q)
            else:
                self.q = q

        S = self.q.parent()

        CombinatorialFreeModule.__init__(self, S, NonNegativeIntegers(),
                                         category=GradedAlgebrasWithBasis(S))

    def _repr_(self):
        return "The quantum divided power algebra over %s with quantum parameter %s" % (self.base_ring(), self.q)

    @cached_method
    def one(self):
        """
        Return the unit of the algebra ``self``,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: from sage.algebras.quantum_divided_power_algebra import UnivariateQuantumDividedPowerAlgebra
            sage: A = UnivariateQuantumDividedPowerAlgebra(ZZ)
            sage: A.one()
            B[0]
        """
        u = NonNegativeIntegers.from_integer(0)
        return self.monomial(u)

    def product_on_basis(self, left, right):
        r"""
        Return the product of the basis elements of ``self`` indexed
        by the nonnegative integers ``left`` and ``right``.

        As per :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis`.

        INPUT:

        - ``left``, ``right`` - non-negative integers determining
          monomials in this algebra

        OUTPUT:

        the product of the two corresponding monomials, as an element
        of ``self``.

        EXAMPLES::

            sage: from sage.algebras.quantum_divided_power_algebra import UnivariateQuantumDividedPowerAlgebra
            sage: A = UnivariateQuantumDividedPowerAlgebra(ZZ)
            sage: B = A.basis()
            sage: B[2]*B[3]
            (q^6+q^5+2*q^4+2*q^3+2*q^2+q+1)*B[5]
        """
        return self.term(left + right, q_binomial(left + right, left, self.q))

    def degree_on_basis(self, t):
        """
        The degree of the basis element of ``self`` indexed by the
        nonnegative integer ``t`` in the graded module ``self``.

        INPUT:

        - ``t`` -- the index of an element of the basis of ``self``,
          i.e. a non-negative integer

        OUTPUT:

        an integer, the degree of the corresponding basis element

        EXAMPLES::

            sage: from sage.algebras.quantum_divided_power_algebra import UnivariateQuantumDividedPowerAlgebra
            sage: A = UnivariateQuantumDividedPowerAlgebra(ZZ)
            sage: A.degree_on_basis(3)
            3
            sage: type(A.degree_on_basis(2))
            <type 'sage.rings.integer.Integer'>
        """
        return Integer(t)

    @cached_method
    def algebra_generators(self):
        r"""
        A family of indices (in this case, nonnegative integers) such
        that the basis elements of ``self`` indexed by these indices
        generate the algebra ``self``.

        As per :meth:`Algebras.ParentMethods.algebra_generators`.

        EXAMPLES::

            sage: from sage.algebras.quantum_divided_power_algebra import UnivariateQuantumDividedPowerAlgebra
            sage: A = UnivariateQuantumDividedPowerAlgebra(ZZ)
            sage: A.algebra_generators()
            Family (Non negative integers)
        """
        return Family(NonNegativeIntegers())
