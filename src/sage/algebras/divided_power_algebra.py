r"""
A minimal implementation of the divided power algebra as a graded Hopf
algebra with basis.

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
from sage.categories.all import GradedHopfAlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.categories.rings import Rings
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.rings.arith import binomial
from sage.rings.integer import Integer


class UnivariateDividedPowerAlgebra(CombinatorialFreeModule):
    r"""
    An example of a graded Hopf algebra with basis: the divided
    power algebra in one variable.

    This class illustrates a minimal implementation of the divided
    power algebra.

    Let `R` be a commutative ring. The divided power algebra `DPA(R)`
    in one variable `t` over `R` is the free `R`-module with basis
    `(t_0, t_1, t_2, t_3, \ldots)`. This free `R`-module `DPA(R)` is
    made into an `R`-algebra by setting

    .. MATH::

        t_i t_j = \binom{i+j}{i} t_{i+j} ,

    and into a coalgebra by setting

    .. MATH::

        \Delta t_n = \sum_{k=0}^{n} t_k \otimes t_{n-k} .

    These two structures combine to form a Hopf algebra structure
    on `DPA(R)`. Its counit sends every `t_n` to `\delta_{n, 0}`;
    its unity is `t_0`; its antipode sends every `t_n` to
    `(-1)^n t_n`. This Hopf algebra is graded, with `t_n` having
    degree `n`; it is furthermore commutative and cocommutative.

    When `R` is a `\QQ`-algebra (that is, the positive integers are
    invertible in `R`), the divided power algebra `DPA(R)` is
    isomorphic to the polynomial ring `R[t]` by the Hopf algebra
    isomorphism which sends every `t_i` to `t^i / i!`. In other
    cases, it can have a structure significantly differing from
    that of the polynomial ring (for example, it fails to be
    Noetherian over a field `R` of positive characteristic).

    The divided power algebra `DPA(R)` is always isomorphic to the
    shuffle algebra
    (:class:`sage.algebras.shuffle_algebra.ShuffleAlgebra`)
    with one generator `x` over `R`. The isomorphism (implemented
    here as the method :meth:`to_shuffle_algebra`, with inverse
    map :meth:`from_shuffle_algebra`) sends `t_i` to the word
    `x x \cdots x` (with `i` factors `x`).

    .. NOTE::

        Due to this being a toy implementation (and essentially a
        particular case of the shuffle algebra), the univariate
        divided power algebra is not available in the global
        namespace for immediate interactive use. Instead, it
        needs to be explicitly imported before using::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(Zmod(9)); A
            The divided power algebra over Ring of integers modulo 9

        If you want to implement an algebra which needs not be
        imported in order to be called, you need to add a
        ``lazy_import`` statement to ``src/sage/algebras/all.py``.
        In the case of the univariate divided power algebra, it
        would look as follows::

            lazy_import('sage.algebras.divided_power_algebra', 'UnivariateDividedPowerAlgebra')

    INPUT:

    - ``R``: base ring (a commutative ring).

    OUTPUT:

    The univariate divided power algebra `DPA(R)` over `R`, as
    a graded Hopf algebra with basis.

    EXAMPLES::

        sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
        sage: A = UnivariateDividedPowerAlgebra(ZZ); A
        The divided power algebra over Integer Ring
        sage: TestSuite(A).run()
        sage: A_bas = A.basis()
        sage: A.one()
        B[0]
        sage: A_bas[2] * A_bas[5]
        21*B[7]
        sage: A_bas[2].coproduct()
        B[0] # B[2] + B[1] # B[1] + B[2] # B[0]
        sage: A_bas[2].antipode()
        B[2]
        sage: A_bas[3].antipode()
        -B[3]
        sage: f = A.to_shuffle_algebra(letter='u')(3*A_bas[2] - 4*A_bas[3])
        sage: f
        3*B[word: uu] - 4*B[word: uuu]
        sage: A.from_shuffle_algebra(letter='u')(f)
        3*B[2] - 4*B[3]
    """
    def __init__(self, R):
        if not R in Rings():
            raise ValueError('R is not a ring')
        GHWBR = GradedHopfAlgebrasWithBasis(R)
        CombinatorialFreeModule.__init__(self, R, NonNegativeIntegers(),
                                         category=GHWBR)

    def _repr_(self):
        return "The divided power algebra over %s" % (self.base_ring())

    @cached_method
    def one(self):
        """
        Return the unit of the algebra ``self``,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(ZZ)
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

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: B = UnivariateDividedPowerAlgebra(ZZ).basis()
            sage: B[3]*B[4]
            35*B[7]
            sage: B = UnivariateDividedPowerAlgebra(Zmod(5)).basis()
            sage: B[3]*B[4]
            0
        """
        return self.term(left + right, self._base(binomial(left + right, left)))

    def coproduct_on_basis(self, t):
        r"""
        Return the coproduct of the basis element of ``self`` indexed
        by the nonnegative integer ``t``.

        EXAMPLES::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(ZZ)
            sage: B = A.basis()
            sage: A.coproduct(B[4])
            B[0] # B[4] + B[1] # B[3] + B[2] # B[2] + B[3] # B[1] + B[4] # B[0]
            sage: A.coproduct(B[0])
            B[0] # B[0]
        """
        AA = self.tensor(self)
        return AA.sum_of_monomials(((k, t - k)
                                    for k in range(t + 1)))

    def counit_on_basis(self, t):
        """
        Return the counit of the basis element of ``self`` indexed
        by the nonnegative integer ``t``.

        EXAMPLES::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(ZZ)
            sage: B = A.basis()
            sage: A.counit(B[3])
            0
            sage: A.counit(B[0])
            1
        """
        if t == 0:
            return self.base_ring().one()

        return self.base_ring().zero()

    def antipode_on_basis(self, t):
        """
        Return the antipode of the basis element of ``self`` indexed
        by the nonnegative integer ``t``.

        EXAMPLES::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(ZZ)
            sage: B = A.basis()
            sage: A.antipode(B[4]+B[5])
            B[4] - B[5]
            sage: A.antipode(2*B[0]-B[1])
            2*B[0] + B[1]
        """
        if t % 2 == 0:
            return self.monomial(t)

        return self.term(t, -1)

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

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(ZZ)
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

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(ZZ)
            sage: A.algebra_generators()
            Family (Non negative integers)
        """
        return Family(NonNegativeIntegers())

    @cached_method
    def to_shuffle_algebra(self, letter="x"):
        r"""
        Return the canonical isomorphism from the univariate divided
        power algebra ``self`` to the shuffle algebra in one
        generator ``letter`` over the base ring of ``self``.

        See :class:`UnivariateDividedPowerAlgebra` for a definition
        of this isomorphism.

        EXAMPLES::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(QQ)
            sage: A_bas = A.basis()
            sage: tosh_x = A.to_shuffle_algebra()
            sage: tosh_x(7*A.one() + 8*A_bas[3] - 9*A_bas[4])
            7*B[word: ] + 8*B[word: xxx] - 9*B[word: xxxx]
            sage: tosh_u = A.to_shuffle_algebra(letter="u")
            sage: tosh_u(7*A.one() + 8*A_bas[3] - 9*A_bas[4])
            7*B[word: ] + 8*B[word: uuu] - 9*B[word: uuuu]
        """
        from sage.algebras.shuffle_algebra import ShuffleAlgebra
        ShA = ShuffleAlgebra(self._base, [letter])
        words = ShA.basis().keys()
        def the_map_on_the_basis(n):
            return ShA.monomial(words(letter * n))
        return self.module_morphism(on_basis=the_map_on_the_basis,
                                    codomain=ShA)

    @cached_method
    def from_shuffle_algebra(self, letter="x"):
        r"""
        Return the canonical isomorphism from the univariate divided
        power algebra ``self`` to the shuffle algebra in one
        generator ``letter`` over the base ring of ``self``.

        See :class:`UnivariateDividedPowerAlgebra` for a definition
        of this isomorphism.

        EXAMPLES::

            sage: from sage.algebras.divided_power_algebra import UnivariateDividedPowerAlgebra
            sage: A = UnivariateDividedPowerAlgebra(QQ)
            sage: frosh_x = A.from_shuffle_algebra()
            sage: ShA_x = ShuffleAlgebra(QQ, 'x')
            sage: frosh_x(3 * ShA_x(Word('')) - 6 * ShA_x(Word('xx')))
            3*B[0] - 6*B[2]
            sage: frosh_u = A.from_shuffle_algebra(letter="u")
            sage: ShA_u = ShuffleAlgebra(QQ, 'u')
            sage: frosh_u(3 * ShA_u(Word('')) - 6 * ShA_u(Word('uu')))
            3*B[0] - 6*B[2]
        """
        from sage.algebras.shuffle_algebra import ShuffleAlgebra
        ShA = ShuffleAlgebra(self._base, [letter])
        def the_map_on_the_basis(word):
            return self.monomial(len(word))
        return ShA.module_morphism(on_basis=the_map_on_the_basis,
                                   codomain=self)


