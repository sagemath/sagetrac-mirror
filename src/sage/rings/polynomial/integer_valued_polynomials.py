# -*- coding: utf-8 -*-
r"""
Integer-valued polynomial rings

AUTHORS:

- Frédéric Chapoton (2013-03): Initial version
"""
# *****************************************************************************
#  Copyright (C) 2013 Frédéric Chapoton
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.arith.misc import (binomial, factorial)
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.rings import Rings
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.family import Family


class IntegerValuedPolynomialRing(CombinatorialFreeModule):
    r"""
    The integer-valued polynomial ring over a base ring.

    Integer-valued polynomial rings are commutative and associative
    algebras, with a basis indexed by non-negative integers.

    The basis used here is given by `B[i] = \binom{i+n}{i}` for `i \in \NN`.

    The product of two monomials
    `B[n_1] \cdot B[n_2]` is given by the sum over ...

    There is a nice formula for the product (ref ?)

    REFERENCES:

    - :wikipedia:`Integer-valued polynomial`

    INPUT:

    - ``R`` -- ring

    EXAMPLES::

        sage: F = IntegerValuedPolynomialRing(QQ); F
        Integer-Valued Polynomial Ring over Rational Field

        sage: F.gen()
        B[1]

        sage: S = IntegerValuedPolynomialRing(ZZ); S
        Integer-Valued Polynomial Ring over Integer Ring
        sage: S.base_ring()
        Integer Ring

        sage: G = IntegerValuedPolynomialRing(S); G
        Integer-Valued Polynomial Ring over Integer-Valued Polynomial
        Ring over Integer Ring
        sage: G.base_ring()
        Integer-Valued Polynomial Ring over Integer Ring

    Integer-valued polynomial rings commute with their base ring::

        sage: K = IntegerValuedPolynomialRing(QQ)
        sage: a = K.gen()
        sage: K.is_commutative()
        True
        sage: L = IntegerValuedPolynomialRing(K)
        sage: c = L.gen()
        sage: L.is_commutative()
        True
        sage: s = a * c^3; s
        B[1]*B[1] + (-6*B[1])*B[2] + 6*B[1]*B[3]
        sage: parent(s)
        Integer-Valued Polynomial Ring over Integer-Valued Polynomial
        Ring over Rational Field

    Integer-valued polynomial rings are commutative::

        sage: c^3 * a == c * a * c * c
        True

    We can also manipulate elements in the basis and coerce elements from our
    base field::

        sage: F = IntegerValuedPolynomialRing(QQ)
        sage: B = F.basis()
        sage: B[2] * B[3]
        3*B[3] - 12*B[4] + 10*B[5]
        sage: 1 - B[2] * B[2] / 2
        B[0] - 1/2*B[2] + 3*B[3] - 3*B[4]
    """
    @staticmethod
    def __classcall_private__(cls, R):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = IntegerValuedPolynomialRing(QQ)
            sage: F2 = IntegerValuedPolynomialRing(QQ)
            sage: F1 is F2
            True
        """
        return super().__classcall__(cls, R)

    def __init__(self, R):
        r"""
        Initialize ``self``.

        INPUT:

        - `R` -- a commutative ring

        EXAMPLES::

            sage: F = IntegerValuedPolynomialRing(QQ); F
            Integer-Valued Polynomial Ring over Rational Field
            sage: TestSuite(F).run()

        TESTS::

            sage: IntegerValuedPolynomialRing(24)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        cat = AlgebrasWithBasis(R).Commutative()
        CombinatorialFreeModule.__init__(self, R, NonNegativeIntegers(),
                                         latex_prefix="",
                                         category=cat)

    def _repr_(self) -> str:
        r"""
        Text representation of this integer-valued polynomial ring.

        EXAMPLES::

            sage: F = IntegerValuedPolynomialRing(QQ)
            sage: F  # indirect doctest
            Integer-Valued Polynomial Ring over Rational Field

            sage: IntegerValuedPolynomialRing(ZZ)
            Integer-Valued Polynomial Ring over Integer Ring
        """
        br = self.base_ring()
        return f"Integer-Valued Polynomial Ring over {br}"

    @cached_method
    def one_basis(self):
        r"""
        Return the number 0, which index of `1` of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = IntegerValuedPolynomialRing(QQ)
            sage: A.one_basis()
            0
            sage: A.one()
            B[0]
        """
        return self.basis().keys()(0)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = IntegerValuedPolynomialRing(QQ)
            sage: A.an_element()
            B[2] + 3*B[3]
        """
        B = self.basis()
        return B[2] + 3 * B[3]

    # mettre ici @cached_method ?
    def product_on_basis(self, n1, n2):
        r"""
        Return the product of basis elements ``n1`` and ``n2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``n1``, ``n2`` -- integers

        EXAMPLES::

            sage: A = IntegerValuedPolynomialRing(QQ)
            sage: A.product_on_basis(0, 1)
            B[1]
        """
        i = ZZ(n1)
        j = ZZ(n2)
        if j < i:
            j, i = i, j

        return self.sum((-1)**k * i.binomial(k) * (i + j - k).binomial(i) *
                        self.basis()[i + j - k] for k in range(i + 1))

    def from_polynomial(self, p):
        """
        Convert a polynomial into the ring of integer-valued polynomials.

        This raises a ``ValueError`` if this is not possible.

        INPUT:

        - `p` -- a polynomial in one variable

        EXAMPLES::

            sage: A = IntegerValuedPolynomialRing(ZZ)
            sage: B = A.basis()
            sage: B[5].polynomial()
            1/120*x^5 + 1/8*x^4 + 17/24*x^3 + 15/8*x^2 + 137/60*x + 1
            sage: A.from_polynomial(_)
            B[5]
            sage: x = polygen(QQ, 'x')
            sage: A.from_polynomial(x)
            -B[0] + B[1]
        """
        B = self.basis()
        x = p.parent().gen()
        remain = p
        result = self.zero()
        while remain != 0:
            N = remain.degree()
            top_coeff = remain.leading_coefficient() * factorial(N)
            try:
                top_coeff = self.base_ring()(top_coeff)
            except TypeError as exc:
                msg = 'not a polynomial with integer'
                msg += f' values: {top_coeff}'
                raise ValueError(msg) from exc
            remain += -top_coeff * binomial(N + x, N)
            result += top_coeff * B[N]
        return result

    def from_h_vector(self, h):
        """
        Convert from some `h`-vector.

        INPUT:

        - ``h`` -- a tuple or vector

        EXAMPLES::

            sage: A = IntegerValuedPolynomialRing(ZZ)
            sage: B = A.basis()
            sage: ex = B[2]+B[4]
            sage: A.from_h_vector(ex.h_vector())
            B[2] + B[4]
        """
        B = self.basis()
        d = len(h) - 1
        m = matrix(QQ, d + 1, d + 1,
                   lambda j, i: (-1)**(d - j) * binomial(d - i, d - j))
        v = vector(QQ, [h[i] for i in range(d + 1)])
        return self.sum(Integer(c) * B[i] for i, c in enumerate(m * v))

    def gen(self):
        r"""
        Return the generator of the algebra.

        EXAMPLES::

            sage: F = IntegerValuedPolynomialRing(ZZ)
            sage: F.gen()
            B[1]
        """
        return self.algebra_generators()[0]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::

            sage: A = IntegerValuedPolynomialRing(ZZ); A
            Integer-Valued Polynomial Ring over Integer Ring
            sage: A.algebra_generators()
            Family (B[1],)
        """
        NonNeg = self.basis().keys()
        return Family([self.monomial(NonNeg(1))])

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        INPUT:

        - ``x`` -- an integer or something convertible

        EXAMPLES::

            sage: R = IntegerValuedPolynomialRing(QQ)
            sage: x = R.gen()
            sage: R(3)
            3*B[0]
            sage: R(x)
            B[1]
        """
        if x in NonNegativeIntegers():
            W = self.basis().keys()
            return self.monomial(W(x))

        P = x.parent()
        if isinstance(P, IntegerValuedPolynomialRing):
            if P is self:
                return x
            if P is not self.base_ring():
                return self.element_class(self, x.monomial_coefficients())

        # ok, not a integer-valued polynomial ring element (or should not be
        # viewed as one).
        if isinstance(x, str):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x, locals=self.gens_dict())
        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self, {})
        return self.from_base_ring_from_one_basis(x)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this integer-valued polynomial ring,

        - anything that coerces to the base ring of this
          integer-valued polynomial ring,

        - any integer-valued polynomial ring on the same variables,
          whose base ring coerces to the base ring of this
          integer-valued polynomial ring.

        EXAMPLES::

            sage: F = IntegerValuedPolynomialRing(GF(7)); F
            Integer-Valued Polynomial Ring over Finite Field of size 7

        Elements of the integer-valued polynomial ring canonically coerce in::

            sage: x = F.gen()
            sage: F.coerce(x*x) # indirect doctest
            6*B[1] + 2*B[2]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            B[0]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to
            Integer-Valued Polynomial Ring over Finite Field of size 7

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*B[0]

        The integer-valued polynomial ring over `\ZZ` on `x` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = IntegerValuedPolynomialRing(ZZ)
            sage: Gx = G.gen()
            sage: z = F.coerce(Gx**2); z
            -B[1] + 2*B[2]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the shuffle
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(x^3+x)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Integer-Valued Polynomial
            Ring over Finite Field of size 7 to Integer-Valued Polynomial
            Ring over Integer Ring
        """
        try:
            R = x.parent()

            # integer-valued polynomial rings in the same variable
            # over any base that coerces in:
            if isinstance(R, IntegerValuedPolynomialRing):
                if self.has_coerce_map_from(R.base_ring()):
                    return self(x)
                raise TypeError("no natural map between bases of"
                                " integer-valued polynomial rings")

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this integer-valued
        # polynomial ring.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        INPUT:

        - ``R`` -- a commutative ring

        The things that coerce into ``self`` are

        - Integer-Valued Polynomial Rings in the same variables over a base
          with a coercion map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = IntegerValuedPolynomialRing(ZZ)
            sage: G = IntegerValuedPolynomialRing(QQ)
            sage: H = IntegerValuedPolynomialRing(ZZ)
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            True
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ,'x'))
            False
        """
        # integer-valued polynomial rings in the same variable over any base
        # that coerces in:
        if isinstance(R, IntegerValuedPolynomialRing):
            return self.base_ring().has_coerce_map_from(R.base_ring())
        return self.base_ring().has_coerce_map_from(R)

    class Element(CombinatorialFreeModule.Element):
        def __call__(self, v):
            """
            Evaluation at some value ``v``

            EXAMPLES::

                 sage: F = IntegerValuedPolynomialRing(ZZ)
                 sage: B = F.gen()
                 sage: f = B**2+4*B+6
                 sage: f(1/3)
                 118/9
            """
            return self.polynomial()(v)

        def shift(self, j=1):
            """
            Shift all indices by `j`.

            INPUT:

            - `j` -- integer (default: 1)

            This corresponds to a summation operator from `0` to `x`.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.gen()
                sage: (B+1).shift()
                B[1] + B[2]
                sage: (B+1).shift(3)
                B[3] + B[4]
            """
            A = self.parent()
            return A.sum(c * A.monomial(i + j) for i, c in self)

        def umbra(self):
            """
            Return the Bernoulli umbra.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.gen()
                sage: (B+1).ombre()
                3/2

            TESTS::

                sage: [(B**n).ombre() for n in range(1, 13)]
                [1/2, 1/6, 0, -1/30, 0, 1/42, 0, -1/30, 0, 5/66, 0, -691/2730]
            """
            return self.shift().derivative_at_minus_one()

        ombre = umbra

        def delta(self):
            r"""
            Return the image by the operator `\Delta`.

            The operator `\Delta` is defined on polynomials by ::

                `f \mapsto f(x+1)-f(x)`.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.basis()
                sage: B[5].delta()
                B[0] + B[1] + B[2] + B[3] + B[4]
            """
            return self.variable_shift() - self

        def variable_shift(self, k=1):
            r"""
            Return the image by the shift of variables.

            On polynomials, the action is the shift
            on variables `x \mapsto x + 1`.

            INPUT:

            - `k` -- integer (default: 1)

            EXAMPLES::

                sage: A = IntegerValuedPolynomialRing(ZZ)
                sage: B = A.basis()
                sage: B[5].variable_shift()
                B[0] + B[1] + B[2] + B[3] + B[4] + B[5]
            """
            if k == 0:
                return self

            A = self.parent()
            resu = A.sum(c * A.monomial(j) for i, c in self
                         for j in range(i + 1))
            if k == 1:
                return resu
            return resu.variable_shift(k - 1)

        def variable_unshift(self, k=1):
            r"""
            Return the image by the unshift of variables.

            On polynomials, the action is the shift
            on variables `x \mapsto x - k`.

            INPUT:

            - `k` -- integer (default: 1)

            EXAMPLES::

                sage: A = IntegerValuedPolynomialRing(ZZ)
                sage: B = A.basis()
                sage: B[5].variable_unshift()
                -B[4] + B[5]
            """
            if k == 0:
                return self

            A = self.parent()
            resu = (A.sum(c * A.monomial(i) for i, c in self) -
                    A.sum(c * A.monomial(i - 1) for i, c in self if i))
            if k == 1:
                return resu
            return resu.variable_unshift(k - 1)

        def derivative_at_minus_one(self):
            """
            Return the derivative at `-1`.

            This is sometimes useful when `-1` is a root.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.gen()
                sage: (B+1).derivative_at_minus_one()
                1
            """
            return QQ.sum(c / QQ(i) for i, c in self if i)

        special_value = derivative_at_minus_one

        def polynomial(self):
            """
            Convert to a standard polynomial in `x`.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.gen()
                sage: (B+1).polynomial()
                x + 2
            """
            x = polygen(QQ, 'x')
            return sum(c * binomial(x + i, i) for i, c in self)

        def to_other_basis(self):
            """
            Convert to the binomial(x,k) basis.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.gen()
                sage: (B+1).to_other_basis()
                B[0] - B[1]
            """
            B = PositiveBasis(self.base_ring())
            base = B.basis()
            return B.sum((-1)**i * c * base[i] for i, c in self)

        def h_vector(self):
            """
            Return the numerator of the generating series of values.

            If ``self`` is an Ehrhart polynomial, this is the `h`-vector.

            EXAMPLES::

                sage: x = polygen(QQ,'x')
                sage: A = IntegerValuedPolynomialRing(ZZ)
                sage: ex = A.from_polynomial((1+x)**3)
                sage: ex.h_vector()
                (0, 1, 4, 1)
            """
            d = max(self.support(), default=-1)
            m = matrix(QQ, d + 1, d + 1,
                       lambda j, i: (-1)**(d - j) * (d - i).binomial(d - j))
            v = vector(QQ, [self.coefficient(i) for i in range(d + 1)])
            return m * v

        def h_polynomial(self):
            """
            Return the `h`-vector as a polynomial.

            EXAMPLES::

                sage: x = polygen(QQ,'x')
                sage: A = IntegerValuedPolynomialRing(ZZ)
                sage: ex = A.from_polynomial((1+x)**3)
                sage: ex.h_polynomial()
                z^3 + 4*z^2 + z
            """
            anneau = PolynomialRing(self.parent().base_ring(), 'z')
            return anneau(list(self.h_vector()))

        def sum(self):
            """
            Return the sum of coefficients.

            This is related to the evaluation at 0.

            EXAMPLES::

                sage: F = IntegerValuedPolynomialRing(ZZ)
                sage: B = F.basis()
                sage: (B[2]*B[4]).sum()
                1
            """
            return sum(c for i, c in self)


# =====     Another basis for the same algebra     =====


class PositiveBasis(CombinatorialFreeModule):
    r"""
    The integer-valued polynomial ring over a base ring.

    Integer-valued polynomial rings are commutative and associative
    algebras, with a basis indexed by integers.

    The basis used here is given by `B[i] = \binom{n}{i}` for `i \in \NN`.

    There is a nice formula for the product (ref ?)
    The product of two monomials
    `B[n_1] \cdot B[n_2]` is given by the sum over ...

    The product of two monomials is a positive linear combination of monomials.

    There is a conversion formula between the two bases::

        `\binom{-x}{i} = (-1)^i \binom{x-1+i}{i}`

    REFERENCES:

    - :wikipedia:`Integer-valued polynomial`

    INPUT:

    - ``R`` -- ring

    EXAMPLES::

        sage: F = PositiveBasis(QQ); F
        Integer-Valued Polynomial Ring over Rational Field

        sage: F.gen()
        B[1]

        sage: S = PositiveBasis(ZZ); S
        Integer-Valued Polynomial Ring over Integer Ring
        sage: S.base_ring()
        Integer Ring

        sage: G = PositiveBasis(S); G
        Integer-Valued Polynomial Ring over Integer-Valued Polynomial
        Ring over Integer Ring
        sage: G.base_ring()
        Integer-Valued Polynomial Ring over Integer Ring

    Integer-valued polynomial rings commute with their base ring::

        sage: K = PositiveBasis(QQ)
        sage: a = K.gen()
        sage: K.is_commutative()
        True
        sage: L = PositiveBasis(K)
        sage: c = L.gen()
        sage: L.is_commutative()
        True
        sage: s = a * c^3; s
        B[1]*B[1] + 6*B[1]*B[2] + 6*B[1]*B[3]
        sage: parent(s)
        Integer-Valued Polynomial Ring over Integer-Valued Polynomial
        Ring over Rational Field

    integer-valued polynomial rings are commutative::

        sage: c^3 * a == c * a * c * c
        True

    We can also manipulate elements in the basis and coerce elements from our
    base field::

        sage: F = PositiveBasis(QQ)
        sage: B = F.basis()
        sage: B[2] * B[3]
        3*B[3] + 12*B[4] + 10*B[5]
        sage: 1 - B[2] * B[2] / 2
        B[0] - 1/2*B[2] - 3*B[3] - 3*B[4]
    """
    @staticmethod
    def __classcall_private__(cls, R):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = PositiveBasis(QQ)
            sage: F2 = PositiveBasis(QQ)
            sage: F1 is F2
            True
        """
        return super().__classcall__(cls, R)

    def __init__(self, R):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = PositiveBasis(QQ); F
            Integer-Valued Polynomial Ring over Rational Field
            sage: TestSuite(F).run()

        TESTS::

            sage: PositiveBasis(24)
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        cat = AlgebrasWithBasis(R).Commutative()
        CombinatorialFreeModule.__init__(self, R, NonNegativeIntegers(),
                                         latex_prefix="",
                                         category=cat)

    def _repr_(self) -> str:
        r"""
        Text representation of this integer-valued polynomial ring.

        EXAMPLES::

            sage: F = PositiveBasis(QQ)
            sage: F  # indirect doctest
            Integer-Valued Polynomial Ring over Rational Field

            sage: PositiveBasis(ZZ)
            Integer-Valued Polynomial Ring over Integer Ring
        """
        br = self.base_ring()
        return f"Integer-Valued Polynomial Ring over {br}"

    @cached_method
    def one_basis(self):
        r"""
        Return the number 0, which index of `1` of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = PositiveBasis(QQ)
            sage: A.one_basis()
            0
            sage: A.one()
            B[0]
        """
        return self.basis().keys()(0)

    def _an_element_(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = PositiveBasis(QQ)
            sage: A.an_element()
            B[2] + 3*B[3]
        """
        B = self.basis()
        return B[2] + 3 * B[3]

    # mettre ici @cached_method ?
    def product_on_basis(self, n1, n2):
        r"""
        Return the product of basis elements ``n1`` and ``n2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``n1``, ``n2`` -- integers

        EXAMPLES::

            sage: A = PositiveBasis(QQ)
            sage: A.product_on_basis(0, 1)
            B[1]
        """
        i = ZZ(n1)
        j = ZZ(n2)
        if j < i:
            j, i = i, j

        return self.sum(binomial(i, k) * binomial(i + j - k, i) *
                        self.basis()[i + j - k]
                        for k in range(i + 1))

    def from_polynomial(self, p):
        """
        Convert a polynomial into the ring of integer-valued polynomials.

        This raises a ``ValueError`` if this is not possible.

        INPUT:

        - `p` -- a polynomial in one variable

        EXAMPLES::

            sage: A = PositiveBasis(ZZ)
            sage: B = A.basis()
            sage: B[5].polynomial()
            1/120*x^5 - 1/12*x^4 + 7/24*x^3 - 5/12*x^2 + 1/5*x
            sage: A.from_polynomial(_)
            B[5]
            sage: x = polygen(QQ, 'x')
            sage: A.from_polynomial(x)
            B[1]
        """
        B = self.basis()
        x = p.parent().gen()
        remain = p
        result = self.zero()
        while remain != 0:
            N = remain.degree()
            top_coeff = remain.leading_coefficient() * factorial(N)
            try:
                top_coeff = self.base_ring()(top_coeff)
            except TypeError as exc:
                raise ValueError('not a polynomial with integer values') from exc
            remain += -top_coeff * binomial(x, N)
            result += top_coeff * B[N]
        return result

    def gen(self):
        r"""
        Return the generator of the algebra.

        EXAMPLES::

            sage: F = PositiveBasis(ZZ)
            sage: F.gen()
            B[1]
        """
        return self.algebra_generators()[0]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::

            sage: A = PositiveBasis(ZZ); A
            Integer-Valued Polynomial Ring over Integer Ring
            sage: A.algebra_generators()
            Family (B[1],)
        """
        NonNeg = self.basis().keys()
        return Family([self.monomial(NonNeg(1))])

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = PositiveBasis(QQ)
            sage: x = R.gen()
            sage: R(3)
            3*B[0]
            sage: R(x)
            B[1]
        """
        if x in NonNegativeIntegers():
            W = self.basis().keys()
            return self.monomial(W(x))

        P = x.parent()
        if isinstance(P, PositiveBasis):
            if P is self:
                return x
            if P is not self.base_ring():
                return self.element_class(self, x.monomial_coefficients())

        # ok, not a integer-valued polynomial ring element (or should not be
        # viewed as one).
        if isinstance(x, str):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x, locals=self.gens_dict())
        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self, {})
        return self.from_base_ring_from_one_basis(x)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this integer-valued polynomial ring,

        - anything that coerces to the base ring of this
          integer-valued polynomial ring,

        - any integer-valued polynomial ring on the same variables,
          whose base ring coerces to the base ring of this
          integer-valued polynomial ring.

        EXAMPLES::

            sage: F = PositiveBasis(GF(7)); F
            Integer-Valued Polynomial Ring over Finite Field of size 7

        Elements of the integer-valued polynomial ring canonically coerce in::

            sage: x = F.gen()
            sage: F.coerce(x*x) # indirect doctest
            B[1] + 2*B[2]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            B[0]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to
            Integer-Valued Polynomial Ring over Finite Field of size 7

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*B[0]

        The integer-valued polynomial ring over `\ZZ` on `x` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = PositiveBasis(ZZ)
            sage: Gx = G.gen()
            sage: z = F.coerce(Gx**2); z
            B[1] + 2*B[2]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the shuffle
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(x^3+x)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Integer-Valued Polynomial
            Ring over Finite Field of size 7 to Integer-Valued Polynomial
            Ring over Integer Ring
        """
        try:
            R = x.parent()

            # integer-valued polynomial rings in the same variable
            # over any base that coerces in:
            if isinstance(R, PositiveBasis):
                if self.has_coerce_map_from(R.base_ring()):
                    return self(x)
                raise TypeError("no natural map between bases of"
                                " integer-valued polynomial rings")

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this integer-valued
        # polynomial ring.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - Integer-Valued Polynomial Rings in the same variables over a base
          with a coercion map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = PositiveBasis(ZZ)
            sage: G = PositiveBasis(QQ)
            sage: H = PositiveBasis(ZZ)
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            True
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
            sage: F.has_coerce_map_from(PolynomialRing(ZZ,'x'))
            False
        """
        # integer-valued polynomial rings in the same variable over any base
        # that coerces in:
        if isinstance(R, PositiveBasis):
            return self.base_ring().has_coerce_map_from(R.base_ring())
        return self.base_ring().has_coerce_map_from(R)

    class Element(CombinatorialFreeModule.Element):
        def __call__(self, v):
            """
            Evaluation at some value ``v``

            EXAMPLES::

                 sage: F = PositiveBasis(ZZ)
                 sage: B = F.gen()
                 sage: f = B**2+4*B+6
                 sage: f(1/3)
                 67/9
            """
            return self.polynomial()(v)

        def shift(self):
            """
            Shift all indices by 1.

            EXAMPLES::

                sage: F = PositiveBasis(ZZ)
                sage: B = F.gen()
                sage: (B+1).shift()
                B[1] + B[2]
            """
            A = self.parent()
            return A.sum(c * A.monomial(i + 1) for i, c in self)

        def polynomial(self):
            """
            Convert to a standard polynomial in `x`.

            EXAMPLES::

                sage: F = PositiveBasis(ZZ)
                sage: B = F.gen()
                sage: (B+1).polynomial()
                x + 1
            """
            x = polygen(QQ, 'x')
            return sum(c * binomial(x, i) for i, c in self)


# def produ(D, j, d, i):
#     """
#     Compute the product of monomials `H^(D)_j` and `H^(d)_i`.

#     This means binomial(x+j,D) and binomial(x+i,d).

#     The result is given by coordinates in the basis `H^(D+d)`.

#     EXAMPLES::

#         sage: produ(5,1,4,1)
#         (0, 1, 20, 60, 40, 5, 0, 0, 0, 0)
#     """
#     return vector(ZZ, [0] * i + [binomial(D + i - j, D - k) *
#                                  binomial(d + j - i, k)
#                                  for k in range(d + D + 1 - i)])


# def produ_ok(D, j, d, i):
#     """
#     EXAMPLES::

#         sage: produ_ok(5,1,4,1)
#         (0, 1, 20, 60, 40, 5, 0, 0, 0, 0)
#         sage: all(produ_ok(7,7,6,i) == produ(7,7,6,i) for i in range(6))
#         True
#         sage: all(produ_ok(7,i,6,2) == produ(7,i,6,2) for i in range(6))
#         True
#     """
#     A = IntegerValuedPolynomialRing(ZZ)
#     x = polygen(QQ, 'x')
#     return A.from_polynomial(binomial(i + x, d) *
#                              binomial(x + j, D)).h_vector()
