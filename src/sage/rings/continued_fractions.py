r"""
Continued fractions

A continued fraction is a representation of a real number in terms of a sequence
of integers denoted `[a_0; a_1, a_2, \ldots]`. The well known decimal expansion
is another way of representing a real number by a sequence of integers. The
value of a continued fraction is defined recursively as:

.. MATH::

    [a_0; a_1, a_2, \ldots] = a_0 + \frac{1}{[a_1; a_2, \ldots]}

In this expansion, all coefficients `a_n` are integers and only the value `a_0`
may be non positive. Note that `a_0` is nothing else but the floor (this remark
provides a way to build the continued fraction expansion from a given real
number).

It is quite remarkable that

- finite expansions correspond to rationals
- ultimately periodic expansions correspond to quadratic numbers (ie numbers of
  the form `a + b \sqrt{D}` with `a` and `b` rationals and `D` square free
  integer)
- two real numbers `x` and `y` have the same tail (up to a shift) if and only if
  there are integers `a,b,c,d` with `|ad - bc| = 1` and such that
  `y = (ax + b) / (cx + d)`.

Moreover, the rational numbers obtained by truncation of the expansion of a real
number gives its so-called best approximations. For more informations on
continued fractions, you may have a look at :wikipedia:`Continued_fraction`.

EXAMPLES:

If you want to create the continued fraction of some real number you may either
use its method continued_fraction (if it exists) or call ``CFF`` (which
stands for real continued fraction) or the function
:func:`continued_fraction`::

    sage: cf = CFF(13/27); cf
    [0; 2, 13]
    sage: 0 + 1/(2 + 1/13)
    13/27
    sage: cf.value()
    13/27

    sage: (22/45).continued_fraction()
    [0; 2, 22]
    sage: 0 + 1/(2 + 1/22)
    22/45

    sage: continued_fraction(-132/17)
    [-8; 4, 4]
    sage: -8 + 1/(4 + 1/4)
    -132/17

It is also possible to create a continued fraction from a list of digits::

    sage: cf = CFF([-3,1,2,3,4,1,2,1])
    sage: cf.value()
    -465/202

For quadratic numbers, the syntax is quite similar and the digits are
represented as a pair of tuples, the preperiod and the period::

    sage: K.<sqrt2> = QuadraticField(2)
    sage: cf = CFF(sqrt2); cf
    [1; (2)*]
    sage: cf.value()
    sqrt2
    sage: cf.preperiod()
    (1,)
    sage: cf.period()
    (2,)

    sage: (3*sqrt2 + 1/2).continued_fraction()
    [4; (1, 2, 1, 7)*]

    sage: cf = CFF([(1,2,3),(1,4)]); cf
    [1; 2, 3, (1, 4)*]
    sage: cf.value()
    -2/23*sqrt2 + 36/23

For an irrational and non quadratic number the continued fraction has, in general, no
particular structure. It is still possible to make computations::

    sage: cf1 = continued_fraction(pi); cf1
    [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
    sage: cf2 = continued_fraction(2*pi/7 + tan(2/5)); cf2
    [1; 3, 8, 3, 1, 33, 3, 1, 2, 1, 1, 1, 5, 2, 17, 3, 2, 1, 33, 1, ...]
    sage: cf1 + cf2
    [4; 2, 6, 13, 7, 2, 8, 2, 5, 1, 1, 5, 3, 9, 2, 14, 1, 1, 36, 2, ...]

On the following we can remark how the tail may change even in the same
quadratic field::

    sage: for i in xrange(20): print CFF(i*sqrt2)
    [0]
    [1; (2)*]
    [2; (1, 4)*]
    [4; (4, 8)*]
    [5; (1, 1, 1, 10)*]
    [7; (14)*]
    ...
    [24; (24, 48)*]
    [25; (2, 5, 6, 5, 2, 50)*]
    [26; (1, 6, 1, 2, 3, 2, 26, 2, 3, 2, 1, 6, 1, 52)*]

Nevertheless, the tail is preserved under invertible integer homographies::

    sage: apply_homography =  lambda m,z: (m[0,0]*z + m[0,1]) / (m[1,0]*z+m[1,1])
    sage: m1 = SL2Z.random_element()
    sage: m2 = SL2Z.random_element()
    sage: a = sqrt2/3
    sage: a.continued_fraction()
    [0; 2, (8, 4)*]
    sage: b = apply_homography(m1, a)
    sage: b.continued_fraction()
    [0; 1, 2, 1, 1, 1, 1, 6, (8, 4)*]
    sage: c = apply_homography(m2, a)
    sage: c.continued_fraction()
    [0; 1, 26, 1, 2, 2, (8, 4)*]
    sage: d = apply_homography(m1**2*m2**3, a)
    sage: d.continued_fraction()
    [0; 1, 2, 1, 1, 1, 1, 5, 2, 1, 1, 1, 1, 5, 26, 1, 2, 1, 26, 1, 2, 1, 26, 1, 2, 2, (8, 4)*]

It is possible to make arithmetic operations on continued fractions::

    sage: c1 = CFF([0,3,3,2,1,4,2]); c1
    [0; 3, 3, 2, 1, 4, 2]
    sage: c2 = CFF([-4,2,1,3]); c2
    [-4; 2, 1, 3]
    sage: c3 = -c1; c3
    [-1; 1, 2, 3, 2, 1, 4, 2]
    sage: c1+c2
    [-4; 1, 2, 628, 2]
    sage: c1+c3
    [0]

    sage: c1/c2
    [-1; 1, 10, 1, 142]
    sage: ~c3
    [-4; 1, 2, 2, 1, 4, 2]

And they can be used to create matrices, polynomials, vectors, etc::

    sage: a = random_matrix(CFF, 4)
    sage: a
    [      [0; 7]    [0; 4, 3]          [0]       [0; 9]]
    [[0; 1, 4, 4]          [1]          [0]    [0; 7, 5]]
    [      [0; 9] [0; 2, 8, 6]       [0; 9]    [0; 3, 3]]
    [   [0; 4, 5]       [0; 8]          [0]          [1]]
    sage: f = a.charpoly()
    sage: f
    [1]*x^4 + ([-3; 1, 2, 1, 15])*x^3 + ... + [-1; 1, 165, 1, 1, 1, 1, 1, 1, 53, 1, 5]
    sage: f(a)
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    sage: vector(CFF, [1/2, 2/3, 3/4, 4/5])
    ([0; 2], [0; 1, 2], [0; 1, 3], [0; 1, 4])

The unary operations (negative and inversion) are quite fast but binary
operations are quite slow. It is then not adviced, if speed is needed, to use
them as the base class in a computation.

.. TODO::

    - Improve numerical approximation (the method
      :meth:`~ContinuedFraction_generic._mpfr_` is quite slow compared to the
      same method for an element of a number field)

    - Make a class for infinite precision real number built from an infinite
      list (ie an infinite word)

    - Make a class for non standard continued fractions of the form `a_0 +
      b_0/(a_1 + b_1/(...))` (the standard continued fractions are when all
      `b_n= 1` while the Hirzebruch-Jung continued fractions are the one for
      which `b_n = -1` for all `n`). See
      :wikipedia:`Generalized_continued_fraction`.

AUTHORS:

- Niles Johnson (2010-08): ``random_element()`` should pass on ``*args`` and
  ``**kwds`` (:trac:`3893`).

- Vincent Delecroix (2013): cleaning, refactorisation, documentation (:trac:`14567`)
"""
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element, FieldElement

from field import Field
from integer import Integer
from infinity import Infinity
from integer_ring import ZZ
from rational_field import QQ
from real_mpfr import RR
from real_mpfi import RealIntervalField

ZZ_0 = Integer(0)
ZZ_1 = Integer(1)
ZZ_2 = Integer(2)


def two_last_convergents(x):
    """
    Given the list ``x`` that consists of numbers, return the two last
    convergents `p_{n-1}, q_{n-1}, p_n, q_n`.

    This function is principally used to compute the value of a ultimately
    periodic continued fraction.

    EXAMPLES::

        sage: from sage.rings.continued_fractions import two_last_convergents
        sage: two_last_convergents([])
        (0, 1, 1, 0)
        sage: two_last_convergents([0])
        (1, 0, 0, 1)
        sage: two_last_convergents([-1,1,3,2])
        (-1, 4, -2, 9)
    """
    p0, p1 = 0, 1
    q0, q1 = 1, 0
    for a in x:
        p0, p1 = p1, a*p1+p0
        q0, q1 = q1, a*q1+q0
    return p0, q0, p1, q1


class ContinuedFraction_generic(FieldElement):
    r"""
    Generic class for (standard) continued fractions.

    A derived class should implements at least:

    - ``quotient(self, n)``: return the ``n``-th quotient of ``self``

    - ``value(self)``: return the value of ``self`` (a real number)

    This generic class provides:

    - computation of convergents :meth:`convergent`, :meth:`pn` and :meth:`qn`

    - comparison :meth:`__cmp__`

    - elementary arithmetic function :meth:`floor`, :meth:`ceil`, :meth:`sign`

    - numerical approximations :meth:`_mpfr_`

    All other methods rely on :meth:`value` (and not on convergents) and may
    fail at execution.
    """
    def __init__(self, parent):
        r"""
        INPUT:

        - ``parent`` -- the parent of ``self``

        TESTS::

            sage: TestSuite(CFF(3)).run()
        """
        FieldElement.__init__(self, parent)
        self._pn = [ZZ_0, ZZ_1]
        self._qn = [ZZ_1, ZZ_0]

    def __abs__(self):
        """
        Return absolute value of self.

        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: abs(a)
            [0; 22, 1, 7, 2]
            sage: QQ(abs(a))
            17/389
        """
        if self.quotient(0) >= 0:
            return self
        return self.__neg__()

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = CFF([(),(1,)])     # the golden ratio
            sage: d = CFF([(-1,),(1,)])
            sage: d < a and a < b and b < c
            True
            sage: d >= a
            False
        """
        i = 0
        while True:
            a = self.quotient(i)
            b = other.quotient(i)
            test = cmp(a,b)
            if test == 1:  # a > b
                return -1 if i % 2 else 1
            if test == -1:  # b > a
                return 1 if i % 2 else -1
            if a == 0 and b == 0 and i:  # rational case
                return 0
            i += 1

    def _mpfr_(self, R):
        r"""
        Return a numerical approximation of ``self`` in the real mpfr ring ``R``.

        It is guaranteed that the result is exact: when the rounding mode of
        ``R`` is 'RNDN' then the result is the nearest binary number of ``R`` to
        ``self``. The other rounding mode are 'RNDD' (toward +infinity), 'RNDU'
        (toward -infinity) and 'RNDZ' (toward zero).

        EXAMPLES::

            sage: continued_fraction(1/2).n()
            0.500000000000000
            sage: continued_fraction([0,4]).n()
            0.250000000000000
            sage: continued_fraction([12,1,3,4,2,2,3,1,2]).n(digits=4)
            12.76

            sage: continued_fraction(12/7).n(digits=13) == (12/7).n(digits=13)
            True
            sage: continued_fraction(-14/333).n(digits=21) == (-14/333).n(digits=21)
            True

            sage: a = (106*pi - 333) / (355 - 113*pi) - 292
            sage: cf = continued_fraction(a); cf
            [0; 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, 2, 1, ...]
            sage: cf.n(digits=3)
            0.635
            sage: cf.n(digits=4)
            0.6346
            sage: cf.n(digits=5)
            0.63459
            sage: cf.n(digits=6)
            0.634591
            sage: cf.n(digits=7)
            0.6345910
            sage: cf.n(digits=8)
            0.63459101

            sage: K.<a> = NumberField(x^3-2, 'a', embedding=1.25)
            sage: b = 504/253*a^2 + 635/253*a + 661/253
            sage: cf = continued_fraction(b); cf
            [8; 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, 2, 1, 3, 4, 1, 1, 2, 14, 3, ...]
            sage: cf.n(digits=3)
            8.94
            sage: cf.n(digits=6)
            8.93715
            sage: cf.n(digits=7)
            8.937154
            sage: cf.n(digits=8)
            8.9371541
            sage: cf.n(digits=9)
            8.93715414
            sage: cf.n(digits=10)
            8.937154138
            sage: cf.n(digits=11)
            8.9371541378

        TESTS::

        We check that the rounding works as expected, at least in the rational
        case::

            sage: for _ in xrange(100):
            ....:     a = QQ.random_element(num_bound=1<<64)
            ....:     cf = continued_fraction(a)
            ....:     for prec in 17,24,53,128,256:
            ....:         for rnd in 'RNDN','RNDD','RNDU','RNDZ':
            ....:             R = RealField(prec=prec, rnd=rnd)
            ....:             assert R(cf) == R(a)
        """
        # 1. integer case
        if self.quotient(1) == 0:
            return R(self.quotient(0))

        # 2. negative numbers
        # TODO: it is possible to deal with negative values. The only problem is
        # that we need to find the good value for N (which involves
        # self.quotient(k) for k=0,1,2)
        if self.quotient(0) < 0:
            rnd = R.rounding_mode()
            if rnd == 'RNDN' or rnd == 'RNDZ':
                return -R(-self)
            elif rnd == 'RNDD':
                r = R(-self)
                s,m,e = r.sign_mantissa_exponent()
                if e < 0:
                    return -(R(m+1) >> (-e))
                return -(R(m+1) << e)
            else:
                r = R(-self)
                s,m,e = r.sign_mantissa_exponent()
                if e < 0:
                    return -(R(m-1) >> (-e))
                return -(R(m-1) << e)

        # 3. positive non integer
        if self.quotient(0) == 0:  # 0 <= self < 1
            N = R.prec() + self.quotient(1).nbits() - 1
            if self.quotient(2) == 0 and self.quotient(1) % (1 << (self.quotient(1).nbits()-1)) == 0:
                # if self is of the form [0; 2^N] then we need the following
                N -= 1
        else:  # self > 1
            N = R.prec() - self.quotient(0).nbits()

        # even/odd convergents are respectively below/above
        k = 0
        p_even = self.pn(2*k)
        p_odd = self.pn(2*k+1)
        q_even = self.qn(2*k)
        q_odd = self.qn(2*k+1)
        m_even = (p_even << N) // q_even      # floor((2^N p_even) / q_even)
        m_odd = (p_odd << N + q_odd - 1) // q_odd  # ceil((2^N p_odd) / q_odd)
        while (m_odd - m_even) > 1:
            k += 1
            p_even = self.pn(2*k)
            p_odd = self.pn(2*k+1)
            q_even = self.qn(2*k)
            q_odd = self.qn(2*k+1)
            m_even = (p_even << N) // q_even
            m_odd = ((p_odd << N) + q_odd - 1) // q_odd

        assert m_odd.nbits() == R.prec() or m_even.nbits() == R.prec()

        if m_even == m_odd:  # no need to worry (we have a decimal number)
            return R(m_even) >> N

        # check ordering
        # m_even/2^N <= p_even/q_even <= self <= p_odd/q_odd <= m_odd/2^N
        assert m_odd == m_even + 1
        assert m_even / (ZZ_1 << N) <= p_even/q_even
        assert p_even / q_even <= p_odd / q_odd
        assert p_odd / q_odd <= m_odd / (ZZ_1 << N)

        rnd = R.rounding_mode()
        if rnd == 'RNDN':  # round to the nearest
            # in order to find the nearest approximation we possibly need to
            # augment our precision on convergents.
            while True:
                assert not(p_odd << (N+1) <= (2*m_odd-1) * q_odd) or not(p_even << (N+1) >= (2*m_even+1) * q_even)
                if p_odd << (N+1) <= (2*m_odd-1) * q_odd:
                    return R(m_even) >> N
                if p_even << (N+1) >= (2*m_even+1) * q_even:
                    return R(m_odd) >> N
                k += 1
                p_even = self.pn(2*k)
                p_odd = self.pn(2*k+1)
                q_even = self.qn(2*k)
                q_odd = self.qn(2*k+1)
        elif rnd == 'RNDU':  # round up (toward +infinity)
            return R(m_odd) >> N
        elif rnd == 'RNDD':  # round down (toward -infinity)
            return R(m_even) >> N
        elif rnd == 'RNDZ':  # round toward zero
            if m_even.sign() == 1:
                return R(m_even) >> N
            else:
                return R(m_odd) >> N
        else:
            raise ValueError("%s unknown rounding mode" % rnd)

    def __float__(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: float(a)
            -0.043701799485861184
            sage: float(-17/389)
            -0.043701799485861184
        """
        return float(self._mpfr_(RR))

    def pn(self, n):
        """
        Return the numerator of the `n`-th partial convergent of ``self``.

        EXAMPLES::

            sage: c = continued_fraction(pi); c
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: c.pn(0)
            3
            sage: c.pn(12)
            80143857
            sage: c.pn(152)
            3943771611212266962743738812600748213157266596588744951727393497446921245353005283
        """
        n = Integer(n)

        p = self._pn
        q = self._qn

        if n < -2:
            raise ValueError("n must be at least -2")

        for k in xrange(len(p), n+3):
            x = self.quotient(k-2)
            if x == 0 and k != 2:
                return p[-1]
            p.append(x*p[k-1] + p[k-2])
            q.append(x*q[k-1] + q[k-2])

        return p[n+2]

    def qn(self, n):
        """
        Return the denominator of the ``n``-th partial convergent of ``self``.

        EXAMPLES::

            sage: c = continued_fraction(pi); c
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: c.qn(0)
            1
            sage: c.qn(12)
            25510582
            sage: c.qn(152)
            1255341492699841451528811722575401081588363886480089431843026103930863337221076748
        """
        self.pn(n)   # ! silent computation of qn
        if len(self._qn) < n+3:
            return self._qn[-1]
        return self._qn[n+2]

    def convergent(self, n):
        """
        Return the ``n``-th partial convergent to self.

        EXAMPLES::

            sage: a = CFF(pi); a
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: a.convergent(3)
            355/113
            sage: a.convergent(15)
            411557987/131002976
        """
        return self.pn(n) / self.qn(n)

    def __getitem__(self, n):
        r"""
        Return the ``n``-th partial quotient of ``self`` or a continued fraction
        associated to a sublist of the partial quotients of ``self``.

        TESTS::

            sage: cf1 = continued_fraction(pi); cf1
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
            sage: cf2 = continued_fraction(QuadraticField(2).gen()); cf2
            [1; (2)*]
            sage: cf3 = continued_fraction(4/17); cf3
            [0; 4, 4]
            sage: cf1[3:17]
            [1; 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2]
            sage: cf2[:10]
            [1; 2, 2, 2, 2, 2, 2, 2, 2, 2]
            sage: cf3[1:16]
            [4; 4]

        Be careful that the truncation of an infinite continued fraction might
        be shorter by one::

            sage: len(CFF(golden_ratio)[:8])
            7
        """
        if isinstance(n, slice):
            if self.length() == Infinity:
                if n.step is None:
                    step = 1
                else:
                    try:
                        step = n.step.__index__()
                    except (AttributeError,ValueError):
                        raise ValueError("step should be a non zero integer")
                if n.start is None:
                    if step < 0:
                        raise ValueError("if step is negative, start can not be None")
                    start = 0
                else:
                    try:
                        start = n.start.__index__()
                    except (AttributeError,ValueError):
                        raise ValueError("start should be an integer")
                if n.stop is None:
                    if step > 0:
                        raise ValueError("infinite list!")
                    stop = -1
                else:
                    try:
                        stop = n.stop.__index__()
                    except (AttributeError, ValueError):
                        raise ValueError("stop should be an integer")
                return self.parent()([self.quotient(i) for i in xrange(start,stop,step)])
            start, stop, step = n.indices(self.length())
            return self.parent()([self.quotient(i) for i in xrange(start,stop,step)])
        try:
            n = n.__index__()
        except (AttributeError, ValueError):
            raise ValueError("n (=%s) should be an integer" % n)
        if n < 0:
            raise ValueError("n (=%s) should be positive" % n)
        q = self.quotient(n)
        if n > 0 and q == 0:
            raise IndexError("index out of range")

    def __iter__(self):
        r"""
        Iterate over the partial quotient of self.

        EXAMPLES::

            sage: cf = continued_fraction(pi)
            sage: i = iter(cf)
            sage: [i.next() for _ in xrange(10)]
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1]
            sage: [i.next() for _ in xrange(10)]
            [3, 1, 14, 2, 1, 1, 2, 2, 2, 2]
            sage: [i.next() for _ in xrange(10)]
            [1, 84, 2, 1, 1, 15, 3, 13, 1, 4]
        """
        yield self.quotient(0)
        i = 1
        while True:
            q = self.quotient(i)
            if not q:
                break
            yield q
            i += 1

    def __int__(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: int(a)
            -1
        """
        return int(self.quotient(0))

    def __long__(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: long(a)
            -1L
        """
        return long(self.quotient(0))

    def sign(self):
        r"""
        Returns the sign of self as an Integer.

        The sign is defined to be ``0`` if ``self`` is ``0``, ``1`` if ``self``
        is positive and ``-1`` if ``self`` is negative.

        EXAMPLES::

            sage: continued_fraction(tan(pi/7)).sign()
            1
            sage: continued_fraction(-34/2115).sign()
            -1
            sage: continued_fraction([0]).sign()
            0
        """
        if self.quotient(0) == 0:
            if self.quotient(1) == 0:
                return ZZ_0
            return ZZ_1
        return self.quotient(0).sign()

    def floor(self):
        r"""
        Return the floor of ``self``.

        EXAMPLES::

            sage: cf = CFF([2,1,2,3])
            sage: cf.floor()
            2
        """
        return self.quotient(0)

    def ceil(self):
        r"""
        Return the ceil of ``self``.

        EXAMPLES::

            sage: cf = CFF([2,1,3,4])
            sage: cf.ceil()
            3
        """
        if self.length() == 1:
            return self.quotient(0)
        return self.quotient(0)+1

    def __neg__(self):
        r"""
        Return the opposite of ``self``.

        EXAMPLES::

            sage: -CFF(pi) == CFF(-pi)
            True
        """
        return self.parent()(-self.value())

    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: ~CFF(pi) == CFF(~pi)
            True
        """
        return self.parent()(~self.value())

    def _add_(self, other):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a+b; c
            [-1; 1, 23, 3, 5]
            sage: c.value()
            -16/389

        Note:

        The algorithm is very naive !
        """
        try:
            return self.parent()(self.value() + other.value())
        except (ValueError,TypeError):
            raise NotImplementedError

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a - b; c
            [-1; 1, 20, 1, 1, 1, 1, 3]
            sage: c.value()
            -18/389
        """
        try:
            return self.parent()(self.value() - other.value())
        except (ValueError,TypeError):
            raise NotImplementedError

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a * b; c
            [-1; 1, 8900, 4, 4]
            sage: c.value(), (-1/389)*(17/389)
            (-17/151321, -17/151321)
        """
        try:
            return self.parent()(self.value() * other.value())
        except (ValueError,TypeError):
            raise NotImplementedError

    def _div_(self, other):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a / b; c
            [-17]
            sage: c.value(), (17/389) / (-1/389)
            (-17, -17)
        """
        try:
            return self.parent()(self.value() / other.value())
        except (ValueError,TypeError):
            raise NotImplementedError

    def __nonzero__(self):
        """
        Return False if self is zero.

        EXAMPLES::

            sage: CFF(0).is_zero()    # indirect doctest
            True
            sage: CFF(1).is_zero()    # indirect doctest
            False
            sage: CFF([(),(1,)]).is_zero()     # indirect doctest
            False
            sage: CFF([(0,),(1,2)]).is_zero()  # indirect doctest
            False
        """
        return bool(self.quotient(0)) or bool(self.quotient(1))

    def is_zero(self):
        r"""
        Test whether ``self`` is zero.

        EXAMPLES::

            sage: CFF(0).is_zero()
            True
            sage: CFF((0,1)).is_zero()
            False
            sage: CFF(-1/2).is_zero()
            False
            sage: CFF(pi).is_zero()
            False
        """
        return self.quotient(0) == ZZ_0 and self.quotient(1) == ZZ_0

    def is_one(self):
        r"""
        Test whether ``self`` is one.

        EXAMPLES::

            sage: CFF(1).is_one()
            True
            sage: CFF(5/4).is_one()
            False
            sage: CFF(0).is_one()
            False
            sage: CFF(pi).is_one()
            False
        """
        return self.quotient(0) == ZZ_1 and self.quotient(1) == ZZ_0

    def additive_order(self):
        """
        Return the additive order of this continued fraction,
        which we defined to be the additive order of its value.

        EXAMPLES::

            sage: CFF(-1).additive_order()
            +Infinity
            sage: CFF(0).additive_order()
            1
        """
        return Infinity if self else ZZ_1

    def multiplicative_order(self):
        """
        Return the multiplicative order of this continued fraction,
        which we defined to be the multiplicative order of its value.

        EXAMPLES::

            sage: CFF(-1).multiplicative_order()
            2
            sage: CFF(1).multiplicative_order()
            1
            sage: CFF(pi).multiplicative_order()
            +Infinity
        """
        if self.is_zero():
            return Infinity
        if self.is_one():
            return ZZ_1
        if (-self).is_one():
            return ZZ_2
        return Infinity


class ContinuedFraction_periodic(ContinuedFraction_generic):
    r"""
    Continued fraction associated with rational or quadratic number.

    A rational number has a finite continued fraction expansion (or ultimately
    0). The one of a quadratic number, ie a number of the form `a + b \sqrt{D}`
    with `a` and `b` rational, is ultimately periodic.

    .. NOTE::

        This class stores a tuple ``_x1`` for the preperiod and a tuple ``_x2``
        for the period. In the purely periodic case ``_x1`` is empty while in
        the rational case ``_x2`` is the tuple ``(0,)``.
    """
    def __init__(self, parent, x1, x2, check=True):
        r"""
        TESTS::

            sage: TestSuite(CFF([(1,2),(3,4,5)])).run()
        """
        ContinuedFraction_generic.__init__(self, parent)
        self._x1 = tuple(x1)
        self._x2 = tuple(x2)

    def period(self):
        r"""
        Return the periodic part of ``self``.

        EXAMPLES::

            sage: K.<sqrt3> = QuadraticField(3)
            sage: cf = continued_fraction(sqrt3); cf
            [1; (1, 2)*]
            sage: cf.period()
            (1, 2)

            sage: for k in xsrange(2,40):
            ....:     if not k.is_square():
            ....:         s = QuadraticField(k).gen()
            ....:         cf = continued_fraction(s)
            ....:         print '%2d %d %s'%(k, len(cf.period()), cf)
             2 1 [1; (2)*]
             3 2 [1; (1, 2)*]
             5 1 [2; (4)*]
             6 2 [2; (2, 4)*]
             7 4 [2; (1, 1, 1, 4)*]
             8 2 [2; (1, 4)*]
            10 1 [3; (6)*]
            11 2 [3; (3, 6)*]
            12 2 [3; (2, 6)*]
            13 5 [3; (1, 1, 1, 1, 6)*]
            14 4 [3; (1, 2, 1, 6)*]
            ...
            35 2 [5; (1, 10)*]
            37 1 [6; (12)*]
            38 2 [6; (6, 12)*]
            39 2 [6; (4, 12)*]
        """
        if self._x2[0] == ZZ_0:
            return ()
        return self._x2

    def preperiod(self):
        r"""
        Return the preperiodic part of ``self``.

        EXAMPLES::

            sage: K.<sqrt3> = QuadraticField(3)
            sage: cf = continued_fraction(sqrt3); cf
            [1; (1, 2)*]
            sage: cf.preperiod()
            (1,)

            sage: cf = continued_fraction(sqrt3/7); cf
            [0; 4, (24, 8)*]
            sage: cf.preperiod()
            (0, 4)
        """
        return self._x1

    def quotient(self, n):
        r"""
        Return the ``n``-th quotient of ``self``.

        EXAMPLES::

            sage: cf = CFF([(12,5),(1,3)])
            sage: [cf.quotient(i) for i in xrange(10)]
            [12, 5, 1, 3, 1, 3, 1, 3, 1, 3]
        """
        n = int(n)
        if n < 0:
            raise ValueError("n (=%d) should be positive" % n)
        if n < len(self._x1):
            return self._x1[n]
        return self._x2[(n-len(self._x1)) % len(self._x2)]

    def length(self):
        r"""
        Returns the number of partial quotients of ``self``.

        EXAMPLES::

            sage: CFF(2/5).length()
            3
            sage: cf = CFF([(0,1),(2,)]); cf
            [0; 1, (2)*]
            sage: cf.length()
            +Infinity
        """
        if len(self._x2) > 1 or self._x2[0] != ZZ_0:
            return Infinity
        return ZZ(len(self._x1))

    def __cmp__(self, other):
        r"""
        EXAMPLES::

            sage: a = CFF([(0,),(1,2,3,1,2,3,1)]); a.n()
            0.694249167819459
            sage: b = CFF([(0,),(1,2,3)]); b.n()
            0.694254176766073
            sage: c = CFF([(0,1),(2,3)]); c.n()
            0.696140478029631
            sage: d = CFF([(0,1,2),(3,)]); d.n()
            0.697224362268005
            sage: a < b and a < c and a < d
            True
            sage: b < c and b < d and c < d
            True
            sage: b == c
            False
            sage: c > c
            False
            sage: b >= d
            False
        """
        if isinstance(other, ContinuedFraction_periodic):
            n = max(len(self._x1) + 2*len(self._x2),
                    len(other._x1) + 2*len(other._x2))
            for i in xrange(n):
                a = self.quotient(i)
                b = other.quotient(i)
                test = cmp(a,b)
                if test == 1:
                    return -1 if i % 2 else 1
                if test == -1:
                    return 1 if i % 2 else -1
            return 0

        return ContinuedFraction_generic.__cmp__(self, other)

    def value(self):
        r"""
        Return the value of ``self`` as a quadratic number (with square free
        discriminant).

        EXAMPLES:

            Some purely periodic examples::

            sage: cf = CFF([(),(2,)]); cf
            [(2)*]
            sage: v = cf.value(); v
            sqrt2 + 1
            sage: v.continued_fraction()
            [(2)*]

            sage: cf = CFF([(),(1,2)]); cf
            [(1, 2)*]
            sage: v = cf.value(); v
            1/2*sqrt3 + 1/2
            sage: v.continued_fraction()
            [(1, 2)*]

            sage: cf = CFF([(),(3,2,1)])
            sage: cf.value().continued_fraction() == cf
            True
            sage: cf = CFF([(),(1,3,1,5)])
            sage: cf.value().continued_fraction() == cf
            True

            Some ultimately periodic but non periodic examples::

            sage: cf = CFF([(1,),(2,)]); cf
            [1; (2)*]
            sage: v = cf.value(); v
            sqrt2
            sage: v.continued_fraction()
            [1; (2)*]

            sage: cf = CFF([(1,3),(1,2)]); cf
            [1; 3, (1, 2)*]
            sage: v = cf.value(); v
            -sqrt3 + 3
            sage: v.continued_fraction()
            [1; 3, (1, 2)*]

            sage: cf = CFF([(-5,18), (1,3,1,5)])
            sage: cf.value().continued_fraction() == cf
            True
            sage: cf = CFF([(-1,),(1,)])
            sage: cf.value().continued_fraction() == cf
            True

        TESTS::

            sage: a1 = ((0,1),(2,3))
            sage: a2 = ((-12,1,1),(2,3,2,4))
            sage: a3 = ((1,),(1,2))
            sage: a4 = ((-2,2),(1,124,13))
            sage: a5 = ((0,),(1,))
            sage: for a in a1,a2,a3,a4,a5:
            ....:     cf = CFF(a)
            ....:     assert cf.value().continued_fraction() == cf
        """
        if self._x1 and self._x1[0] < 0:
            return -(-self).value()

        if len(self._x2) == 1 and not self._x2[0]:
            return self._rational_()

        # determine the equation for the purely periodic cont. frac. determined
        # by self._x2
        p0,q0,p1,q1 = two_last_convergents(self._x2)

        # now x is one of the root of the equation
        #   q1 x^2 + (q0 - p1) x - p0 = 0
        from sage.rings.number_field.number_field import QuadraticField
        from sage.misc.functional import squarefree_part
        D = (q0-p1)**2 + 4*q1*p0
        DD = squarefree_part(D)
        Q = QuadraticField(DD, 'sqrt%d' % DD)
        x = ((p1 - q0) + (D/DD).sqrt() * Q.gen()) / (2*q1)

        # we add the preperiod
        p0,q0,p1,q1 = two_last_convergents(self._x1)
        return (p1*x + p0) / (q1*x + q0)

    def _repr_(self):
        r"""
        TESTS::

            sage: a = continued_fraction(pi.n()); a
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: a.rename('continued fraction of pi')
            sage: a
            continued fraction of pi

            sage: CFF([(0,1),(2,)])
            [0; 1, (2)*]
            sage: CFF([(),(1,3)])
            [(1, 3)*]
        """
        if len(self._x2) == 1 and not self._x2[0]:  # rational
            if len(self._x1) == 1:
                return '[%d]' % self._x1[0]
            return '[%d; ' % self._x1[0] + ', '.join(str(a) for a in self._x1[1:]) + ']'

        period = '(' + ', '.join(str(a) for a in self._x2) + ')*'
        if not self._x1:  # purely periodic case
            return '[' + period + ']'

        if len(self._x1) == 1:
            return '[%d; ' % self._x1[0] + period + ']'
        return '[%d; ' % self._x1[0] + ', '.join(str(a) for a in self._x1[1:]) + ', ' + period + ']'

    def __reduce__(self):
        r"""
        Pickling support.

        EXAMPLES::

            sage: a = CFF([1,2,3])
            sage: loads(dumps(a)) == a
            True

            sage: a = CFF([(-1,2),(3,1,4)])
            sage: loads(dumps(a)) == a
            True
        """
        return self.parent(),((self._x1,self._x2),)

    def convergents(self):
        """
        Return the list of partial convergents of ``self``.

        EXAMPLES::

            sage: a = CFF(23/157); a
            [0; 6, 1, 4, 1, 3]
            sage: a.convergents()
            [0, 1/6, 1/7, 5/34, 6/41, 23/157]
        """
        return [self.pn(n) / self.qn(n) for n in xrange(len(self))]

    def __len__(self):
        """
        Return the number of terms in this continued fraction.

        EXAMPLES::

            sage: len(continued_fraction([1,2,3,4,5]) )
            5
        """
        if len(self._x2) == 1 and self._x2[0] == ZZ_0:
            return len(self._x1)
        raise ValueError("the length is infinite")

    def _rational_(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: a._rational_()
            -17/389
            sage: QQ(a)
            -17/389
        """
        if len(self._x2) > 1 or self._x2[0] != 0:
            raise ValueError("this is not a rational!")
        n = len(self)
        return self.pn(n-1) / self.qn(n-1)

    def _latex_(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: latex(a)
            -1+ \frac{\displaystyle 1}{\displaystyle 1+ \frac{\displaystyle 1}{\displaystyle 21+ \frac{\displaystyle 1}{\displaystyle 1+ \frac{\displaystyle 1}{\displaystyle 7+ \frac{\displaystyle 1}{\displaystyle 2}}}}}
        """
        if self._x2 != (0,):
            raise NotImplementedError("latex not implemented for non rational continued fractions")
        v = self._x1
        if len(v) == 0:
            return '0'
        s = str(v[0])
        for i in range(1,len(v)):
            s += '+ \\frac{\\displaystyle 1}{\\displaystyle %s' % v[i]
        s += '}'*(len(v)-1)
        return s

    def __invert__(self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES::

            sage: a = CFF(13/25)
            sage: ~a == CFF(25/13)
            True
            sage: a * ~a
            [1]

            sage: a = CFF(-17/253)
            sage: ~a == CFF(-253/17)
            True
            sage: a * ~a
            [1]

            sage: K.<sqrt5> = QuadraticField(5)
            sage: a1 = (sqrt5+1)/2
            sage: c1 = a1.continued_fraction(); c1
            [(1)*]
            sage: ~c1
            [0; (1)*]
            sage: ~c1 * c1
            [1]
            sage: c2 = (sqrt5/3 + 1/7).continued_fraction(); c2
            [0; 1, (7, 1, 17, ..., 1, 2)*]
            sage: ~c2
            [1; (7, 1, 17, ..., 1, 2)*]
            sage: c2 * ~c2
            [1]
        """
        if not self:
            raise ZeroDivisionError("Rational division by 0")
        if self._x1:
            if self._x1[0] < 0:
                return -(-self).__invert__()
            if self._x1[0] == 0:
                return self.parent()([self._x1[1:],self._x2])
        return self.parent()([(0,) + self._x1, self._x2])

    def __neg__(self):
        """
        Return additive inverse of self.

        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1; 1, 21, 1, 7, 2]
            sage: -a
            [0; 22, 1, 7, 2]
            sage: -(-a)
            [-1; 1, 21, 1, 7, 2]
            sage: QQ(-a)
            17/389

            sage: -CFF([2])
            [-2]
            sage: -CFF([1,2])
            [-2; 2]

            sage: a = CFF([(1,2),(1,)]); a
            [1; 2, (1)*]
            sage: -a
            [-2; (1)*]
            sage: -(-a)
            [1; 2, (1)*]

            sage: a = CFF([(),(1,2,3)]); a
            [(1, 2, 3)*]
            sage: -a
            [-2; 1, 1, (3, 1, 2)*]
            sage: -(-a)
            [(1, 2, 3)*]

            sage: -CFF([0])
            [0]
            sage: -CFF([1])
            [-1]
            sage: -CFF([-1])
            [1]
        """
        x1 = self._x1

        # integer case
        if len(self._x2) == 1 and self._x2[0] == 0 and len(self._x1) == 1:
            return self.parent()(-self._x1[0])

        # to make all other cases work, we need 3 elements in x1
        if len(x1) < 3:
            x1 = self._x1 + self._x2
        if len(x1) < 3:
            x1 += self._x2

        # we are yet sure that there are at least two elements...
        if x1[0] == 0 and x1[1] == 0:
            return self.parent().zero()

        if len(x1) < 3:
            x1 += self._x2

        if x1[1] == 1:
            return self.parent()([(-x1[0]-1, x1[2]+1) + x1[3:], self._x2])
        return self.parent()([(-x1[0]-1, ZZ_1, x1[1]-1) + x1[2:], self._x2])


class ContinuedFraction_irrational(ContinuedFraction_generic):
    r"""
    Continued fraction of an exact irrational number.

    EXAMPLES::

        sage: continued_fraction(pi)
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        sage: continued_fraction(e)
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, ...]
    """
    def __init__(self, parent, x):
        r"""
        INPUT:

        - ``parent`` -- the parent of that continued fraction

        - ``x`` -- the real number from which we want the continued fraction

        TESTS::

            sage: TestSuite(CFF(pi)).run()
        """
        ContinuedFraction_generic.__init__(self, parent)
        self._x0 = x

        self._xa = RealIntervalField(53)(self._x0)   # an approximation of the
                                                     # last element of the orbit
                                                     # under the Gauss map
        self._quotients = []

    def length(self):
        r"""
        Return infinity

        EXAMPLES::

            sage: CFF(pi).length()
            +Infinity
        """
        return Infinity

    def __len__(self):
        r"""
        TESTS::

            sage: len(CFF(pi))
            Traceback (most recent call last):
            ...
            ValueError: the length is infinite!
        """
        raise ValueError("the length is infinite!")

    def __reduce__(self):
        r"""
        Pickling support.

        TESTS::

            sage: cf = continued_fraction(pi)
            sage: loads(dumps(cf)) == cf
            True
        """
        return self.parent(),(self.value(),)

    def __cmp__(self, other):
        r"""
        Comparison.

        EXAMPLES::

            sage: CFF(pi) > CFF(e)
            True
            sage: CFF(pi) > CFF(e+4)
            False
        """
        try:
            # The following is crazy and prevent us from using cmp(self.value(),
            # other.value()). On sage-5.10.beta2:
            #     sage: cmp(pi, 4)
            #     -1
            #     sage: cmp(pi, pi+4)
            #     1
            if self.value() == other.value():
                return 0
            if self.value() - other.value() > 0:
                return 1
            return -1
        except StandardError:
            return ContinuedFraction_generic.__cmp__(self, other)

    def _repr_(self):
        r"""
        String representation.

        EXAMPLES::

            sage: continued_fraction(pi) # indirect doctest
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        """
        return '[%d; ' % self.quotient(0) + ', '.join(str(self.quotient(i)) for i in xrange(1,20)) + ", ...]"

    def quotient(self, n):
        r"""
        Returns the ``n``-th quotient of ``self``.

        EXAMPLES::

            sage: cf = continued_fraction(pi)
            sage: cf.quotient(27)
            13
            sage: cf.quotient(2552)
            152
            sage: cf.quotient(10000)   # long time
            5

        The algorithm is not efficient with element of the symbolic ring and,
        if possible, one can always prefer number fields elements. The reason is
        that, given a symbolic element ``x``, there is no automatic way to
        evaluate in ``RIF`` an expression of the form ``(a*x+b)/(c*x+d)`` where
        both the numerator and the denominator are extremely small::

            sage: a1 = pi
            sage: c1 = continued_fraction(a1)
            sage: p0 = c1.pn(12); q0 = c1.qn(12)
            sage: p1 = c1.pn(13); q1 = c1.qn(13)
            sage: num = (q0*a1 - p0); num.n()
            1.49011611938477e-8
            sage: den = (q1*a1 - p1); den.n()
            -2.98023223876953e-8
            sage: a1 = -num/den
            sage: RIF(a1)
            [-infinity .. +infinity]

        The same computation with an element of a number field instead of
        ``pi`` gives a very satisfactory answer::

            sage: K.<a2> = NumberField(x^3 - 2, embedding=1.25)
            sage: c2 = continued_fraction(a2)
            sage: p0 = c2.pn(111); q0 = c2.qn(111)
            sage: p1 = c2.pn(112); q1 = c2.qn(112)
            sage: num = (q0*a2 - p0); num.n()
            -4.56719261665907e46
            sage: den = (q1*a2 - p1); den.n()
            -3.65375409332726e47
            sage: a2 = -num/den
            sage: b2 = RIF(a2); b2
            1.002685823312715?
            sage: b2.absolute_diameter()
            8.88178419700125e-16

        The consequence is that the precision needed with ``c1`` grows when we
        compute larger and larger partial quotients::

            sage: c1.quotient(100)
            2
            sage: c1._xa.parent()
            Real Interval Field with 353 bits of precision
            sage: c1.quotient(200)
            3
            sage: c1._xa.parent()
            Real Interval Field with 753 bits of precision
            sage: c1.quotient(300)
            5
            sage: c1._xa.parent()
            Real Interval Field with 1053 bits of precision

            sage: c2.quotient(200)
            6
            sage: c2._xa.parent()
            Real Interval Field with 53 bits of precision
            sage: c2.quotient(500)
            1
            sage: c2._xa.parent()
            Real Interval Field with 53 bits of precision
            sage: c2.quotient(1000)
            1
            sage: c2._xa.parent()
            Real Interval Field with 53 bits of precision
        """
        x = self._xa
        for k in xrange(len(self._quotients), n+1):
            if x.lower().is_infinity() or x.upper().is_infinity() or x.lower().floor() != x.upper().floor():
                orbit = lambda z: -(self.qn(k-2)*z-self.pn(k-2))/(self.qn(k-1)*z-self.pn(k-1))
                x = x.parent()(orbit(self._x0))

                # It may happen that the above line fails to give an
                # approximation with the expected number of digits (see in the
                # examples). In that case, we augment the precision.
                while x.lower().is_infinity() or x.upper().is_infinity() or x.lower().floor() != x.upper().floor():
                    self._prec = x.parent().prec() + 100
                    x = RealIntervalField(self._prec)(orbit(self._x0))

            self._quotients.append(x.unique_floor())
            x = (x-x.unique_floor())
            if not x:
                raise ValueError("the number is rational")
            x = ~x

        self._xa = x
        return self._quotients[n]

    def value(self):
        r"""
        Return the value of ``self`` (the number from which it was built).

        EXAMPLES::

            sage: cf = continued_fraction(e)
            sage: cf.value()
            e
        """
        return self._x0


def check_and_reduce_pair(x1,x2=None):
    r"""
    There are often two ways to represent a given continued fraction. This
    function makes it canonical.

    In the very special case of the number `0` we return the pair
    ``((0,),(0,))``.

    TESTS::

        sage: from sage.rings.continued_fractions import check_and_reduce_pair
        sage: check_and_reduce_pair([])
        ((0,), (0,))
        sage: check_and_reduce_pair([-1,1])
        ((0,), (0,))
        sage: check_and_reduce_pair([1,1,1])
        ((1, 2), (0,))
        sage: check_and_reduce_pair([1,3],[2,3])
        ((1,), (3, 2))
        sage: check_and_reduce_pair([1,2,3],[2,3,2,3,2,3])
        ((1,), (2, 3))
        sage: check_and_reduce_pair([0,0],[0,0,0])
        ((0,), (0,))
        sage: check_and_reduce_pair([1,2,0],[0])
        ((1, 2), (0,))
    """
    y1 = map(ZZ,x1)

    if x2 is None or (not x2) or all(a == 0 for a in x2):
        y2 = [0]
        if not y1:
            y1 = [ZZ_0]
        elif len(y1) > 1 and y1[-1] == 1:
            y1.pop(-1)
            y1[-1] += 1

    else:
        y2 = map(ZZ,x2)
        if any(b <= ZZ_0 for b in y2):
            raise ValueError("the elements of the period can not be negative")

    # add possibly some element of x1 into the period
    while y1 and y1[-1] == y2[-1]:
        y1.pop(-1)
        y2.insert(0,y2.pop(-1))

    # some special cases to treat
    if len(y2) == 1 and y2[0] == 0:
        if not y1:
            y1 = [ZZ_0]
        elif len(y1) > 1 and y1[-1] == 1:
            y1.pop(-1)
            y1[-1] += 1

    # check that y2 is not a pure power (in a very naive way!!)
    n2 = len(y2)
    for i in xrange(1,(n2+2)/2):
        if n2 % i == 0 and y2[:-i] == y2[i:]:
            y2 = y2[:i]
            break

    # check that at then end y1 has no zeros in it
    for i in xrange(1,len(y1)):
        if y1[i] <= 0:
            raise ValueError("all quotient except the first must be positive")

    return tuple(y1),tuple(y2)


class ContinuedFractionField(UniqueRepresentation, Field):
    """
    A common parent for continued fractions.

    The continued fraction is a field isomorphic to `\RR`. Nevertheless, it is
    not a good idea to use continued fraction to make calculus: operations are
    very slow.

    .. SEEALSO::

        :func:`continued_fraction`

    EXAMPLES::

        sage: CFF
        Real field with infinite precision (continued fractions)
        sage: CFF([0,1,3,2])
        [0; 1, 3, 2]
        sage: CFF(133/25)
        [5; 3, 8]
        sage: CFF(cos(pi/12))
        [0; 1, 28, 2, 1, 7, 21, 1, 8, 1, 3, 10, 1, 2, 3, 1, 8, 2, 1, 2, ...]

        sage: CFF.category()
        Category of fields
    """
    def __init__(self):
        r"""
        TESTS::

            sage: TestSuite(CFF).run()
        """
        Field.__init__(self, self)

        # build the main element classes
        self.element_class = self._element_class_periodic = self.__make_element_class__(ContinuedFraction_periodic)
        self._element_class_number = self.__make_element_class__(ContinuedFraction_irrational)

        self._populate_coercion_lists_(convert_method_name='continued_fraction')

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(CFF)
            \Bold{CFF}
        """
        return "\\Bold{CFF}"

    def _repr_(self):
        """
        EXAMPLES::

            sage: CFF
            Real field with infinite precision (continued fractions)
        """
        return "Real field with infinite precision (continued fractions)"

    def an_element(self):
        r"""
        Returns a continued fraction.

        EXAMPLES::

            sage: CFF.an_element()
            [-1; 2, 3]
        """
        return self([-1,2,3])

    def _element_constructor_(self, x, bits=None, nterms=None, check=True, value=None):
        """
        Construct a continued fraction from the given data.

        INPUT:

            - `x` -- a number

            - ``bits`` -- integer (optional) the number of bits of the
              input number to use when computing with continued fraction.

            - ``nterms`` -- integer (optional)

            - ``check`` -- boolean (optional) -- whether or not we check the

        EXAMPLES::

            sage: CFF(1.5)
            [1; 2]
            sage: CFF(23/17)
            [1; 2, 1, 5]

            sage: continued_fraction(e)
            [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, ...]
            sage: continued_fraction(pi)
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2,
            ...]

            sage: CFF([1,2,3])
            [1; 2, 3]
            sage: CFF(15/17)
            [0; 1, 7, 2]

            sage: CFF(pi, bits=20)
            [3; 7]
            sage: CFF(pi, bits=80)
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 3]
            sage: CFF(pi, bits=100)
            [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, 2, 1, 1, 15, 3]
            sage: CFF(pi, nterms=3)
            [3; 7, 15]
            sage: CFF(pi, nterms=10)
            [3; 7, 15, 1, 292, 1, 1, 1, 3]
        """
        if isinstance(x, Element) and x.parent() is self:
            return x

        if isinstance(x, (list, tuple)):  # a digit expansion
            if len(x) == 2 and isinstance(x[0], (list,tuple)) and isinstance(x[1], (list,tuple)):
                x1 = tuple(ZZ(a) for a in x[0])
                x2 = tuple(ZZ(a) for a in x[1])
                x1, x2 = check_and_reduce_pair(x1, x2)
            else:
                x1, x2 = check_and_reduce_pair(x)
            return self._element_class_periodic(self, x1, x2)

        else:  # a number
            if x in QQ:  # rational
                return QQ(x).continued_fraction()
            if bits is None and nterms is None and isinstance(x, Element):  # an irrational
                from sage.symbolic.ring import SR
                if x.parent() is SR:
                    # we try a conversion to the algebraic real field (in order
                    # to make more accurate computations and detect purely
                    # periodic expansion).
                    from sage.rings.qqbar import AA
                    try:
                        x = AA(x)
                    except StandardError:
                        pass
                    else:
                        K,a,_ = x.as_number_field_element(minimal=True)
                        # actually, the method as_number_field_element is not
                        # well implemented as it is not initialized with a
                        # coerce embedding...
                        try:
                            return a.continued_fraction()
                        except AttributeError:
                            pass

                if x.parent() is SR or x.parent().is_exact():
                    # TODO: even if ``x`` belongs to the symbolic ring, we don't
                    # know if it is an exact real number (like ``pi`` or ``e``)
                    # or if it is a floating point approximation...
                    from sage.rings.real_mpfi import RIF
                    try:
                        RIF(x)
                    except StandardError:
                        raise ValueError("the number %s does not seem to be a real number" % x)
                    return self._element_class_number(self, x)

            # now work with RIF
            if bits is None:
                try:
                    bits = x.parent().prec()
                except AttributeError:
                    bits = 53
            x = RealIntervalField(bits)(x)
            cf = []
            while True:
                try:
                    a = x.unique_floor()
                except ValueError:
                    x1, x2 = check_and_reduce_pair(cf)
                    cf = self._element_class_periodic(self, x1, x2)
                    break
                cf.append(a)
                x = ~(x - a)

            if nterms is None:
                return cf
            return self([cf.quotient(i) for i in xrange(nterms)])

    def is_field(self, proof=True):
        """
        Return True.

        EXAMPLES::

            sage: CFF.is_field()
            True
        """
        return True

    def is_exact(self):
        r"""
        Return True.

        EXAMPLES::

            sage: CFF.is_exact()
            True
        """
        return True

    def is_finite(self):
        """
        Return False, since the continued fraction field is not finite.

        EXAMPLES::

            sage: CFF.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Return 0, since the continued fraction field has characteristic 0.

        EXAMPLES::

            sage: c = CFF.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return ZZ_0

    def order(self):
        """
        EXAMPLES::

            sage: CFF.order()
            +Infinity
        """
        return Infinity

    def random_element(self, first_digit=0, preperiod=True, period=False,
                       x=1, y=10, distribution=None):
        """
        Return a somewhat random continued fraction (the result is either
        finite or ultimately periodic).

        INPUT:

        - ``first_digit`` -- an optional value for the first digit (default
          ``0``)

        - ``preperiod`` -- boolean -- whether or not have a preperiod (default
          is ``True``)

        - ``period`` -- boolean -- whether or not have a period (default is
          ``False``)

        - ``x``, ``y`` -- optional bounds for the digits (default are ``1`` and
          ``10``). Note that these options are sent to ``ZZ.random_element``.

        - ``distribution`` -- an optional distribution for the digits

        EXAMPLES::

            sage: CFF.random_element()
            [0; 4, 7]
            sage: CFF.random_element()
            [0; 7, 6]

            sage: CFF.random_element(preperiod=False, period=True)
            [0; (6)*]
            sage: CFF.random_element(preperiod=False, period=True)
            [0; (1, 3, 4, 5, 8)*]
        """
        from sage.misc.prandom import random
        x1 = [ZZ(first_digit)]
        x2 = []
        if preperiod:
            if random() > 0.1:
                x1.append(ZZ.random_element(x,y,distribution))
                while random() > 0.5:  # from here they may be arbitrarily many
                                       # elements but the mean length is 2
                    x1.append(ZZ.random_element(x,y,distribution))
        if period:
            x2.append(ZZ.random_element(x,y,distribution))
            while random() > 0.2:
                x2.append(ZZ.random_element(x,y,distribution))

        x1,x2 = check_and_reduce_pair(x1,x2)
        return self._element_class_periodic(self, x1, x2)

CFF = ContinuedFractionField()


def continued_fraction_list(x, type="std", partial_convergents=False, bits=None, nterms=None):
    r"""
    Returns the (finite) continued fraction of ``x`` as a list.

    The continued fraction expansion of ``x`` are the coefficients `a_i` in

    .. MATH::

        x = a_0 + 1/(a_1 + 1/(...))

    with `a_0` integer and `a_1`, `...` positive integers. The Hirzebruch-Jung
    continued fraction is the one for which the `+` signs are replaced with `-`
    signs

    .. MATH::

        x = a_0 - 1/(a_1 - 1/(...))

    .. SEEALSO::

        :func:`continued_fraction`

    INPUT:

    - ``x`` -- exact rational or floating-point number. The number to
      compute the continued fraction of.

    - ``type`` -- either "std" (default) for standard continued fractions or
      "hj" for Hirzebruch-Jung ones.

    - ``partial_convergents`` -- boolean. Whether to return the partial convergents.

    - ``bits`` -- integer. the precision of the real interval field
      that is used internally.

    - ``nterms`` -- integer. The upper bound on the number of terms in
      the continued fraction expansion to return.

    OUTPUT:

    A lits of integers, the coefficients in the continued fraction expansion of
    ``x``. If ``partial_convergents`` is set to ``True``, then return a pair
    containing the coefficient list and the partial convergents list is
    returned.

     EXAMPLES::

        sage: continued_fraction_list(45/19)
        [2, 2, 1, 2, 2]
        sage: 2 + 1/(2 + 1/(1 + 1/(2 + 1/2)))
        45/19

        sage: continued_fraction_list(45/19,type="hj")
        [3, 2, 3, 2, 3]
        sage: 3 - 1/(2 - 1/(3 - 1/(2 - 1/3)))
        45/19

    Specifying ``bits`` or ``nterms`` modify the length of the output::

        sage: continued_fraction_list(e, bits=20)
        [2, 1, 2, 1, 1, 4, 2]
        sage: continued_fraction_list(sqrt(2)+sqrt(3), bits=30)
        [3, 6, 1, 5, 7, 2]
        sage: continued_fraction_list(pi, bits=53)
        [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]

        sage: continued_fraction_list(log(3/2), nterms=15)
        [0, 2, 2, 6, 1, 11, 2, 1, 2, 2, 1, 4, 3, 1, 1]
        sage: continued_fraction_list(tan(sqrt(pi)), nterms=20)
        [-5, 9, 4, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1, 1, 1, 2, 4, 3, 1, 63]

    When the continued fraction is infinite (ie ``x`` is an irrational number)
    and the parameters ``bits`` and ``nterms`` are not specified then a warning
    is raised::

        sage: continued_fraction_list(sqrt(2))
        doctest:...: UserWarning: the continued fraction of sqrt(2) seems infinite, return only the first 20 terms
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: continued_fraction_list(sqrt(4/19))
        doctest:...: UserWarning: the continued fraction of 2*sqrt(1/19) seems infinite, return only the first 20 terms
        [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 16]

    An examples with the list of partial convergents::

        sage: continued_fraction_list(RR(pi), partial_convergents=True)
        ([3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3],
         [(3, 1),
          (22, 7),
          (333, 106),
          (355, 113),
          (103993, 33102),
          (104348, 33215),
          (208341, 66317),
          (312689, 99532),
          (833719, 265381),
          (1146408, 364913),
          (4272943, 1360120),
          (5419351, 1725033),
          (80143857, 25510582),
          (245850922, 78256779)])

    TESTS::

        sage: continued_fraction_list(1 + 10^-10, nterms=3)
        [1, 10000000000]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, nterms=3)
        [1, 100000000000000000000, 2688]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, nterms=5)
        [1, 100000000000000000000, 2688, 8, 1]
        sage: continued_fraction_list(1 + 10^-20 - e^-100, nterms=5)
        [1, 100000000000000000000, 2688, 8, 1]
    """
    if type == "hj":
        if bits is not None:
            x = RealIntervalField(bits)(x).simplest_rational()
        else:
            try:
                x = QQ(x)
            except (TypeError,ValueError):
                x = QQ(x.n())
        l = x.continued_fraction_list(type="hj")
        ## The C-code in sage.rings.rational is much more faster than the pure
        ## Python below
        ## v = []
        ## while True:
        ##    div, mod = divmod(x.numerator(), x.denominator())
        ##    if mod == 0:
        ##        v.append(div)
        ##        break
        ##    v.append(div+1)
        ##    if nterms is not None and len(v) >= nterms:
        ##        break
        ##    x = 1/(div+1-x)
        ## return v
        if nterms is None:
            return l
        return l[:nterms]
    if type != "std":
        raise ValueError("type must be either \"std\" or \"hj\"")
    if bits is not None:
        x = RealIntervalField(bits)(x)
    cf = CFF(x)
    if nterms:
        limit = min(cf.length(), nterms)
    elif cf.length() != Infinity:
        limit = cf.length()
    else:
        import warnings
        warnings.warn("the continued fraction of %s seems infinite, return only the first 20 terms" % x)
        limit = 20
    if partial_convergents:
        return [cf.quotient(i) for i in xrange(limit)], [(cf.pn(i),cf.qn(i)) for i in xrange(limit)]
    return [cf.quotient(i) for i in xrange(limit)]


def continued_fraction(x, bits=None, nterms=None):
    r"""
    Return the continued fraction of ``x``.

    INPUT:

        - `x` -- a number or a list of digits or two list of digits (preperiod
          and period)

        - ``bits`` -- None (default) or a positive integer

        - ``nterms`` -- None (default) or a positive integer

        - ``bits`` -- None (default) or a positive integer. If not ``None`` then
          use an approximation of ``x`` with ``bits`` precision.

        - ``nterms`` -- None (default) or a positive integer. If not ``None``
          then return the continued fraction that consists only on the first
          ``nterms`` of the continued fraction of ``x``.

    EXAMPLES:

    A continued fraction may be initialized by a number or by its digit
    expansion::

        sage: continued_fraction(12/571)
        [0; 47, 1, 1, 2, 2]
        sage: continued_fraction([3,2,1,4])
        [3; 2, 1, 4]

        sage: c = continued_fraction(golden_ratio); c
        [(1)*]
        sage: c.convergent(12)
        377/233
        sage: fibonacci(14)/fibonacci(13)
        377/233

        sage: c = continued_fraction(pi); c
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        sage: a = c.convergent(3); a
        355/113
        sage: a.n()
        3.14159292035398
        sage: pi.n()
        3.14159265358979

    When possible, it is adviced to use anything else but the symbolic ring.
    Here we present an other way of dealing with the golden mean and the cube
    root of 2::

        sage: K.<sqrt5> = NumberField(x^2-5, embedding=2.23)
        sage: my_golden_ratio = (1 + sqrt5)/2
        sage: cf = continued_fraction((1+sqrt5)/2); cf
        [(1)*]
        sage: cf.convergent(12)
        377/233
        sage: cf.period()
        (1,)

        sage: K.<a> = NumberField(x^3-2, embedding=1.25)
        sage: cf = continued_fraction(a); cf
        [1; 3, 1, 5, 1, 1, 4, 1, 1, 8, 1, 14, 1, 10, 2, 1, 4, 12, 2, 3, ...]

    It is possible to pass from quadratic numbers to their ultimately periodic
    continued fraction expansion::

        sage: c = continued_fraction(my_golden_ratio); c
        [(1)*]
        sage: c.preperiod()
        ()
        sage: c.period()
        (1,)

        sage: c = continued_fraction(2/3+sqrt5/5); c
        [1; 8, (1, 3, 1, 1, 3, 9)*]
        sage: c.preperiod()
        (1, 8)
        sage: c.period()
        (1, 3, 1, 1, 3, 9)

        sage: cf = continued_fraction([(1,1),(2,8)]); cf
        [1; 1, (2, 8)*]
        sage: cf.value()
        2/11*sqrt5 + 14/11

    Some well known mathematical constants have continued fractions with a lot
    of structure::

        sage: continued_fraction(e)
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, ...]
        sage: continued_fraction(tan(1))
        [1; 1, 1, 3, 1, 5, 1, 7, 1, 9, 1, 11, 1, 13, 1, 15, 1, 17, 1, 19, ...]
        sage: continued_fraction(tanh(1))
        [0; 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, ...]

    But others do not::

        sage: continued_fraction(pi)
        [3; 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, ...]
        sage: continued_fraction(euler_gamma)
        [0; 1, 1, 2, 1, 2, 1, 4, 3, 13, 5, 1, 1, 8, 1, 2, 4, 1, 1, 40, ...]

    Note that initial rounding can result in incorrect trailing digits::

        sage: continued_fraction(RealField(39)(e))
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 2]
        sage: continued_fraction(RealIntervalField(39)(e))
        [2; 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10]
    """
    return CFF(x, bits=bits, nterms=nterms)


def Hirzebruch_Jung_continued_fraction_list(x, bits=None, nterms=None):
    r"""
    Return the Hirzebruch-Jung continued fraction of ``x`` as a list.

    This function is deprecated since :trac:`14567`. See
    :func:`continued_fraction_list` and the documentation therein.

    INPUT:

    - ``x`` -- exact rational or something that can be numerically
      evaluated. The number to compute the continued fraction of.

    - ``bits`` -- integer (default: the precision of ``x``). the
      precision of the real interval field that is used
      internally. This is only used if ``x`` is not an exact fraction.

    - ``nterms`` -- integer (default: None). The upper bound on the
      number of terms in the continued fraction expansion to return.
      A lits of integers, the coefficients in the Hirzebruch-Jung continued
      fraction expansion of ``x``.

    EXAMPLES::

        sage: Hirzebruch_Jung_continued_fraction_list(17/11)
        doctest:...: DeprecationWarning: Hirzebruch_Jung_continued_fraction_list(x) is replaced by
            continued_fraction_list(x,type="hj")
        or for rationals
            x.continued_fraction_list(type="hj")
        See http://trac.sagemath.org/14567 for details.
        [2, 3, 2, 2, 2, 2]
    """
    from sage.misc.superseded import deprecation
    deprecation(14567, 'Hirzebruch_Jung_continued_fraction_list(x) is replaced by\n\tcontinued_fraction_list(x,type="hj")\nor for rationals\n\tx.continued_fraction_list(type="hj")')
    return continued_fraction_list(x, type="hj", bits=bits, nterms=nterms)

# Unpickling support is needed as the file sage.rings.contfrac is renamed
# sage.rings.continued_fractions in #14567
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.rings.contfrac', 'ContinuedFractionField_class',ContinuedFractionField)
