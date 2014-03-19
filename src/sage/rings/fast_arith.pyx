"""
Basic arithmetic with c-integers.
"""

#*****************************************************************************
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


###################################################################
# We define the following functions in this file, both
# for int (up to bound = 2**31 - 1) and longlong (up to 2**63 - 1).
# The function definitions are identical except for the types.
# Some of their input can be at most sqrt(bound), since
# it is necessary to multiply numbers and reduce the product
# modulo n, where n is at most bound.
#
#   * abs_int -- absolute value of integer
#   * sign_int -- sign of integer
#   * c_gcd_int -- gcd of two ints
#   * gcd_int -- python export of c_gcd_int
#   * c_xgcd_int -- extended gcd of two ints
#   * c_inverse_mod_int -- inverse of an int modulo another int
#   * c_rational_recon_int -- rational reconstruction of ints
#   * rational_recon_int -- export of rational reconstruction for ints
#
#  The long long functions are the same, except they end in _longlong.
#
###################################################################

# The int definitions

include "sage/ext/gmp.pxi"
include "sage/ext/stdsage.pxi"
include "sage/libs/pari/decl.pxi"

cdef extern from "pari/pari.h":
    cdef long NEXT_PRIME_VIADIFF(long, unsigned char*)

from sage.libs.pari.all import pari
from sage.rings.integer cimport Integer

from libc.stdint cimport uint_fast32_t, uint_fast64_t

cpdef prime_range(start, stop=None, algorithm="pari_primes", bint py_ints=False):
    r"""
    List of all primes between start and stop-1, inclusive.  If the
    second argument is omitted, returns the primes up to the first
    argument.

    This function is closely related to (and can use) the primes iterator.
    Use algorithm "pari_primes" when both start and stop are not too large,
    since in all cases this function makes a table of primes up to
    stop. If both are large, use algorithm "pari_isprime" instead.

    Algorithm "pari_primes" is faster for most input, but crashes for larger input.
    Algorithm "pari_isprime" is slower but will work for much larger input.

    INPUT:

        - ``start`` -- lower bound

        - ``stop`` -- upper bound

        - ``algorithm`` -- string, one of:

             - "pari_primes": Uses PARI's primes function.  Generates all primes up to stop.
                              Depends on PARI's primepi function.

             - "pari_isprime": Uses a mod 2 wheel and PARI's isprime function by calling
                             the primes iterator.

        - ``py_ints`` -- boolean (default False), return Python ints rather than Sage Integers (faster)


    EXAMPLES:
        sage: prime_range(10)
        [2, 3, 5, 7]
        sage: prime_range(7)
        [2, 3, 5]
        sage: prime_range(2000,2020)
        [2003, 2011, 2017]
        sage: prime_range(2,2)
        []
        sage: prime_range(2,3)
        [2]
        sage: prime_range(5,10)
        [5, 7]
        sage: prime_range(-100,10,"pari_isprime")
        [2, 3, 5, 7]
        sage: prime_range(2,2,algorithm="pari_isprime")
        []
        sage: prime_range(10**16,10**16+100,"pari_isprime")
        [10000000000000061, 10000000000000069, 10000000000000079, 10000000000000099]
        sage: prime_range(10**30,10**30+100,"pari_isprime")
        [1000000000000000000000000000057, 1000000000000000000000000000099]
        sage: type(prime_range(8)[0])
        <type 'sage.rings.integer.Integer'>
        sage: type(prime_range(8,algorithm="pari_isprime")[0])
        <type 'sage.rings.integer.Integer'>

    TESTS:
        sage: len(prime_range(25000,2500000))
        180310
        sage: prime_range(2500000)[-1].is_prime()
        True

    AUTHORS:
      - William Stein (original version)
      - Craig Citro (rewrote for massive speedup)
      - Kevin Stueve (added primes iterator option) 2010-10-16
      - Robert Bradshaw (speedup using Pari prime table, py_ints option)
    """
    cdef Integer z
    cdef long c_start, c_stop, p
    cdef unsigned char* pari_prime_ptr
    if algorithm == "pari_primes":
        if stop is None:
            # In this case, "start" is really stop
            c_start = 1
            c_stop = start
        else:
            c_start = start
            c_stop = stop
            if c_stop <= c_start:
                return []
            if c_start < 1:
                c_start = 1
        if maxprime() < c_stop:
            pari.init_primes(c_stop)
        pari_prime_ptr = <unsigned char*>diffptr
        p = 0
        res = []
        while p < c_start:
            NEXT_PRIME_VIADIFF(p, pari_prime_ptr)
        while p < c_stop:
            if py_ints:
                res.append(p)
            else:
                z = <Integer>PY_NEW(Integer)
                mpz_set_ui(z.value, p)
                res.append(z)
            NEXT_PRIME_VIADIFF(p, pari_prime_ptr)

    elif algorithm == "pari_isprime":
        from sage.rings.arith import primes
        res = list(primes(start, stop))
    else:
        raise ValueError("algorithm argument must be either ``pari_primes`` or ``pari_isprime``")
    return res

cdef class QualityBounder(object):
    """
    Given projective rational coordinates, returns a lower and upper
    bound on the quality of the projective point.

    INPUT:

        - ``coordinates`` - projective coordinates

    EXAMPLES::

        sage: from sage.rings.fast_arith import QualityBounder
        sage: qb = QualityBounder()
        sage: qb(1, 8, 9)
        (log(9)/log(6), log(9)/log(6))
        sage: for i in range(100):
        ....:     x = 1
        ....:     while not x.is_squarefree() or x == 1:
        ....:         x = ZZ.random_element(-10**5, 10**5)
        ....:     if qb(x, 1) != (1, 1):
        ....:         print x
        sage: qb(ProjectiveSpace(2, QQ)(3, 4, 5))
        (log(5)/log(30), log(5)/log(30))
    """

    cdef uint_fast32_t *_primes
    cdef uint_fast32_t _len, _prime_bound

    def __cinit__(self, uint_fast32_t prime_bound=10**5):
        """
        Allocates memory for the c tables used for trial division.
        """

        if prime_bound > 1:
            self._prime_bound = prime_bound

            from sage.functions.prime_pi import prime_pi
            self._len = prime_pi(prime_bound-1)
        else:
            self._prime_bound = 1u
            self._len = 0u

        # use c tables
        self._primes = <uint_fast32_t *>sage_malloc(sizeof(uint_fast32_t)*self._len)

    def __dealloc__(self):
        """
        Deallocates memory for the c tables used for trial division.
        """
        sage_free(self._primes)

    def __init__(self, uint_fast32_t prime_bound=10**5):
        """
        Constructs an object for bounding the quality of projective
        points over the rationals.

        INPUT:

            - ``prime_bound`` - upper bound used in trial division

        EXAMPLES::

            sage: from sage.rings.fast_arith import QualityBounder
            sage: qb1 = QualityBounder(prime_bound=10**5)
            sage: qb2 = QualityBounder(prime_bound=10**6)
            sage: map(n, qb1(614474125456, 24966161193627, 25580635319083))
            [0.374709057573204, 0.698836005185383]
            sage: map(n, qb2(614474125456, 24966161193627, 25580635319083))
            [0.536657184562553, 0.536657184562553]
        """
        tmp = prime_range(prime_bound, py_ints=True)

        cdef uint_fast32_t i
        # populate prime table
        for i in range(self._len):
            self._primes[i] = tmp.pop()

    @staticmethod
    def _fix_up_coordinates(*coordinates):
        """
        Fixes up projective coordinates to relatively prime
        positive integral coordinates.

        INPUT:

            - ``coordinates`` - projective coordinates

        EXAMPLES::

            sage: from sage.rings.fast_arith import QualityBounder
            sage: QualityBounder._fix_up_coordinates(2/3, 4/5, 1)
            [10, 12, 15]
            sage: QualityBounder._fix_up_coordinates(2/3, 1/2, 7/3, -19/6)
            [4, 3, 14, 19]
            sage: QualityBounder._fix_up_coordinates(ProjectiveSpace(3, QQ)(1/2, 1/3, 1/5, 1/7))
            [105, 70, 42, 30]
            sage: QualityBounder._fix_up_coordinates(1/2, 0)
            [1, 0]
        """
        if len(coordinates) == 1:
            try:
                coordinates = tuple(coordinates[0])
            except TypeError:
                coordinates = (coordinates[0],)

        # fixup coordinates to be integral
        from sage.rings.rational_field import QQ
        try:
            coordinates = map(QQ, coordinates)
        except TypeError:
            raise NotImplementedError(
                    "quality bounds are only implemented over the rationals")

        for q in coordinates:
            for i in range(len(coordinates)):
                coordinates[i] *= q.denominator()

        coordinates = map(Integer, coordinates)
        coordinates = map(abs, coordinates)

        if max(coordinates) == 0:
            raise ValueError("projective points must have a non-zero coordinate")

        # need to divide out by the gcd
        from sage.rings.integer import GCD_list
        cdef Integer gcd = GCD_list(coordinates)
        cdef Integer tmp

        for tmp in coordinates:
            # use mpz directly to avoid having to recast
            # integral rationals as integers
            mpz_divexact(tmp.value, tmp.value, gcd.value)

        return coordinates

    def _rad_bounds(self, *coordinates):
        """
        Underlying method that gets bounds on the radical (and hence
        the quality).

        INPUT:

            - ``coordinates`` -- projective coordinates, assumed to be
              relatively prime non-negative Sage Integers
              (not all zero)

        WARNING:

            This will modify the values of the input integers.

        EXAMPLES::

            sage: from sage.rings.fast_arith import QualityBounder
            sage: qb = QualityBounder()
            sage: a = next_prime(2*10**5)**3*2**8
            sage: b = next_prime(4*10**5)**2*3**4
            sage: c = a+b
            sage: a,b,c
            (2048092161382406912, 12960583206561, 2048105121965613473)
            sage: qb._rad_bounds(a, b, c)
            (399793800000, 15730863037389643841037540121480095666228392706)
            sage: a,b,c
            (8000360005400027, 160007200081, 3073742197051)
        """
        cdef uint_fast32_t i, prime, not_exact
        cdef uint_fast64_t prime_bound_squared
        cdef Integer rad1, rad2, tmp, tmp2

        rad1 = <Integer>PY_NEW(Integer)
        rad2 = <Integer>PY_NEW(Integer)
        not_exact = len(coordinates)
        prime_bound_squared = <uint_fast64_t>self._prime_bound
        prime_bound_squared *= prime_bound_squared

        sig_on()

        mpz_set_ui(rad1.value, 1ul)

        # trial divide primes
        for tmp in coordinates:
            if not mpz_cmp_ui(tmp.value, 1ul):
                not_exact -= 1u
                continue
            for i in range(self._len):
                if mpz_divisible_ui_p(tmp.value, self._primes[i]):
                    prime = self._primes[i]
                    if not mpz_divisible_ui_p(rad1.value, prime):
                        mpz_mul_ui(rad1.value, rad1.value, prime)
                    for tmp2 in coordinates:
                        while mpz_divisible_ui_p(tmp2.value, prime):
                            mpz_divexact_ui(tmp2.value, tmp2.value, prime)
                    if not mpz_cmp_ui(tmp.value, 1ul):
                        not_exact -= 1u
                        break

        for tmp in coordinates:
            if not mpz_cmp_ui(tmp.value, 1ul):
                continue
            if mpz_cmp_ui(tmp.value, prime_bound_squared) < 0:
                # must be a prime, via sieving

                mpz_mul(rad1.value, rad1.value, tmp.value)
                # store value, since it will be changed to 1 in for loop
                mpz_set(rad2.value, tmp.value)
                for tmp2 in coordinates:
                    while mpz_divisible_p(tmp2.value, rad2.value):
                        mpz_divexact(tmp2.value, tmp2.value, rad2.value)
                not_exact -= 1u

        if not_exact:
            # lower bound = prod p * prime_bound
            mpz_mul_ui(rad2.value, rad1.value, self._prime_bound)

            # upper bound = prod p * remainder
            for tmp in coordinates:
                mpz_mul(rad1.value, rad1.value, tmp.value)

        else:
            mpz_set(rad2.value, rad1.value)

        sig_off()

        return rad2, rad1

    def __call__(self, *coordinates):
        """
        TESTS::

            sage: from sage.rings.fast_arith import QualityBounder
            sage: qb = QualityBounder(2)
            sage: for p in Permutations((1, 4, 8, 16)):
            ....:     assert qb(*p) == (log(16, 512), 4)
            sage: qb = QualityBounder(3)
            sage: for p in Permutations((1, 4, 8, 16)):
            ....:     assert qb(*p) == (4, 4)
            sage: qb(1, 8, 9)
            (log(9)/log(18), log(9)/log(6))
            sage: QualityBounder(4)(1, 8, 9)
            (log(9)/log(6), log(9)/log(6))
        """
        coordinates = QualityBounder._fix_up_coordinates(*coordinates)

        if min(coordinates) == 0:
            return 0, 0

        height = Integer(max(coordinates)) # make explicit copy
        rad1, rad2 = self._rad_bounds(*coordinates)

        from sage.functions.log import log
        # deal with annoying edge cases where the radical is 1
        try:
            ret1 = log(height, rad2)
        except (ValueError, ZeroDivisionError):
            ret1 = Integer(1)/Integer(3)
        try:
            ret2 = log(height, rad1)
        except (ValueError, ZeroDivisionError):
            from sage.rings.infinity import infinity
            ret2 = infinity

        return ret1, ret2

cdef class arith_int:
    cdef public int abs_int(self, int x) except -1:
        if x < 0:
            return -x
        return x

    cdef public int sign_int(self, int n) except -2:
        if n < 0:
            return -1
        return 1

    cdef public int c_gcd_int(self, int a, int b) except -1:
        cdef int c
        if a==0:
            return self.abs_int(b)
        if b==0:
            return self.abs_int(a)
        if a<0: a=-a
        if b<0: b=-b
        while(b):
            c = a % b
            a = b
            b = c
        return a


    def gcd_int(self, int a, int b):
        return self.c_gcd_int(a,b)


    cdef public int c_xgcd_int(self, int a, int b, int* ss, int* tt) except -1:
        cdef int psign, qsign, p, q, r, s, c, quot, new_r, new_s

        if a == 0:
            ss[0] = 0
            tt[0] = self.sign_int(b)
            return self.abs_int(b)

        if b == 0:
            ss[0] = self.sign_int(a)
            tt[0] = 0
            return self.abs_int(a)

        psign = 1; qsign = 1

        if a<0: a = -a; psign = -1
        if b<0: b = -b; qsign = -1

        p = 1; q = 0; r = 0; s = 1
        while (b):
            c = a % b; quot = a/b
            a = b; b = c
            new_r = p - quot*r
            new_s = q - quot*s
            p = r; q = s
            r = new_r; s = new_s

        ss[0] = p*psign
        tt[0] = q*qsign

        return a

    def xgcd_int(self, int a, int b):
        cdef int g, s, t
        g = self.c_xgcd_int(a,b, &s, &t)
        return (g,s,t)

    cdef public int c_inverse_mod_int(self, int a, int m) except -1:
        if a == 1 or m<=1: return a%m   # common special case
        cdef int g, s, t
        g = self.c_xgcd_int(a,m, &s, &t)
        if g != 1:
            raise ArithmeticError, "The inverse of %s modulo %s is not defined."%(a,m)
        s = s % m
        if s < 0:
            s = s + m
        return s


    def inverse_mod_int(self, int a, int m):
        return self.c_inverse_mod_int(a, m)

    cdef int c_rational_recon_int(self, int a, int m, int* n, int* d) except -1:
        cdef int u, v, u0, u1, u2, v0, v1, v2, q, t0, t1, t2, x, y
        cdef float bnd

        if m>46340:
            raise OverflowError, "The modulus m(=%s) should be at most 46340"%m
            return -1

        a = a % m

        if a==0 or m == 0:
            n[0] = 0
            d[0] = 1
            return 0

        if m<0: m = -m
        if a<0: a = m - a
        if a==1:
            n[0] = 1
            d[0] = 1
            return 0

        u = m
        v = a
        bnd = sqrt(m/2.0)
        u0=1; u1=0; u2=u
        v0=0; v1=1; v2=v
        while self.abs_int(v2) > bnd:
            q = u2/v2   # floor is implicit
            t0=u0-q*v0; t1=u1-q*v1; t2=u2-q*v2
            u0=v0; u1=v1; u2=v2
            v0=t0; v1=t1; v2=t2;

        x = self.abs_int(v1); y = v2
        if v1<0:  y = -1*y
        if x<=bnd and self.c_gcd_int(x,y)==1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_int(self, int a, int m):
        """
        Rational reconstruction of a modulo m.
        """
        cdef int n, d
        self.c_rational_recon_int(a, m, &n, &d)
        return (n,d)


# The long long versions are next.
cdef class arith_llong:

    cdef public long long abs_longlong(self, long long x) except -1:
        if x < 0:
            return -x
        return x

    cdef public long long sign_longlong(self, long long n) except -2:
        if n < 0:
            return -1
        return 1

    cdef public long long c_gcd_longlong(self, long long a, long long b) except -1:
        cdef long long c
        if a==0:
            return self.abs_longlong(b)
        if b==0:
            return self.abs_longlong(a)
        if a<0: a=-a
        if b<0: b=-b
        while(b):
            c = a % b
            a = b
            b = c
        return a


    def gcd_longlong(self, long long a, long long b):
        return self.c_gcd_longlong(a,b)


    cdef public long long c_xgcd_longlong(self, long long a, long long b,
                                          long long *ss,
                                          long long *tt) except -1:
        cdef long long psign, qsign, p, q, r, s, c, quot, new_r, new_s


        if a == 0:
            ss[0] = 0
            tt[0] = self.sign_longlong(b)
            return self.abs_longlong(b)

        if b == 0:
            ss[0] = self.sign_longlong(a)
            tt[0] = 0
            return self.abs_longlong(a)

        psign = 1; qsign = 1

        if a<0: a = -a; psign = -1
        if b<0: b = -b; qsign = -1

        p = 1; q = 0; r = 0; s = 1
        while (b):
            c = a % b; quot = a/b
            a = b; b = c
            new_r = p - quot*r
            new_s = q - quot*s
            p = r; q = s
            r = new_r; s = new_s

        ss[0] = p*psign
        tt[0] = q*qsign


        return a

    cdef public long long c_inverse_mod_longlong(self, long long a, long long m) except -1:
        cdef long long g, s, t
        g = self.c_xgcd_longlong(a,m, &s, &t)
        if g != 1:
            raise ArithmeticError("The inverse of %s modulo %s is not defined."%(a,m))
        s = s % m
        if s < 0:
            s = s + m
        return s

    def inverse_mod_longlong(self, long long a, long long m):
        return self.c_inverse_mod_longlong(a, m)

    cdef long long c_rational_recon_longlong(self, long long a, long long m,
                                             long long *n, long long *d) except -1:
        cdef long long u, v, u0, u1, u2, v0, v1, v2, q, t0, t1, t2, x, y
        cdef float bnd

        if m > 2147483647:
            raise OverflowError, "The modulus m(=%s) must be at most 2147483647"%m
            return -1

        a = a % m

        if a==0 or m == 0:
            n[0] = 0
            d[0] = 1
            return 0

        if m<0: m = -m
        if a<0: a = m - a
        if a==1:
            n[0] = 1
            d[0] = 1
            return 0

        u = m
        v = a
        bnd = sqrt(m/2.0)
        u0=1; u1=0; u2=u
        v0=0; v1=1; v2=v
        while self.abs_longlong(v2) > bnd:
            q = u2/v2   # floor is implicit
            t0=u0-q*v0; t1=u1-q*v1; t2=u2-q*v2
            u0=v0; u1=v1; u2=v2
            v0=t0; v1=t1; v2=t2;

        x = self.abs_longlong(v1); y = v2
        if v1<0:  y = -1*y
        if x<=bnd and self.gcd_longlong(x,y)==1:
            n[0] = y
            d[0] = x
            return 0

        n[0] = 0
        d[0] = 0
        return 0

    def rational_recon_longlong(self, long long a, long long m):
        """
        Rational reconstruction of a modulo m.
        """
        cdef long long n, d
        self.c_rational_recon_longlong(a, m, &n, &d)
        return (n,d)
