# -*- coding: utf-8 -*-
r"""
A ``PowComputer`` for relative extensions

This module provides helper classes for the various kinds of relative `p`-adic
extensions. You should never have to access these directly, unless you are
working on linkages or other low-level `p`-adics code within the Sage library.

AUTHORS:

- David Roe, Julian Rüth (2017-06-11): initial version

"""
#*****************************************************************************
#       Copyright (C) 2017 David Roe <roed.math@gmail.com>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from cysignals.memory cimport sig_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from sage.libs.gmp.mpz cimport mpz_init, mpz_clear, mpz_pow_ui

from cpython.object cimport Py_EQ, Py_NE
from sage.structure.richcmp cimport richcmp_not_equal
from sage.rings.integer cimport Integer
from sage.rings.all import ZZ
from sage.misc.cachefunc import cached_method

cdef class PowComputer_relative(PowComputer_class):
    r"""
    Base class for a ``PowComputer`` for use in `p`-adics implemented by Sage
    Polynomials.

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
        sage: R.<a> = ZqFM(25)
        sage: S.<x> = R[]
        sage: f = x^3 - 5*x - 5*a
        sage: PC = PowComputer_relative_maker(3, 20, 20, 60, False, f, 'fixed-mod')

    TESTS::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative
        sage: isinstance(PC, PowComputer_relative)

    """
    def __cinit__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly, shift_seed):
        r"""
        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: PC = PowComputer_relative_maker(3, 20, 20, 60, False, f, 'fixed-mod')

        """
        self.__allocated = 4

    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly, shift_seed):
        r"""
        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod') # indirect doctest
            sage: TestSuite(PC).run()

        """
        PowComputer_class.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)
        self.e = poly.degree()
        self.f = 1

        self.modulus = poly

        self.tmp_cconv_out = poly.parent()()
        self.tmp_clist = poly.parent()()
        self.tmp_ccmp_a = poly.parent()()
        self.tmp_ccmp_b = poly.parent()()

        self.base_ring = poly.base_ring()
        self.poly_ring = poly.parent()

    def __dealloc__(self):
        r"""
        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]
            sage: f = x^3 - 5*x - 5*a
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod')
            sage: del PC

        """

    def __reduce__(self):
        r"""
        Return a picklable representation of this ``PowComputer``.

        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]
            sage: f = x^3 - 5*x - 5*a
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod') # indirect doctest
            sage: loads(dumps(PC)) == PC

        """
        return PowComputer_relative_maker, (self.prime, self.cache_limit, self.prec_cap, self.ram_prec_cap, self.in_field, self.polynomial(), self._prec_type)

    def _repr_(self):
        r"""
        Return a string representation of this ``PowComputer``.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]
            sage: f = x^3 - 5*x - 5*a
            sage: PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod') # indirect doctest
            Relative PowComputer for modulus x^3 - 5*x - 5*a

        """
        return "Relative PowComputer for modulus %s"%(self.modulus,)

    cdef unsigned long capdiv(self, unsigned long n):
        r"""
        Return `\lceil n/e \rceil`.
        """
        if self.e == 1: return n
        if n == 0: return 0
        return (n - 1)/self.e + 1

    def polynomial(self, n=None, var='x'):
        r"""
        Return the modulus of the `p`-adic extension that is handled by this
        ``PowComputer``.

        EXAMPLES::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_maker
            sage: R.<a> = ZqFM(25)
            sage: S.<x> = R[]; f = x^3 - 5*x - 5*a
            sage: PC = PowComputer_relative_maker(5, 20, 20, 60, False, f, 'fixed-mod') # indirect doctest
            sage: PC.polynomial() is f
            True

        """
        return self.modulus


cdef class PowComputer_relative_eis(PowComputer_relative):
    r"""
    A ``PowComputer`` for a relative extension defined by an Eisenstein polynomial

    For a description of inputs see :func:`PowComputer_relative_maker`.

    EXAMPLES::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_eis
        sage: R.<x> = ZZ[]
        sage: f = x^3 - 25*x + 5
        sage: PC = PowComputer_relative_eis(5, 20, 20, 60, False, f)

    TESTS::

        sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_eis
        sage: isinstance(PC, PowComputer_relative_eis)

    """
    def __init__(self, Integer prime, long cache_limit, long prec_cap, long ram_prec_cap, bint in_field, poly, shift_seed):
        r"""
        TESTS::

            sage: from sage.rings.padics.pow_computer_relative import PowComputer_relative_eis
            sage: R.<x> = ZZ[]
            sage: f = x^3 - 25*x + 5
            sage: A = PowComputer_relative_eis(5, 20, 20, 60, False, f)
            sage: TestSuite(A).run()

        """
        PowComputer_relative.__init__(self, prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)

        self.e = self.modulus.degree()
        self.f = 1
        self._shift_seed = shift_seed

    cdef Polynomial_generic_dense invert(self, Polynomial_generic_dense a, long prec):
        r"""
        Return the inverse of ``a``.

        INPUT:

        - ``a`` -- a `p`-adic element, represented as a reduced
          Polynomial in ``poly_ring``

        - ```prec`` -- a ``long``, the required precision

        OUTPUT:

        A polynomial ``b`` such that ``a*b`` is one modulo `π^\mathrm{prec}`

        EXAMPLES::

            sage: TODO: invert something

        """
        # TODO: At the moment, we use an xgcd. We should use Newton iterations
        # instead to make this faster.
        K = self.base_ring.fraction_field()
        Qpmodulus = self.modulus.change_ring(K)
        one, _, b = Qpmodulus.xgcd(a)
        if not one.is_one():
            raise ValueError("element has no inverse")
        return b.change_ring(self.base_ring)

    @cached_method
    def px_pow(self, r):
        r"""
        Return `p/π^r` where π is the uniformizer and `p` is the uniformizer of
        the base ring (not necessarily an integer.)

        INPUT:

        - ``r`` -- a non-negative integer

        OUTPUT:

        A reduced polynomial in π

        EXAMPLES::

            sage: TODO

        """
        if r < 0:
            raise ValueError("r must be non-negative")
        elif r == 0:
            return self.poly_ring(self.base_ring.uniformizer())
        elif r < self.e:
            # Split modulus in half:
            # modulus = p*c - π^r*d, where c and d are integral polynomials, and c has unit const term.
            # Then p/π^r = d/c.
            R = self.modulus.parent()
            c = R(self._shift_seed.list()[:r])
            d = -R(self.modulus.list()[r:])
            c_inverse = self.invert(c, self.ram_prec_cap)
            return (d*c_inverse) % self.modulus
        elif r == self.e:
            return self.invert(self._shift_seed, self.ram_prec_cap)
        else:
            raise NotImplementedError

    @cached_method
    def pxe_pow(self, r):
        r"""
        Return the ``r``-th power of the unit `p/π^e` where `e` is the relative
        ramification index and `p` is the uniformizer of the base ring (not necessarily an integer.)

        INPUT:

        - ``r`` -- a non-negative integer

        OUTPUT:

        A reduced polynomial in π

        EXAMPLES::

            sage: TODO

        """
        if r < 0:
            raise ValueError("r must be non-negative")
        elif r == 0:
            return self.poly_ring.one()
        elif r == 1:
            return self.px_pow(self.e)
        elif r%2:
            return (self.pxe_pow(r-1) * self.pxe_pow(1)) % self.modulus
        else:
            return (self.pxe_pow(r//2)*self.pxe_pow(r//2)) % self.modulus

    @cached_method
    def uniformizer_pow(self, r):
        r"""
        Return the ``r``-th power of the uniformizer.

        INPUT:

        - ``r`` -- a non-negative integer

        OUTPUT:

        A reduced polynomial in π

        EXAMPLES::

            sage: TODO

        """
        if r < 0:
            raise ValueError("r must be non-negative")
        elif r == 0:
            return self.poly_ring.one()
        elif r < self.e:
            return self.poly_ring.one() << r
        elif r%2:
            return (self.uniformizer_pow(r-1) << 1) % self.modulus
        else:
            return (self.uniformizer_pow(r//2) * self.uniformizer_pow(r//2)) % self.modulus

def PowComputer_relative_maker(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed, prec_type):
    r"""
    Create a ``PowComputer``.

    INPUT:

    - ``prime`` -- a uniformizer in the base ring

    - ``cache_limit`` -- a non-negative integer, controlling the caching. The
      ``PowComputer`` caches frequently used things like powers of ``prime``.
      This parameter, e.g., controls up to which power these are cached.

    - ``prec_cap`` -- the power of ``prime`` modulo which elements of largest
      precision are defined

    - ``ram_prec_cap`` -- approximately ``e*prec_cap``, where ``e`` is
      the relative ramification degree of the extension.  For a ramified
      extension this is what Sage calls the precision cap of the ring.  In
      fact, it's possible to have rings with precision cap not a multiple of
      `e`, in which case the actual relationship between ``ram_prec_cap`` and
      ``prec_cap`` is that ``prec_cap = ceil(n/e)``

    - ``in_field`` -- a boolean; whether the associated ring is actually a
      field

    - ``poly`` -- the polynomial defining the extension

    - `prec_type`` -- one of ``"capped-rel"``, ``"capped-abs"`` or
      ``"fixed-mod"``, ``"floating-point"``, the precision type of the ring

    .. NOTE::

        Because of the way templates work, this function imports the class of
        its return value from the appropriate element files.  This means that
        the returned PowComputer will have the appropriate compile-time-type
        for Cython.

    EXAMPLES::

        sage: TODO

    """
    PC = PowComputer_relative_eis(prime, cache_limit, prec_cap, ram_prec_cap, in_field, poly, shift_seed)
    PC._prec_type = prec_type
    return PC
