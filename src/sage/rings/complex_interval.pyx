"""
Arbitrary Precision Complex Intervals

This is a simple complex interval package, using intervals which are
axis-aligned rectangles in the complex plane.  It has very few special
functions, and it does not use any special tricks to keep the size of
the intervals down.

AUTHORS:

These authors wrote ``complex_number.pyx``:

- William Stein (2006-01-26): complete rewrite
- Joel B. Mohler (2006-12-16): naive rewrite into pyrex
- William Stein(2007-01): rewrite of Mohler's rewrite

Then ``complex_number.pyx`` was copied to ``complex_interval.pyx`` and
heavily modified:

- Carl Witty (2007-10-24): rewrite to become a complex interval package

- Travis Scrimshaw (2012-10-18): Added documentation to get full coverage.

.. TODO::

    Implement :class:`ComplexIntervalFieldElement` multiplicative
    order similar to :class:`ComplexNumber` multiplicative
    order with ``_set_multiplicative_order(n)`` and
    :meth:`ComplexNumber.multiplicative_order()` methods.
"""

#################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math
import operator

from sage.structure.element cimport FieldElement, RingElement, Element, ModuleElement
from complex_number cimport ComplexNumber

import complex_interval_field
from complex_field import ComplexField
import sage.misc.misc
import integer
import infinity
import real_mpfi
import real_mpfr
cimport real_mpfr

include "sage/ext/stdsage.pxi"

cdef double LOG_TEN_TWO_PLUS_EPSILON = 3.321928094887363 # a small overestimate of log(10,2)

def is_ComplexIntervalFieldElement(x):
    """
    Check if ``x`` is a :class:`ComplexIntervalFieldElement`.

    EXAMPLES::

        sage: from sage.rings.complex_interval import is_ComplexIntervalFieldElement as is_CIFE
        sage: is_CIFE(CIF(2))
        True
        sage: is_CIFE(CC(2))
        False
    """
    return isinstance(x, ComplexIntervalFieldElement)

cdef class ComplexIntervalFieldElement(sage.structure.element.FieldElement):
    """
    A complex interval.

    EXAMPLES::

        sage: I = CIF.gen()
        sage: b = 1.5 + 2.5*I
        sage: TestSuite(b).run()
    """
    cdef ComplexIntervalFieldElement _new(self):
        """
        Quickly creates a new initialized complex interval with the
        same parent as ``self``.
        """
        cdef ComplexIntervalFieldElement x
        x = PY_NEW(ComplexIntervalFieldElement)
        x._parent = self._parent
        x._prec = self._prec
        mpfi_init2(x.__re, self._prec)
        mpfi_init2(x.__im, self._prec)
        return x

    def __init__(self, parent, real, imag=None):
        """
        Initialize a complex number (interval).

        EXAMPLES::

            sage: CIF(1.5, 2.5)
            1.5000000000000000? + 2.5000000000000000?*I
            sage: CIF((1.5, 2.5))
            1.5000000000000000? + 2.5000000000000000?*I
            sage: CIF(1.5 + 2.5*I)
            1.5000000000000000? + 2.5000000000000000?*I
        """
        cdef real_mpfi.RealIntervalFieldElement rr, ii
        self._parent = parent
        self._prec = self._parent._prec
        mpfi_init2(self.__re, self._prec)
        mpfi_init2(self.__im, self._prec)

        if imag is None:
            if PY_TYPE_CHECK(real, ComplexNumber):
                real, imag = (<ComplexNumber>real).real(), (<ComplexNumber>real).imag()
            elif PY_TYPE_CHECK(real, ComplexIntervalFieldElement):
                real, imag = (<ComplexIntervalFieldElement>real).real(), (<ComplexIntervalFieldElement>real).imag()
            elif isinstance(real, sage.libs.pari.all.pari_gen):
                real, imag = real.real(), real.imag()
            elif isinstance(real, list) or isinstance(real, tuple):
                re, imag = real
                real = re
            elif isinstance(real, complex):
                real, imag = real.real, real.imag
            else:
                imag = 0
        try:
            R = parent._real_field()
            rr = R(real)
            ii = R(imag)
            mpfi_set(self.__re, rr.value)
            mpfi_set(self.__im, ii.value)
        except TypeError:
            raise TypeError, "unable to coerce to a ComplexIntervalFieldElement"


    def  __dealloc__(self):
        mpfi_clear(self.__re)
        mpfi_clear(self.__im)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CIF(1.5) # indirect doctest
            1.5000000000000000?
            sage: CIF(1.5, 2.5) # indirect doctest
            1.5000000000000000? + 2.5000000000000000?*I
        """
        return self.str(10)

    def __hash__(self):
        """
        Return the hash value of ``self``.

        EXAMPLES::

            sage: hash(CIF(1.5)) # indirect doctest
            1517890078            # 32-bit
            -3314089385045448162  # 64-bit
            sage: hash(CIF(1.5, 2.5)) # indirect doctest
            -1103102080           # 32-bit
            3834538979630251904   # 64-bit
        """
        return hash(self.str())

    def __getitem__(self, i):
        """
        Returns either the real or imaginary component of ``self`` depending
        on the choice of ``i``: real (``i=0``), imaginary (``i=1``)

        INPUTS:

        - ``i`` - 0 or 1

          - ``0`` -- will return the real component of ``self``
          - ``1`` -- will return the imaginary component of ``self``

        EXAMPLES::

            sage: z = CIF(1.5, 2.5)
            sage: z[0]
            1.5000000000000000?
            sage: z[1]
            2.5000000000000000?
        """
        if i == 0:
            return self.real()
        elif i == 1:
            return self.imag()
        raise IndexError, "i must be between 0 and 1."

    def __reduce__( self ):
        """
        Pickling support.

        TESTS::

            sage: a = CIF(1 + I)
            sage: loads(dumps(a)) == a
            True
        """
        # TODO: This is potentially slow -- make a 1 version that
        # is native and much faster -- doesn't use .real()/.imag()
        return (make_ComplexIntervalFieldElement0, (self._parent, self.real(), self.imag()))

    def str(self, base=10, style=None):
        """
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: CIF(1.5).str()
            '1.5000000000000000?'
            sage: CIF(1.5, 2.5).str()
            '1.5000000000000000? + 2.5000000000000000?*I'
            sage: CIF(1.5, -2.5).str()
            '1.5000000000000000? - 2.5000000000000000?*I'
            sage: CIF(0, -2.5).str()
            '-2.5000000000000000?*I'
            sage: CIF(1.5).str(base=3)
            '1.1111111111111111111111111111111112?'
            sage: CIF(1, pi).str(style='brackets')
            '[1.0000000000000000 .. 1.0000000000000000] + [3.1415926535897931 .. 3.1415926535897936]*I'

        .. SEEALSO::

            - :meth:`RealIntervalFieldElement.str`
        """
        s = ""
        if not self.real().is_zero():
            s = self.real().str(base=base, style=style)
        if not self.imag().is_zero():
            y  =  self.imag()
            if s!="":
                if y < 0:
                    s = s+" - "
                    y = -y
                else:
                    s = s+" + "
            s = s+"%s*I"%y.str(base=base, style=style)
        if len(s) == 0:
            s = "0"
        return s

    def plot(self, pointsize=10, **kwds):
        r"""
        Plot a complex interval as a rectangle.

        EXAMPLES::

            sage: sum(plot(CIF(RIF(1/k, 1/k), RIF(-k, k))) for k in [1..10])
            Graphics object consisting of 20 graphics primitives

        Exact and nearly exact points are still visible::

            sage: plot(CIF(pi, 1), color='red') + plot(CIF(1, e), color='purple') + plot(CIF(-1, -1))
            Graphics object consisting of 6 graphics primitives

        A demonstration that `z \mapsto z^2` acts chaotically on `|z|=1`::

            sage: z = CIF(0, 2*pi/1000).exp()
            sage: g = Graphics()
            sage: for i in range(40):
            ...       z = z^2
            ...       g += z.plot(color=(1./(40-i), 0, 1))
            ...
            sage: g
            Graphics object consisting of 80 graphics primitives
        """
        from sage.plot.polygon import polygon2d
        x, y = self.real(), self.imag()
        x0, y0 = x.lower(), y.lower()
        x1, y1 = x.upper(), y.upper()
        g = polygon2d([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)],
                thickness=pointsize/4, **kwds)
        # Nearly empty polygons don't show up.
        g += self.center().plot(pointsize= pointsize, **kwds)
        return g

    def _latex_(self):
        """
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: latex(CIF(1.5, -2.5)) # indirect doctest
            1.5000000000000000? - 2.5000000000000000?i
            sage: latex(CIF(0, 3e200)) # indirect doctest
            3.0000000000000000? \times 10^{200}i
        """
        import re
        s = self.str().replace('*I', 'i')
        return re.sub(r"e(-?\d+)", r" \\times 10^{\1}", s)

    def bisection(self):
        """
        Returns the bisection of ``self`` into four intervals whose union is
        ``self`` and intersection is :meth:`center()`.

        EXAMPLES::

            sage: z = CIF(RIF(2, 3), RIF(-5, -4))
            sage: z.bisection()
            (3.? - 5.?*I, 3.? - 5.?*I, 3.? - 5.?*I, 3.? - 5.?*I)
            sage: for z in z.bisection():
            ...       print z.real().endpoints(), z.imag().endpoints()
            (2.00000000000000, 2.50000000000000) (-5.00000000000000, -4.50000000000000)
            (2.50000000000000, 3.00000000000000) (-5.00000000000000, -4.50000000000000)
            (2.00000000000000, 2.50000000000000) (-4.50000000000000, -4.00000000000000)
            (2.50000000000000, 3.00000000000000) (-4.50000000000000, -4.00000000000000)

            sage: z = CIF(RIF(sqrt(2), sqrt(3)), RIF(e, pi))
            sage: a, b, c, d = z.bisection()
            sage: a.intersection(b).intersection(c).intersection(d) == CIF(z.center())
            True

            sage: zz = a.union(b).union(c).union(c)
            sage: zz.real().endpoints() == z.real().endpoints()
            True
            sage: zz.imag().endpoints() == z.imag().endpoints()
            True
        """
        cdef ComplexIntervalFieldElement a00 = self._new()
        mpfr_set(&a00.__re.left, &self.__re.left, GMP_RNDN)
        mpfi_mid(&a00.__re.right, self.__re)
        mpfr_set(&a00.__im.left, &self.__im.left, GMP_RNDN)
        mpfi_mid(&a00.__im.right, self.__im)

        cdef ComplexIntervalFieldElement a01 = self._new()
        mpfr_set(&a01.__re.left, &a00.__re.right, GMP_RNDN)
        mpfr_set(&a01.__re.right, &self.__re.right, GMP_RNDN)
        mpfi_set(a01.__im, a00.__im)

        cdef ComplexIntervalFieldElement a10 = self._new()
        mpfi_set(a10.__re, a00.__re)
        mpfi_mid(&a10.__im.left, self.__im)
        mpfr_set(&a10.__im.right, &self.__im.right, GMP_RNDN)

        cdef ComplexIntervalFieldElement a11 = self._new()
        mpfi_set(a11.__re, a01.__re)
        mpfi_set(a11.__im, a10.__im)

        return a00, a01, a10, a11

    def is_exact(self):
        """
        Returns whether this complex interval is exact (i.e. contains exactly
        one complex value).

        EXAMPLES::

            sage: CIF(3).is_exact()
            True
            sage: CIF(0, 2).is_exact()
            True
            sage: CIF(-4, 0).sqrt().is_exact()
            True
            sage: CIF(-5, 0).sqrt().is_exact()
            False
            sage: CIF(0, 2*pi).is_exact()
            False
            sage: CIF(e).is_exact()
            False
            sage: CIF(1e100).is_exact()
            True
            sage: (CIF(1e100) + 1).is_exact()
            False
        """
        return mpfr_equal_p(&self.__re.left, &self.__re.right) and \
               mpfr_equal_p(&self.__im.left, &self.__im.right)

    def diameter(self):
        """
        Returns a somewhat-arbitrarily defined "diameter" for this interval.

        The diameter of an interval is the maximum of the diameter of the real
        and imaginary components, where diameter on a real interval is defined
        as absolute diameter if the interval contains zero, and relative
        diameter otherwise.

        EXAMPLES::

            sage: CIF(RIF(-1, 1), RIF(13, 17)).diameter()
            2.00000000000000
            sage: CIF(RIF(-0.1, 0.1), RIF(13, 17)).diameter()
            0.266666666666667
            sage: CIF(RIF(-1, 1), 15).diameter()
            2.00000000000000
        """
        cdef real_mpfr.RealNumber diam
        diam = real_mpfr.RealNumber(self._parent._real_field()._middle_field(), None)
        cdef mpfr_t tmp
        mpfr_init2(tmp, self.prec())
        mpfi_diam(diam.value, self.__re)
        mpfi_diam(tmp, self.__im)
        mpfr_max(diam.value, diam.value, tmp, GMP_RNDU)
        mpfr_clear(tmp)
        return diam

    def overlaps(self, ComplexIntervalFieldElement other):
        """
        Returns ``True`` if ``self`` and other are intervals with at least
        one value in common.

        EXAMPLES::

            sage: CIF(0).overlaps(CIF(RIF(0, 1), RIF(-1, 0)))
            True
            sage: CIF(1).overlaps(CIF(1, 1))
            False
        """
        return mpfr_greaterequal_p(&self.__re.right, &other.__re.left) \
           and mpfr_greaterequal_p(&other.__re.right, &self.__re.left) \
           and mpfr_greaterequal_p(&self.__im.right, &other.__im.left) \
           and mpfr_greaterequal_p(&other.__im.right, &self.__im.left)

    def intersection(self, other):
        """
        Returns the intersection of the two complex intervals ``self`` and
        ``other``.

        EXAMPLES::

            sage: CIF(RIF(1, 3), RIF(1, 3)).intersection(CIF(RIF(2, 4), RIF(2, 4))).str(style='brackets')
            '[2.0000000000000000 .. 3.0000000000000000] + [2.0000000000000000 .. 3.0000000000000000]*I'
            sage: CIF(RIF(1, 2), RIF(1, 3)).intersection(CIF(RIF(3, 4), RIF(2, 4)))
            Traceback (most recent call last):
            ...
            ValueError: intersection of non-overlapping intervals
        """

        cdef ComplexIntervalFieldElement x = self._new()
        cdef ComplexIntervalFieldElement other_intv
        if PY_TYPE_CHECK(other, ComplexIntervalFieldElement):
            other_intv = other
        else:
            # Let type errors from _coerce_ propagate...
            other_intv = self._parent(other)

        mpfi_intersect(x.__re, self.__re, other_intv.__re)
        mpfi_intersect(x.__im, self.__im, other_intv.__im)
        if mpfr_less_p(&x.__re.right, &x.__re.left) \
           or mpfr_less_p(&x.__im.right, &x.__im.left):
            raise ValueError, "intersection of non-overlapping intervals"

        return x

    def union(self, other):
        """
        Returns the smallest complex interval including the
        two complex intervals ``self`` and ``other``.

        EXAMPLES::

            sage: CIF(0).union(CIF(5, 5)).str(style='brackets')
            '[0.00000000000000000 .. 5.0000000000000000] + [0.00000000000000000 .. 5.0000000000000000]*I'
        """
        cdef ComplexIntervalFieldElement x = self._new()
        cdef ComplexIntervalFieldElement other_intv
        if PY_TYPE_CHECK(other, ComplexIntervalFieldElement):
            other_intv = other
        else:
            # Let type errors from _coerce_ propagate...
            other_intv = self._parent(other)

        mpfi_union(x.__re, self.__re, other_intv.__re)
        mpfi_union(x.__im, self.__im, other_intv.__im)
        return x

    def center(self):
        """
        Returns the closest floating-point approximation to the center
        of the interval.

        EXAMPLES::

            sage: CIF(RIF(1, 2), RIF(3, 4)).center()
            1.50000000000000 + 3.50000000000000*I
        """
        cdef complex_number.ComplexNumber center
        center = complex_number.ComplexNumber(self._parent._middle_field(), None)
        mpfi_mid(center.__re, self.__re)
        mpfi_mid(center.__im, self.__im)

        return center

    def __contains__(self, other):
        """
        Test whether ``other`` is totally contained in ``self``.

        EXAMPLES::

            sage: CIF(1, 1) in CIF(RIF(1, 2), RIF(1, 2))
            True
        """
        # This could be more efficient (and support more types for "other").
        return (other.real() in self.real()) and (other.imag() in self.imag())

    def contains_zero(self):
        """
        Returns ``True`` if ``self`` is an interval containing zero.

        EXAMPLES::

            sage: CIF(0).contains_zero()
            True
            sage: CIF(RIF(-1, 1), 1).contains_zero()
            False
        """
        return mpfi_has_zero(self.__re) and mpfi_has_zero(self.__im)

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add ``self`` and ``right``.

        EXAMPLES::

            sage: CIF(2,-3)._add_(CIF(1,-2))
            3 - 5*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()
        mpfi_add(x.__re, self.__re, (<ComplexIntervalFieldElement>right).__re)
        mpfi_add(x.__im, self.__im, (<ComplexIntervalFieldElement>right).__im)
        return x

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract ``self`` by ``right``.

        EXAMPLES::

            sage: CIF(2,-3)._sub_(CIF(1,-2))
            1 - 1*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()
        mpfi_sub(x.__re, self.__re, (<ComplexIntervalFieldElement>right).__re)
        mpfi_sub(x.__im, self.__im, (<ComplexIntervalFieldElement>right).__im)
        return x

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiply ``self`` and ``right``.

        EXAMPLES::

            sage: CIF(2,-3)._mul_(CIF(1,-2))
            -4 - 7*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()
        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)
        mpfi_mul(t0, self.__re, (<ComplexIntervalFieldElement>right).__re)
        mpfi_mul(t1, self.__im, (<ComplexIntervalFieldElement>right).__im)
        mpfi_sub(x.__re, t0, t1)
        mpfi_mul(t0, self.__re, (<ComplexIntervalFieldElement>right).__im)
        mpfi_mul(t1, self.__im, (<ComplexIntervalFieldElement>right).__re)
        mpfi_add(x.__im, t0, t1)
        mpfi_clear(t0)
        mpfi_clear(t1)
        return x

    def norm(self):
        """
        Returns the norm of this complex number.

        If `c = a + bi` is a complex number, then the norm of `c` is defined as
        the product of `c` and its complex conjugate:

        .. MATH::

            \text{norm}(c)
            =
            \text{norm}(a + bi)
            =
            c \cdot \overline{c}
            =
            a^2 + b^2.

        The norm of a complex number is different from its absolute value.
        The absolute value of a complex number is defined to be the square
        root of its norm. A typical use of the complex norm is in the
        integral domain `\ZZ[i]` of Gaussian integers, where the norm of
        each Gaussian integer `c = a + bi` is defined as its complex norm.

        .. SEEALSO::

            - :meth:`sage.rings.complex_double.ComplexDoubleElement.norm`

        EXAMPLES::

            sage: CIF(2, 1).norm()
            5
            sage: CIF(1, -2).norm()
            5
        """
        return self.norm_c()

    cdef real_mpfi.RealIntervalFieldElement norm_c(ComplexIntervalFieldElement self):
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)

        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)

        mpfi_sqr(t0, self.__re)
        mpfi_sqr(t1, self.__im)

        mpfi_add(x.value, t0, t1)

        mpfi_clear(t0)
        mpfi_clear(t1)
        return x

    cdef real_mpfi.RealIntervalFieldElement abs_c(ComplexIntervalFieldElement self):
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)

        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)

        mpfi_sqr(t0, self.__re)
        mpfi_sqr(t1, self.__im)

        mpfi_add(x.value, t0, t1)
        mpfi_sqrt(x.value, x.value)

        mpfi_clear(t0)
        mpfi_clear(t1)
        return x

    cpdef RingElement _div_(self, RingElement right):
        """
        Divide ``self`` by ``right``.

        EXAMPLES::

            sage: CIF(2,-3)._div_(CIF(1,-2))
            1.600000000000000? + 0.200000000000000?*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()
        cdef mpfi_t a, b, t0, t1, right_nm
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)
        mpfi_init2(a, self._prec)
        mpfi_init2(b, self._prec)
        mpfi_init2(right_nm, self._prec)

        mpfi_sqr(t0, (<ComplexIntervalFieldElement>right).__re)
        mpfi_sqr(t1, (<ComplexIntervalFieldElement>right).__im)
        mpfi_add(right_nm, t0, t1)

        mpfi_div(a, (<ComplexIntervalFieldElement>right).__re, right_nm)
        mpfi_div(b, (<ComplexIntervalFieldElement>right).__im, right_nm)

        ## Do this: x.__re =  a * self.__re + b * self.__im
        mpfi_mul(t0, a, self.__re)
        mpfi_mul(t1, b, self.__im)
        mpfi_add(x.__re, t0, t1)

        ## Do this: x.__im =  a * self.__im - b * self.__re
        mpfi_mul(t0, a, self.__im)
        mpfi_mul(t1, b, self.__re)
        mpfi_sub(x.__im, t0, t1)
        mpfi_clear(t0)
        mpfi_clear(t1)
        mpfi_clear(a)
        mpfi_clear(b)
        mpfi_clear(right_nm)
        return x

    def __rdiv__(self, left):
        """
        Divide ``left`` by ``self``.

        EXAMPLES::

            sage: CIF(2,-3).__rdiv__(CIF(1,-2))
            0.6153846153846154? - 0.0769230769230769?*I
        """
        return ComplexIntervalFieldElement(self._parent, left)/self

    def __pow__(self, right, modulus):
        r"""
        Compute `x^y`.

        If `y` is an integer, uses multiplication;
        otherwise, uses the standard definition `\exp(\log(x) \cdot y)`.

        .. WARNING::

            If the interval `x` crosses the negative real axis, then we use a
            non-standard definition of `\log()` (see the docstring for
            :meth:`argument()` for more details). This means that we will not
            select the principal value of the power, for part of the input
            interval (and that we violate the interval guarantees).

        EXAMPLES::

            sage: C.<i> = ComplexIntervalField(20)
            sage: a = i^2; a
            -1
            sage: a.parent()
            Complex Interval Field with 20 bits of precision
            sage: a = (1+i)^7; a
            8 - 8*I
            sage: (1+i)^(1+i)
            0.27396? + 0.58370?*I
            sage: a.parent()
            Complex Interval Field with 20 bits of precision
            sage: (2+i)^(-39)
            1.688?e-14 + 1.628?e-14*I

        If the interval crosses the negative real axis, then we don't use the
        standard branch cut (and we violate the interval guarantees)::

            sage: (CIF(-7, RIF(-1, 1)) ^ CIF(0.3)).str(style='brackets')
            '[0.99109735947126309 .. 1.1179269966896264] + [1.4042388462787560 .. 1.4984624123369835]*I'
            sage: CIF(-7, -1) ^ CIF(0.3)
            1.117926996689626? - 1.408500714575360?*I
        """
        if isinstance(right, (int, long, integer.Integer)):
            return sage.rings.ring_element.RingElement.__pow__(self, right)
        return (self.log() * self.parent()(right)).exp()

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: t = CIF((1, 1.1), 2.5); t
            1.1? + 2.5000000000000000?*I
            sage: magma(t) # optional - magma # indirect doctest
            1.05000000000000 + 2.50000000000000*$.1
            sage: t = ComplexIntervalField(100)((1, 4/3), 2.5); t
            2.? + 2.5000000000000000000000000000000?*I
            sage: magma(t) # optional - magma
            1.16666666666666666666666666670 + 2.50000000000000000000000000000*$.1
        """
        return "%s![%s, %s]" % (self.parent()._magma_init_(magma), self.center().real(), self.center().imag())

    def _interface_init_(self, I=None):
        """
        Raise a ``TypeError``.

        This function would return the string representation of ``self``
        that makes sense as a default representation of a complex
        interval in other computer algebra systems. But, most other
        computer algebra systems do not support interval arithmetic,
        so instead we just raise a ``TypeError``.

        Define the appropriate ``_cas_init_`` function if there is a
        computer algebra system you would like to support.

        EXAMPLES::

            sage: n = CIF(1.3939494594)
            sage: n._interface_init_()
            Traceback (most recent call last):
            ...
            TypeError

        Here a conversion to Maxima happens, which results in a ``TypeError``::

            sage: a = CIF(2.3)
            sage: maxima(a)
            Traceback (most recent call last):
            ...
            TypeError
        """
        raise TypeError


    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(CIF(RIF(e, pi), RIF(sqrt(2), sqrt(3))), verify=True)
            # Verified
            CIF(RIF(RR(2.7182818284590451), RR(3.1415926535897936)), RIF(RR(1.4142135623730949), RR(1.7320508075688774)))
            sage: sage_input(ComplexIntervalField(64)(2)^I, preparse=False, verify=True)
            # Verified
            RIF64 = RealIntervalField(64)
            RR64 = RealField(64)
            ComplexIntervalField(64)(RIF64(RR64('0.769238901363972126565'), RR64('0.769238901363972126619')), RIF64(RR64('0.638961276313634801076'), RR64('0.638961276313634801184')))
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: ComplexIntervalField(15)(3+I).log()._sage_input_(sib, False)
            {call: {call: {atomic:ComplexIntervalField}({atomic:15})}({call: {call: {atomic:RealIntervalField}({atomic:15})}({call: {call: {atomic:RealField}({atomic:15})}({atomic:1.15125})}, {call: {call: {atomic:RealField}({atomic:15})}({atomic:1.15137})})}, {call: {call: {atomic:RealIntervalField}({atomic:15})}({call: {call: {atomic:RealField}({atomic:15})}({atomic:0.321655})}, {call: {call: {atomic:RealField}({atomic:15})}({atomic:0.321777})})})}
        """
        # Interval printing could often be much prettier,
        # but I'm feeling lazy :)
        return sib(self.parent())(sib(self.real()), sib(self.imag()))

    def prec(self):
        """
        Return precision of this complex number.

        EXAMPLES::

            sage: i = ComplexIntervalField(2000).0
            sage: i.prec()
            2000
        """
        return self._parent.prec()

    def real(self):
        """
        Return real part of ``self``.

        EXAMPLES::

            sage: i = ComplexIntervalField(100).0
            sage: z = 2 + 3*i
            sage: x = z.real(); x
            2
            sage: x.parent()
            Real Interval Field with 100 bits of precision
        """
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)
        mpfi_set(x.value, self.__re)
        return x

    def imag(self):
        """
        Return imaginary part of ``self``.

        EXAMPLES::

            sage: i = ComplexIntervalField(100).0
            sage: z = 2 + 3*i
            sage: x = z.imag(); x
            3
            sage: x.parent()
            Real Interval Field with 100 bits of precision
        """
        cdef real_mpfi.RealIntervalFieldElement x
        x = real_mpfi.RealIntervalFieldElement(self._parent._real_field(), None)
        mpfi_set(x.value, self.__im)
        return x

    def __neg__(self):
        """
        Return the negation of ``self``.

        EXAMPLES::

            sage: CIF(1.5, 2.5).__neg__()
            -1.5000000000000000? - 2.5000000000000000?*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()
        mpfi_neg(x.__re, self.__re)
        mpfi_neg(x.__im, self.__im)
        return x

    def __pos__(self):
        """
        Return the "positive" of ``self``, which is just ``self``.

        EXAMPLES::

            sage: CIF(1.5, 2.5).__pos__()
            1.5000000000000000? + 2.5000000000000000?*I
        """
        return self

    def __abs__(self):
        """
        Return the absolute value of ``self``.

        EXAMPLES::

            sage: CIF(1.5, 2.5).__abs__()
            2.915475947422650?
        """
        return self.abs_c()

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        EXAMPLES::

            sage: I = CIF.0
            sage: a = ~(5+I) # indirect doctest
            sage: a * (5+I)
            1.000000000000000? + 0.?e-16*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()

        cdef mpfi_t t0, t1
        mpfi_init2(t0, self._prec)
        mpfi_init2(t1, self._prec)

        mpfi_sqr(t0, self.__re)
        mpfi_sqr(t1, self.__im)

        mpfi_add(t0, t0, t1)         # now t0 is the norm
        mpfi_div(x.__re, self.__re, t0)   #     x.__re = self.__re/norm

        mpfi_neg(t1, self.__im)
        mpfi_div(x.__im, t1, t0)  #     x.__im = -self.__im/norm

        mpfi_clear(t0)
        mpfi_clear(t1)

        return x

    def __int__(self):
        """
        Convert ``self`` to an ``int``.

        EXAMPLES::

            sage: int(CIF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex interval to int
        """
        raise TypeError, "can't convert complex interval to int"

    def __long__(self):
        """
        Convert ``self`` to a ``lon``.

        EXAMPLES::

            sage: long(CIF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex interval to long
        """
        raise TypeError, "can't convert complex interval to long"

    def __float__(self):
        """
        Convert ``self`` to a ``float``.

        EXAMPLES::

            sage: float(CIF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex interval to float
        """
        raise TypeError, "can't convert complex interval to float"

    def __complex__(self):
        """
        Convert ``self`` to a ``complex``.

        EXAMPLES::

            sage: complex(CIF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex interval to complex
        """
        raise TypeError, "can't convert complex interval to complex"

    def __nonzero__(self):
        """
        Return ``True`` if ``self`` is not known to be exactly zero.

        EXAMPLES::

            sage: CIF(RIF(0, 0), RIF(0, 0)).__nonzero__()
            False
            sage: bool(CIF(RIF(0, 0), RIF(0, 0)))
            False
            sage: CIF(RIF(1), RIF(0)).__nonzero__()
            True
            sage: CIF(RIF(0), RIF(1)).__nonzero__()
            True
            sage: CIF(RIF(1, 2), RIF(0)).__nonzero__()
            True
            sage: CIF(RIF(-1, 1), RIF(-1, 1)).__nonzero__()
            True
        """
        return self.real().__nonzero__() or self.imag().__nonzero__()

    def __richcmp__(left, right, int op):
        r"""
        As with the real interval fields this never returns false positives.
        Thus, `a == b` is ``True`` iff both `a` and `b` represent the same
        one-point interval. Likewise `a != b` is ``True`` iff `x != y` for all
        `x \in a, y \in b`.

        EXAMPLES::

            sage: CIF(0) == CIF(0)
            True
            sage: CIF(0) == CIF(1)
            False
            sage: CIF.gen() == CIF.gen()
            True
            sage: CIF(0) == CIF.gen()
            False
            sage: CIF(0) != CIF(1)
            True
            sage: -CIF(-3).sqrt() != CIF(-3).sqrt()
            True

        These intervals overlap, but contain unequal points::

            sage: CIF(3).sqrt() == CIF(3).sqrt()
            False
            sage: CIF(3).sqrt() != CIF(3).sqrt()
            False

        In the future, complex interval elements may be unordered,
        but or backwards compatibility we order them lexicographically::

            sage: CDF(-1) < -CDF.gen() < CDF.gen() < CDF(1)
            True
            sage: CDF(1) >= CDF(1) >= CDF.gen() >= CDF.gen() >= 0 >= -CDF.gen() >= CDF(-1)
            True
        """
        return (<Element>left)._richcmp(right, op)

    cdef _richcmp_c_impl(left, Element right, int op):
        cdef ComplexIntervalFieldElement lt, rt
        lt = left
        rt = right
        if op == 2: #==
            # intervals a == b iff a<=b and b <= a
            # (this gives a result with two comparisons, where the
            # obvious approach would use three)
            return mpfr_lessequal_p(&lt.__re.right, &rt.__re.left) \
                and mpfr_lessequal_p(&rt.__re.right, &lt.__re.left) \
                and mpfr_lessequal_p(&lt.__im.right, &rt.__im.left) \
                and mpfr_lessequal_p(&rt.__im.right, &lt.__im.left)
        elif op == 3: #!=
            return mpfr_less_p(&lt.__re.right, &rt.__re.left) \
                or mpfr_less_p(&rt.__re.right, &lt.__re.left) \
                or mpfr_less_p(&lt.__im.right, &rt.__im.left) \
                or mpfr_less_p(&rt.__im.right, &lt.__im.left)
        else:
            # Eventually we probably want to disable comparison of complex
            # intervals, just like python complexes will be unordered.
            ## raise TypeError, "no ordering relation is defined for complex numbers"
            diff = left - right
            real_diff = diff.real()
            imag_diff = diff.imag()
            if op == 0: #<
                return real_diff < 0 or (real_diff == 0 and imag_diff < 0)
            elif op == 1: #<=
                return real_diff < 0 or (real_diff == 0 and imag_diff <= 0)
            elif op == 4: #>
                return real_diff > 0 or (real_diff == 0 and imag_diff > 0)
            elif op == 5: #>=
                return real_diff > 0 or (real_diff == 0 and imag_diff >= 0)

    def __cmp__(left, right):
        """
        Intervals are compared lexicographically on the 4-tuple:
        ``(x.real().lower(), x.real().upper(),
        x.imag().lower(), x.imag().upper())``

        EXAMPLES::

            sage: a = CIF(RIF(0,1), RIF(0,1))
            sage: b = CIF(RIF(0,1), RIF(0,2))
            sage: c = CIF(RIF(0,2), RIF(0,2))
            sage: cmp(a, b)
            -1
            sage: cmp(b, c)
            -1
            sage: cmp(a, c)
            -1
            sage: cmp(a, a)
            0
            sage: cmp(b, a)
            1
        """
        return (<Element>left)._cmp(right)


    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        """
        Intervals are compared lexicographically on the 4-tuple:
        ``(x.real().lower(), x.real().upper(),
        x.imag().lower(), x.imag().upper())``

        TESTS::

            sage: tests = []
            sage: for rl in (0, 1):
            ...       for ru in (rl, rl + 1):
            ...           for il in (0, 1):
            ...               for iu in (il, il + 1):
            ...                   tests.append((CIF(RIF(rl, ru), RIF(il, iu)), (rl, ru, il, iu)))
            sage: for (i1, t1) in tests:
            ...       for (i2, t2) in tests:
            ...           assert(cmp(i1, i2) == cmp(t1, t2))
        """
        cdef int a, b
        a = mpfi_nan_p(left.__re)
        b = mpfi_nan_p((<ComplexIntervalFieldElement>right).__re)
        if a != b:
            return -1

        cdef int i
        i = mpfr_cmp(&left.__re.left, &(<ComplexIntervalFieldElement>right).__re.left)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&left.__re.right, &(<ComplexIntervalFieldElement>right).__re.right)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&left.__im.left, &(<ComplexIntervalFieldElement>right).__im.left)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(&left.__im.right, &(<ComplexIntervalFieldElement>right).__im.right)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        return 0

    ########################################################################
    # Transcendental (and other) functions
    ########################################################################

    def argument(self):
        r"""
        The argument (angle) of the complex number, normalized
        so that `-\pi < \theta.lower() \leq \pi`.

        We raise a ``ValueError`` if the interval strictly contains 0,
        or if the interval contains only 0.

        .. WARNING::

            We do not always use the standard branch cut for
            argument!  If the interval crosses the negative real axis,
            then the argument will be an interval whose lower bound is
            less than `\pi` and whose upper bound is more than `\pi`; in
            effect, we move the branch cut away from the interval.

        EXAMPLES::

            sage: i = CIF.0
            sage: (i^2).argument()
            3.141592653589794?
            sage: (1+i).argument()
            0.785398163397449?
            sage: i.argument()
            1.570796326794897?
            sage: (-i).argument()
            -1.570796326794897?
            sage: (RR('-0.001') - i).argument()
            -1.571796326461564?
            sage: CIF(2).argument()
            0
            sage: CIF(-2).argument()
            3.141592653589794?

        Here we see that if the interval crosses the negative real
        axis, then the argument can exceed `\pi`, and we
        we violate the standard interval guarantees in the process::

            sage: CIF(-2, RIF(-0.1, 0.1)).argument().str(style='brackets')
            '[3.0916342578678501 .. 3.1915510493117365]'
            sage: CIF(-2, -0.1).argument()
            -3.091634257867851?
        """
        if mpfi_has_zero(self.__re) and mpfi_has_zero(self.__im):

            if mpfi_is_zero(self.__re) and mpfi_is_zero(self.__im):
                raise ValueError, "Can't take the argument of complex zero"
            if not mpfi_is_nonpos(self.__re) and not mpfi_is_nonneg(self.__re) \
               and not mpfi_is_nonpos(self.__im) and not mpfi_is_nonneg(self.__im):
                raise ValueError, "Can't take the argument of interval strictly containing zero"

            # Now if we exclude zero from the interval, we know that the
            # argument of the remaining points is bounded.  Check which
            # axes the interval extends along (we can deduce information
            # about the quadrants from information about the axes).

            which_axes = [False, False, False, False]
            if not mpfi_is_nonpos(self.__re):
                which_axes[0] = True
            if not mpfi_is_nonpos(self.__im):
                which_axes[1] = True
            if not mpfi_is_nonneg(self.__re):
                which_axes[2] = True
            if not mpfi_is_nonneg(self.__im):
                which_axes[3] = True

            lower = None
            for i in range(-1, 3):
                if which_axes[i % 4] and not which_axes[(i - 1) % 4]:
                    if lower is not None:
                        raise ValueError, "Can't take the argument of line-segment interval strictly containing zero"
                    lower = i

            for i in range(lower, lower+4):
                if which_axes[i % 4] and not which_axes[(i + 1) % 4]:
                    upper = i
                    break

            fld = self.parent()._real_field()
            return fld.pi() * fld(lower, upper) * fld(0.5)

        else:

            # OK, we know that the interval is bounded away from zero
            # in either the real or the imaginary direction (or both).
            # We'll handle the "bounded away in the imaginary direction"
            # case first.

            fld = self.parent()._real_field()

            if mpfi_is_strictly_pos(self.__im):
                return (-self.real() / self.imag()).arctan() + fld.pi()/2
            if mpfi_is_strictly_neg(self.__im):
                return (-self.real() / self.imag()).arctan() - fld.pi()/2

            if mpfi_is_strictly_pos(self.__re):
                return (self.imag() / self.real()).arctan()

            # The only remaining case is that self.__re is strictly
            # negative and self.__im contains 0.  In that case, we
            # return an interval containing pi.

            return (self.imag() / self.real()).arctan() + fld.pi()

    def arg(self):
        """
        Same as :meth:`argument()`.

        EXAMPLES::

            sage: i = CIF.0
            sage: (i^2).arg()
            3.141592653589794?
        """
        return self.argument()

    def crosses_log_branch_cut(self):
        """
        Returns ``True`` if this interval crosses the standard branch cut
        for :meth:`log()` (and hence for exponentiation) and for argument.
        (Recall that this branch cut is infinitesimally below the
        negative portion of the real axis.)

        EXAMPLES::

            sage: z = CIF(1.5, 2.5) - CIF(0, 2.50000000000000001); z
            1.5000000000000000? + -1.?e-15*I
            sage: z.crosses_log_branch_cut()
            False
            sage: CIF(-2, RIF(-0.1, 0.1)).crosses_log_branch_cut()
            True
        """

        if mpfi_is_nonneg(self.__re):
            return False
        if mpfi_is_nonneg(self.__im):
            return False
        if mpfi_is_neg(self.__im):
            return False
        return True

    def conjugate(self):
        """
        Return the complex conjugate of this complex number.

        EXAMPLES::

            sage: i = CIF.0
            sage: (1+i).conjugate()
            1 - 1*I
        """
        cdef ComplexIntervalFieldElement x
        x = self._new()

        mpfi_set(x.__re, self.__re)
        mpfi_neg(x.__im, self.__im)
        return x

    def exp(self):
        r"""
        Compute `e^z` or `\exp(z)` where `z` is the complex number ``self``.

        EXAMPLES::

            sage: i = ComplexIntervalField(300).0
            sage: z = 1 + i
            sage: z.exp()
            1.46869393991588515713896759732660426132695673662900872279767567631093696585951213872272450? + 2.28735528717884239120817190670050180895558625666835568093865811410364716018934540926734485?*I
        """
        mag = self.real().exp()
        theta = self.imag()
        re = theta.cos() * mag
        im = theta.sin() * mag
        return ComplexIntervalFieldElement(self._parent, re, im)

    def log(self,base=None):
        """
        Complex logarithm of `z`.

        .. WARNING::

            This does always not use the standard branch cut for complex log!
            See the docstring for :meth:`argument()` to see what we do instead.

        EXAMPLES::

            sage: a = CIF(RIF(3, 4), RIF(13, 14))
            sage: a.log().str(style='brackets')
            '[2.5908917751460420 .. 2.6782931373360067] + [1.2722973952087170 .. 1.3597029935721503]*I'
            sage: a.log().exp().str(style='brackets')
            '[2.7954667135098274 .. 4.2819545928390213] + [12.751682453911920 .. 14.237018048974635]*I'
            sage: a in a.log().exp()
            True

        If the interval crosses the negative real axis, then we don't
        use the standard branch cut (and we violate the interval guarantees)::

            sage: CIF(-3, RIF(-1/4, 1/4)).log().str(style='brackets')
            '[1.0986122886681095 .. 1.1020725100903968] + [3.0584514217013518 .. 3.2247338854782349]*I'
            sage: CIF(-3, -1/4).log()
            1.102072510090397? - 3.058451421701352?*I

        Usually if an interval contains zero, we raise an exception::

            sage: CIF(RIF(-1,1),RIF(-1,1)).log()
            Traceback (most recent call last):
            ...
            ValueError: Can't take the argument of interval strictly containing zero

        But we allow the exact input zero::

            sage: CIF(0).log()
            [-infinity .. -infinity]

        If a base is passed from another function, we can accommodate this::

            sage: CIF(-1,1).log(2)
            0.500000000000000? + 3.399270106370396?*I
        """
        if self == 0:
            from real_mpfi import RIF
            return RIF(0).log()
        theta = self.argument()
        rho = abs(self)
        if base is None or base == 'e':
            return ComplexIntervalFieldElement(self._parent, rho.log(), theta)
        else:
            from real_mpfr import RealNumber, RealField
            return ComplexIntervalFieldElement(self._parent, rho.log()/RealNumber(RealField(self.prec()),base).log(), theta/RealNumber(RealField(self.prec()),base).log())

    def sqrt(self, bint all=False, **kwds):
        """
        The square root function.

        .. WARNING::

            We approximate the standard branch cut along the negative real
            axis, with ``sqrt(-r^2) = i*r`` for positive real ``r``; but if
            the interval crosses the negative real axis, we pick the root with
            positive imaginary component for the entire interval.

        INPUT:

        - ``all`` -- bool (default: ``False``); if ``True``, return a list
          of all square roots.

        EXAMPLES::

            sage: CIF(-1).sqrt()^2
            -1
            sage: sqrt(CIF(2))
            1.414213562373095?
            sage: sqrt(CIF(-1))
            1*I
            sage: sqrt(CIF(2-I))^2
            2.00000000000000? - 1.00000000000000?*I
            sage: CC(-2-I).sqrt()^2
            -2.00000000000000 - 1.00000000000000*I

        Here, we select a non-principal root for part of the interval, and
        violate the standard interval guarantees::

            sage: CIF(-5, RIF(-1, 1)).sqrt().str(style='brackets')
            '[-0.22250788030178321 .. 0.22250788030178296] + [2.2251857651053086 .. 2.2581008643532262]*I'
            sage: CIF(-5, -1).sqrt()
            0.222507880301783? - 2.247111425095870?*I
        """
        if self.is_zero():
            return [self] if all else self
        if mpfi_is_zero(self.__im) and not mpfi_has_zero(self.__re):
            if mpfr_sgn(&self.__re.left) > 0:
                x = ComplexIntervalFieldElement(self._parent, self.real().sqrt(), 0)
            else:
                x = ComplexIntervalFieldElement(self._parent, 0, (-self.real()).sqrt())
        else:
            theta = self.argument()/2
            rho = abs(self).sqrt()
            x = ComplexIntervalFieldElement(self._parent, rho*theta.cos(), rho*theta.sin())
        if all:
            return [x, -x]
        else:
            return x


    def is_square(self):
        r"""
        This function always returns ``True`` as `\CC` is algebraically closed.

        EXAMPLES::

            sage: CIF(2, 1).is_square()
            True
        """
        return True

#     def algdep(self, n, **kwds):
#         """
#         Returns a polynomial of degree at most $n$ which is approximately
#         satisfied by this complex number.  Note that the returned polynomial
#         need not be irreducible, and indeed usually won't be if $z$ is a good
#         approximation to an algebraic number of degree less than $n$.

#         ALGORITHM: Uses the PARI C-library algdep command.

#         INPUT: Type algdep? at the top level prompt. All additional
#         parameters are passed onto the top-level algdep command.

#         EXAMPLES::
#
#             sage: C = ComplexIntervalField()
#             sage: z = (1/2)*(1 + sqrt(3.0) *C.0); z
#             0.500000000000000 + 0.866025403784439*I
#             sage: p = z.algdep(5); p
#             x^5 + x^2
#             sage: p.factor()
#             (x + 1) * x^2 * (x^2 - x + 1)
#             sage: z^2 - z + 1
#             0.000000000000000111022302462516
#         """
#         import sage.rings.arith
#         return sage.rings.arith.algdep(self,n, **kwds)

#     def algebraic_dependancy( self, n ):
#         return self.algdep( n )

def make_ComplexIntervalFieldElement0( fld, re, im ):
    """
    Construct a :class:`ComplexIntervalFieldElement` for pickling.

    TESTS::

        sage: a = CIF(1 + I)
        sage: loads(dumps(a)) == a # indirect doctest
        True
    """
    x = ComplexIntervalFieldElement( fld, re, im )
    return x



def create_ComplexIntervalFieldElement(s_real, s_imag=None, int pad=0, min_prec=53):
    r"""
    Return the complex number defined by the strings ``s_real`` and ``s_imag``
    as an element of ``ComplexIntervalField(prec=n)``, where `n` potentially
    has slightly more (controlled by pad) bits than given by `s`.

    INPUT:

    - ``s_real`` -- a string that defines a real number (or something whose
      string representation defines a number)

    - ``s_imag`` -- a string that defines a real number (or something whose
      string representation defines a number)

    - ``pad`` -- an integer at least 0.

    - ``min_prec`` -- number will have at least this many bits of precision,
      no matter what.

    EXAMPLES::

        sage: ComplexIntervalFieldElement('2.3')
        2.300000000000000?
        sage: ComplexIntervalFieldElement('2.3','1.1')
        2.300000000000000? + 1.100000000000000?*I
        sage: ComplexIntervalFieldElement(10)
        10
        sage: ComplexIntervalFieldElement(10,10)
        10 + 10*I
        sage: ComplexIntervalFieldElement(1.000000000000000000000000000,2)
        1 + 2*I
        sage: ComplexIntervalFieldElement(1,2.000000000000000000000)
        1 + 2*I
        sage: ComplexIntervalFieldElement(1.234567890123456789012345, 5.4321098654321987654321)
        1.234567890123456789012350? + 5.432109865432198765432000?*I

    TESTS:

    Make sure we've rounded up ``log(10,2)`` enough to guarantee
    sufficient precision (:trac:`10164`).  This is a little tricky
    because at the time of writing, we don't support intervals long
    enough to trip the error.  However, at least we can make sure that
    we either do it correctly or fail noisily::

        sage: c_CIFE = sage.rings.complex_interval.create_ComplexIntervalFieldElement
        sage: for kp in range(2,6):
        ...       s = '1.' + '0'*10**kp + '1'
        ...       try:
        ...           assert c_CIFE(s,0).real()-1 != 0
        ...           assert c_CIFE(0,s).imag()-1 != 0
        ...       except TypeError:
        ...           pass

    """
    if s_imag is None:
        s_imag = 0

    if not isinstance(s_real, str):
        s_real = str(s_real).strip()
    if not isinstance(s_imag, str):
        s_imag = str(s_imag).strip()
    #if base == 10:
    bits = max(int(LOG_TEN_TWO_PLUS_EPSILON*len(s_real)),
               int(LOG_TEN_TWO_PLUS_EPSILON*len(s_imag)))
    #else:
    #    bits = max(int(math.log(base,2)*len(s_imag)),int(math.log(base,2)*len(s_imag)))

    C = complex_interval_field.ComplexIntervalField(prec=max(bits+pad, min_prec))
    return ComplexIntervalFieldElement(C, s_real, s_imag)
