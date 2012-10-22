"""
p-Adic Generic

A generic superclass for all p-adic parents.

AUTHORS:

- David Roe

- Genya Zaytman: documentation

- David Harvey: doctests

- Julian Rueth (2012-09-05): added (x)gcd for polynomials

"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                     2007 William Stein <wstein@gmail.com>
#                     2012 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import infinity
from local_generic import LocalGeneric
from sage.rings.ring import PrincipalIdealDomain
from sage.rings.integer import Integer
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.padics.padic_printing import pAdicPrinter
from sage.misc.cachefunc import cached_method

class pAdicGeneric(PrincipalIdealDomain, LocalGeneric):
    def __init__(self, base, p, prec, print_mode, names, element_class):
        """
        Initialization.

        INPUTS::

            - base -- Base ring.
            - p -- prime
            - print_mode -- dictionary of print options
            - names -- how to print the uniformizer
            - element_class -- the class for elements of this ring

        EXAMPLES::

            sage: R = Zp(17) #indirect doctest
        """
        LocalGeneric.__init__(self, base, prec, names, element_class)
        self._printer = pAdicPrinter(self, print_mode)

    def _modified_print_mode(self, print_mode):
        """
        Returns a dictionary of print options, starting with self's
        print options but modified by the options in the dictionary
        print_mode.

        INPUTS::

            - print_mode -- dictionary with keys in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']

        EXAMPLES::

            sage: R = Zp(5)
            sage: R._modified_print_mode({'mode': 'bars'})['ram_name']
            '5'
        """
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
            if not print_mode.has_key(option):
                print_mode[option] = self._printer.dict()[option]
        return print_mode

    def ngens(self):
        """
        Returns the number of generators of self.

        We conventionally define this as 1: for base rings, we take a
        uniformizer as the generator; for extension rings, we take a
        root of the minimal polynomial defining the extension.

        EXAMPLES::

            sage: Zp(5).ngens()
            1
            sage: Zq(25,names='a').ngens()
            1
        """
        return 1

    def gens(self):
        """
        Returns a list of generators.

        EXAMPLES::

            sage: R = Zp(5); R.gens()
            [5 + O(5^21)]
            sage: Zq(25,names='a').gens()
            [a + O(5^20)]
            sage: S.<x> = ZZ[]; f = x^5 + 25*x -5; W.<w> = R.ext(f); W.gens()
            [w + O(w^101)]
        """
        return [self.gen()]

    def _gcd_univariate_polynomial_fixed(self, f, g):
        """
        This method computes the greatest common divisor of the polynomials
        ``f`` and ``g``, treating them as if they were defined over a
        fixed-modulus ring.

        This is a helper method for
        :meth:`sage.rings.padics.generic_nodes.pAdicFixedModRingGeneric._gcd_univariate_polynomial`.

        This method handles some special cases and then calls
        :meth:`__xgcd_univariate_polynomial_fixed` which should be consulted
        for further details.

        INPUT:

            - ``f``, ``g`` -- two polynomials defined over ``self``.

        OUTPUT:

        A polynomial defined over ``self``, and the precision to which this
        result is accurate.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        TESTS:

        We check that this method handles correctly the special cases that it
        is supposed to take from :meth:`__xgcd_univariate_polynomial_fixed`.

        ``f`` or ``g`` are zero::

            sage: S = ZpFM(3,3)
            sage: R.<t> = S[]
            sage: f = 6*t^2 + 3*t + 3
            sage: S._gcd_univariate_polynomial_fixed(f,R.zero())
            ((3 + O(3^3))*t^2 + (2*3 + 3^2 + O(3^3))*t + (2*3 + 3^2 + O(3^3)), 3)
            sage: S._gcd_univariate_polynomial_fixed(R.zero(),f)
            ((3 + O(3^3))*t^2 + (2*3 + 3^2 + O(3^3))*t + (2*3 + 3^2 + O(3^3)), 3)

        ``f`` and ``g`` are constant::

            sage: S._gcd_univariate_polynomial_fixed(R(1),R(2))
            ((1 + O(3^3)), 3)
            sage: S._gcd_univariate_polynomial_fixed(R(3),R(9))
            ((3 + O(3^3)), 3)
            sage: S._gcd_univariate_polynomial_fixed(R(0),R(3))
            ((3 + O(3^3)), 3)

        ``f`` is constant, ``g`` is not::

            sage: S._gcd_univariate_polynomial_fixed(R(1),f)
            ((1 + O(3^3)), 3)
            sage: S._gcd_univariate_polynomial_fixed(R(3),f)
            ((3 + O(3^3)), 3)
            sage: S._gcd_univariate_polynomial_fixed(R(9),f)
            ((3 + O(3^3)), 3)

        """
        if f.parent() is not g.parent():
            raise ValueError("arguments must have the same parent")
        # we refuse to handle polynomials with leading zero coefficients - it
        # is not clear what a gcd should be in such a case. (this can only
        # happen for polynomials over capped relative rings.)
        if len(list(f)) != f.degree() + 1 or len(list(g)) != g.degree() + 1:
            raise ValueError("Cannot compute the gcd of polynomials with leading zero coefficients.")
        # the caller must have taken care of coefficients with negative
        # valuation
        if any([c.valuation()<0 for c in list(f)]+[c.valuation()<0 for c in list(g)]):
            raise ValueError("Cannot compute the gcd of polynomials with coefficients with negative valuation.")

        precision_cap = f.parent().base_ring().precision_cap()

        # constant polynomials are treated by the base ring
        if f.degree()<=0 and g.degree()<=0:
            return f.parent()([f[0].gcd(g[0])]), precision_cap

        # if f is zero, we can assume that is an exact zero, and get g as their
        # gcd (which is a result of maximal degree).
        if f.is_zero():
            return g*g.leading_coefficient().unit_part().inverse_of_unit(), precision_cap
        # same for g
        if g.is_zero():
            return f*f.leading_coefficient().unit_part().inverse_of_unit(), precision_cap

        # if exactly one of the polynomials is constant but non-zero, then the
        # gcd is simply the gcd of their contents.
        if f.degree()==0 or g.degree()==0:
            if f.parent().base_ring().is_field():
                return f.parent().one(), precision_cap
            return f.parent()([f.content().gen().gcd(g.content().gen())]), precision_cap

        return tuple(self.__xgcd_univariate_polynomial_fixed(f,g)[0:2])

    def _xgcd_univariate_polynomial_fixed(self, f, g):
        """
        This method computes the extended greatest common divisor of the
        polynomials ``f`` and ``g``, treating them as if they were defined over
        a fixed-modulus ring.

        This is a helper method for
        :meth:`sage.rings.padics.generic_nodes.pAdicFixedModRingGeneric._xgcd_univariate_polynomial`.

        This method handles some special cases and then calls
        :meth:`__xgcd_univariate_polynomial_fixed` which should be consulted
        for further details.

        INPUT:

            - ``f``, ``g`` -- two polynomials defined over ``self``.

        OUTPUT:

            A polynomial defined over ``self``, and the precision to which this
            result is accurate.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        TESTS:

        We check that this method handles correctly the special cases that it
        is supposed to take from :meth:`__xgcd_univariate_polynomial_fixed`.

        ``f`` or ``g`` are zero::

            sage: S = ZpFM(3,2)
            sage: R.<t> = S[]
            sage: f = 6*t^2 + 3*t + 3
            sage: S._xgcd_univariate_polynomial_fixed(f,R.zero())
            ((3 + O(3^2))*t^2 + (2*3 + O(3^2))*t + (2*3 + O(3^2)), 2, 2 + 3 + O(3^2), 0)
            sage: S._xgcd_univariate_polynomial_fixed(R.zero(),f)
            ((3 + O(3^2))*t^2 + (2*3 + O(3^2))*t + (2*3 + O(3^2)), 2, 0, 2 + 3 + O(3^2))

        ``f`` and ``g`` are constant::

            sage: S._xgcd_univariate_polynomial_fixed(R(1),R(2))
            ((1 + O(3^2)), 2, 0, (2 + 3 + O(3^2)))
            sage: S._xgcd_univariate_polynomial_fixed(R(3),R(9))
            ((3 + O(3^2)), 2, (1 + O(3^2)), 0)
            sage: S._xgcd_univariate_polynomial_fixed(R(0),R(3))
            ((3 + O(3^2)), 2, 0, (1 + O(3^2)))

        ``f`` is constant and ``g`` is not::

            sage: S._xgcd_univariate_polynomial_fixed(R(1),f)
            ((1 + O(3^2)), 2, 1 + O(3^2), 0)
            sage: S._xgcd_univariate_polynomial_fixed(R(3),f)
            ((3 + O(3^2)), 2, 1 + O(3^2), 0)
            sage: S._xgcd_univariate_polynomial_fixed(R(9),f)
            ((3 + O(3^2))*t^2 + (2*3 + O(3^2))*t + (2*3 + O(3^2)), 2, 0, 2 + 3 + O(3^2))

        """
        if f.parent() is not g.parent():
            raise ValueError("arguments must have the same parent")
        # we refuse to handle polynomials with leading zero coefficients - it
        # is not clear what a gcd should be in such a case. (this can only
        # happen for polynomials over capped relative rings.)
        if len(list(f)) != f.degree() + 1 or len(list(g)) != g.degree() + 1:
            raise ValueError("Cannot compute the gcd of polynomials with leading zero coefficients.")
        # the caller must have taken care of coefficients with negative
        # valuation
        if any([c.valuation()<0 for c in list(f)]+[c.valuation()<0 for c in list(g)]):
            raise ValueError("Cannot compute the gcd of polynomials with coefficients with negative valuation.")

        precision_cap = f.parent().base_ring().precision_cap()

        # constant polynomials are treated by the base ring
        if f.degree()<=0 and g.degree()<=0:
            xgcd = f[0].xgcd(g[0])
            return f.parent()(xgcd[0]), precision_cap, f.parent()(xgcd[1]), f.parent()(xgcd[2])

        # if f is zero, we can assume that is an exact zero, and get g as their
        # xgcd (which is a result of maximal degree).
        if f.is_zero():
            factor = g.leading_coefficient().unit_part().inverse_of_unit()
            return g*factor, precision_cap, f.parent().zero(), factor
        # same for g
        if g.is_zero():
            factor = f.leading_coefficient().unit_part().inverse_of_unit()
            return f*factor, precision_cap, factor, g.parent().zero()

        # if exactly one of the polynomials is constant but non-zero, then that
        # is the xgcd
        if f.degree()==0:
            factor = f.leading_coefficient().unit_part().inverse_of_unit()
            return f*factor, precision_cap, factor, g.parent().zero()
        if g.degree()==0:
            factor = g.leading_coefficient().unit_part().inverse_of_unit()
            return g*factor, precision_cap, f.parent().zero(), factor

        xgcd = self.__xgcd_univariate_polynomial_fixed(f,g)
        return tuple(self.__xgcd_univariate_polynomial_fixed(f,g)[2:])

    def __xgcd_shift(self, c, shift, precision_cap):
        """
        Helper method for :meth:`__xgcd_univariate_polynomial_fixed` which
        performs a left shift of ``shift`` digits on ``c``. Reports a
        ``PrecisionError`` if the shift exceeds the ``precision_cap``.

        INPUT:

        - ``c`` -- an element of ``self``

        - ``shift`` -- an integer

        - ``precision_cap`` -- a positive integer

        OUTPUT:

        Returns ``c<<shift`` or raises a ``ValueError`` if ``-shift`` exceeds
        ``precision_cap``.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        TESTS::

            sage: R = Zp(3,20,'fixed-mod')
            sage: R._pAdicGeneric__xgcd_shift(R(3^19), -19, 20)
            1 + O(3^20)
            sage: R._pAdicGeneric__xgcd_shift(R(3),19, 20)
            O(3^20)
            sage: R._pAdicGeneric__xgcd_shift(R(0,20), -20, 20)
            O(3^20)
            sage: R._pAdicGeneric__xgcd_shift(R(0,20), -21, 20)
            Traceback (most recent call last):
            ...
            PrecisionError: Insufficient precision to compute greatest common divisor - computation requires at least 1 additional digit(s).

        """
        if -shift > precision_cap:
            raise PrecisionError("Insufficient precision to compute greatest common divisor - computation requires at least %s additional digit(s)."%(-shift-precision_cap))
        assert self.__xgcd_valuation(c,precision_cap) >= -shift
        return c<<shift

    def __xgcd_valuation(self, x, precision_cap):
        """
        Helper method for :meth:`__xgcd_univariate_polynomial_fixed` which
        returns the valuation fo ``x`` taking into account the current precision ``precision_cap``.

        INPUT:

        - ``x`` -- an element of ``self``

        - ``precision_cap`` -- a positive integer

        OUTPUT:

        Returns the minimum of ``x.valuation()`` and ``precision_cap``.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        TESTS::

            sage: R = Zp(3,20,'fixed-mod')
            sage: R._pAdicGeneric__xgcd_valuation(R(3^19), 20)
            19
            sage: R._pAdicGeneric__xgcd_valuation(R(3^19), 2)
            2

        """
        ret = x.valuation()
        if ret is infinity:
            return precision_cap
        if ret > precision_cap:
            return precision_cap
        assert ret >= 0
        return ret

    def __xgcd_change_precision(self, old_precision_cap, new_precision_cap):
        """
        Helper method for :meth:`__xgcd_univariate_polynomial_fixed` to report
        a ``PrecisionError`` if a precision change happens to a negative
        precision.

        INPUT:

        - ``old_precision_cap`` -- a positive integer

        - ``new_precision_cap`` -- an integer less or equal to ``old_precision_cap``

        OUTPUT:

        Returns ``new_precision_cap`` or throws a ``PrecisionError`` if
        ``new_precision_cap`` is not positive.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        TESTS::

            sage: R = Zp(3)
            sage: R._pAdicGeneric__xgcd_change_precision(2,1)
            1
            sage: R._pAdicGeneric__xgcd_change_precision(2,-1)
            Traceback (most recent call last):
            ...
            PrecisionError: Insufficient precision to compute greatest common divisor - computation requires at least 2 additional digit(s).

        """
        if new_precision_cap <= 0:
            raise PrecisionError("Insufficient precision to compute greatest common divisor - computation requires at least %s additional digit(s)."%(-new_precision_cap+1))
        assert new_precision_cap <= old_precision_cap
        return new_precision_cap

    def __xgcd_univariate_polynomial_fixed(self, f, g, use_substitution=None):
        """
        Helper method for :meth:`_gcd_univariate_polynomial_fixed` and
        :meth:`_xgcd_univariate_polynomial_fixed` to compute the gcd and xgcd
        of ``f`` and ``g``, treating them as if they were defined over a
        fixed-modulus ring.

        The greatest common divisor of input polynomials which are only given
        up to some precision is not unique. Such polynomials represent all
        polynomials over the p-adic field which they approximate. As an
        example, over the 3-adics, the polynomial `f=g=(1 + O(3))t + O(3)`
        represents the polynomial t as well as t + 3. One could claim that the
        greatest common divisor of `f` and `g` could be any of t, t + 3, and 1.
        This method returns, however, amongst all approximations of the
        greatest common divisor of two polynomials represented by ``f`` and
        ``g``, one with maximal degree. In the example, it returns `(1 + O(3))t
        + O(3)` which has the maximal possible degree 1.

        INPUT:

        - ``f``, ``g`` -- univariate polynomials defined over ``self`` whose
          coefficients have nonnegative valuation

        - ``use_substition`` -- (default: None) boolean value which controls
          whether or not to apply an initial substition; if ``None`` then the
          implementation tries to pick the better choice in the given case.

        OUTPUT:

        A tuple ``gcd, prec_gcd, xgcd, prec_xgcd, xgcd_u, xgcd_v``:

        ``gcd`` is a greatest common divisor of ``f`` and ``g`` which is
        accurate to ``prec_gcd`` digits, i.e., ``gcd`` divides both ``f`` and
        ``g`` if they are reduced to precision ``prec``, and no polynomial of
        higher degree divides ``f`` and ``g``.

        ``xgcd`` is a polynomial which is a multiple of ``gcd`` by some element
        of nonnegative valuation. The polynomials ``xgcd_u`` and ``xgcd_v``
        satisfy a Bezout identity ``xgcd = xgcd_u * f + xgcd_v * g`` when
        reducing all polynomials to precision ``prec_xgcd``.

        If ``self`` is a field, then ``xgcd`` and ``gcd`` are equal, and
        ``prec_gcd`` and ``prec_xgcd`` are equal.

        ALGORITHM:

            The implemented algorithm echelonizes the sylvester matrix of the
            polynomials. The coefficients of the greatest common divisor are
            given by the last non-zero row of the echelonized matrix.

            Since the input is only given with finite precision, this
            echelonization essentially takes place over a finite ring
            `\mathbb{Z}/n\mathbb{Z}`. Special care has to be taken if during
            the echelonization a pivot of positive valuation has to be used.
            For example, mod 9, `3\cdot 4=3\cdot 1` so if the pivot is 3, then
            the echelonization process can produce two different matrices: one
            that is obtained by eliminating with the row as it is, and another
            one obtained by first multiplying the row by 4.

            The implementation treats the coefficients of the input polynomials
            as fixed-mod padic numbers. It keeps track of their precision in a
            variable ``precision_cap`` and does not rely on the precision
            handling they might provide themselves.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        EXAMPLES:

        The example mentioned above::

            sage: S = ZpFM(3,1)
            sage: R.<t> = S[]
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(t,t)
            ((1 + O(3))*t, 1, (1 + O(3))*t, 1, (1 + O(3)), 0)

        Gcd computation might result in a loss of precision, especially if the
        leading coefficients of the input polynomials have positive valuation::

            sage: S = ZpFM(3,8)
            sage: R.<t> = S[]
            sage: f = 3553*t^3 + 5039*t^2 + 5347*t + 1859; f
            (1 + 2*3 + 3^2 + 2*3^3 + 3^4 + 2*3^5 + 3^6 + 3^7 + O(3^8))*t^3 + (2 + 2*3 + 3^2 + 2*3^4 + 2*3^5 + 2*3^7 + O(3^8))*t^2 + (1 + 3^5 + 3^6 + 2*3^7 + O(3^8))*t + (2 + 3 + 2*3^2 + 2*3^3 + 3^4 + 3^5 + 2*3^6 + O(3^8))
            sage: g = 1068*t^3 + 1190*t^2 + 6090*t + 5446; g
            (2*3 + 3^2 + 3^4 + 3^5 + 3^6 + O(3^8))*t^3 + (2 + 2*3^3 + 2*3^4 + 3^5 + 3^6 + O(3^8))*t^2 + (2*3 + 3^2 + 3^5 + 2*3^6 + 2*3^7 + O(3^8))*t + (1 + 2*3^2 + 3^4 + 3^5 + 3^6 + 2*3^7 + O(3^8))
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,g)
            ((1 + O(3^8)), 8, (3^2 + O(3^8)), 8, (3 + 2*3^2 + 2*3^3 + 3^4 + 2*3^5 + 3^7 + O(3^8))*t^2 + (1 + 3 + 3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^7 + O(3^8))*t + (1 + 2*3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 2*3^7 + O(3^8)), (1 + 3^3 + 2*3^5 + 2*3^6 + 3^7 + O(3^8))*t^2 + (3^3 + 3^6 + O(3^8))*t + (1 + 3 + 3^2 + 3^3 + 2*3^5 + 2*3^7 + O(3^8)))

        In the following series of examples, the degree of the result decreases
        with increasing precision::

            sage: R.<t> = ZZ[]
            sage: f = (t + 8)^2 * (t - 19)
            sage: S = ZpFM(3,2)
            sage: g = f.change_ring(S)
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(g,g.derivative())
            ((1 + O(3^2))*t^2 + (1 + O(3^2))*t + (1 + O(3^2)), 1, (3 + O(3^2))*t^2 + (3 + O(3^2))*t + (3 + O(3^2)), 2, (2*3 + O(3^2))*t + (2*3 + O(3^2)), (1 + O(3^2))*t^2)

        The meaning of this result is the following: ``g`` divides ``f`` and
        ``f.derivative()``, at least if we reduce them mod `p`. At the same
        time there is no polynomial of higher degree dividing both unreduced
        polynomials ``f`` and ``f.derivative()``.

        This comes as no surprise since ``f`` is ``(t - 1)^3`` mod `p`.
        Increasing the precision of the input polynomial, ``f`` is not of this
        form anymore::

            sage: S = ZpFM(3,15)
            sage: g = f.change_ring(S)
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(g,g.derivative())
            ((1 + O(3^15))*t + (2 + 2*3 + O(3^15)), 10, (3^5 + O(3^15))*t + (2*3^5 + 2*3^6 + O(3^15)), 15, (3 + 3^2 + 3^6 + 2*3^7 + 2*3^9 + 2*3^10 + 2*3^11 + 3^13 + O(3^15))*t + (3^5 + 2*3^7 + 2*3^8 + 2*3^9 + 3^11 + 3^12 + 3^13 + O(3^15)), (2 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 3^5 + 2*3^7 + 2*3^11 + 3^12 + 2*3^13 + 3^14 + O(3^15))*t^2 + (1 + 3 + 2*3^4 + 3^7 + 2*3^8 + 3^9 + 3^10 + 2*3^11 + 2*3^12 + 2*3^13 + O(3^15))*t)

        But, by construction, the gcd is of course at least of degree one,
        independent of the precision::

            sage: S = ZpFM(3,30)
            sage: g = f.change_ring(S)
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(g,g.derivative())
            ((1 + O(3^30))*t + (2 + 2*3 + O(3^30)), 25, (3^5 + O(3^30))*t + (2*3^5 + 2*3^6 + O(3^30)), 30, (3 + 3^2 + 3^6 + 2*3^7 + 2*3^9 + 2*3^10 + 2*3^13 + 3^15 + 3^16 + 3^17 + 3^18 + 2*3^19 + 2*3^20 + 3^21 + 3^22 + 3^23 + 2*3^24 + 2*3^28 + O(3^30))*t + (3^5 + 2*3^7 + 2*3^8 + 2*3^9 + 3^11 + 3^12 + 3^13 + 3^15 + 2*3^16 + 3^17 + 3^18 + 3^23 + 2*3^25 + 2*3^26 + 2*3^27 + 3^29 + O(3^30)), (2 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 3^5 + 2*3^7 + 2*3^10 + 2*3^11 + 2*3^13 + 3^14 + 3^15 + 3^16 + 3^17 + 3^20 + 3^21 + 3^22 + 2*3^24 + 2*3^25 + 2*3^26 + 2*3^28 + O(3^30))*t^2 + (1 + 3 + 2*3^4 + 3^7 + 2*3^8 + 3^9 + 2*3^10 + 3^11 + 2*3^15 + 2*3^16 + 2*3^17 + 3^18 + 2*3^19 + 3^20 + 3^21 + 2*3^23 + 3^24 + 3^27 + 2*3^28 + 2*3^29 + O(3^30))*t)

        In some cases the input precision might not be sufficient to provide a
        result::

            sage: R.<t> = ZZ[]
            sage: f = 18*t^2 + 20*t + 6
            sage: g = 9*t^4 + 7*t^3 + 16*t^2 + 23*t + 6
            sage: S = ZpFM(3,3)
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f.change_ring(S), g.change_ring(S))
            Traceback (most recent call last):
            ...
            PrecisionError: Insufficient precision to compute greatest common divisor - computation requires at least 2 additional digit(s).

        Increasing the precision helps in such cases::

            sage: S = ZpFM(3,6)
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f.change_ring(S), g.change_ring(S))
            ((1 + O(3^6))*t + (3 + 2*3^3 + 3^5 + O(3^6)), 3, (1 + O(3^6))*t + (3 + 2*3^3 + 3^5 + O(3^6)), 3, (3^5 + O(3^6))*t^3 + (2*3^4 + 2*3^5 + O(3^6))*t^2 + (3^2 + 3^3 + 2*3^4 + 2*3^5 + O(3^6))*t + (2 + 3 + 2*3^2 + 3^3 + O(3^6)), (3^5 + O(3^6))*t)

        Again, a further increase of the precision can decrease the degree of the greatest common divisor::

            sage: S = ZpFM(3,7)
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f.change_ring(S), g.change_ring(S))
            ((1 + O(3^7)), 7, (3^3 + O(3^7)), 4, (3^2 + 2*3^3 + 2*3^5 + O(3^7))*t^3 + (1 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^7))*t^2 + (1 + 2*3^2 + 3^3 + 2*3^5 + 2*3^6 + O(3^7))*t + (2 + 3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^7)), (3^2 + 3^3 + 3^4 + 3^5 + 3^6 + O(3^7))*t + (1 + 3 + 2*3^2 + 3^3 + 3^4 + 3^5 + 2*3^6 + O(3^7)))

        A more involved example::

            sage: S = ZpFM(3,20)
            sage: R.<t> = S[]
            sage: a0,a1,a2 = 12345678,90123456,78901234
            sage: f = (t - a0) * (t - a1)^2 * (t - a2)^3
            sage: g,prec,_,_,_,_ = S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,f.derivative())
            sage: h = ((t - a1) * (t - a2)^2).map_coefficients(lambda c:c.add_bigoh(prec))
            sage: g == h
            True

        TESTS:

        The following pairs of polynomials are by construction prime mod p.
        Therefore, their greatest common divisor must be 1::

            sage: S = ZpFM(3,30)
            sage: R.<t> = S[]
            sage: n = 30
            sage: factors_1 = [ t - a*3 + 1 for a in range(n*3) ]
            sage: factors_2 = [ t - a*3 + 2 for a in range(n*3) ]
            sage: polynomials_1 = [ prod(factors_1[i:i+3]) for i in range(n) ]
            sage: polynomials_2 = [ prod(factors_2[i:i+3]) for i in range(n) ]
            sage: from itertools import product
            sage: polynomials = list(product(polynomials_1, polynomials_2))

            sage: all([ S._xgcd_univariate_polynomial_fixed(p1,p2)[0].is_one() for (p1,p2) in polynomials ]) # long time
            True

            sage: for p1,p2 in polynomials: # long time
            ...     g,prec,_,_,_,_ = S._xgcd_univariate_polynomial_fixed(p1,p2*p1)
            ...     assert g == p1.map_coefficients(lambda c:c.add_bigoh(prec))

        The same tests for polynomials whose leading coefficient has positive
        valuation::

            sage: S = ZpFM(3,30)
            sage: R.<t> = S[]
            sage: factors_1 = [ 3*t - a*9 + 1 for a in range(n*3) ]
            sage: factors_2 = [ 3*t - a*9 + 2 for a in range(n*3) ]
            sage: polynomials_1 = [ prod(factors_1[i:i+3]) for i in range(n) ]
            sage: polynomials_2 = [ prod(factors_2[i:i+3]) for i in range(n) ]
            sage: from itertools import product
            sage: polynomials = list(product(polynomials_1, polynomials_2))

            sage: all([ S._xgcd_univariate_polynomial_fixed(p1,p2)[0].is_one() for (p1,p2) in polynomials ]) # long time
            True

            sage: for p1,p2 in polynomials: # long time
            ...     g,prec,_,_,_,_ = S._xgcd_univariate_polynomial_fixed(p1,p2*p1)
            ...     assert g == p1.map_coefficients(lambda c:c.add_bigoh(prec))
            sage: all([ S._xgcd_univariate_polynomial_fixed(p1,p2)[0].is_one() for (p1,p2) in polynomials ]) # long time
            True

        Some cases that came up during development of the method::

            sage: S = ZpFM(3,8)
            sage: R.<t> = S[]
            sage: f = 2151*t^2 + 723*t + 1296
            sage: g = 6273*t^2 + 2934*t + 4311
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,g)
            ((3 + O(3^8)), 8, (3^2 + O(3^8)), 7, (3 + 2*3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^8))*t + (2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 3^6 + 2*3^7 + O(3^8)), (3 + 3^6 + O(3^8))*t + (2 + 3 + 2*3^2 + 2*3^3 + 2*3^6 + 3^7 + O(3^8)))

        ::

            sage: S = ZpFM(3,8)
            sage: R.<t> = S[]
            sage: f = 1575*t^2 + 4436*t + 5750
            sage: g = 504*t^2 + 2369*t + 6494
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,g)
            ((1 + O(3^8)), 8, (3 + O(3^8)), 6, (2*3^2 + 3^3 + 2*3^5 + 3^6 + O(3^8))*t + (2 + 3 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + O(3^8)), (2*3^2 + 2*3^3 + 2*3^6 + 2*3^7 + O(3^8))*t + (1 + 2*3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 3^6 + 3^7 + O(3^8)))

        A case where a initial substitution leads to lower precision loss than
        computing the echelonization of the Sylvester matrix at once::

            sage: S = ZpFM(3,20)
            sage: R.<t> = S[]
            sage: f = 2118930912*t^2 + 2699481150*t + 398851877
            sage: g = 2507761323*t^4 + 1524377763*t^3 + 2528732334*t^2 + 2318123756*t + 3016937729
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,g,use_substitution=True)[1]
            17
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,g,use_substitution=False)[1]
            12
            sage: S._pAdicGeneric__xgcd_univariate_polynomial_fixed(f,g)[1] # the heuristic picks use_substitution = True in this case
            17

        """
        polys = [f,g]

        # During the echelonization process, pivots of positive valuation kill
        # a lot of precision. Doing an initial substitution if a leading
        # coefficient has positive valuation might kill less precision.

        # In such a substitution, there are two kinds of precision loss possible:
        # - when we shift a coefficient down, e.g., 3+O(3^2) -> 1+O(3); this is
        #   kept track of in precision_loss_down
        # - when we shift a coefficient up, e.g., 1+O(3) -> O(3), we cannot say
        #   anything about the coefficients that were shifted out; this is kept
        #   track of in precision_loss_up

        # We explicitely keep track of the precision, so that this method also
        # applies to fixed mod elements
        precision_cap = min(f.parent().base_ring().precision_cap(), min([c.precision_absolute() for c in f]), min([c.precision_absolute() for c in g]))
        if precision_cap == 0:
            raise PrecisionError("Insufficient precision to compute greatest common divisor - computation requires at least 1 additional digit.")
        assert precision_cap >= 0

        # A rough estimate of the loss of precision if we do not apply any
        # substitution initially
        if self.__xgcd_valuation(f.leading_coefficient(),precision_cap) == self.__xgcd_valuation(g.leading_coefficient(),precision_cap):
            estimated_precision_loss = 0
        else:
            lower_poly = f
            if self.__xgcd_valuation(f.leading_coefficient(),precision_cap) > self.__xgcd_valuation(g.leading_coefficient(),precision_cap):
                lower_poly = g
            estimated_precision_loss = self.__xgcd_valuation(lower_poly.leading_coefficient(),precision_cap) - min([self.__xgcd_valuation(c,precision_cap) for c in list(lower_poly)])
            estimated_precision_loss *= max([f.degree(),g.degree()])

        # We estimate the loss of precision if we did the following:
        # Do a substitution of the form t -> t/p^N and multiply
        # everything by p^M such that the leading coefficient is a unit
        estimated_precision_loss_substitution = infinity
        try:
            if all([poly.leading_coefficient().valuation() > min([c.valuation() for c in poly]) for poly in polys]):
                N = 0
                for poly in polys:
                    leading_subst_power = (self.__xgcd_valuation(poly.leading_coefficient(),precision_cap)/poly.degree()).ceil()
                    other_subst_powers = [ ( (self.__xgcd_valuation(poly.leading_coefficient(),precision_cap) - self.__xgcd_valuation(poly[i],precision_cap)) / (poly.degree()-i) ).ceil() for i in range(1, poly.degree()) ]
                    N = max(N, leading_subst_power, *other_subst_powers)
                M = [ poly.degree()*N - self.__xgcd_valuation(poly.leading_coefficient(),precision_cap) for poly in polys ]
            else:
                N = 0
                M = [0,0]

            precision_loss_up, precision_loss_down = (0, 0)
            shifted_polys = []
            for i,poly in enumerate(polys):
                new_coeffs = []
                for j,x in enumerate(poly):
                    shift = M[i]-N*j
                    new_coeffs.append( self.__xgcd_shift(x, shift, precision_cap) )
                    precision_loss_down = max(precision_loss_down, -shift)
                    precision_loss_up = max(precision_loss_up, shift)
                shifted_polys.append(poly.parent()(new_coeffs))

            estimated_precision_loss_substitution = precision_loss_up+precision_loss_down
        except PrecisionError:
            # if we ran into precision errors during this estimation, then
            # estimated_precision_loss_substitution is infinite and we won't
            # substitute
            pass

        # don't rely on the automatic choice if use_substition parameter is
        # present
        if use_substitution is not None:
            estimated_precision_loss_substitution = -1 if use_substitution else infinity

        # the precision loss appears to be better if we do such a substitution
        if estimated_precision_loss_substitution <= estimated_precision_loss:
            polys = shifted_polys
            # the final result can at most be accurate to this precision
            max_precision = precision_cap - precision_loss_up
            precision_cap = self.__xgcd_change_precision( precision_cap, precision_cap - precision_loss_down )
            assert not all([poly.leading_coefficient().valuation() > min([c.valuation() for c in poly]) for poly in polys])
        else:
            max_precision = precision_cap
            N = 0
            M = [0,0]

        # Echelonize the Sylvester matrix
        from copy import copy
        S = copy(polys[0].sylvester_matrix(polys[1]))

        s,l,precision_cap = self.__xgcd_echelon(S, precision_cap)
        assert len(l)==S.nrows()
        assert not s[0].is_zero(precision_cap)
        assert s[0].parent()==l[0].parent()
        assert all([c.precision_absolute()>=precision_cap for c in s])
        assert all([c.precision_absolute()>=precision_cap for c in l])

        leading_valuation = max(f.leading_coefficient().valuation(),g.leading_coefficient().valuation())
        if precision_cap <= leading_valuation:
            raise PrecisionError("Insufficient precision to compute greatest common divisor - computation requires at least %s additional digit(s)."%(leading_valuation-precision_cap+1))

        # undo the substitution
        # again we might lose precision when shifting down
        precision_loss_down = 0
        for i,x in enumerate(s):
            shift = (len(s)-i-1)*N - min(*M)
            precision_loss_down = max(precision_loss_down, -shift)
            if i==0 and shift + self.__xgcd_valuation(x,precision_cap) >= precision_cap:
                raise PrecisionError("Insufficient precision to compute greatest common divisor - computation requires at least %s additional digit(s)."%(shift + self.__xgcd_valuation(x,precision_cap) - precision_cap + 1))
            s[i] = self.__xgcd_shift(x, shift, precision_cap)
        for i,x in enumerate(l):
            shift = (g.degree()-i-1 if i<g.degree() else len(l)-i-1)*N + (M[0] if i<g.degree() else M[1]) - min(*M)
            precision_loss_down = max(precision_loss_down, -shift)
            l[i] = self.__xgcd_shift(x, shift, precision_cap)
        precision_cap = self.__xgcd_change_precision(precision_cap, min( max_precision, precision_cap - precision_loss_down ) )

        assert all([c.precision_absolute()>=precision_cap for c in s])
        assert all([c.precision_absolute()>=precision_cap for c in l])

        if self.__xgcd_valuation(s[0],precision_cap) >= precision_cap:
            raise PrecisionError("Insufficient precision to compute greatest common divisor - computation requires at least %s additional digit(s)."%(self.__xgcd_valuation(s[0],precision_cap) - precision_cap + 1))

        # construct the xgcd polynomials
        s.reverse()
        xgcd = f.parent()(s).change_ring(s[0].parent())
        assert len(xgcd.list())==xgcd.degree()+1
        xgcd_polys = [ l[:g.degree()], l[g.degree():] ]
        for poly in xgcd_polys: poly.reverse()
        xgcd_polys = [ f.parent()(poly).change_ring(l[0].parent()) for poly in xgcd_polys ]
        assert all([c.precision_absolute()>=precision_cap for poly in xgcd_polys for c in poly])
        # normalize the leading coefficient
        xgcd_polys = [ xgcd.leading_coefficient().unit_part().inverse_of_unit().lift_to_precision(precision_cap) * poly for poly in xgcd_polys ]
        assert all([c.precision_absolute()>=precision_cap for poly in xgcd_polys for c in poly])
        xgcd *= xgcd.leading_coefficient().unit_part().inverse_of_unit().lift_to_precision(precision_cap)
        assert [ c.is_equal_to(d, precision_cap) for c,d in zip(xgcd_polys[0]*f + xgcd_polys[1]*g, xgcd) ], "%s * %s + %s * %s == %s != %s (mod p^%s)"%(xgcd_polys[0],f,xgcd_polys[1],g,xgcd_polys[0]*f + xgcd_polys[1]*g,xgcd,precision_cap)

        # the content of the gcd is the gcd of the contents of the input polynomials
        contents = f.content(), g.content(), xgcd.content()
        contents = [ c.gen() if hasattr(c,'gen') else c for c in contents ]
        gcd_valuation = self.__xgcd_valuation(contents[0].gcd(contents[1]),precision_cap)
        shift = gcd_valuation - self.__xgcd_valuation(contents[2],precision_cap)
        assert shift<=0
        gcd = xgcd.map_coefficients(lambda c: self.__xgcd_shift(c, shift, precision_cap))
        gcd_precision_cap = precision_cap + shift

        # if the gcd is constant, then it can be lifted to any precision
        if gcd.degree()==0:
            gcd_precision_cap = f.parent().base_ring().precision_cap()
            gcd = gcd.map_coefficients(lambda c: c.lift_to_precision(gcd_precision_cap))

        assert gcd.degree() == xgcd.degree()

        return gcd,gcd_precision_cap,xgcd,precision_cap,xgcd_polys[0],xgcd_polys[1]

    def __xgcd_echelon(self, S, precision_cap):
        """
        Helper method for :meth:`__xgcd_univariate_polynomial_fixed` which
        performs echelonization on a sylvester matrix ``S``.

        INPUT:

            - ``S`` -- a square matrix over ``self``

            - ``precision_cap`` -- an integer, precision to which the entries
              of ``S`` are accurate

        OUTPUT:

        A tuple ``s, l, precision_cap``:

            - ``s`` -- a list of coefficients of the last nonzero row in the
              echelonization of the matrix ``S``, starting from the first
              nonzero entry

            - ``l`` -- a list which as a vector multiplied from the left on
              ``S`` gives the last non-zero row of the echelonzation of ``S``

            - ``precision_cap`` -- an integer, the precision to which the
              entries of ``s`` and ``l`` are accurate

        ALGORITHM:

            Uses Gauss elimination with row pivoting. We cannot use the
            echelonizations implemented for general matrices since we have to
            keep track of precision loss that results from certain operations.

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        TESTS::

            sage: R.<t> = ZZ[]
            sage: f = 18*t^2 + 20*t + 6
            sage: g = 9*t^4 + 7*t^3 + 16*t^2 + 23*t + 6
            sage: Z = ZpFM(3,20)
            sage: S = f.change_ring(Z).sylvester_matrix(g.change_ring(Z))
            sage: Z._pAdicGeneric__xgcd_echelon(S,4)
            ([1 + O(3^20), 3 + 2*3^3 + 2*3^4 + 2*3^7 + 2*3^8 + 2*3^11 + 2*3^12 + 2*3^15 + 2*3^16 + 2*3^19 + O(3^20)], [O(3^20), O(3^20), O(3^20), 2 + 3 + 2*3^2 + 2*3^3 + 3^4 + 3^5 + 2*3^6 + 2*3^7 + 3^8 + 3^9 + 2*3^10 + 2*3^11 + 3^12 + 3^13 + 2*3^14 + 2*3^15 + 3^16 + 3^17 + 2*3^18 + 2*3^19 + O(3^20), O(3^20), O(3^20)], 2)
            sage: Z._pAdicGeneric__xgcd_echelon(S,20)
            ([3^3 + 2*3^4 + 3^5 + 3^6 + 3^8 + 3^13 + 3^14 + 3^15 + 3^17 + 3^18 + 2*3^19 + O(3^20)], [2*3^2 + 2*3^3 + 2*3^5 + 3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 3^10 + 2*3^11 + 3^12 + 2*3^15 + 3^16 + 3^17 + 2*3^19 + O(3^20), 3^2 + 2*3^3 + 2*3^4 + 3^5 + 2*3^6 + 2*3^7 + 3^10 + 2*3^11 + 3^12 + 3^13 + 3^14 + 3^17 + 2*3^18 + 2*3^19 + O(3^20), 3^4 + 3^6 + 3^7 + 3^9 + 2*3^11 + 2*3^12 + 3^15 + 2*3^16 + 2*3^18 + 3^19 + O(3^20), 2*3^3 + 3^4 + 3^5 + 2*3^6 + 3^8 + 3^11 + 2*3^12 + 2*3^13 + 3^14 + 2*3^16 + O(3^20), 3^3 + 2*3^6 + 2*3^7 + 2*3^9 + 3^11 + 3^13 + 3^15 + 3^16 + 2*3^18 + O(3^20), 1 + O(3^20)], 17)

        """
        from copy import copy
        # L keeps track of the steps performed in the echelonization; similar
        # to the L in an LU decomposition
        L = copy(S.parent().one())

        for c in range(S.ncols()):
            # search for a pivot of minimal valuation in the column c
            pivot = c
            for r in range(c,S.nrows()):
                if self.__xgcd_valuation(S[r,c],precision_cap) < self.__xgcd_valuation(S[pivot,c],precision_cap):
                    pivot = r

            S.swap_rows(pivot,c)
            L.swap_rows(pivot,c)

            # if the entire column is zero, then everything to the lower right
            # of c,c has to be zero as well since the Sylvester matrix encodes
            # linear combinations of two polynomials: if it is not possible to
            # produce a linear combination for a polynomial of degree d, then
            # one cannot produce a polynomial of degree less than d
            if S[c,c].is_zero(precision_cap):
                assert c!=0
                for rr in range(c,S.nrows()):
                    for cc in range(c,S.ncols()):
                        # however, this can happen when a loss of precision
                        # turned a non-zero entry into a zero
                        if not S[rr,cc].is_zero(precision_cap):
                            raise PrecisionError("Insufficient precision to compute greatest common divisor")
                # a previous row must have been the last nonzero row
                if S[c-1,c-1].is_zero(precision_cap):
                    raise PrecisionError("Insufficient precision to compute greatest common divisor")
                return S[c-1].list()[c-1:],L[c-1].list(),precision_cap

            # check if this is the last nonzero row
            is_last_nonzero_row = True
            for rr in range(c+1,S.nrows()):
                for cc in range(c,S.nrows()):
                    if not S[rr,cc].is_zero(precision_cap):
                        is_last_nonzero_row = False
            if is_last_nonzero_row:
                return S[c].list()[c:],L[c].list(),precision_cap

            # standard row elimination
            L.rescale_row(c,S[c,c].unit_part().inverse_of_unit().lift_to_precision(precision_cap))
            S.rescale_row(c,S[c,c].unit_part().inverse_of_unit().lift_to_precision(precision_cap))

            assert S[c,c].is_equal_to(S.base_ring().one()<<S[c,c].valuation(),precision_cap),"%s != pi^%s (mod pi^%s)"%(S[c,c],S[c,c].valuation(),precision_cap)

            # during the echelonization process we lose precision if the
            # valuation of the pivot exceeds the valuation of an entry in the
            # pivot's row
            precision_loss_up = self.__xgcd_valuation(S[c,c],precision_cap) - min([self.__xgcd_valuation(S[c,cc],precision_cap) for cc in range(c+1,S.ncols())])
            if precision_loss_up > 0:
                precision_cap = self.__xgcd_change_precision(precision_cap, precision_cap - precision_loss_up)
            for r in range(c+1,S.nrows()):
                if S[r,c].is_zero(precision_cap):
                    continue
                shift = self.__xgcd_valuation(S[r,c],precision_cap)-self.__xgcd_valuation(S[c,c],precision_cap)
                assert shift>=0
                L.add_multiple_of_row(r,c,-S[r,c].unit_part().lift_to_precision(precision_cap)<<shift)
                S.add_multiple_of_row(r,c,-S[r,c].unit_part().lift_to_precision(precision_cap)<<shift)
                assert S[r,c].is_zero(precision_cap),"%s != 0 (mod pi^%s)"%(S[r,c]._x,precision_cap)

        # must return before leaving the loop
        assert False

    def __cmp__(self, other):
        """
        Returns 0 if self == other, and 1 or -1 otherwise.

        We consider two p-adic rings or fields to be equal if they are
        equal mathematically, and also have the same precision cap and
        printing parameters.

        EXAMPLES::

            sage: R = Qp(7)
            sage: S = Qp(7,print_mode='val-unit')
            sage: R == S
            False
            sage: S = Qp(7,type='capped-rel')
            sage: R == S
            True
            sage: R is S
            True
        """
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        if self.prime() < other.prime():
            return -1
        elif self.prime() > other.prime():
            return 1
        try:
            if self.halting_parameter() < other.halting_parameter():
                return -1
            elif self.halting_parameter() > other.halting_parameter():
                return 1
        except AttributeError:
            pass
        if self.precision_cap() < other.precision_cap():
            return -1
        elif self.precision_cap() > other.precision_cap():
            return 1
        return self._printer.cmp_modes(other._printer)

    #def ngens(self):
    #    return 1

    #def gen(self, n = 0):
    #    if n != 0:
    #        raise IndexError, "only one generator"
    #    return self(self.prime())

    def print_mode(self):
        r"""
        Returns the current print mode as a string.

        INPUT::

            self -- a p-adic field

        OUTPUT::

            string -- self's print mode

        EXAMPLES::

            sage: R = Qp(7,5, 'capped-rel')
            sage: R.print_mode()
            'series'
        """
        return self._printer._print_mode()

    #def element_class(self):
    #    return self._element_class

    def characteristic(self):
        r"""
        Returns the characteristic of self, which is always 0.

        INPUT::

            self -- a p-adic parent

        OUTPUT::

            integer -- self's characteristic, i.e., 0

        EXAMPLES::

            sage: R = Zp(3, 10,'fixed-mod'); R.characteristic()
            0
        """
        return Integer(0)

    def prime(self):
        """
        Returns the prime, ie the characteristic of the residue field.

        INPUT::

            self -- a p-adic parent

        OUTPUT::

            integer -- the characteristic of the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: R.prime()
            3
        """
        return self.prime_pow._prime()

    def uniformizer_pow(self, n):
        """
        Returns p^n, as an element of self.

        If n is infinity, returns 0.

        EXAMPLES::

            sage: R = Zp(3, 5, 'fixed-mod')
            sage: R.uniformizer_pow(3)
            3^3 + O(3^5)
            sage: R.uniformizer_pow(infinity)
            O(3^5)
        """
        if n is infinity:
            return self(0)
        return self(self.prime_pow.pow_Integer_Integer(n))

    def _unram_print(self):
        """
        For printing.  Will be None if the unramified subextension of self is of degree 1 over Z_p or Q_p.

        EXAMPLES::

            sage: Zp(5)._unram_print()
        """
        return None

    def residue_characteristic(self):
        """
        Returns the prime, i.e., the characteristic of the residue field.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            integer -- the characteristic of the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: R.residue_characteristic()
            3
        """
        return self.prime()

    def residue_class_field(self):
        """
        Returns the residue class field.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_class_field()
            sage: k
            Finite Field of size 3
        """
        from sage.rings.finite_rings.constructor import GF
        return GF(self.prime())

    def residue_field(self):
        """
        Returns the residue class field.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            the residue field

        EXAMPLES::

            sage: R = Zp(3,5,'fixed-mod')
            sage: k = R.residue_field()
            sage: k
            Finite Field of size 3
        """
        return self.residue_class_field()

    def residue_system(self):
        """
        Returns a list of elements representing all the residue classes.

        INPUT::

            self -- a p-adic ring

        OUTPUT::

            list of elements -- a list of elements representing all the residue classes

        EXAMPLES::

            sage: R = Zp(3, 5,'fixed-mod')
            sage: R.residue_system()
            [O(3^5), 1 + O(3^5), 2 + O(3^5)]
        """
        return [self(i) for i in self.residue_class_field()]

    def teichmuller(self, x, prec = None):
        r"""
        Returns the teichmuller representative of x.

        INPUT::

            - self -- a p-adic ring
            - x -- something that can be cast into self

        OUTPUT::

            - element -- the teichmuller lift of x

        EXAMPLES::

            sage: R = Zp(5, 10, 'capped-rel', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Qp(5, 10,'capped-rel','series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5, 10, 'capped-abs', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5, 10, 'fixed-mod', 'series')
            sage: R.teichmuller(2)
            2 + 5 + 2*5^2 + 5^3 + 3*5^4 + 4*5^5 + 2*5^6 + 3*5^7 + 3*5^9 + O(5^10)
            sage: R = Zp(5,5)
            sage: S.<x> = R[]
            sage: f = x^5 + 75*x^3 - 15*x^2 +125*x - 5
            sage: W.<w> = R.ext(f)
            sage: y = W.teichmuller(3); y
            3 + 3*w^5 + w^7 + 2*w^9 + 2*w^10 + 4*w^11 + w^12 + 2*w^13 + 3*w^15 + 2*w^16 + 3*w^17 + w^18 + 3*w^19 + 3*w^20 + 2*w^21 + 2*w^22 + 3*w^23 + 4*w^24 + O(w^25)
            sage: y^5 == y
            True
            sage: g = x^3 + 3*x + 3
            sage: A.<a> = R.ext(g)
            sage: b = A.teichmuller(1 + 2*a - a^2); b
            (4*a^2 + 2*a + 1) + 2*a*5 + (3*a^2 + 1)*5^2 + (a + 4)*5^3 + (a^2 + a + 1)*5^4 + O(5^5)
            sage: b^125 == b
            True

        AUTHORS:

        - Initial version: David Roe
        - Quadratic time version: Kiran Kedlaya <kedlaya@math.mit.edu> (3/27/07)
        """
        if prec is None:
            prec = self.precision_cap()
        else:
            prec = min(Integer(prec), self.precision_cap())
        ans = self(x, prec)
        ans._teichmuller_set()
        return ans

    def teichmuller_system(self):
        r"""
        Returns a set of teichmuller representatives for the invertible elements of $\ZZ / p\ZZ$.

        INPUT:

        - self -- a p-adic ring

        OUTPUT:

        - list of elements -- a list of teichmuller representatives for the invertible elements of $\ZZ / p\ZZ$

        EXAMPLES::

            sage: R = Zp(3, 5,'fixed-mod', 'terse')
            sage: R.teichmuller_system()
            [1 + O(3^5), 242 + O(3^5)]

        NOTES:

        Should this return 0 as well?
        """
        return [self.teichmuller(i.lift()) for i in self.residue_class_field() if i != 0]

#     def different(self):
#         raise NotImplementedError

#     def automorphisms(self):
#         r"""
#         Returns the group of automorphisms of $\Z_p$, i.e. the trivial group.
#         """
#         raise NotImplementedError

#     def galois_group(self):
#         r"""
#         Returns the Galois group of $\Z_p$, i.e. the trivial group.
#         """
#         raise NotImplementedError

#     def hasGNB(self):
#         r"""
#         Returns whether or not $\Z_p$ has a Gauss Normal Basis.
#         """
#         raise NotImplementedError

    def extension(self, modulus, prec = None, names = None, print_mode = None, halt = None, **kwds):
        """
        Create an extension of this p-adic ring.

        EXAMPLES::

            sage: k = Qp(5)
            sage: R.<x> = k[]
            sage: l.<w> = k.extension(x^2-5); l
            Eisenstein Extension of 5-adic Field with capped relative precision 20 in w defined by (1 + O(5^20))*x^2 + (O(5^21))*x + (4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19 + 4*5^20 + O(5^21))

            sage: F = list(Qp(19)['x'](cyclotomic_polynomial(5)).factor())[0][0]
            sage: L = Qp(19).extension(F, names='a')
            sage: L
            Unramified Extension of 19-adic Field with capped relative precision 20 in a defined by (1 + O(19^20))*x^2 + (5 + 2*19 + 10*19^2 + 14*19^3 + 7*19^4 + 13*19^5 + 5*19^6 + 12*19^7 + 8*19^8 + 4*19^9 + 14*19^10 + 6*19^11 + 5*19^12 + 13*19^13 + 16*19^14 + 4*19^15 + 17*19^16 + 8*19^18 + 4*19^19 + O(19^20))*x + (1 + O(19^20))
        """
        from sage.rings.padics.factory import ExtensionFactory
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'print_mode': print_mode}
        else:
            if not isinstance(print_mode, dict):
                print_mode = dict(print_mode)
            for option in ['mode', 'pos', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
                if print_mode.has_key(option):
                    print_mode["print_" + option] = print_mode[option]
                    del print_mode[option]
                elif not print_mode.has_key("print_" + option):
                    if kwds.has_key("print_" + option):
                        print_mode["print_" + option] = kwds["print_" + option]
                    else:
                        print_mode["print_" + option] = self._printer.dict()[option]
            for option in ['ram_name', 'unram_name', 'var_name']:
                if not print_mode.has_key(option):
                    if kwds.has_key(option):
                        print_mode[option] = kwds[option]
                    else:
                        print_mode[option] = self._printer.dict()[option]
        return ExtensionFactory(base=self, premodulus=modulus, prec=prec, halt=halt, names=names, check = True, **print_mode)

    @cached_method
    def _log_unit_part_p(self):
        """
        Compute the logarithm of the unit-part of `p`.

        If `\pi` is the uniformizer in this ring, then we can uniquely write
        `p=\pi^e u` where `u` is a `\pi`-adic unit. This method computes the
        logarithm of `u`.

        This is a helper method for
        :meth:`sage.rings.padics.padic_generic_element.pAdicGenericElement.log`.

        TESTS::

            sage: R = Qp(3,5)
            sage: R._log_unit_part_p()
            O(3^5)

            sage: S.<x> = ZZ[]
            sage: W.<pi> = R.extension(x^3-3)
            sage: W._log_unit_part_p()
            O(pi^15)

            sage: W.<pi> = R.extension(x^3-3*x-3)
            sage: W._log_unit_part_p()
            2 + pi + 2*pi^2 + pi^4 + pi^5 + 2*pi^7 + 2*pi^8 + pi^9 + 2*pi^10 + pi^11 + pi^12 + 2*pi^14 + O(pi^15)

        """
        return self(self.prime()).unit_part().log()

    @cached_method
    def _exp_p(self):
        """
        Compute the exponential of `p`.

        This is a helper method for
        :meth:`sage.rings.padics.padic_generic_element.pAdicGenericElement.exp`.

        TESTS::

            sage: R = Qp(3, 5)
            sage: R._exp_p()
            1 + 3 + 3^2 + 2*3^3 + 2*3^4 + O(3^5)

            sage: S.<x> = ZZ[]
            sage: W.<pi> = R.extension(x^3-3)
            sage: W._exp_p()
            1 + pi^3 + pi^6 + 2*pi^9 + 2*pi^12 + O(pi^15)
            sage: R._exp_p() == W._exp_p()
            True

            sage: W.<pi> = R.extension(x^3-3*x-3)
            sage: W._exp_p()
            1 + pi^3 + 2*pi^4 + pi^5 + pi^7 + pi^9 + pi^10 + 2*pi^11 + pi^12 + pi^13 + 2*pi^14 + O(pi^15)
            sage: R._exp_p() == W._exp_p()
            True

        """
        p = self.prime()
        if p == 2:
            # the exponential of 2 does not exists, so we compute the
            # exponential of 4 instead.
            p = 4
        return self(p)._exp(self.precision_cap())

def local_print_mode(obj, print_options, pos = None, ram_name = None):
    r"""
    Context manager for safely temporarily changing the print_mode
    of a p-adic ring/field.

    EXAMPLES::

        sage: R = Zp(5)
        sage: R(45)
        4*5 + 5^2 + O(5^21)
        sage: with local_print_mode(R, 'val-unit'):
        ...       print R(45)
        ...
        5 * 9 + O(5^21)

    NOTES::

        For more documentation see localvars in parent_gens.pyx
    """
    if isinstance(print_options, str):
        print_options = {'mode': print_options}
    elif not isinstance(print_options, dict):
        raise TypeError, "print_options must be a dictionary or a string"
    if pos is not None:
        print_options['pos'] = pos
    if ram_name is not None:
        print_options['ram_name'] = ram_name
    for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
        if not print_options.has_key(option):
            print_options[option] = obj._printer.dict()[option]
    return pAdicPrinter(obj, print_options)
