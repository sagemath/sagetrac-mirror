"""
`p`-Adic Generic Nodes

This file contains a bunch of intermediate classes for the `p`-adic
parents, allowing a function to be implemented at the right level of
generality.

AUTHORS:

- David Roe

- Julian Rueth (2012-09-05): (x)gcd for polynomials
"""
#*****************************************************************************
#       Copyright (C) 2009 David Roe <roed@math.harvard.edu>
#                     2012 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.padics.local_generic import LocalGeneric
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.ring import EuclideanDomain, Field
from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

class CappedAbsoluteGeneric(LocalGeneric):
    def is_capped_absolute(self):
        """
        Returns whether this `p`-adic ring bounds precision in a
        capped absolute fashion.

        The absolute precision of an element is the power of `p` modulo
        which that element is defined.  In a capped absolute ring, the
        absolute precision of elements are bounded by a constant
        depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_absolute()
            True
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_absolute()
            False
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpCA(5)._prec_type()
            'capped-abs'
        """
        return 'capped-abs'

class CappedRelativeGeneric(LocalGeneric):
    def is_capped_relative(self):
        """
        Returns whether this `p`-adic ring bounds precision in a capped
        relative fashion.

        The relative precision of an element is the power of p modulo
        which the unit part of that element is defined.  In a capped
        relative ring, the relative precision of elements are bounded
        by a constant depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_relative()
            False
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_relative()
            True
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: Zp(5)._prec_type()
            'capped-rel'
        """
        return 'capped-rel'

class FixedModGeneric(LocalGeneric):
    def is_fixed_mod(self):
        """
        Returns whether this `p`-adic ring bounds precision in a fixed
        modulus fashion.

        The absolute precision of an element is the power of p modulo
        which that element is defined.  In a fixed modulus ring, the
        absolute precision of every element is defined to be the
        precision cap of the parent.  This means that some operations,
        such as division by `p`, don't return a well defined answer.

        EXAMPLES::

            sage: R = ZpFM(5,15)
            sage: R.is_fixed_mod()
            True
            sage: R(5^7,absprec=9)
            5^7 + O(5^15)
            sage: S = ZpCA(5, 15)
            sage: S.is_fixed_mod()
            False
            sage: S(5^7,absprec=9)
            5^7 + O(5^9)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpFM(5)._prec_type()
            'fixed-mod'
        """
        return 'fixed-mod'

class CappedRelativeRingGeneric(CappedRelativeGeneric):
    pass
class CappedRelativeFieldGeneric(CappedRelativeGeneric):#, sage.rings.ring.Field):
    pass


def is_pAdicRing(R):
    """
    Returns ``True`` if and only if ``R`` is a `p`-adic ring (not a
    field).

    EXAMPLES::

        sage: is_pAdicRing(Zp(5))
        True
        sage: is_pAdicRing(RR)
        False
    """
    return isinstance(R, pAdicRingGeneric)

class pAdicRingGeneric(pAdicGeneric, EuclideanDomain):
    def is_field(self, proof = True):
        """
        Returns whether this ring is actually a field, ie ``False``.

        EXAMPLES::

            sage: Zp(5).is_field()
            False
        """
        return False


    def krull_dimension(self):
        r"""
        Returns the Krull dimension of self, i.e. 1

        INPUT:

        - self -- a `p`-adic ring

        OUTPUT:

        - the Krull dimension of self.  Since self is a `p`-adic ring,
          this is 1.

        EXAMPLES::

            sage: Zp(5).krull_dimension()
            1
        """
        return 1

def is_pAdicField(R):
    """
    Returns ``True`` if and only if ``R`` is a `p`-adic field.

    EXAMPLES::

        sage: is_pAdicField(Zp(17))
        False
        sage: is_pAdicField(Qp(17))
        True
    """
    return isinstance(R, pAdicFieldGeneric)

class pAdicFieldGeneric(pAdicGeneric, Field):
    pass

    #def class_field(self, group=None, map=None, generators=None):
    #    raise NotImplementedError

    #def composite(self, subfield1, subfield2):
    #    raise NotImplementedError

    #def norm_equation(self):
    #    raise NotImplementedError

    #def norm_group(self):
    #    raise NotImplementedError

    #def norm_group_discriminant(self, group=None, map=None, generators=None):
    #    raise NotImplementedError

    #def number_of_extensions(self, degree, discriminant=None, e=None, f=None):
    #    raise NotImplementedError

    #def list_of_extensions(self, degree, discriminant=None, e=None, f=None):
    #    raise NotImplementedError

    #def subfield(self, list):
    #    raise NotImplementedError

    #def subfield_lattice(self):
    #    raise NotImplementedError

    #def subfields_of_degree(self, n):
    #    raise NotImplementedError

class pAdicFixedModRingGeneric(pAdicRingGeneric, FixedModGeneric):
    def _gcd_univariate_polynomial(self, f, g):
        """
        Compute a greatest common divisor of the polynomials ``f`` and ``g``.

        This is a helper method for
        :meth:`sage.rings.polynomial.polynomial_element.Polynomial.gcd`. Its
        implementation relies on
        :meth:`sage.rings.padics.padic_generic.pAdicGeneric._gcd_univariate_polynomial_fixed`
        which should be consulted for further details and examples.

        INPUT:

            - ``f``, ``g`` -- two polynomials defined over ``self``.

        OUTPUT:

        A polynomial defined over ``self``, and the precision to which this
        result is accurate (see
        :meth:`sage.rings.padics.padic_generic.pAdicGeneric.__xgcd_univariate_polynomial_fixed`
        for the precise meaning of this accuracy.)

        AUTHORS:

        - Julian Rueth (2012-09-05): initial version

        EXAMPLES::

            sage: R.<t> = ZpFM(3,20)[]
            sage: (t + 1).gcd( (t - 1) * (t + 1) )
            ((1 + O(3^20))*t + (1 + O(3^20)), 20)
            sage: (t^3).gcd( t^5 )
            ((1 + O(3^20))*t^3, 20)

        Also works over extensions::

            sage: K = ZpFM(3,20)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension( a^2 - 3 ) # Eisenstein extension
            sage: R.<t> = L[]
            sage: a0,a1,a2 = 12345678+a,90123456-a,78901234*a
            sage: f = (t - a0) * (t - a1)^2 * (t - a2)^3
            sage: g,prec = f.gcd(f.derivative())
            sage: h = ((t - a1) * (t - a2)^2) # the correct result
            sage: g.degree() == h.degree()
            True
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(g),list(h))]) # g and h are equal mod a^prec
            True

            sage: R.<a> = K[]
            sage: L.<a> = K.extension( a^2 - 2 ) # unramified extension
            sage: R.<t> = L[]
            sage: a0,a1,a2 = 12345678+a,90123456-a,78901234*a
            sage: f = (t - a0) * (t - a1)^2 * (t - a2)^3
            sage: g,prec = f.gcd(f.derivative())
            sage: h = ((t - a1) * (t - a2)^2) # the correct result
            sage: g.degree() == h.degree()
            True
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(g),list(h))]) # g and h are equal mod p^prec
            True

        TESTS:

        Check that the examples from :trac:`13439` work::

            sage: R.<t> = ZpFM(3,3)[]
            sage: f = 3*t + 7
            sage: g = 5*t + 9
            sage: f.gcd(f*g)
            ((3 + O(3^3))*t + (1 + 2*3 + O(3^3)), 2)

            sage: R.<t> = ZpFM(3,20)[]
            sage: f = 729*490473657*t + 257392844
            sage: g = 225227399*t - 59049*8669753175
            sage: h,prec = f.gcd(f*g)
            sage: h.degree() == f.degree()
            True
            sage: h *= f.leading_coefficient().unit_part()
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(f),list(h))]) # f and h are equal mod p^prec
            True

        """
        return self._gcd_univariate_polynomial_fixed(f, g)

    def _xgcd_univariate_polynomial(self, f, g):
        """
        Compute an extended greatest common divisor of the polynomials ``f``
        and ``g``.

        This is a helper method for
        :meth:`sage.rings.polynomial.polynomial_element.Polynomial.xgcd`. Its
        implementation relies on
        :meth:`sage.rings.padics.padic_generic.pAdicGeneric._xgcd_univariate_polynomial_fixed`
        which should be consulted for further details and examples.

        INPUT:

            - ``f``, ``g`` -- two polynomials defined over ``self``.

        OUTPUT:

        A tuple ``r,prec,s,t`` which satisfies ``r = s*f + t*g`` when reduced
        to precision ``prec``.  (see
        :meth:`sage.rings.padics.padic_generic.pAdicGeneric.__xgcd_univariate_polynomial_fixed`
        for the precise meaning of the precision ``prec``.)

        AUTHORS:

        - Julian Rueth (2012-10-22): initial version

        EXAMPLES::

            sage: R.<t> = ZpFM(3,20)[]
            sage: (t + 1).xgcd( (t - 1) * (t + 1) )
            ((1 + O(3^20))*t + (1 + O(3^20)), 20, (1 + O(3^20)), 0)
            sage: (t^3).xgcd( t^5 )
            ((1 + O(3^20))*t^3, 20, (1 + O(3^20)), 0)

        Also works over extensions::

            sage: K = ZpFM(3,20)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension( a^2 - 3 ) # Eisenstein extension
            sage: R.<t> = L[]
            sage: a0,a1,a2 = 12345678+a,90123456-a,78901234*a
            sage: f = (t - a0) * (t - a1)^2 * (t - a2)^3
            sage: g,prec,u,v = f.xgcd(f.derivative())
            sage: h = (t - a1) * (t - a2)^2 # the correct result
            sage: g.degree() == h.degree()
            True
            sage: h *= g.leading_coefficient()
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(g),list(h))]) # g and h are equal mod a^prec
            True
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(g),list(u*f+v*f.derivative()))]) # the equation g = u*f + v*f.derivative() is satisfied mod a^prec
            True

            sage: R.<a> = K[]
            sage: L.<a> = K.extension( a^2 - 2 ) # unramified extension
            sage: R.<t> = L[]
            sage: a0,a1,a2 = 12345678+a,90123456-a,78901234*a
            sage: f = (t - a0) * (t - a1)^2 * (t - a2)^3
            sage: g,prec,u,v = f.xgcd(f.derivative())
            sage: h = (t - a1) * (t - a2)^2 # the correct result
            sage: g.degree() == h.degree()
            True
            sage: h *= g.leading_coefficient()
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(g),list(h))]) # g and h are equal mod a^prec
            True
            sage: all([(c1-c2).is_zero(prec) for c1,c2 in zip(list(g),list(u*f+v*f.derivative()))]) # the equation g = u*f + v*f.derivative() is satisfied mod a^prec
            True

        TESTS:

        Check that the examples from :trac:`13439` work::

            sage: R.<t> = ZpFM(3,3)[]
            sage: f = 3*t + 7
            sage: g = 5*t + 9
            sage: f.xgcd(f*g)
            ((3 + O(3^3))*t + (1 + 2*3 + O(3^3)), 2, (1 + O(3^3)), 0)

            sage: R.<t> = ZpFM(3,20)[]
            sage: f = 729*490473657*t + 257392844
            sage: g = 225227399*t - 59049*8669753175
            sage: h,_,_,_ = f.xgcd(f*g)
            sage: h.degree() == f.degree()
            True

        """
        return self._xgcd_univariate_polynomial_fixed(f, g)

class pAdicCappedAbsoluteRingGeneric(pAdicRingGeneric, CappedAbsoluteGeneric):
    pass
class pAdicCappedRelativeRingGeneric(pAdicRingGeneric, CappedRelativeRingGeneric):
    pass
class pAdicCappedRelativeFieldGeneric(pAdicFieldGeneric, CappedRelativeFieldGeneric):
    pass

class pAdicRingBaseGeneric(pAdicBaseGeneric, pAdicRingGeneric):
    def construction(self):
        """
        Returns the functorial construction of self, namely,
        completion of the rational numbers with respect a given prime.

        Also preserves other information that makes this field unique
        (e.g. precision, rounding, print mode).

        EXAMPLE::

            sage: K = Zp(17, 8, print_mode='val-unit', print_sep='&')
            sage: c, L = K.construction(); L
            Integer Ring
            sage: c(L)
            17-adic Ring with capped relative precision 8
            sage: K == c(L)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(self.prime(),
                                  self.precision_cap(),
                                  {'print_mode':self._printer.dict(), 'type':self._prec_type(), 'names':self._names}),
                ZZ)

    def random_element(self, algorithm='default'):
        r"""
        Returns a random element of self, optionally using the
        algorithm argument to decide how it generates the
        element. Algorithms currently implemented:

        - default: Choose `a_i`, `i >= 0`, randomly between `0` and
          `p-1` until a nonzero choice is made. Then continue choosing
          `a_i` randomly between `0` and `p-1` until we reach
          precision_cap, and return `\sum a_i p^i`.

        EXAMPLES::

            sage: Zp(5,6).random_element()
            3 + 3*5 + 2*5^2 + 3*5^3 + 2*5^4 + 5^5 + O(5^6)
            sage: ZpCA(5,6).random_element()
            4*5^2 + 5^3 + O(5^6)
            sage: ZpFM(5,6).random_element()
            2 + 4*5^2 + 2*5^4 + 5^5 + O(5^6)
        """
        if (algorithm == 'default'):
            if self.is_capped_relative():
                i = 0
                a_i = ZZ.random_element(self.prime())
                while a_i.is_zero():
                    i += 1
                    a_i = ZZ.random_element(self.prime())
                return self((self.prime()**i)*(a_i + self.prime()*ZZ.random_element(self.prime_pow.pow_Integer_Integer(self.precision_cap()-1))))
            else:
                return self(ZZ.random_element(self.prime_pow.pow_Integer_Integer(self.precision_cap())))
        else:
            raise NotImplementedError, "Don't know %s algorithm"%algorithm

    #def unit_group(self):
    #    raise NotImplementedError

    #def unit_group_gens(self):
    #    raise NotImplementedError

    #def principal_unit_group(self):
    #    raise NotImplementedError

class pAdicFieldBaseGeneric(pAdicBaseGeneric, pAdicFieldGeneric):
    def composite(self, subfield1, subfield2):
        r"""
        Returns the composite of two subfields of self, i.e., the
        largest subfield containing both

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``subfield1`` -- a subfield
        - ``subfield2`` -- a subfield

        OUTPUT:

        - the composite of subfield1 and subfield2

        EXAMPLES::

            sage: K = Qp(17); K.composite(K, K) is K
            True
        """
        #should be overridden for extension fields
        if (subfield1 is self) and (subfield2 is self):
            return self
        raise ValueError, "Arguments must be subfields of self."

    def subfields_of_degree(self, n):
        r"""
        Returns the number of subfields of self of degree `n`

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``n`` -- an integer

        OUTPUT:

        - integer -- the number of subfields of degree ``n`` over self.base_ring()

        EXAMPLES::

            sage: K = Qp(17)
            sage: K.subfields_of_degree(1)
            1
        """
        if n == 1:
            return 1
        else:
            return 0

    def subfield(self, list):
        r"""
        Returns the subfield generated by the elements in list

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``list`` -- a list of elements of ``self``

        OUTPUT:

        - the subfield of ``self`` generated by the elements of list

        EXAMPLES::

            sage: K = Qp(17); K.subfield([K(17), K(1827)]) is K
            True
        """
        for x in list:
            if not self.__contains__(x):
                raise TypeError, "Members of the list of generators must be elements of self."
        return self

    def construction(self):
        """
        Returns the functorial construction of ``self``, namely,
        completion of the rational numbers with respect a given prime.

        Also preserves other information that makes this field unique
        (e.g. precision, rounding, print mode).

        EXAMPLE::

            sage: K = Qp(17, 8, print_mode='val-unit', print_sep='&')
            sage: c, L = K.construction(); L
            Rational Field
            sage: c(L)
            17-adic Field with capped relative precision 8
            sage: K == c(L)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(self.prime(),
                                  self.precision_cap(),
                                  {'print_mode':self._printer.dict(), 'type':self._prec_type(), 'names':self._names}),
                QQ)
