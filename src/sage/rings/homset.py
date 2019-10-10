"""
Space of homomorphisms between two rings
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from sage.categories.homset import HomsetWithBase
from sage.categories.rings import Rings
_Rings = Rings()

from . import morphism
from . import quotient_ring

def is_RingHomset(H):
    """
    Return ``True`` if ``H`` is a space of homomorphisms between two rings.

    EXAMPLES::

        sage: from sage.rings.homset import is_RingHomset as is_RH
        sage: is_RH(Hom(ZZ, QQ))
        True
        sage: is_RH(ZZ)
        False
        sage: is_RH(Hom(RR, CC))
        True
        sage: is_RH(Hom(FreeModule(ZZ,1), FreeModule(QQ,1)))
        False
    """
    return isinstance(H, RingHomset_generic)


def RingHomset(R, S, category = None):
    """
    Construct a space of homomorphisms between the rings ``R`` and ``S``.

    For more on homsets, see :func:`Hom()`.

    EXAMPLES::

        sage: Hom(ZZ, QQ) # indirect doctest
        Set of Homomorphisms from Integer Ring to Rational Field

    """
    if quotient_ring.is_QuotientRing(R):
        return RingHomset_quo_ring(R, S, category = category)
    return RingHomset_generic(R, S, category = category)


class RingHomset_generic(HomsetWithBase):
    """
    A generic space of homomorphisms between two rings.

    EXAMPLES::

        sage: Hom(ZZ, QQ)
        Set of Homomorphisms from Integer Ring to Rational Field
        sage: QQ.Hom(ZZ)
        Set of Homomorphisms from Rational Field to Integer Ring
    """
    def __init__(self, R, S, category = None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Hom(ZZ, QQ)
            Set of Homomorphisms from Integer Ring to Rational Field
        """
        if category is None:
            category = _Rings
        HomsetWithBase.__init__(self, R, S, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Hom(ZZ, QQ) # indirect doctest
            Set of Homomorphisms from Integer Ring to Rational Field
        """
        return "Set of Homomorphisms from %s to %s"%(self.domain(), self.codomain())

    def has_coerce_map_from(self, x):
        """
        The default for coercion maps between ring homomorphism spaces is
        very restrictive (until more implementation work is done).

        Currently this checks if the domains and the codomains are equal.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: H2 = Hom(QQ, ZZ)
            sage: H.has_coerce_map_from(H2)
            False
        """
        return (x.domain() == self.domain() and x.codomain() == self.codomain())

    def _coerce_impl(self, x):
        """
        Check to see if we can coerce ``x`` into a homomorphism with the
        correct rings.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: phi = H([1])
            sage: H2 = Hom(QQ, QQ)
            sage: phi2 = H2(phi); phi2 # indirect doctest
            Ring endomorphism of Rational Field
              Defn: 1 |--> 1
            sage: H(phi2) # indirect doctest
            Ring morphism:
              From: Integer Ring
              To:   Rational Field
              Defn: 1 |--> 1
        """
        from sage.categories.map import Map
        if not (isinstance(x, Map) and x.category_for().is_subcategory(Rings())):
            raise TypeError
        if x.parent() is self:
            return x
        # Case 1: the parent fits
        if x.parent() == self:
            if isinstance(x, morphism.RingHomomorphism_im_gens):
                return morphism.RingHomomorphism_im_gens(self, x.im_gens())
            elif isinstance(x, morphism.RingHomomorphism_cover):
                return morphism.RingHomomorphism_cover(self)
            elif isinstance(x, morphism.RingHomomorphism_from_base):
                return morphism.RingHomomorphism_from_base(self, x.underlying_map())
        # Case 2: unique extension via fraction field
        try:
            if isinstance(x, morphism.RingHomomorphism_im_gens) and x.domain().fraction_field().has_coerce_map_from(self.domain()):
                return morphism.RingHomomorphism_im_gens(self, x.im_gens())
        except Exception:
            pass
        # Case 3: the homomorphism can be extended by coercion
        try:
            return x.extend_codomain(self.codomain()).extend_domain(self.domain())
        except Exception:
            pass
        # Last resort, case 4: the homomorphism is induced from the base ring
        if self.domain()==self.domain().base() or self.codomain()==self.codomain().base():
            raise TypeError
        try:
            x = self.domain().base().Hom(self.codomain().base())(x)
            return morphism.RingHomomorphism_from_base(self, x)
        except Exception:
            raise TypeError

    def __call__(self, im_gens, check=True, base_map=None):
        """
        Create a homomorphism.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: H([1])
            Ring morphism:
              From: Integer Ring
              To:   Rational Field
              Defn: 1 |--> 1

        You can provide a morphism on the base::

            sage: k = GF(9)
            sage: z2 = k.gen()
            sage: cc = k.frobenius_endomorphism()
            sage: R.<x> = k[]
            sage: H = Hom(R, R)
            sage: phi = H([x^2], base_map=cc); phi
            Ring endomorphism of Univariate Polynomial Ring in x over Finite Field in z2 of size 3^2
              Defn: x |--> x^2
                    with map of base ring
            sage: phi(z2*x) == z2^3 * x^2
            True

            sage: R.<x> = ZZ[]
            sage: K.<a> = GF(7^2)
            sage: L.<u> = K.extension(x^3 - 3)
            sage: phi = L.hom([u^7], base_map=K.frobenius_endomorphism())
            sage: phi(u) == u^7
            True
            sage: phi(a) == a^7
            True

        TESTS::

            sage: H = Hom(ZZ, QQ)
            sage: H == loads(dumps(H))
            True
        """
        from sage.categories.map import Map
        if isinstance(im_gens, Map):
            if base_map is not None:
                raise ValueError("Cannot specify base_map when providing a map")
            return self._coerce_impl(im_gens)
        else:
            return morphism.RingHomomorphism_im_gens(self, im_gens, base_map=base_map, check=check)

    def natural_map(self):
        """
        Returns the natural map from the domain to the codomain.

        The natural map is the coercion map from the domain ring to the
        codomain ring.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: H.natural_map()
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        f = self.codomain().coerce_map_from(self.domain())
        if f is None:
            raise TypeError("natural coercion morphism from %s to %s not defined"%(self.domain(), self.codomain()))
        return f

    def zero(self):
        r"""
        Return the zero element of this homset.

        EXAMPLES:

        Since a ring homomorphism maps 1 to 1, there can only be a zero
        morphism when mapping to the trivial ring::

            sage: Hom(ZZ, Zmod(1)).zero()
            Ring morphism:
              From: Integer Ring
              To:   Ring of integers modulo 1
              Defn: 1 |--> 0
            sage: Hom(ZZ, Zmod(2)).zero()
            Traceback (most recent call last):
            ...
            ValueError: homset has no zero element

        """
        if not self.codomain().is_zero():
            raise ValueError("homset has no zero element")
        # there is only one map in this homset
        return self.an_element()


class RingHomset_quo_ring(RingHomset_generic):
    """
    Space of ring homomorphisms where the domain is a (formal) quotient
    ring.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)
        sage: phi = S.hom([b,a]); phi
        Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
          Defn: a |--> b
                b |--> a
        sage: phi(a)
        b
        sage: phi(b)
        a

    TESTS:

    We test pickling of a homset from a quotient.

    ::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)
        sage: H = S.Hom(R)
        sage: H == loads(dumps(H))
        True

    We test pickling of actual homomorphisms in a quotient::

        sage: phi = S.hom([b,a])
        sage: phi == loads(dumps(phi))
        True
    """
    def __call__(self, im_gens, base_map=None, check=True):
        """
        Create a homomorphism.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quotient(x^2 + y^2)
            sage: H = S.Hom(R)
            sage: phi = H([b,a]); phi
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: a |--> b
                    b |--> a

        """
        if isinstance(im_gens, morphism.RingHomomorphism_from_quotient):
            return morphism.RingHomomorphism_from_quotient(self, im_gens._phi())
        try:
            pi = self.domain().cover()
            phi = pi.domain().hom(im_gens, base_map=base_map, check=check)
            return morphism.RingHomomorphism_from_quotient(self, phi)
        except (NotImplementedError, ValueError):
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError("images do not define a valid homomorphism")

    def _coerce_impl(self, x):
        """
        Check to see if we can coerce ``x`` into a homomorphism with the
        correct rings.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quotient(x^2 + y^2)
            sage: H = S.Hom(R)
            sage: phi = H([b,a])
            sage: R2.<x,y> = PolynomialRing(ZZ, 2)
            sage: H2 = Hom(R2, S)
            sage: H2(phi) # indirect doctest
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Integer Ring
              To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
              Defn:   Coercion map:
                      From: Multivariate Polynomial Ring in x, y over Integer Ring
                      To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
                    then
                      Ring morphism:
                      From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                      Defn: a |--> b
                            b |--> a
                    then
                      Coercion map:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)

        """
        if not isinstance(x, morphism.RingHomomorphism_from_quotient):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return morphism.RingHomomorphism_from_quotient(self, x._phi())
        raise TypeError
