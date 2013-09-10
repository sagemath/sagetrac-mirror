"""
The discrete valuation of a `p`-adic ring

This file makes it possible to use `p`-adics, integers, and rationals in the
general discrete valuation framework.

AUTHORS:

- Julian Rueth (2013-03-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from discrete_valuation import DiscreteValuation
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.functional import coerce

def pAdicValuation(domain, p = None):
    r"""
    Create a ``p``-adic valuation on ``domain``.

    INPUT:

    - ``domain`` -- a ring

    - ``p`` -- a prime or ``None`` (default: ``None``), if ``None``, tries
      to automatically determine the appropriate prime for ``domain``

    EXAMPLES::

        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

    ``p`` may be ``None`` for `p`-adic rings::

        sage: pAdicValuation(Qp(2))
        2-adic valuation
        sage: pAdicValuation(ZZ)
        Traceback (most recent call last):
        ...
        ValueError: prime must be specified for this ring

    """
    from sage.rings.all import ZZ, QQ
    from sage.rings.padics.padic_generic import pAdicGeneric
    from sage.misc.functional import coerce
    if domain is ZZ or domain is QQ:
        if p is None:
            raise ValueError("prime must be specified for this ring")
        return pAdicValuation_int(domain, p)
    elif isinstance(domain, pAdicGeneric):
        if p is None:
            p = domain.prime()
        if p != domain.prime():
            raise NotImplementedError("p-adic valuation not implemented for this ring")
        return pAdicValuation_padic(domain)
    else:
        raise NotImplementedError("p-adic valuation not implemented for this ring")

class pAdicValuation_base(UniqueRepresentation, DiscreteValuation):
    r"""
    Common base class for `p`-adic valuations.

    INPUT:

    - ``domain`` -- an integral domain whose elements provide a method
      ``valuation`` which takes ``prime`` as a parameter

    - ``prime`` -- a prime

    - ``check`` -- a boolean (default: ``True``) whether to perform checks on
      the parameters

    EXAMPLES::

        sage: pAdicValuation(ZZ, 3)
        3-adic valuation

        sage: pAdicValuation(QQ, 5)
        5-adic valuation

     For `p`-adic rings, ``prime`` has to match the prime of the ring.

        sage: v = pAdicValuation(Zp(3), 2); v
        Traceback (most recent call last):
        ...
        NotImplementedError: p-adic valuation not implemented for this ring

    TESTS::

        sage: TestSuite(pAdicValuation(ZZ, 3)).run()
        sage: TestSuite(pAdicValuation(QQ, 5)).run()
        sage: TestSuite(pAdicValuation(Zp(5), 5)).run()

    """
    @staticmethod
    def __classcall__(cls, domain, prime, check=True):
        r"""
        Normalize such that ``check`` is not part of the key for the
        :class:`sage.structure.factory.UniqueRepresentation` and the ``prime``
        is an integer.

        TESTS::

            sage: from sage.rings.padics.padic_valuation import pAdicValuation_int
            sage: pAdicValuation_int(ZZ, 3, check=False) is pAdicValuation_int(ZZ, 3, check=True)
            True

        """
        from sage.rings.all import ZZ
        prime = ZZ(prime)
        if check:
            if not prime.is_prime():
                raise ValueError("prime is not prime")
        return super(pAdicValuation_base, cls).__classcall__(cls, domain, prime)

    def __init__(self, domain, prime):
        """
        Initialization.

        TESTS::

            sage: type(pAdicValuation(ZZ, 3))
            <class 'sage.rings.padics.padic_valuation.pAdicValuation_int'>

        """
        DiscreteValuation.__init__(self, domain)
        self._prime = prime

    def prime(self):
        """
        Return the prime `p` of this `p`-adic valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).prime()
            3

        """
        return self._prime

    def _repr_(self):
        """
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)._repr_()
            '3-adic valuation'

        """
        return "%s-adic valuation"%(self.prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element in the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3).reduce(4)
            1

        """
        x = coerce(self.domain(), x)
        return self.residue_field()(x)

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: xbar = v.reduce(4)
            sage: v.lift(xbar)
            1

        """
        x = coerce(self.residue_field(), x)
        return self.domain()(x)

    def uniformizer(self):
        """
        Return a uniformizer of this valuation, i.e., `p` as an element of the domain.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.uniformizer()
            3

        """
        return self.domain()(self._prime)

    def residue_field(self):
        """
        Return the residue field of the ring of integers of this valuation.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.residue_field()
            Finite Field of size 3

        """
        from sage.rings.all import GF
        return GF(self._prime,names=('u',))

    def shift(self, c, v):
        """
        Multiply ``c`` by a ``v``th power of the uniformizer.

        INPUT:

        - ``c`` -- an element of the domain of this valuation

        - ``v`` -- an integer

        OUTPUT:

        If the resulting element has negative valuation, then it will be
        brought to the fraction field, otherwise it is an element of the
        domain.

        EXAMPLES::

            sage: v = pAdicValuation(ZZ, 3)
            sage: v.shift(2,3)
            54
            sage: v.shift(2,-3)
            2/27

        """
        from sage.rings.all import ZZ
        c = coerce(self.domain(), c)
        v = ZZ(v)

        c = self.domain().fraction_field()(c)
        ret = c*self.uniformizer()**v

        if self(c) >= -v:
            ret = self.domain()(ret)
        return ret

class pAdicValuation_padic(pAdicValuation_base):
    """
    The `p`-adic valuation of a `p`-adic ring.

    INPUT:

    - ``domain`` -- a `p`-adic ring

    EXAMPLES::

        sage: pAdicValuation(Qp(2)) #indirect doctest
        2-adic valuation

    """
    @staticmethod
    def __classcall__(cls, domain):
        r"""
        Infer the prime from the ``domain``.

        TESTS::

            sage: from sage.rings.padics.padic_valuation import pAdicValuation_padic
            sage: K = Qp(3)
            sage: pAdicValuation_padic(K) is pAdicValuation_padic(K)
            True

        """
        return super(pAdicValuation_padic, cls).__classcall__(cls, domain, domain.prime(), check=False)

    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT::

        - ``x`` --  an element in the domain of this valuation

        EXAMPLES::

            sage: K = Qp(3)
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^2 - 3)
            sage: v = pAdicValuation(L)
            sage: v(a)
            1/2

        """
        from sage.misc.functional import coerce
        x = coerce(self.domain(), x)

        return x.ordp(self.prime())

    def reduce(self, x):
        """
        Reduce ``x`` modulo the ideal of elements of positive valuation.

        INPUT:

        - ``x`` -- an element of the domain of this valuation

        OUTPUT:

        An element of the :meth:`residue_field`.

        EXAMPLES::

            sage: R = Zp(3)
            sage: pAdicValuation(Zp(3)).reduce(R(4))
            1

        """
        from sage.misc.functional import coerce
        x = coerce(self.domain(), x)

        return self.residue_field()(x.residue())

    def lift(self, x):
        """
        Lift ``x`` from the residue field to the domain of this valuation.

        INPUT:

        - ``x`` -- an element of the :meth:`residue_field`

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = pAdicValuation(R)
            sage: xbar = v.reduce(R(4))
            sage: v.lift(xbar)
            1 + O(3^20)

        """
        from sage.misc.functional import coerce
        x = coerce(self.residue_field(), x)

        return self.domain()(x).lift_to_precision()

    def uniformizer(self):
        """
        Return a uniformizer of this valuation, i.e., `p` as an element of the
        domain.

        EXAMPLES::

            sage: v = pAdicValuation(Zp(3))
            sage: v.uniformizer()
            3 + O(3^21)

        """
        return self.domain().uniformizer()

    def residue_field(self):
        """
        Return the residue field of the ring of integers of this valuation.

        EXAMPLES::

            sage: v = pAdicValuation(Zp(3))
            sage: v.residue_field()
            Finite Field of size 3

        """
        return self.domain().residue_field()

    def shift(self, c, v):
        """
        Multiply ``c`` by a ``v``th power of the uniformizer.

        INPUT:

        - ``c`` -- an element of the domain of this valuation

        - ``v`` -- an integer

        OUTPUT:

        If the resulting element has negative valuation, then it will be brought
        to the fraction field, otherwise it is an element of the domain.

        EXAMPLES::

            sage: R = Zp(3)
            sage: v = pAdicValuation(Zp(3))
            sage: v.shift(R(2),3)
            2*3^3 + O(3^23)

        """
        from sage.rings.all import ZZ
        c = coerce(self.domain(), c)
        v = ZZ(v)

        return c<<v

class pAdicValuation_int(pAdicValuation_base):
    def _call_(self, x):
        """
        Evaluate this valuation at ``x``.

        INPUT::

        - ``x`` --  an element in the domain of this valuation

        EXAMPLES::

            sage: pAdicValuation(ZZ, 3)(9)
            2

        """
        from sage.misc.functional import coerce
        x = coerce(self.domain(), x)

        return x.valuation(self.prime())
