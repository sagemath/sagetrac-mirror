# coding: UTF-8
r"""
Sharmir Secret Sharing

Implements the original versions of perfectly secure secret sharing 
as proposed by Shamir in [Shamir1979]_. Note that this code is for educational 
purposes only and demonstrate the basic algorithms. 

AUTHORS:

- Thomas Loruenser (2013): initial version

REFERENCES:

.. [Shamir1979] Shamir, A. (1979). How to share a secret. 
   Communications of the ACM, 22(11), 612–613. :doi:`10.1145/359168.359176`
"""
###############################################################################
# Copyright 2013, Thomas Loruenser <thomas.loruenser@ait.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

from sage.structure.sage_object import SageObject

class ShamirSS(SageObject):
    r"""
    Shamir secret sharing.
   
    This class implements the original version of perfectly secure secret sharing 
    as proposed by Shamir in [Shamir1979]_. It is a very basic implementation
    intended for educational purposes only.

    INPUT:

    - ``k``  --  (default: ``3``) the threshold for reconstruction.
    - ``n``  --  (default: ``7``) the number of shares.
    - ``order`` --  (default: ``2^8``) field order for data to share.

    EXAMPLES::

        sage: from sage.crypto.smc.shamir_ss import ShamirSS

    Generate shares::

        sage: k = 3; n = 7
        sage: sss = ShamirSS(n,k)
        sage: secret = 42
        sage: shares = sss.share(secret)

    Reconstruct secret (Lagrange interpolation)::

        sage: secret == sss.reconstruct(shares)
        True

    Reconstruct with error shares (Belekamp-Welsch decoder)::

        sage: shares[0] = (shares[0][0], shares[0][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        True

    Input vectors::

        sage: sss = ShamirSS()
        sage: secret = [42, 43, 44, 45]
        sage: shares = sss.share(secret)
        sage: secret == sss.reconstruct(shares)
        True

    TESTS:

    More random input::

        sage: secret = randint(0,255)
        sage: shares = sss.share(secret)
        sage: secret == sss.reconstruct(shares, decoder='lg')
        True
        sage: secret == sss.reconstruct(shares, decoder='bw')
        True

    Secret must be integer representation in field::

        sage: secret = 333
        sage: shares = sss.share(secret)
        Traceback (most recent call last):
        ...
        TypeError: secret must be within 0 and field order.

    Working in prime fields::

        sage: from sage.rings.arith import random_prime

        sage: order = random_prime(10**10)
        sage: secret = randint(0, order-1)
        sage: sss = ShamirSS(7, 3, order)
        sage: shares = sss.share(secret)
        sage: secret == sss.reconstruct(shares)
        True

        sage: shares[0] = (shares[0][0], shares[0][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        True

        sage: shares[1] = (shares[1][0], shares[1][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        True

        sage: shares[-1] = (shares[-1][0], shares[-1][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        False

    Working in extension fields::

        sage: order = 2**randint(3,15)
        sage: secret = randint(0, order)
        sage: sss = ShamirSS(7, 3, order)
        sage: shares = sss.share(secret)
        sage: secret == sss.reconstruct(shares[:-2])
        True

        sage: shares[0] = (shares[0][0], shares[0][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        True

        sage: shares[1] = (shares[1][0], shares[1][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        True

        sage: shares[-1] = (shares[-1][0], shares[-1][1]+1)
        sage: secret == sss.reconstruct(shares, decoder='bw')
        False
    """
    def __init__(self, n=7, k=3, order=2**8):
        r"""
        Sharmir secret sharing.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS
            sage: sss = ShamirSS()
            sage: secret = 42
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares)
            True
        """
        self._k = k  # threshold
        self._n = n  # number shares
        self._order = order  # order of field

        from sage.rings.finite_rings.constructor import FiniteField
        self._F = FiniteField(self._order, 'a')
        if not self._F.is_prime_field() and not hasattr(self._F, 'fetch_int'):
            raise TypeError("field order not supported")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        self._P = PolynomialRing(self._F, 'x')

    ### begin module private api

    def _latex_(self):
        r"""
        Return Latex representation of self.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS
            sage: sss=ShamirSS()
            sage: latex(sss)
            `(7,3)`-Shamir secret sharing over the field `\Bold{F}_{2^{8}}`
        """
        from sage.misc.latex import latex
        return "`({},{})`-Shamir secret sharing over the field `{}`".format(
            self._n, self._k, latex(self._F))


    def _to_GF(self, x):
        r""" 
        Convert integer representation to finite field

        INPUT:

        - ``x`` --  the integer representation to be converted.

        OUTPUT:

        The finite field representation.

        EXAMPLES::
        
            sage: from sage.crypto.smc.shamir_ss import ShamirSS
            sage: sss = ShamirSS()
            sage: sss._to_GF(42)
            a^5 + a^3 + a
            sage: sss._to_GF(255)
            a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1
        """
        # input checking
        if x > self._F.order():
            raise TypeError("secret must be within 0 and field order.")

        # convert to field type
        if self._F.is_prime_field():
            return self._F(x)
        else:
            return self._F.fetch_int(x)


    def _to_Int(self, x):
        r""" 
        Convert field representation to integer

        INPUT:

        - ``x`` -- the field element to be converted.

        OUTPUT:

        The integer representation of the field element.

        EXAMPLES::
        
            sage: from sage.crypto.smc.shamir_ss import ShamirSS
            sage: sss = ShamirSS()
            sage: secret = 42
            sage: test = sss._to_Int(sss._to_GF(42))
            sage: test == secret
            True
        """
        from sage.rings.all import Integer
        if self._F.is_prime_field():
            return Integer(x)
        else:
            return x.integer_representation()


    def _rec_berlekamp_welsh(self, points):
        r"""
        Reconstruct with Berlekamp-Welsh decoding.

        INPUT:

        - ``points`` -- shares as list of tuples.

        OUTPUT:

        Reconstructed secret.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS

        Decoding with errors::

            sage: sss = ShamirSS()
            sage: secret = 42
            sage: shares = sss.share(secret)
            sage: shares[0] = (shares[0][0], shares[0][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

            sage: k = 4; n = 10
            sage: sss = ShamirSS(n,k)
            sage: secret = 84
            sage: shares = sss.share(secret)
            sage: shares[0] = (shares[0][0], shares[0][1]+1)
            sage: shares[1] = (shares[1][0], shares[1][1]+1)
            sage: shares[-1] = (shares[-1][0], shares[-1][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

        """
        from berlekamp_welsh import berlekamp_welsh
        polycoeffs =  berlekamp_welsh(self._k-1, points).coeffs()
        return polycoeffs


    def _rec_lagrange(self, points):
        r""" 
        Reconstruct with Lagrange interpolation.

        INPUT:

        - ``points`` --  shares as list of tuples.

        OUTPUT:

        Reconstructed secret.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS

        Erasure decoding::

            sage: sss = ShamirSS()
            sage: secret = 42
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares)
            True

            sage: k = 4; n = 10
            sage: sss = ShamirSS(n,k)
            sage: secret = 42
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares)
            True
        """
        polycoeffs = self._P.lagrange_polynomial(points).coeffs()
        if len(polycoeffs) != self._k:
            raise ValueError("lagrange polynomial degree mismatch.")
        return polycoeffs


    def _repr_(self):
        r"""
        Return String representation of self.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS
            sage: sss=ShamirSS()
            sage: print(sss)
            (7,3)-Shamir secret sharing over Finite Field in a of size 2^8
        """
        return "({},{})-Shamir secret sharing over {}".format(self._n, self._k, 
                                                              self._F)
    ### begin public api

    def reconstruct(self, shares, decoder='lg'):
        r"""
        Reconstruct shares.

        INPUT:

        - ``shares`` -- a list of shares ((x,y)-tuples of integer) or list of it.
        - ``decoder`` -- (default: ``'lg'``) decoder used to reconstruct secret. Must
            be one of the supported types ``'lg'`` or ``'bw'``.

        OUTPUT:

        The reconstructed secret or list of secrets.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS

        Simple interface::

            sage: sss = ShamirSS()
            sage: secret = 42
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares)
            True

        Decoding with errors::

            sage: shares[0] = (shares[0][0], shares[0][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

            sage: k = 4; n = 10
            sage: sss = ShamirSS(n,k)
            sage: secret = randint(0, 255)
            sage: shares = sss.share(secret)
            sage: shares[0] = (shares[0][0], shares[0][1]+1)
            sage: shares[1] = (shares[1][0], shares[1][1]+1)
            sage: shares[-1] = (shares[-1][0], shares[-1][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

        TESTS:

        Working in prime fields::

            sage: from sage.rings.arith import random_prime
            sage: order = random_prime(10**10)
            sage: secret = randint(0, order-1)
            sage: sss = ShamirSS(7, 3, order)
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares)
            True

            sage: shares[0] = (shares[0][0], shares[0][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

            sage: shares[1] = (shares[1][0], shares[1][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

            sage: shares[-1] = (shares[-1][0], shares[-1][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            False

        Working in extension fields::

            sage: order = 2**randint(3, 15)
            sage: secret = randint(0, order-1)
            sage: sss = ShamirSS(7, 3, order)
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares[:-2])
            True

            sage: shares[0] = (shares[0][0], shares[0][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

            sage: shares[1] = (shares[1][0], shares[1][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            True

            sage: shares[-1] = (shares[-1][0], shares[-1][1]+1)
            sage: secret == sss.reconstruct(shares, decoder='bw')
            False
        """
        # make shares iterable
        if type(shares[0]) == tuple:
            shares = [shares]

        # set decoder
        if decoder == 'lg':
            decode = self._rec_lagrange
        elif decoder == 'bw':
            decode = self._rec_berlekamp_welsh
        else:
            raise ValueError("unknown decoder.")

        # reconstruct secret
        secret = []
        for element in shares:
            # convert to field
            points = [(self._to_GF(x), self._to_GF(y)) for x,y in element]
            # call decoder
            secret.append(self._to_Int(decode(points)[0]))
        if len(secret) == 1:
            secret = secret[0]
        return secret


    def share(self, secret):
        r"""
        Generate shares.

        A polynomial of degree `k-1` is generated at random with the secret
        being the constant coefficient. It is then evaluated at points starting
        from `1`.

        INPUT:

        - ``secret`` -- the secret to be shared as integer or list of integer.

        OUTPUT:

        The shares or a list of shares, if list input.

        EXAMPLES::

            sage: from sage.crypto.smc.shamir_ss import ShamirSS

        Simple interface::

            sage: sss = ShamirSS()
            sage: secret = 42
            sage: shares = sss.share(secret)
            sage: [i+1 == share[0]  for i, share in enumerate(shares)]
            [True, True, True, True, True, True, True]

        Input vector::

            sage: sss = ShamirSS()
            sage: secret = [42, 43, 44, 45]
            sage: shares = sss.share(secret)
            sage: secret == sss.reconstruct(shares)
            True
        """
        # make input iterable
        from sage.rings.integer import Integer
        if not type(secret) == list:
            secret = [secret]
        
        # generate shares
        shares = []
        for s in secret:
            # random polynomial with s as constant coefficient
            ssp = self._to_GF(s)
            for i in range(1, self._k):
                ssp += self._F.random_element() * self._P.gen()**i

            # evaluate polynomial at different points (shares)
            shares.append([(i, self._to_Int(ssp(self._to_GF(i)))) for i in range(1, self._n+1)])
        
        if len(shares) == 1:
            shares = shares[0]
        return shares


# vim: set fileencoding=UTF-8 filetype=python :
