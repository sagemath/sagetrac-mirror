r"""
Rivest, Shamir, Adleman public-key encryption scheme.

The Rivest, Shamir, Adleman public-key encryption scheme. The Rivest-Shamir-Adleman (RSA) scheme has since that time reigned supreme as the most widely accepted and implemented general-purpose approach to public-key encryption. See also the `Wikipedia article <http://en.wikipedia.org/wiki/RSA>`_ on this scheme.

REFERENCES:

.. William Stallings. *Cryptography and Network Security, Priciples and Practices*, Fourth Edition. Prentice Hall, 16 November 2005.

.. Anoop MS. *Public Key Cryptography, Applications Algorithms and Mathematical Explanations*. Tata Elxsi Ltd, India anoopms@tataelxsi.co.in.

.. Minh Van Nguyen. *Number Theory and the RSA Public Key Cryptosystem*. nguyenminh2@gmail.com, July 2009.

AUTHORS:

-Ajeesh R (2011-07)- initial procedural version released as public domain software.

"""

###########################################################################
# Copyright (c) 2011
# AJEESH R <iamajeeshr@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################

from sage.crypto.cryptosystem import PublicKeyCryptosystem
from sage.rings.arith import gcd
from sage.rings.arith import xgcd
from sage.rings.arith import is_prime
from sage.rings.arith import power_mod
from sage.rings.arith import euler_phi
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod import Mod as mod
from sage.rings.integer_ring import ZZ

class RSACryptosystem(PublicKeyCryptosystem):
    r"""
    The Rivest, Shamir, Adleman public-key encryption scheme.

    The RSA encryption and decryption algorithms as described in
    :func:`encrypt() <RSACryptosystem.encrypt>` and
    :func:`decrypt() <RSACryptosystem.decrypt>`, respectively, make use of the

    EXAMPLES:

    The following is an encryption/decryption example.::

        sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
        sage: rc = RSACryptosystem(); rc
        The Rivest, Shamir, Adleman public-key encryption scheme.
        sage: p = 31; q = 61;
        sage: pubkey = rc.public_key(p, q); pubkey
        (4951760154835678088235319297L, 17)
        sage: prikey = rc.private_key(p, q); prikey
        (4951760154835678088235319297L, 4077920125612805357425763753)
        sage: P = 72697676798779827668
        sage: C = rc.encrypt(P, pubkey); C
        2467165704948727396791981601
        sage: M = rc.decrypt(C, prikey); M
        72697676798779827668
        sage: M == P
        True

    Generate a pair of public/private keys. Use the public key to
    encrypt a plaintext. Then decrypt the resulting ciphertext using the
    private key. Finally, compare the decrypted message with the original
    plaintext.::

        sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
        sage: rc = RSACryptosystem()
        sage: P = 72697676798779827668
        sage: C = rc.encrypt(P, pubkey)
        sage: M = rc.decrypt(C, prikey)
        sage: M == P
        True
    """

    def __init__(self):
        """
        Construct the RSA public-key encryption scheme.

        OUTPUT:

        - A ``RSACryptosystem`` object representing the RSA public-key encryption scheme.

        See the class docstring of ``RSACryptosystem`` for detailed documentation.

        EXAMPLES::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: rc == loads(dumps(rc))
            True
        """
        # no internal data for now; nothing to initialize
        pass

    def __eq__(self, other):
        """
        Compare this ``RSACryptosystem`` object with ``other``.

        INPUT:

        - ``other`` -- a ``RSACryptosystem`` object.

        OUTPUT:

        - ``True`` if both ``self`` and ``other`` are ``RSACryptosystem``
          objects. ``False`` otherwise.

        Two objects are ``RSACryptosystem`` objects if their string
        representations are the same.

        EXAMPLES::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc1 = RSACryptosystem()
            sage: rc2 = RSACryptosystem()
            sage: rc1 == rc2
            True
        """
        if self.__repr__() == other.__repr__():
            return True
        else:
            return False

    def __repr__(self):
        """
        A string representation of this RSA public-key encryption scheme.

        OUTPUT:

        - A string representation of this RSA public-key encryption scheme.

        EXAMPLES::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: RSACryptosystem()
            The Rivest, Shamir, Adleman public-key encryption scheme.
        """
        return "The Rivest, Shamir, Adleman public-key encryption scheme."

    def decrypt(self, C, K):
        r"""
        Apply the RSA public-key encryption scheme to decrypt the ciphertext ``C``
        using the private key ``K``.

        INPUT:

        - ``C`` -- a ciphertext resulting from encrypting a plaintext using
          the RSA public-key encryption algorithm.

        - ``K`` -- a private key `(n, d)` where `n` is the product of two Mersenne
          primes `p` and `q`, `d` is computed using the extended euclidean algorithm.

        OUTPUT:

        - The plaintext resulting from decrypting the ciphertext ``C`` using
          the RSA public-key decryption algorithm.

        ALGORITHM:

        The RSA public-key decryption algorithm is described as follows:

        #. Let `C` be the ciphertext `C = (M^e) mod n`.
        #. Let `(n, d)` be the private key whose corresponding
           public key is `(n, e)`.
        #. The plaintext is `M = (C^d) mod n`.

        EXAMPLES:

        The following is a decryption example. Here we decrypt a string of numbers.::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: p = 31; q = 61;
            sage: C = 2467165704948727396791981601
            sage: K = rc.private_key(p, q); K
            (4951760154835678088235319297L, 4077920125612805357425763753)
            sage: P = rc.decrypt(C, K); P
            72697676798779827668

        TESTS:

        The private key `K = (n, d)` must be such that `p` and `q` are
        distinct prime numbers.::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: C = 2467165704948727396791981601
            sage: K = rc.private_key(31, 31);
            Traceback (most recent call last):
            ...
            ValueError: p and q must be distinct primes.
            sage: K = rc.private_key(31, 61)
            sage: P =  rc.decrypt(C, K); P
            72697676798779827668
        """
        # private key
        (n, d) = K

        # perform the decryption
        return power_mod(C, d, n)

    def encrypt(self, P, K):
        r"""
        Apply the RSA public-key encryption scheme to encrypt the plaintext ``P`` using
        the public key ``K``.

        INPUT:

        - ``P`` -- a non-empty string of plaintext. The string ``""`` is
          an empty string, whereas ``" "`` is a string consisting of one
          white space character. The plaintext can be a string of numbers or
          a string of ASCII characters.

        - ``K`` -- a public key, which is the product of two Mersenne primes and e, 0 < e < euler_phi(n)
          such that also satisfy the requirement `\gcd(e, euler_phi(n))=1`.

        OUTPUT:

        - The ciphertext resulting from encrypting ``P`` using the public
          key ``K``.

        ALGORITHM:

        The RSA public-key encryption algorithm is described as follows:

        #. Let `(n, e)` be a public key, where `n = pq` is the product of two
           distinct Mersenne primes `p` and `q` and `e` is a random positive 
           integer that is co-prime to euler_phi(n).
        #. Let `P = a string` be the message (plaintext).
        #. The ciphertext is `C = (P^e) mod n`.


        EXAMPLES:

        The following is an encryption example. Here, we encrypt a string of number.::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: p = 101; q = 103;
            sage: P = 4
            sage: C = rc.encrypt(P, K); C
            4932417489295796679844298689

        Now encrypt another string of numbers. The result is random; no seed is
        provided to the encryption function so the function generates its
        own random seed.::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: p = 31; q = 61;
            sage: P = 72697676798779827668
            sage: K = rc.public_key(p, q);
            sage: C = rc.encrypt(P, K);
            sage: C
            2467165704948727396791981601

        TESTS:

        The plaintext cannot be an empty string. ::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: rc.encrypt("", 3)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext cannot be an empty string.
        """
        # public key
        (n, e) = K

        # sanity check
        if P == "":
            raise ValueError("The plaintext cannot be an empty string.")

        else:
            return power_mod(P, e, n)

    def private_key(self, p, q):
        r"""
        Return the RSA private key corresponding to the
        distinct Mersenne primes ``p`` and ``q``.

        INPUT:

        - ``p`` -- a prime number.

        - ``q`` -- a prime number.

        OUTPUT:

        - The RSA private key `(n, d)` where
          `de is congruent to (1 mod euler_phi(n))`.

        Both ``p`` and ``q`` must be distinct primes. Let `p` be a
        positive prime. Then `p` is a Mersenne prime if `p` is equal to `(2^p)-1`.

        EXAMPLES:

        Obtain two distinct primes and compute the RSA
        private key corresponding to two Mersenne primes generated from the distinct primes.::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: rc.private_key(19, 23)
            (4398037598209L, 1173935455331)

        Choose two distinct primes, compute the RSA
        private key corresponding to two Mersenne primes generated from the distinct primes , and test that the
        resulting private key `(n, d)` satisfies
        `\mod(de, euler_phi(n)) = 1`.::

            sage: p = 31
            sage: q = 61
            sage: p = (2 ^ p)-1; q = (2 ^ q)-1
            sage: is_prime(p); is_prime(q)
            True
            True

        TESTS:

        Both of the input ``p`` and ``q`` must be distinct primes. ::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: rc.private_key(78307, 78307)
            Traceback (most recent call last):
            ...
            ValueError: p and q must be distinct primes.
        """
        if p == q:
            raise ValueError("p and q must be distinct primes.")
        if is_prime(p) and is_prime(q):
            p = (2**p)-1
            q = (2**q)-1

            n = p * q

            e = 3
            while 1:
                if gcd(euler_phi(n), e) == 1:
                    break
                else:
                    e = e + 2

            # Generate a number d such that de = 1 (mod euler_phi(n))
            bezout = xgcd(e, euler_phi(n))
            d = Integer(mod(bezout[1], euler_phi(n)))
            while mod(d * e, euler_phi(n)) != 1:
                d = Integer(mod(bezout[1], euler_phi(n)))
        else:
            raise ValueError("p and q must be distinct primes.")
        return (n, d)

    def public_key(self, p, q):
        r"""
        Return RSA public key corresponding to the
        distinct Mersenne primes ``p`` and ``q``.

        INPUT:

        - ``p`` -- a prime number.

        - ``q`` -- a prime number.

        OUTPUT:

        - The RSA public key `n = pq` and `e`.

        Both ``p`` and ``q`` must be distinct primes. Let `p` be a
        positive prime. Then `p` is a Mersenne prime if `p` is equal to `(2^p)-1`.

        EXAMPLES:

        Obtain two distinct Mersenne primes and compute the RSA
        public key corresponding to those two Mersenne primes::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: rc.public_key(31, 61)
            (4951760154835678088235319297L, 17)

        Choose two distinct primes, compute the RSA
        public key corresponding to two Mersenne primes generated from the distinct primes.::

            sage: p = 31
            sage: q = 61
            sage: p = (2 ^ p)-1; q = (2 ^ q)-1
            sage: is_prime(p); is_prime(q)
            True
            True

        TESTS:

        The input ``p`` and ``q`` must be distinct primes. ::

            sage: from sage.crypto.public_key.rsa_cryptosystem import RSACryptosystem
            sage: rc = RSACryptosystem()
            sage: rc.public_key(3, 3)
            Traceback (most recent call last):
            ...
            ValueError: p and q must be distinct primes.
        """
        if p == q:
            raise ValueError("p and q must be distinct primes.")
        if is_prime(p) and is_prime(q):
            p = (2**p)-1
            q = (2**q)-1

            n = p * q
            # Generate a number e so that gcd(euler_phi(n), e) = 1
            e = 3
            while 1:
                if gcd(euler_phi(n), e) == 1: 
                    break
                else:
                    e = e + 2
        else:
            raise ValueError("p and q must be distinct primes.")
        return (n, e)
