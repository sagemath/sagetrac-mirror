r"""
Rivest, Shamir, Adleman (RSA) public-key encryption scheme.

The Rivest, Shamir, Adleman public-key encryption scheme.
The Rivest-Shamir-Adleman (RSA) scheme is a widely used approach to
public-key encryption.

REFERENCES:

 * David Ireland. *RSA Algorithm*. http://www.di-mgt.com.au/rsa_alg.html. 19 March 2015.

 * William Stallings. *Cryptography and Network Security, Priciples and Practices*, Fourth Edition. Prentice Hall, 16 November 2005.

 * Brian Raiter. *Prime Number Hide-and-Seek: How the RSA Cipher Works*. http://www.muppetlabs.com/~breadbox/txt/rsa.html. Accessed 20 March 2015.

 * Thomas Pornin. *Should RSA public exponent be only in {3, 5, 17, 257 or 65537} due to security considerations?*. http://security.stackexchange.com/questions/2335/should-rsa-public-exponent-be-only-in-3-5-17-257-or-65537-due-to-security-c

 * `Wikipedia: RSA Operation <http://en.wikipedia.org/wiki/RSA_(cryptosystem)#Operation>`_


AUTHORS:

 * Peter Story (2015-03) Major rewrite with the goal of a pedagogical implementation.
 * Ajeesh R (2011-07) Initial procedural version released as public domain software.

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
from sage.rings.arith import euler_phi
from sage.rings.arith import gcd
from sage.rings.arith import is_prime
from sage.rings.arith import power_mod
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod import Mod as mod
from sage.rings.integer_ring import ZZ
from sage.sets.primes import Primes

class RSACryptosystem(PublicKeyCryptosystem):
    r"""
    The Rivest, Shamir, Adleman (RSA) public-key encryption scheme.

    Note: Throughout the documentation, the Euler Phi function is represented
    as `\phi`.

    EXAMPLES:

    The following is an encryption/decryption example.::

        sage: rc = RSACryptosystem(); rc
        The Rivest, Shamir, Adleman public-key encryption scheme.
        sage: p = 101; q = 103
        sage: public, private = rc.calculate_keys(p, q); public, private
        ((10403, 257), (10403, 5993))
        sage: P = 9001
        sage: C = rc.encrypt(P, public); C
        6022
        sage: M = rc.decrypt(C, private); M
        9001
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

            sage: rc = RSACryptosystem()
            sage: rc == loads(dumps(rc))
            True
        """
        # In order to implement PublicKeyCryptosystem, it is important to pass: 
        #   self, plaintext_space, ciphertext_space, 
        #   key_space, block_length=1, period=None
        # Currently unsure of what the appropriate values are.
        #PublicKeyCryptosystem.__init__(self, S, S, S)
        pass

    def __eq__(self, other):
        """
        Compare this ``RSACryptosystem`` object with ``other``.

        INPUT:

        - ``other`` -- a ``RSACryptosystem`` object.

        OUTPUT:

        - ``True`` if both ``self`` and ``other`` are ``RSACryptosystem``
          objects. ``False`` otherwise.

        Since ``RSACryptosystem`` objects do not have state, any two
        ``RSACryptosystem`` objects are equivalent. 

        EXAMPLES::

            sage: rc1 = RSACryptosystem()
            sage: rc2 = RSACryptosystem()
            sage: rc1 == rc2
            True

            sage: rc1 = RSACryptosystem()
            sage: rc_str = "The Rivest, Shamir, Adleman public-key encryption scheme."
            sage: rc1 == rc_str
            False
        """
        if (self.__class__ == other.__class__):
            return True
        else:
            return False

    def __repr__(self):
        """
        A string representation of this RSA public-key encryption scheme.

        OUTPUT:

        - A string representation of this RSA public-key encryption scheme.

        EXAMPLES::

            sage: RSACryptosystem()
            The Rivest, Shamir, Adleman public-key encryption scheme.
        """
        return "The Rivest, Shamir, Adleman public-key encryption scheme."

    def decrypt(self, C, K):
        r"""
        Apply the RSA public-key encryption scheme to decrypt the ciphertext ``C``
        using the private key ``K``.

        INPUT:

        - ``C`` -- an integer ciphertext resulting from encrypting a plaintext
          using the RSA public-key encryption algorithm.

        - ``K`` -- a private key, stored as a tuple ``(n, d)`` where ``n``
          is the modulus and ``d`` is inverse of the exponent mod `(\phi(n))`.

        OUTPUT:

        - The plaintext resulting from decrypting the ciphertext ``C`` using
          the RSA public-key decryption algorithm.

        ALGORITHM:

        The RSA public-key decryption algorithm is described as follows:

        #. Let `C` be the ciphertext `C = (M^e) \text{ mod } (n)`.
        #. Let `(n, d)` be the private key whose corresponding
           public key is `(n, e)`.
        #. The plaintext is `M = C^d = (M^e)^d \text{ mod } (n)`.

        EXAMPLES:

        An example of decryption.::

            sage: rc = RSACryptosystem()
            sage: p = 101; q = 103
            sage: public, private = rc.calculate_keys(p, q)
            sage: P = 42
            sage: C = rc.encrypt(P, public); C
            6270
            sage: M = rc.decrypt(C, private); M; P == M
            42
            True

        TESTS:

        Setup for tests.::

            sage: rc = RSACryptosystem()
            sage: p = 101; q = 103; n = p * q
            sage: public, private = rc.calculate_keys(p, q)

        Decryption should recover the original message.::

            sage: P = 42
            sage: C = rc.encrypt(P, public)
            sage: rc.decrypt(C, private) == P
            True
            sage: P = ZZ.random_element(2, n)
            sage: C = rc.encrypt(P, public)
            sage: rc.decrypt(C, private) == P
            True

        Ciphertext outside the range ``[2, n-1]`` should be rejected.::

            sage: C = -5; rc.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The ciphertext must be an integer greater than 1 and less than n.

            sage: C = 0; rc.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The ciphertext must be an integer greater than 1 and less than n.

            sage: C = 1; rc.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The ciphertext must be an integer greater than 1 and less than n.

            sage: C = n; rc.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The ciphertext must be an integer greater than 1 and less than n.

            sage: C = n + 5; rc.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The ciphertext must be an integer greater than 1 and less than n.


        Ciphertext in the range ``[2, n-1]`` should be accepted.::

            sage: C = 2; M = rc.decrypt(C, private)
            sage: C = n - 1; M = rc.decrypt(C, private)
            sage: C = ZZ.random_element(2, n); M = rc.decrypt(C, private)

        """
        # private key
        (n, d) = K

        # sanity check
        if (C < 2) or (C >= n):
            raise ValueError("The ciphertext must be an integer greater than 1 and less than n.")

        # perform the decryption
        return power_mod(C, d, n)

    def encrypt(self, P, K):
        r"""
        Apply the RSA public-key encryption scheme to encrypt the plaintext ``P`` using
        the public key ``K``.

        INPUT:

        - ``P`` -- a number in the range ``[2, n-1]``. Note that a secure
          implementation of the cryptosystem would add padding to this
          plaintext, in order to protect against a number of attacks.

        - ``K`` -- a public key, stored as a tuple ``(n, e)`` where ``n``
          is the modulus and ``e`` is the exponent.

        OUTPUT:

        - The ciphertext resulting from encrypting ``P`` using the public
          key ``K``.

        ALGORITHM:

        The RSA public-key encryption algorithm is described as follows:

        #. Let `(n, e)` be a public key, where `n = pq` is the product of two
           distinct primes `p` and `q` and exponent `e` is a positive 
           integer that is co-prime to `\phi(n)`.
        #. Let `P = a string` be the message (plaintext).
        #. The ciphertext is `C = (P^e) \text{ mod } (n)`.


        EXAMPLES:

        An example of encryption. 
        ``P = 42`` is our plaintext, and ``C = 6270`` is the resulting ciphertext.::

            sage: rc = RSACryptosystem()
            sage: p = 101; q = 103
            sage: public, private = rc.calculate_keys(p, q)
            sage: P = 42
            sage: C = rc.encrypt(P, public); C
            6270


        If the system allowed us to encrypt values greater than ``n-1``, our
        message could be encrypted to the same value a different messsage.
        ``P_1`` is our first plaintext, and ``P_2`` is a plaintext greater than 
        ``n - 1``, so it isn't encypted to a distinct value.::

            sage: rc = RSACryptosystem()
            sage: p = 101; q = 103; n = p * q
            sage: public, private = rc.calculate_keys(p, q)
            sage: e = public[1]
            sage: P_1 = 42
            sage: P_2 = 42 + n
            sage: C_1 = power_mod(P_1, e, n); C_1
            6270
            sage: C_2 = power_mod(P_2, e, n); C_2
            6270
            sage: C_1 == C_2
            True


        TESTS:

        Setup for tests.::

            sage: rc = RSACryptosystem()
            sage: p = 101; q = 103; n = p * q
            sage: public, private = rc.calculate_keys(p, q)

        Plaintext outside the range ``[2, n-1]`` should be rejected.::

            sage: P = -5; rc.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than n.

            sage: P = 0; rc.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than n.

            sage: P = 1; rc.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than n.

            sage: P = n; rc.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than n.

            sage: P = n + 5; rc.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than n.


        Plaintext in the range ``[2, n-1]`` should be accepted.::

            sage: P = 2; C = rc.encrypt(P, public)
            sage: P = n - 1; C = rc.encrypt(P, public)
            sage: P = ZZ.random_element(2, n); C = rc.encrypt(P, public)

        """
        # public key
        n, e = K

        # sanity check
        if (P < 2) or (P >= n):
            raise ValueError("The plaintext must be an integer greater than 1 and less than n.")

        # Encrypt and return the message
        return power_mod(P, e, n)

    def choose_exponent(self, phi_n):
        """
        Choose an exponent relatively prime to the Euler
        Characteristic of ``n``.

        Using a small number is not a security problem. For pedagogical
        purposes, we include two different ways of finding ``e``.

        INPUT:

        - ``phi_n`` -- the Euler Characteristic of ``n``

        OUTPUT:

        - A number relatively prime to ``phi_n``


        EXAMPLES:

        Get a number relatively prime to ``phi_n``::

            sage: rc = RSACryptosystem()
            sage: phi_n = 11020
            sage: e = rc.choose_exponent(phi_n)
            sage: gcd(phi_n, e) == 1
            True

        TESTS:

        The returned value must be relatively prime to ``phi_n``::

            sage: rc = RSACryptosystem()
            sage: phi_n = ZZ.random_element(50, 100)
            sage: e = rc.choose_exponent(phi_n)
            sage: gcd(phi_n, e) == 1
            True


        Exercise the second method of finding a relatively prime exponent.::

            sage: rc = RSACryptosystem()
            sage: e_options = [3, 5, 17, 257, 65537]
            sage: phi_n = prod(e_options)
            sage: e = rc.choose_exponent(phi_n)
            sage: gcd(phi_n, e) == 1
            True
            sage: e not in e_options
            True

        """
        # Several recommendations for e, according to the SSL specification
        e_options = [3, 5, 17, 257, 65537]
        # Seeking an e such that phi_n and e are coprime
        for e_option in e_options:
            # Since each e_option is a prime, we can check efficiently
            # check coprimeness; the only way they wouldn't be coprime would
            # be if phi_n was a multiple of the e_option
            if ((phi_n % e_option) != 0): # They are coprime
                return e_option

        # It is unlikely that none of the e_options will work, but just in case
        # we provide the option to seek additional values.
        P = Primes()
        p = P.first()
        while (phi_n % p == 0):
            p = P.next(p)
        return p

    def calculate_keys(self, p, q, enforce_primes=True):
        r"""
        Return an RSA public and private key pair corresponding to the
        distinct primes ``p`` and ``q``.

        INPUT:

        - ``p`` -- a prime number.

        - ``q`` -- a prime number.

        OUTPUT:

        - The RSA public and private keys as a tuple ``((n, e), (n, d))``.


        EXAMPLES:

        Compute an RSA public and private key pair corresponding to the 
        supplied primes. The ``d`` in the private key should be the 
        multiplicative inverse of `e \text{ mod } \phi(n)`::

            sage: rc = RSACryptosystem()
            sage: public, private = rc.calculate_keys(19, 23); public, private
            ((437, 5), (437, 317))
            sage: n = public[0]; e = public[1]; d = private[1]
            sage: (e * d) % euler_phi(n)
            1

        For teaching purposes, you can supply two numbers that are not primes.
        This is also helpful if you are using very large primes, to avoid the
        expensive computation of checking primeness.
        In this example, we use prime ``17`` and Carmichael number
        ``561 = 3 * 11 * 17`` to generate the keys. Note that decrypting
        doesn't recover the original plaintext.::

            sage: rc = RSACryptosystem()
            sage: public, private = rc.calculate_keys(17, 561, False); public, private
            ((9537, 3), (9537, 2987))
            sage: P = 1234
            sage: C = rc.encrypt(P, public); C
            5794
            sage: M = rc.decrypt(C, private); M
            5722

        In this example we use ``p`` and ``q`` which are relatively prime, but
        ``p`` is not a prime. In this case, only certain plaintexts can be 
        recovered.::

            sage: rc = RSACryptosystem()
            sage: public, private = rc.calculate_keys(4, 19, False); public, private
            ((76, 5), (76, 11))
            sage: P = 10
            sage: M = rc.decrypt(rc.encrypt(P, public), private); M
            48
            sage: P = 24
            sage: M = rc.decrypt(rc.encrypt(P, public), private); M
            24

        However, it is also possible to be build a cryptosystem using a 
        composite which works for all plaintexts.::

            sage: rc = RSACryptosystem()
            sage: public, private = rc.calculate_keys(5, 6, False)
            sage: n = public[0]
            sage: fail_count = 0
            sage: for P in range(2, n):
            ....:     if P != rc.decrypt(rc.encrypt(P, public), private):
            ....:         fail_count += 1
            sage: fail_count
            0

        TESTS:

        Both of the inputs ``p`` and ``q`` must be distinct.::

            sage: rc = RSACryptosystem()
            sage: rc.calculate_keys(23, 23)
            Traceback (most recent call last):
            ...
            ValueError: p and q must be distinct.


        Both of the inputs ``p`` and ``q`` must be prime.::

            sage: rc = RSACryptosystem()
            sage: rc.calculate_keys(13, 21)
            Traceback (most recent call last):
            ...
            ValueError: p and q must be primes.

        Allow ``p`` or ``q`` to be composite if ``enforce_primes`` is set
        to ``False``.::

            sage: rc = RSACryptosystem()
            sage: rc.calculate_keys(13, 21, False)
            ((273, 17), (273, 113))

        The ``d`` in the private key should be the multiplicative inverse
        of `e \text{ mod } (\phi(n))`::

            sage: rc = RSACryptosystem()
            sage: p = next_prime(ZZ.random_element(50, 100)); q = next_prime(ZZ.random_element(150, 200))
            sage: public, private = rc.calculate_keys(p, q)
            sage: n = public[0]; e = public[1]; d = private[1]
            sage: (e * d) % euler_phi(n)
            1

        """
        if p == q:
            raise ValueError("p and q must be distinct.")
        # Unsure of the following comment: 
        # Unfortunately, we cannot enforce the primeness of ``p`` and ``q``, 
        # since primeness cannot be efficiently verified for large primes.::
        if enforce_primes and ((not is_prime(p)) or (not is_prime(q))):
            raise ValueError("p and q must be primes.")

        n = p * q
        # Compute the Euler Characteristic of ``n``
        # We can compute this efficiently because we know ``p`` and ``q``
        phi_n = (p - 1) * (q - 1)

        # Chooses an ``e`` relative prime to ``phi_n``
        e = self.choose_exponent(phi_n)

        # Generate a number ``d`` such that ``de = 1 mod (\phi(n))``
        # ``d`` serves as our decryption key, contained in the private key; 
        # it is the multiplicative inverse of ``e``

        # Get the multiplicative inverse of ``e mod (\phi(n))``
        d = power_mod(e, -1, phi_n)

        private_key = (n, d)
        public_key = (n, e)
        return (public_key, private_key)
