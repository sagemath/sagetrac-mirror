r"""
ElGamal public-key encryption scheme

ElGamal is an asymmetric, non-deterministic public-key encryption system.
Its security is equivalent to Diffie-Hellman, since both rely on the hardness
of the discrete logarithm problem. 
Also see the `Wikipedia article <http://en.wikipedia.org/wiki/ElGamal_encryption>`_ on this scheme.

REFERENCES:

- Alfred Menezes, et al. *Handbook of Applied Cryptography; Chapter 8, Public-Key Encryption*. http://cacr.uwaterloo.ca/hac/about/chap8.pdf. Copyright 1997.
- Thomas Pornin. *When to use RSA and when ElGamal asymmetric encryption*. http://crypto.stackexchange.com/questions/1677/when-to-use-rsa-and-when-elgamal-asymmetric-encryption. Accessed 11 March 2015.

AUTHORS:

- Peter Story (2015-03) - initial release under the GPL

"""

###########################################################################
# Copyright (c) 2015
# Peter Story
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

from sage.misc.randstate import set_random_seed
from sage.rings.arith import euler_phi
from sage.rings.arith import gcd
from sage.rings.arith import is_prime
from sage.rings.arith import power_mod
from sage.rings.arith import primitive_root
from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod import Mod as mod
from sage.rings.integer_ring import ZZ
from sage.sets.primes import Primes

class ElGamalCryptosystem:
    r"""
    The ElGamal public-key encryption scheme.

    EXAMPLES:

    The following is an encryption/decryption example.::

        sage: ec = ElGamalCryptosystem(); ec
        The ElGamal public-key encryption scheme.
        sage: p = 101
        sage: public, private = ec.calculate_keys(p)
        sage: P = 42
        sage: C = ec.encrypt(P, public)
        sage: M = ec.decrypt(C, private); M
        42
        sage: M == P
        True
    """

    def __eq__(self, other):
        """
        Compare this ``ElGamalCryptosystem`` object with ``other``.

        INPUT:

        - ``other`` - an object.

        OUTPUT:

        - ``True`` if both ``self`` and ``other`` are ``ElGamalCryptosystem``
          objects. ``False`` otherwise.

        Since ``ElGamalCryptosystem`` objects do not have state, any two
        ``ElGamalCryptosystem`` objects are equivalent. 

        EXAMPLES::

            sage: ec1 = ElGamalCryptosystem()
            sage: ec2 = ElGamalCryptosystem()
            sage: ec1 == ec2
            True
        """
        if (self.__class__ == other.__class__):
            return True
        else:
            return False

    def __repr__(self):
        """
        A string representation of this ElGamal public-key encryption scheme.

        OUTPUT:

        - A string representation of the ElGamal public-key encryption scheme.

        EXAMPLES::

            sage: ElGamalCryptosystem()
            The ElGamal public-key encryption scheme.
        """
        return "The ElGamal public-key encryption scheme."

    def decrypt(self, C, K):
        r"""
        Apply the ElGamal public-key encryption scheme to decrypt the ciphertext ``C``
        using the private key ``K``.

        INPUT:

        - ``C`` - an integer ciphertext and header, of the form 
          ``(gamma, delta)`` resulting from encrypting a plaintext using the 
          ElGamal public-key encryption algorithm.

        - ``K`` - a private key, stored as a tuple ``(p, a)``.

        OUTPUT:

        - The plaintext resulting from decrypting the ciphertext ``C`` using
          the ElGamal public-key decryption algorithm.

        ALGORITHM:
        The ElGamal public-key decryption algorithm is described as follows:

        #. Let `C` be the ciphertext and header, of the form `(\gamma, \delta)`. 
           Remember that `\gamma = g^k \text { mod } (p)` and 
           `\delta = (g^a)^k * M \text { mod } (p)`, where `M` is the plaintext
           message. 
        #. Let `(p, a)` be the private key. Remember that ``a``
           is the power used to compute `g^a` in the corresponding public key
           `(p, g, g^a)`.
        #. The plaintext is recovered by 
           `M = \gamma^{\phi(p) - a} * \delta \text{ mod } (p)`. This works 
           because it is equivalent to
           `\left(g^k\right)^{\phi(p) - a} \left( g^{k} \right)^a M = \left(g^{\phi(p)}\right)^k M`
           which by Euler's Theorem is equivalent to 
           `1^k M \equiv M`

        EXAMPLES:

        An example of decryption.::

            sage: ec = ElGamalCryptosystem()
            sage: p = 101
            sage: public, private = ec.calculate_keys(p)
            sage: P = 42
            sage: C = ec.encrypt(P, public)
            sage: M = ec.decrypt(C, private); M; P == M
            42
            True


        TESTS:

        Setup for tests.::

            sage: ec = ElGamalCryptosystem()
            sage: p = 101
            sage: public, private = ec.calculate_keys(p)


        Decryption should recover the original message.::

            sage: P = 42
            sage: C = ec.encrypt(P, public)
            sage: M = ec.decrypt(C, private); P == M
            True
            sage: P = ZZ.random_element(2, p)
            sage: C = ec.encrypt(P, public)
            sage: ec.decrypt(C, private) == P
            True

        A header (gamma) outside the range ``[2, p-1]`` should be rejected.::

            sage: P = 42
            sage: C = ec.encrypt(P, public)
            sage: C = (-1, C[1])
            sage: M = ec.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The header must be an integer greater than 0 and less than p.

        Ciphertext (delta) outside the range ``[2, p-1]`` should be rejected.::

            sage: P = 42
            sage: C = ec.encrypt(P, public)
            sage: C = (C[0], -1)
            sage: M = ec.decrypt(C, private)
            Traceback (most recent call last):
            ...
            ValueError: The ciphertext must be an integer greater than 0 and less than p.

        """
        # private key
        p, a = K

        # Header and Ciphertext, respectively
        gamma, delta = C

        # Computer the Euler characteristic of p
        euler_phi = p - 1

        # sanity check
        if not ( (gamma > 0) and (gamma < p) and (gamma in ZZ) ):
            raise ValueError("The header must be an integer greater than 0 and less than p.")
        if not ( (delta > 0) and (delta < p) and (delta in ZZ) ):
            raise ValueError("The ciphertext must be an integer greater than 0 and less than p.")

        # perform the decryption
        return (power_mod(gamma, euler_phi - a, p) * delta) % p

    def encrypt(self, P, K):
        r"""
        Apply the ElGamal public-key encryption scheme to encrypt the plaintext ``P`` using
        the public key ``K``.

        INPUT:

        - ``P`` - the plaintext, a number in the range ``[2, p-1]``. 

        - ``K`` - a public key, stored as a tuple ``(p, g, g^a)`` where ``n``
          is the modulus and ``e`` is the exponent.

        OUTPUT:

        - The ciphertext and header of the form ``(gamma, delta)`` resulting
          from encrypting ``P`` using the public key ``K``. 

        ALGORITHM:

        The ElGamal public-key encryption algorithm is described as follows:

        #. Let `(p, g, g^a)` be a public key, where `g` is a generator for the
           group `\mathds{Z}_p`.
        #. Let `M` be the plaintext message, such that `0 < M < p`.
        #. First, choose an integer `k`, such that `0 < k < p` . 
        #. Next, compute `\gamma = g^k \text{ mod } (p)`, which is known as
           the header. 
        #. Finally, compute the ciphertext `\delta = (g^a)^k * m \text{ mod } (p)`. 
        #. The ciphertext is packaged with the header as `C = (\gamma, \delta)`.


        EXAMPLES:
        An example of encryption. 
        ``P = 42`` is our plaintext, and `C = (\gamma, \delta)` is the
        resulting ciphertext.::

            sage: ec = ElGamalCryptosystem()
            sage: p = 101
            sage: public, private = ec.calculate_keys(p)
            sage: P = 42
            sage: C = ec.encrypt(P, public)

        Note that encryption is not deterministic, because `k` is chosen
        randomly each time a message is encrypted.::

            sage: ec = ElGamalCryptosystem()
            sage: p = 5915587277
            sage: public, private = ec.calculate_keys(p)
            sage: P = 42
            sage: ec.encrypt(P, public) != ec.encrypt(P, public)
            True


        TESTS:

        Setup for tests.::

            sage: ec = ElGamalCryptosystem()
            sage: p = 101
            sage: public, private = ec.calculate_keys(p)

        Plaintext outside the range ``[2, p-1]`` should be rejected.::

            sage: P = 0
            sage: C = ec.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than p.

        Non-integer plaintext should be rejected.::

            sage: P = "cat"
            sage: C = ec.encrypt(P, public)
            Traceback (most recent call last):
            ...
            ValueError: The plaintext must be an integer greater than 1 and less than p.

        """
        # public key
        p, g, g_to_a = K

        # sanity check
        if (P not in ZZ) or (P < 2) or (P >= p):
            raise ValueError("The plaintext must be an integer greater than 1 and less than p.")

        # Choose a temporary number
        set_random_seed()
        k = ZZ.random_element(2, p)

        # Calculate gamma
        gamma = power_mod(g, k, p)

        # Calculate delta
        delta = (power_mod(g_to_a, k, p) * P) % p

        # Return the header and the encrypted message
        return (gamma, delta)

    def calculate_keys(self, p, enforce_primes=True):
        r"""
        Return an ElGamal public and private key pair corresponding to the
        prime ``p``. 

        INPUT:

        - ``p`` - a prime number.

        - ``enforce_primes`` - whether to require that ``p`` be prime.

        OUTPUT:

        - The ElGamal public and private keys as a tuple ``((p, g, g^a), (p, a))``.


        EXAMPLES:

        Compute an ElGamal public and private key pair corresponding to the 
        supplied prime.::

            sage: ec = ElGamalCryptosystem()
            sage: p = 19
            sage: public, private = ec.calculate_keys(p)

        For teaching purposes, you can supply a number which is not a prime.
        This is also helpful if you are using very large primes, to avoid the
        expensive computation of checking primeness.
        Frequently, generating the keys will fail because it will be impossible
        to find a generator for the group.::

            sage: ec = ElGamalCryptosystem()
            sage: n = 8
            sage: public, private = ec.calculate_keys(n, False)
            Traceback (most recent call last):
            ...
            ValueError: no primitive root

        However, some composite numbers do have generators. In this case,
        key generation will succeed.
        Unfortunately, without using a prime, decryption will (usually) fail.::

            sage: ec = ElGamalCryptosystem()
            sage: n = 250
            sage: public, private = ec.calculate_keys(n, False)
            sage: P = 42
            sage: P != ec.decrypt(ec.encrypt(P, public), private)
            True

        TESTS:

        Setup for tests.::
            sage: ec = ElGamalCryptosystem()

        Key generation should succeed when prime ``p`` is supplied.::

            sage: ec.calculate_keys(3)
            ((3, 2, 1), (3, 2))

        Input ``p`` must be an integer.::

            sage: public, private = ec.calculate_keys("my_bad_key")
            Traceback (most recent call last):
            ...
            ValueError: p must be an integer.

        Input ``p`` must be a prime.::

            sage: public, private = ec.calculate_keys(21)
            Traceback (most recent call last):
            ...
            ValueError: p must be a prime.

        Allow ``p`` to be composite if ``enforce_primes`` is set
        to ``False``.::

            sage: public, private = ec.calculate_keys(250, False)

        Of course, composite ``p`` are not guaranteed to work, since
        composites are not guaranteed to have generators.::

            sage: public, private = ec.calculate_keys(8, False)
            Traceback (most recent call last):
            ...
            ValueError: no primitive root

        """
        if p not in ZZ:
            raise ValueError("p must be an integer.")
        if enforce_primes and (not is_prime(p)):
            raise ValueError("p must be a prime.")

        # We are guaranteed a primitive root because ``p`` is prime
        g = primitive_root(p)
        # Choose the private key ``a`` such that ``1 < a < p``
        set_random_seed()
        a = ZZ.random_element(2, p)
        # Calculate ``g^a mod (p)``
        g_to_a = power_mod(g, a, p)

        public_key = (p, g, g_to_a)
        private_key = (p, a)
        return (public_key, private_key)
