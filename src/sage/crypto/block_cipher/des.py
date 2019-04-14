r"""
DES

The Data Encryption Standard.

This file implements the Data Encryption Standard and the corresponding key
schedule as described in [U.S1999]_.

Note that since this is Sage the implementation is neither meant to be
particular secure nor super fast but to be easy to use for research purposes.

EXAMPLES:

Encrypt a message::

    sage: from sage.crypto.block_cipher.des import DES
    sage: des = DES()
    sage: des.encrypt(P=0x01A1D6D039776742, K=0x7CA110454A1A6E57).hex()
    '690f5b0d9a26939b'

And decrypt it again::

    sage: des.decrypt(C=0x690F5B0D9A26939B, K=0x7CA110454A1A6E57).hex()
    '1a1d6d039776742'

Have a look at the used round keys::

    sage: from sage.crypto.block_cipher.des import DES_KS
    sage: ks = DES_KS()
    sage: [k.hex() for k in ks(0x1F08260D1AC2465E)]
    ['103041bfb90e',
     '808540f07bf',
      ...
     '221000f2dd97']

AUTHORS:

- Lukas Stennes (2019-03-29): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Lukas Stennes <lukas.stennes@rub.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF


class DES(SageObject):
    r"""
    This class implements DES described in [U.S1999]_.

    EXAMPLES:

    You can invoke DES encryption/decryption either by calling DES with an
    appropriate falg::

        sage: from sage.crypto.block_cipher.des import DES
        sage: des = DES()
        sage: P = 0x95F8A5E5DD31D900
        sage: K = 0x0
        sage: des(des(P, K, 'encrypt'), K, 'decrypt') == P
        True

    Or by calling encryption/decryption methods directly::

        sage: C = des.encrypt(P, K)
        sage: P == des.decrypt(C, K)
        True

    The number of rounds can be reduced easily::

        sage: des = DES(rounds=15)
        sage: des(des(P, K, 'encrypt'), K, 'decrypt') == P
        True

    You can use hex (i.e. integers) or a list-like bit representation for the
    inputs. If the input is an integer the output will be too. If it is
    list-like the output will be a bit vector::

        sage: P = ZZ(0).digits(2,padto=64)
        sage: K = ZZ(0).digits(2,padto=56)
        sage: list(des(des(P, K, 'encrypt'), K, 'decrypt')) == P
        True
        sage: P = ZZ(0).digits(2,padto=64)
        sage: K = 0x0
        sage: list(des(des(P, K, 'encrypt'), K, 'decrypt')) == P
        True

    .. SEEALSO::

        :class:`DES_KS`
        :mod:`sage.crypto.sboxes`

    .. automethod:: __init__
    .. automethod:: __call__
    """

    def __init__(self, rounds=None, keySchedule=None):
        r"""
        Construct an instance of DES.

        INPUT:

        - ``rounds``  -- integer (default: ``None``); the number of rounds. If
          ``None`` the number of rounds of the key schedule is used.

        - ``keySchedule`` -- (default: ``None``); the key schedule that will be
          used for encryption and decryption. If ``None`` the default DES key
          schedule is used.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: DES() # indirect doctest
            DES block cipher with 16 rounds and the following key schedule:
            Original DES key schedule with 16 rounds

        Reducing the number of rounds is simple. But increasing it is not
        possible::

            sage: DES(rounds=11) # indirect doctest
            DES block cipher with 11 rounds and the following key schedule:
            Original DES key schedule with 16 rounds
            sage: DES(rounds=42) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: number of rounds must be less or equal to the number
            of rounds of the key schedule

        You can use arbitrary key schedules. Since it is the only one
        implemented here the original key schedule is used for demonstration::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES(keySchedule=DES_KS(11)) # indirect doctest
            DES block cipher with 11 rounds and the following key schedule:
            Original DES key schedule with 11 rounds
        """
        if keySchedule is None:
            self._keySchedule = DES_KS()
        else:
            self._keySchedule = keySchedule
        if rounds is None:
            self._rounds = self._keySchedule._rounds
        elif rounds <= self._keySchedule._rounds:
            self._rounds = rounds
        else:
            raise ValueError('number of rounds must be less or equal to the '
                             'number of rounds of the key schedule')
        self._blocksize = 64

    def __call__(self, B, K, algorithm='encrypt'):
        r"""
        Apply DES encryption or decryption on ``B`` using the key ``K``.
        The flag ``algorithm`` controls what action is to be performed on
        ``B``.

        INPUT:

        - ``B`` -- integer or bit list-like; the plaintext or ciphertext

        - ``K`` -- integer or bit list-like; the key

        - ``algorithm`` -- string (default: ``'encrypt'``); a flag to signify
          whether encryption or decryption is to be applied to ``B``. The
          encryption flag is ``'encrypt'`` and the decryption flag is
          ``'decrypt'``

        OUTPUT:

        - The plaintext or ciphertext corresponding to ``B``, obtained using
          the key ``K``. If ``B`` is an integer the output will be too. If
          ``B`` is list-like the output will be a bit vector.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: P = 0x480D39006EE762F2
            sage: K = 0x025816164629B007
            sage: des(P, K, 'encrypt').hex()
            'a1f9915541020b56'
        """
        if algorithm == 'encrypt':
            return self.encrypt(B, K)
        elif algorithm == 'decrypt':
            return self.decrypt(B, K)
        else:
            raise ValueError('Algorithm must be \'encrypt\' or \'decrypt\' and'
                             ' not \'%s\'' % algorithm)

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        DES objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des1 = DES()
            sage: des2 = DES()
            sage: des3 = DES(rounds=11)
            sage: des1 == des2
            True
            sage: des1 == des3
            False
            sage: des2 == 42
            False
        """
        if not isinstance(other, DES):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: DES() # indirect doctest
            DES block cipher with 16 rounds and the following key schedule:
            Original DES key schedule with 16 rounds
        """
        return('DES block cipher with %s rounds and the following key '
               'schedule:\n%s' % (self._rounds, self._keySchedule.__repr__()))

    def encrypt(self, P, K):
        r"""
        Return the ciphertext corresponding to the plaintext ``P``,
        using DES encryption with key ``K``.

        INPUT:

        - ``P`` -- integer or bit list-like; the plaintext that will be
          encrypted.

        - ``K`` -- integer or bit list-like; the key

        OUTPUT:

        - The ciphertext corresponding to ``P``, obtained using the key
          ``K``. If ``P`` is an integer the output will be too. If ``P`` is
          list-like the output will be a bit vector.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: K = 0x133457799BBCDFF1
            sage: P = 0x0123456789ABCDEF
            sage: C = 0x85E813540F0AB405
            sage: des.encrypt(P, K) == C
            True
            sage: K = 0x0101010101010101
            sage: P = 0x95F8A5E5DD31D900
            sage: C = 0x8000000000000000
            sage: des.encrypt(P, K) == C
            True
        """
        state, inputType = _convert_to_vector(P, 64)
        K, _ = _convert_to_vector(K, self._keySchedule._keysize)
        roundKeys = self._keySchedule(K)
        state = self._ip(state)
        L, R = state[0:32], state[32:64]
        for k in roundKeys[:self._rounds]:
            L, R = R, L + self._f(R, k)
        state = vector(GF(2), 64, list(R)+list(L))
        state = self._inv_ip(state)
        return state if inputType == 'vector' else ZZ(list(state)[::-1], 2)

    def decrypt(self, C, K):
        r"""
        Return the plaintext corresponding to the ciphertext ``C``,
        using DES decryption with key ``K``.

        INPUT:

        - ``C`` -- integer or bit list-like; the ciphertext that will be
          decrypted

        - ``K`` -- integer or bit list-like; the key

        OUTPUT:

        - The plaintext corresponding to ``C``, obtained using the key
          ``K``. If ``C`` is an integer the output will be too. If ``C`` is
          list-like the output will be a bit vector.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: K = 0x7CA110454A1A6E57
            sage: C = 0x690F5B0D9A26939B
            sage: des.decrypt(C, K).hex()
            '1a1d6d039776742'
        """
        state, inputType = _convert_to_vector(C, 64)
        K, _ = _convert_to_vector(K, self._keySchedule._keysize)
        roundKeys = self._keySchedule(K)
        state = self._ip(state)
        L, R = state[0:32], state[32:64]
        for k in roundKeys[:self._rounds][::-1]:
            L, R = R, L + self._f(R, k)
        state = vector(GF(2), 64, list(R)+list(L))
        state = self._inv_ip(state)
        return state if inputType == 'vector' else ZZ(list(state)[::-1], 2)

    def _ip(self, B):
        r"""
        Return the initial permutation of ``B``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 64, [0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1])
            sage: des._ip(B)
            (1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0)
        """
        IP = [58, 50, 42, 34, 26, 18, 10, 2,
              60, 52, 44, 36, 28, 20, 12, 4,
              62, 54, 46, 38, 30, 22, 14, 6,
              64, 56, 48, 40, 32, 24, 16, 8,
              57, 49, 41, 33, 25, 17,  9, 1,
              59, 51, 43, 35, 27, 19, 11, 3,
              61, 53, 45, 37, 29, 21, 13, 5,
              63, 55, 47, 39, 31, 23, 15, 7]
        return vector(GF(2), 64, [B[i-1] for i in IP])

    def _f(self, R, K):
        r"""
        Apply the cipher function to ``R`` and ``K``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: R = vector(GF(2), 32, [1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0])
            sage: K = vector(GF(2), 48, [0,0,0,1,1,0,1,1,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1,0])
            sage: des._f(R, K)
            (0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1)
        """
        return self._permutaion(self._sboxes(self._expand(R)+K))

    def _expand(self, R):
        r"""
        Apply the expansion function to ``R``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: R = vector(GF(2), 32, [1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0])
            sage: des._expand(R)
            (0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
        """
        E = [32,  1,  2,  3,  4,  5,
              4,  5,  6,  7,  8,  9,
              8,  9, 10, 11, 12, 13,
             12, 13, 14, 15, 16, 17,
             16, 17, 18, 19, 20, 21,
             20, 21, 22, 23, 24, 25,
             24, 25, 26, 27, 28, 29,
             28, 29, 30, 31, 32,  1]
        return vector(GF(2), 48, [R[i-1] for i in E])

    def _sboxes(self, B):
        r"""
        Apply the Sboxes to ``R``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 48, [0,1,1,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,1])
            sage: des._sboxes(B)
            (0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1)

        .. SEEALSO::

            :mod:`sage.crypto.sboxes`
        """
        from sage.crypto.sboxes import DES_S1_1, DES_S1_2, DES_S1_3, DES_S1_4
        from sage.crypto.sboxes import DES_S2_1, DES_S2_2, DES_S2_3, DES_S2_4
        from sage.crypto.sboxes import DES_S3_1, DES_S3_2, DES_S3_3, DES_S3_4
        from sage.crypto.sboxes import DES_S4_1, DES_S4_2, DES_S4_3, DES_S4_4
        from sage.crypto.sboxes import DES_S5_1, DES_S5_2, DES_S5_3, DES_S5_4
        from sage.crypto.sboxes import DES_S6_1, DES_S6_2, DES_S6_3, DES_S6_4
        from sage.crypto.sboxes import DES_S7_1, DES_S7_2, DES_S7_3, DES_S7_4
        from sage.crypto.sboxes import DES_S8_1, DES_S8_2, DES_S8_3, DES_S8_4
        from itertools import chain
        sbox = [[DES_S1_1, DES_S1_2, DES_S1_3, DES_S1_4],
                [DES_S2_1, DES_S2_2, DES_S2_3, DES_S2_4],
                [DES_S3_1, DES_S3_2, DES_S3_3, DES_S3_4],
                [DES_S4_1, DES_S4_2, DES_S4_3, DES_S4_4],
                [DES_S5_1, DES_S5_2, DES_S5_3, DES_S5_4],
                [DES_S6_1, DES_S6_2, DES_S6_3, DES_S6_4],
                [DES_S7_1, DES_S7_2, DES_S7_3, DES_S7_4],
                [DES_S8_1, DES_S8_2, DES_S8_3, DES_S8_4]]
        B = [B[i:i+6] for i in range(0, 48, 6)]
        B = list(chain.from_iterable([sbox[i][2*int(b[0])+int(b[5])](b[1:5])
                                      for i, b in enumerate(B)]))
        return vector(GF(2), 32, B)

    def _permutaion(self, B):
        r"""
        Apply the permutation function to ``R``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 32, [0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1])
            sage: des._permutaion(B)
            (0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1)
        """
        P = [16,  7, 20, 21,
             29, 12, 28, 17,
              1, 15, 23, 26,
              5, 18, 31, 10,
              2,  8, 24, 14,
             32, 27,  3,  9,
             19, 13, 30,  6,
             22, 11,  4, 25,]
        return vector(GF(2), 32, [B[i-1] for i in P])

    def _inv_ip(self, B):
        r"""
        Apply the inverse permutation function to ``R``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: des = DES()
            sage: B = vector(GF(2), 64, [0,0,0,0,1,0,1,0,0,1,0,0,1,1,0,0,1,1,0,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,0,0])
            sage: des._inv_ip(B)
            (1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1)
        """
        invIP = [40, 8, 48, 16, 56, 24, 64, 32,
                 39, 7, 47, 15, 55, 23, 63, 31,
                 38, 6, 46, 14, 54, 22, 62, 30,
                 37, 5, 45, 13, 53, 21, 61, 29,
                 36, 4, 44, 12, 52, 20, 60, 28,
                 35, 3, 43, 11, 51, 19, 59, 27,
                 34, 2, 42, 10, 50, 18, 58, 26,
                 33, 1, 41,  9, 49, 17, 57, 25]
        return vector(GF(2), 64, [B[i-1] for i in invIP])


class DES_KS(SageObject):
    r"""
    This class implements the DES key schedules described in [U.S1999]_.

    EXAMPLES:

    Initialise the key schedule with a `master\_key` to use it as an iterable::

        sage: from sage.crypto.block_cipher.des import DES_KS
        sage: ks = DES_KS(master_key=0)
        sage: ks[0] == 0x0
        True
        sage: ks[15] == 0x0
        True

    Or omit the `master\_key` and pass a key when calling the key schedule::

        sage: ks = DES_KS()
        sage: K = ks(0x584023641ABA6176)
        sage: K[0] == 0xD0A2E52FA124
        True
        sage: K[15] == 0x42B42AF81183
        True

    .. SEEALSO::

        :class:`DES`

    .. automethod:: __init__
    .. automethod:: __call__
    """

    def __init__(self, rounds=16, master_key=None):
        r"""
        Construct an instance of DES_KS.

        INPUT:

        - ``rounds`` -- integer (default: ``16``); the number of rounds
          ``self`` can create keys for

        - ``master_key`` -- integer of bit list-like; the key that will be used

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS()
            Original DES key schedule with 16 rounds

        .. NOTE::

            If you want to use a DES_KS object as an iterable you have to
            pass a ``master_key`` value on initialisation. Otherwise you can
            omit ``master_key`` and pass a key when you call the object.
        """
        self._rounds = rounds
        self._master_key = master_key
        self._keysize = 64

    def __call__(self, K):
        r"""
        Return all round keys in a list.

        INPUT:

        - ``K`` -- integer or bit list-like; the key

        OUTPUT:

        - A list containing the round keys

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1])
            sage: ks = DES_KS(16, K)
            sage: [k for k in ks]
            [(0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0),
             (0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1),
             ...
             (1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1)]
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,1,1,0,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,0])
            sage: ks = DES_KS()
            sage: ks(K)
            [(0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0),
             (0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1),
             ...
             (1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1)]
            sage: K = 0x133457799bbcdff1
            sage: ks = DES_KS(master_key=K)
            sage: [k.hex() for k in ks]
            ['1b02effc7072',
             '79aed9dbc9e5',
             ...
             'cb3d8b0e17f5']

        .. NOTE::

            If you want to use a DES_KS object as an iterable you have to
            pass a ``master_key`` value on initialisation. Otherwise you can
            omit ``master_key`` and pass a key when you call the object.
        """
        K, inputType = _convert_to_vector(K, self._keysize)
        roundKeys = []
        # ensure that K is a 64 bit vector
        if not any(K[56:]):
            # delete msbs and insert 'parity' bits
            K = list(K)[:56]
            for i in range(7, 64, 8):
                K.insert(i, 0)
            K = vector(GF(2), 64, K)
        C, D = self._pc1(K)
        for i in range(16):
            C, D = self._left_shift(C, i), self._left_shift(D, i)
            roundKeys.append(self._pc2(list(C)+list(D)))
        return roundKeys if inputType == 'vector' else [ZZ(list(k)[::-1], 2)
                                                        for k in roundKeys]

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        DES_KS objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS() == DES_KS() # indirect doctest
            True
            sage: DES_KS() == DES_KS(11) # indirect doctest
            False
        """
        if not isinstance(other, DES_KS):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this DES_KS.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: DES_KS() # indirect doctest
            Original DES key schedule with 16 rounds
        """
        return ('Original DES key schedule with %s rounds' % (self._rounds))

    def __getitem__(self, r):
        r"""
        Computes the sub key for round ``r`` derived from initial master key.

        The key schedule object has to have been initialised with the
        `master_key` argument.

        INPUT:

        - ``r`` integer; the round for which the sub key is computed

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS(master_key=0x0)
            sage: ks[0] ==  0x0 # indirect doctest
            True
            sage: ks[15] ==  0x0 # indirect doctest
            True
        """
        if self._master_key is None:
            raise ValueError('Key not set during initialisation')
        return self(self._master_key)[r]

    def __iter__(self):
        r"""
        Iterate over the DES round keys, derived from `master_key`.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: K = [k for k in DES_KS(master_key=0x0)]
            sage: K[0] == 0x0 # indirect doctest
            True
            sage: K[15] == 0x0 # indirect doctest
            True
       """
        if self._master_key is None:
            raise ValueError('Key not set during initialisation')
        return iter(self(self._master_key))

    def _pc1(self, K):
        r"""
        Return Permuted Choice 1 of ``K``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: K = vector(GF(2),[0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,0,0,0,1])
            sage: C, D = ks._pc1(K)
            sage: C
            (1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1)
            sage: D
            (0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1)
        """
        PC1_C = [57, 49, 41, 33, 25, 17,  9,
                  1, 58, 50, 42, 34, 26, 18,
                 10,  2, 59, 51, 43, 35, 27,
                 19, 11,  3, 60, 52, 44, 46]
        PC1_D = [63, 55, 47, 39, 31, 23, 15,
                  7, 62, 54, 46, 38, 30, 22,
                 14,  6, 61, 53, 45, 37, 29,
                 21, 13,  5, 28, 20, 12,  4]
        C = vector(GF(2), 28, [K[i-1] for i in PC1_C])
        D = vector(GF(2), 28, [K[i-1] for i in PC1_D])
        return C, D

    def _pc2(self, K):
        r"""
        Return Permuted Choice 2 of ``K``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: K = vector(GF(2),[1,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1,0])
            sage: ks._pc2(K)
            (0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0)
        """
        PC2 = [14, 17, 11, 24,  1,  5,
                3, 28, 15,  6, 21, 10,
               23, 19, 12,  4, 26,  8,
               16,  7, 27, 20, 13,  2,
               41, 52, 31, 37, 47, 55,
               30, 40, 51, 45, 33, 48,
               44, 49, 39, 56, 34, 53,
               46, 42, 50, 36, 29, 32]
        return vector(GF(2), 48, [K[i-1] for i in PC2])

    def _left_shift(self, half, i):
        r"""
        Shift ``half`` one or two positions to the left depending on the
        iteration number ``i``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES_KS
            sage: ks = DES_KS()
            sage: bits = vector(GF(2), 6, [1,0,1,0,1,0])
            sage: ks._left_shift(bits, 0)
            (0, 1, 0, 1, 0, 1)
            sage: bits
            (1, 0, 1, 0, 1, 0)
            sage: ks._left_shift(bits, 2)
            (1, 0, 1, 0, 1, 0)
        """
        amount = [1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]
        return vector(GF(2),
                      list(half[amount[i]:]) + list(half[0:amount[i]]))


def _convert_to_vector(I, L):
    r"""
    Convert ``I`` to a bit vector of length ``L``.

    INPUT:

    - ``I`` -- integer or bit list-like

    - ``L`` -- integer; the desired bit length of the ouput

    OUTPUT: a tuple of

    - the ``L``-bit vector representation of ``I``

    - a flag indicating the input type. Either ``'integer'`` or ``'vector'``.

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import _convert_to_vector
        sage: _convert_to_vector(0x1F, 8)
        ((0, 0, 0, 1, 1, 1, 1, 1), 'integer')
        sage: v = vector(GF(2), 4, [1,0,1,0])
        sage: _convert_to_vector(v, 4)
        ((1, 0, 1, 0), 'vector')
    """
    try:
        state = vector(GF(2), L, ZZ(I).digits(2, padto=L)[::-1])
        return state, 'integer'
    except TypeError:
        # ignore the error and try list-like types
        pass
    state = vector(GF(2), L, ZZ(list(I), 2).digits(2, padto=L))
    return state, 'vector'
