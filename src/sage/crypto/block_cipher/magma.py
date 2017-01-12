r"""
Magma Encryption and Decryption Algorithms

Algorithm of the Magma Cipher (GOST 34.12-2015, n=64), ECB mode. Note that
this algorithm is for educational purposes only. It is a
 version of the Magma designed to help beginners understand the
basic structure of algorithm and get a plain text/cipher text values for different number of encryption round.

AUTHORS:

- Ekaterina Maro (2017-01): initial version
"""

###########################################################################
# Copyright (c) 2017
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################

from sage.structure.sage_object import SageObject
from sage.monoids.string_monoid import BinaryStrings
from sage.crypto.mq import SBox
S0=SBox(12, 4, 6, 2, 10, 5, 11, 9, 14, 8, 13, 7, 0, 3, 15, 1)
S1=SBox(6, 8, 2, 3, 9, 10, 5, 12, 1, 14, 4, 7, 11, 13, 0, 15)
S2=SBox(11, 3, 5, 8, 2, 15, 10, 13, 14, 1, 7, 4, 12, 9, 6, 0)
S3=SBox(12, 8, 2, 1, 13, 4, 15, 6, 7, 0, 10, 5, 3, 14, 9, 11)
S4=SBox(7, 15, 5, 10, 8, 1, 6, 13, 0, 9, 3, 14, 11, 4, 2, 12)
S5=SBox(5, 13, 15, 6, 9, 2, 12, 10, 11, 7, 8, 1, 4, 3, 14, 0)
S6=SBox(8, 14, 2, 5, 6, 9, 1, 12, 15, 4, 11, 0, 13, 10, 3, 7)
S7=SBox(1, 7, 14, 13, 0, 5, 8, 3, 4, 15, 10, 6, 9, 12, 11, 2)

class MagmaCipher(SageObject):
    r"""
        This class implements the Magma Encryption Encryption Algorithm (ECB mode)
        described in [http://www.tc26.ru/en/standard/gost/GOST_R_34_12_2015_ENG.pdf].
        For educational purposales only and is not use for practical purposes.
        """
    def __init__(self):
        r"""
        A  variant of the Magma Encryption Algorithm GOST 34.12-2015, n=64, ECB mode.  
        INPUT:

        - ``B`` -- a integer value, where the number of bits is positive and may be different.

        - ``K`` -- a secret key; this may be a not multyple 256 binary string.

        - ``Z`` -- a number of round; this must be a not bigger 32 value.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.magma import MagmaCipher
            sage: magma = MagmaCipher(); magma
            Magma block cipher with 256-bit keys
            sage: K=0xffeeddccbbaa99887766554433221100f0f1f2f3f4f5f6f7f8f9fafbfcfdfeff
            sage: B=0xfedcba9876543210
            sage: magma.encrypt(B, K, 32)
            '4ee901e5c2d8ca3d'
            sage: C=0x4ee901e5c2d8ca3d
            sage: magma.decrypt(C, K, 32)
            'fedcba9876543210'
            """
        Blength = 64


    def __call__(self, B, K, Z, algorithm="encrypt"):
        """"Apply Magma encryption or decryption on the binary string ``B``
        using the key ``K`` with "Z" rounds.  The flag ``algorithm`` controls what action is
        to be performed on ``B``.

        INPUT:

        - ``B`` -- a integer value, where the number of bits is positive and may be different.

        - ``K`` -- a secret key; this may be a not multyple 256 binary string.

        - ``algorithm`` -- (default: ``"encrypt"``) a string; a flag to signify
          whether encryption or decryption is to be applied to the binary
          string ``B``. The encryption flag is ``"encrypt"`` and the decryption
          flag is ``"decrypt"``.

        OUTPUT:

        - The ciphertext (respectively plaintext) corresponding to the
          binary string ``S``.
        """

        # encrypt each 64-bit block
        if algorithm == "encrypt":
            block = B
            key = K
            number_round=Z
            enc=self.encrypt(block, key, number_round)
            return hex(int(enc,16))
        elif algorithm == "decrypt":
            block = B
            key = K
            number_round=Z
            denc=self.decrypt(block, key, number_round)
            return hex(int(denc,16))
        else:
            raise ValueError("Algorithm must be either 'encrypt' or 'decrypt'")

    def __call__(self, B, K, Z, algorithm="decrypt"):
        r"""
        Return an plaintext corresponding to the ciphertext ``B``,
        using Magma decryptionalgorithm with key ``K`` in "Z"-round of decryption ."""
        # decrypt each 64-bit block

        if algorithm == "encrypt":
            block = B
            key = K
            number_round=Z
            enc=self.encrypt(block, key, number_round)
            return hex(int(enc,16))
        elif algorithm == "decrypt":
            block = B
            key = K
            number_round=Z
            denc=self.decrypt(block, key, number_round)
            return hex(int(denc,16))
        else:
            raise ValueError("Algorithm must be either 'encrypt' or 'decrypt'")

    def encrypt(self, P, K, Z):
        N = P.nbits() // 64 # the number of 64-bit blocks
        if P.nbits() % 64!=0:
            N=N+1
        C=0
        S=""
        for j in range (N):
            plain_text=P>>((N-1-j)*64)
            plain_text_l=(plain_text&0xffffffff00000000)>>32
            plain_text_r=plain_text&0xffffffff
            for i in range (Z):
                plain_text_buf=plain_text_r
                if i>=0 and i<24:
                    plain_text_r=(plain_text_r+(K>>(256-(i%8+1)*32)))&0xffffffff #addition with round key
                    key_round=(K>>(256-(i%8+1)*32))&0xffffffff

                if i>=24:
                    plain_text_r=(plain_text_r+(K>>(i%8*32)))&0xffffffff #addition with inverse round key

    #substitution layer
                plain_text_S=0x0
                plain_text_Shift=0x0
                plain_text_S=S7(plain_text_r>>28&0xf)
                plain_text_S=(plain_text_S<<4)^S6(plain_text_r>>24&0xf)
                plain_text_S=(plain_text_S<<4)^S5(plain_text_r>>20&0xf)
                plain_text_S=(plain_text_S<<4)^S4(plain_text_r>>16&0xf)
                plain_text_S=(plain_text_S<<4)^S3(plain_text_r>>12&0xf)
                plain_text_S=(plain_text_S<<4)^S2(plain_text_r>>8&0xf)
                plain_text_S=(plain_text_S<<4)^S1(plain_text_r>>4&0xf)
                plain_text_S=(plain_text_S<<4)^S0(plain_text_r&0xf)

                plain_text_Shift=plain_text_S&0xffffffff
                plain_text_Shift=(plain_text_S<<11)&0xffffffff
                plain_text_Shift=plain_text_Shift^(plain_text_S>>21)
                plain_text_r=plain_text_Shift^plain_text_l

                if i!=Z-1:
                    plain_text_l=plain_text_buf

                if i==Z-1:
                    cipher_text=(plain_text_r<<32)^plain_text_buf
                    C=cipher_text
                    #print  "cipher text is ", hex(cipher_text)
                    S = "".join([S, str(hex(C))])
        return S

    def decrypt(self, P, K, Z):
        N = P.nbits() // 64 # the number of 64-bit blocks
        if P.nbits() % 64!=0:
            N=N+1
        #print N
        C=0
        S=""
        for j in range (N):
            cipher_text=P>>((N-1-j)*64)
            cipher_text_l=(cipher_text&0xffffffff00000000)>>32
            cipher_text_r=cipher_text&0xffffffff
            for i in range (Z):
                cipher_text_buf=cipher_text_r
                if i>=0 and i<8:
                    cipher_text_r=(cipher_text_r+(K>>(256-(i%8+1)*32)))&0xffffffff #addition with round key
                    key_round=(K>>(256-(i%8+1)*32))&0xffffffff

                if i>=8:
                    cipher_text_r=(cipher_text_r+(K>>(i%8*32)))&0xffffffff #addition with round key


        #substitution layer
                cipher_text_S=0x0
                cipher_text_Shift=0x0
                cipher_text_S=S7(cipher_text_r>>28&0xf)
                cipher_text_S=(cipher_text_S<<4)^S6(cipher_text_r>>24&0xf)
                cipher_text_S=(cipher_text_S<<4)^S5(cipher_text_r>>20&0xf)
                cipher_text_S=(cipher_text_S<<4)^S4(cipher_text_r>>16&0xf)
                cipher_text_S=(cipher_text_S<<4)^S3(cipher_text_r>>12&0xf)
                cipher_text_S=(cipher_text_S<<4)^S2(cipher_text_r>>8&0xf)
                cipher_text_S=(cipher_text_S<<4)^S1(cipher_text_r>>4&0xf)
                cipher_text_S=(cipher_text_S<<4)^S0(cipher_text_r&0xf)

                cipher_text_Shift=cipher_text_S&0xffffffff
                cipher_text_Shift=(cipher_text_S<<11)&0xffffffff
                cipher_text_Shift=cipher_text_Shift^(cipher_text_S>>21)
                cipher_text_r=cipher_text_Shift^cipher_text_l

                if i!=Z-1:
                    cipher_text_l=cipher_text_buf

                if i==Z-1:
                    plain_text=(cipher_text_r<<32)^cipher_text_buf
                    C=plain_text
                    #print  "plain text is ", hex(plain_text)
                    S = "".join([S, str(hex(C))])
        return S

    def __repr__(self):
        r"""
        A string representation of this Magma Encryption.
        """
        return "Magma block cipher with 256-bit keys"
