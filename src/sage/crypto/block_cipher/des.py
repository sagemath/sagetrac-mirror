r"""
DES

An implementation of the Data Encryption Standard (DES).  This only  allows
the encryption of a single 64-bit plaintext with a 64-bit key.

It thus provides a foundation for the user to build up triple-DES,
DES in standard NIST block modes (ECB, CBC etc), and for padding.  None
of these are provided in this implementation.

AUTHORS:

- Alasdair McAndrew (2010-04): initial version
"""

###########################################################################
# Copyright (c) 2009 Alasdair McAndrew <amca01@gmail.com>
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

from sage.structure.sage_object import SageObject
from string import join
from sage.rings.integer import Integer

# First some helper functions, not part pf the DES class, but useful


def hex2bin(str):
    """
    turns a hexadecimal string into binary

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import hex2bin, bin2hex
        sage: hexkey='0123456789abcdef'
        sage: a = hex2bin(hexkey); a
        '0000000100100011010001010110011110001001101010111100110111101111'
        sage: bin2hex(a) == hexkey
        True
    """
    tmp = str
    if tmp[:2] != '0x':
        tmp = '0x' + tmp
    if len(tmp) != 18:
        raise ValueError("Invalid hexadecimal string size - must be 16")
    tmp = bin(Integer(tmp))
    if tmp[:2] == '0b':
        tmp = tmp[2:]
    return tmp.zfill(64)


def bin2hex(str):
    """
    turns a binary string into hexadecimal

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import hex2bin, bin2hex
        sage: binkey='0000000100100011010001010110011110001001101010111100110111101111'
        sage: a = bin2hex(binkey); a
        '0123456789abcdef'
        sage: hex2bin(a) == binkey
        True
    """
    tmp = str
    if tmp[:2] != '0b':
        tmp = '0b' + tmp
    if len(tmp) != 66:
        raise ValueError("Invalid binary string size - must be 64")
    tmp = hex(Integer(tmp))
    if tmp[:2] == '0x':
        tmp = tmp[2:]
    return tmp.zfill(16)


def rivest_test(hexstr):
    """
    Applies Rivest's test for DES implementations

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import rivest_test
        sage: hexkey='9474b8e8c73bca7d'
        sage: rivest_test(hexkey)
        '1b1a2ddb4c642438'
    """
    X = {0: hex2bin(hexstr)}
    for i in range(16):
        di = DES(X[i])
        if i % 2 == 0:
            X[i + 1] = di.encrypt(X[i])
        else:
            X[i + 1] = di.decrypt(X[i])
    return bin2hex(X[16])


class DES(SageObject):
    r"""
    This class implements the Data Encryption Standard (DES) asw
    described in NIST standard FIPS PUB 46-3 and available at
    <http://csrc.nist.gov/publications/fips/fips46-3/fips46-3.pdf>

    Input of both plaintext and key are 64 bit binary strings,
    as is the output.  For ease of use, there are two "helper"
    functions hex2bin and bin2hex which translate hexadecimal
    to and from binary.

    EXAMPLES::

        sage: from sage.crypto.block_cipher.des import DES, hex2bin, bin2hex

        sage: hexkey='0123456789abcdef'
        sage: key=hex2bin(hexkey)
        sage: D=DES(key)
        sage: hexpl='aaaabbbbccccdddd'
        sage: pl=hex2bin(hexpl)
        sage: ct=D.encrypt(pl)
        sage: bin2hex(ct)
        '6f32c06fb6ac37bf'
        sage: D.decrypt(ct)==pl
        True

    There is also a test developed by Ronald Rivest in
    "TESTING IMPLEMENTATIONS OF DES" and available at
    <http://people.csail.mit.edu/rivest/Destest.txt>::

        sage: from sage.crypto.block_cipher.des import rivest_test
        sage: rivest_test('9474B8E8C73BCA7D')
        '1b1a2ddb4c642438'

    which is the correct value to obtain so that according to
    the article: "your implementation does not have any of the
    36,568 possible single-fault errors described herein."
    """

    # Permutations and S-boxes

    # PC1 obtains 56 bits from the original 64 bits of the key

    __PC1 = [56, 48, 40, 32, 24, 16, 8,
             0, 57, 49, 41, 33, 25, 17,
             9, 1, 58, 50, 42, 34, 26,
             18, 10, 2, 59, 51, 43, 35,
             62, 54, 46, 38, 30, 22, 14,
             6, 61, 53, 45, 37, 29, 21,
             13, 5, 60, 52, 44, 36, 28,
             20, 12, 4, 27, 19, 11, 3]

     # PC2 is the compression permutation which reduces the key to 48 bits

    __PC2 = [13, 16, 10, 23, 0, 4,
             2, 27, 14, 5, 20, 9,
             22, 18, 11, 3, 25, 7,
             15, 6, 26, 19, 12, 1,
             40, 51, 30, 36, 46, 54,
             29, 39, 50, 44, 32, 47,
             43, 48, 38, 55, 33, 52,
             45, 41, 49, 35, 28, 31]

    # initial permutation IP

    __IP = [57, 49, 41, 33, 25, 17, 9, 1,
            59, 51, 43, 35, 27, 19, 11, 3,
            61, 53, 45, 37, 29, 21, 13, 5,
            63, 55, 47, 39, 31, 23, 15, 7,
            56, 48, 40, 32, 24, 16, 8, 0,
            58, 50, 42, 34, 26, 18, 10, 2,
            60, 52, 44, 36, 28, 20, 12, 4,
            62, 54, 46, 38, 30, 22, 14, 6]

     # Expansion table for turning 32 bit blocks into 48 bits

    __E = [31, 0, 1, 2, 3, 4,
           3, 4, 5, 6, 7, 8,
           7, 8, 9, 10, 11, 12,
           11, 12, 13, 14, 15, 16,
           15, 16, 17, 18, 19, 20,
           19, 20, 21, 22, 23, 24,
           23, 24, 25, 26, 27, 28,
           27, 28, 29, 30, 31, 0]

    # S1
    __sbox = [[14, 4, 13, 1, 2, 15, 11, 8, 3, 10, 6, 12, 5, 9, 0, 7,
               0, 15, 7, 4, 14, 2, 13, 1, 10, 6, 12, 11, 9, 5, 3, 8,
               4, 1, 14, 8, 13, 6, 2, 11, 15, 12, 9, 7, 3, 10, 5, 0,
               15, 12, 8, 2, 4, 9, 1, 7, 5, 11, 3, 14, 10, 0, 6, 13],
              # S2
              [15, 1, 8, 14, 6, 11, 3, 4, 9, 7, 2, 13, 12, 0, 5, 10,
               3, 13, 4, 7, 15, 2, 8, 14, 12, 0, 1, 10, 6, 9, 11, 5,
               0, 14, 7, 11, 10, 4, 13, 1, 5, 8, 12, 6, 9, 3, 2, 15,
               13, 8, 10, 1, 3, 15, 4, 2, 11, 6, 7, 12, 0, 5, 14, 9],
              # S3
              [10, 0, 9, 14, 6, 3, 15, 5, 1, 13, 12, 7, 11, 4, 2, 8,
               13, 7, 0, 9, 3, 4, 6, 10, 2, 8, 5, 14, 12, 11, 15, 1,
               13, 6, 4, 9, 8, 15, 3, 0, 11, 1, 2, 12, 5, 10, 14, 7,
               1, 10, 13, 0, 6, 9, 8, 7, 4, 15, 14, 3, 11, 5, 2, 12],
              # S4
              [7, 13, 14, 3, 0, 6, 9, 10, 1, 2, 8, 5, 11, 12, 4, 15,
               13, 8, 11, 5, 6, 15, 0, 3, 4, 7, 2, 12, 1, 10, 14, 9,
               10, 6, 9, 0, 12, 11, 7, 13, 15, 1, 3, 14, 5, 2, 8, 4,
               3, 15, 0, 6, 10, 1, 13, 8, 9, 4, 5, 11, 12, 7, 2, 14],
              # S5
              [2, 12, 4, 1, 7, 10, 11, 6, 8, 5, 3, 15, 13, 0, 14, 9,
               14, 11, 2, 12, 4, 7, 13, 1, 5, 0, 15, 10, 3, 9, 8, 6,
               4, 2, 1, 11, 10, 13, 7, 8, 15, 9, 12, 5, 6, 3, 0, 14,
               11, 8, 12, 7, 1, 14, 2, 13, 6, 15, 0, 9, 10, 4, 5, 3],
              # S6
              [12, 1, 10, 15, 9, 2, 6, 8, 0, 13, 3, 4, 14, 7, 5, 11,
               10, 15, 4, 2, 7, 12, 9, 5, 6, 1, 13, 14, 0, 11, 3, 8,
               9, 14, 15, 5, 2, 8, 12, 3, 7, 0, 4, 10, 1, 13, 11, 6,
               4, 3, 2, 12, 9, 5, 15, 10, 11, 14, 1, 7, 6, 0, 8, 13],
              # S7
              [4, 11, 2, 14, 15, 0, 8, 13, 3, 12, 9, 7, 5, 10, 6, 1,
               13, 0, 11, 7, 4, 9, 1, 10, 14, 3, 5, 12, 2, 15, 8, 6,
               1, 4, 11, 13, 12, 3, 7, 14, 10, 15, 6, 8, 0, 5, 9, 2,
               6, 11, 13, 8, 1, 4, 10, 7, 9, 5, 0, 15, 14, 2, 3, 12],
              # S8
              [13, 2, 8, 4, 6, 15, 11, 1, 10, 9, 3, 14, 5, 0, 12, 7,
               1, 15, 13, 8, 10, 3, 7, 4, 12, 5, 6, 11, 0, 14, 9, 2,
               7, 11, 4, 1, 9, 12, 14, 2, 0, 6, 10, 13, 15, 3, 5, 8,
               2, 1, 14, 7, 4, 10, 8, 13, 15, 12, 9, 0, 3, 5, 6, 11]]

    # 32-bit permutation function P used on the output of the S-boxes
    __P = [15, 6, 19, 20, 28, 11, 27, 16,
           0, 14, 22, 25, 4, 17, 30, 9,
           1, 7, 23, 13, 31, 26, 2, 8,
           18, 12, 29, 5, 21, 10, 3, 24]

    # final permutation IP**-1
    __FP = [39, 7, 47, 15, 55, 23, 63, 31,
            38, 6, 46, 14, 54, 22, 62, 30,
            37, 5, 45, 13, 53, 21, 61, 29,
            36, 4, 44, 12, 52, 20, 60, 28,
            35, 3, 43, 11, 51, 19, 59, 27,
            34, 2, 42, 10, 50, 18, 58, 26,
            33, 1, 41, 9, 49, 17, 57, 25,
            32, 0, 40, 8, 48, 16, 56, 24]

    def __init__(self, key):
        """
        Sanity checking of arguments.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, hex2bin
            sage: hexkey='0123456789abcdef'
            sage: key=hex2bin(hexkey)
            sage: DES(key)
            DES cryptosystem with 64 bit plaintext and key

        TESTS::

            sage: DES('012')
            Traceback (most recent call last):
            ...
            ValueError: invalid key - must be binary

            sage: DES('010101010101')
            Traceback (most recent call last):
            ...
            ValueError: invalid DES key size. Key must be exactly 64 bits long
        """
        if set(map(Integer, key)) != set([0, 1]):
            raise ValueError("invalid key - must be binary")
        if len(key) != 64:
            raise ValueError("invalid DES key size. Key must be exactly 64 bits long")
        self.__key = key

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, hex2bin
            sage: hexkey='0123456789abcdef'
            sage: key=hex2bin(hexkey)
            sage: DES(key)
            DES cryptosystem with 64 bit plaintext and key
        """
        return "DES cryptosystem with 64 bit plaintext and key"

    def keyschedule(self):
        """
        Obtains all subkeys from an initial 64 bit key k

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, hex2bin
            sage: hexkey='0123456789abcdef'
            sage: key=hex2bin(hexkey)
            sage: D = DES(key)
            sage: D.keyschedule()[0]
            '000010110000001001100111100110110100100110100101'
        """
        tmp = join([self.__key[i] for i in self.__PC1], '')  # Permute the key bits according to PC1
        ks = []
        shifts = [1, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1]
        for i in shifts:
            c = tmp[:28]
            d = tmp[28:]
            cs = c[i:] + c[:i]
            ds = d[i:] + d[:i]
            tmp = cs + ds
            h = join([tmp[i] for i in self.__PC2], '')
            ks.append(h)
            # print hex(Integer('0b'+h)).zfill(12)
        return ks

    def xor(self, str1, str2):
        """
        returns the string made by bit-wise xor of the bits of str1 and str2

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, hex2bin
            sage: hexkey='0123456789abcdef'
            sage: key=hex2bin(hexkey)
            sage: D = DES(key)
            sage: D.xor('1001','1100')
            '0101'
        """
        from operator import __xor__
        return join([str(__xor__(Integer(x), Integer(y)))
                     for (x, y) in zip(str1, str2)], '')

    def apply_sbox(self, str, n):
        """
        Applies sbox n to six bit string str

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, hex2bin
            sage: hexkey='0123456789abcdef'
            sage: key=hex2bin(hexkey)
            sage: D = DES(key)
            sage: D.apply_sbox('100110', 4)
            '1011'
        """
        s = map(Integer, str)
        row = 2 * s[0] + s[5]
        col = 8 * s[1] + 4 * s[2] + 2 * s[3] + s[4]
        return bin(self.__sbox[n][16 * row + col])[2:].zfill(4)

    def _mixright(self, R0, keys, m):
        """
        This is the function which mixes the 32 bit right half with subkey m

        EXAMPLES::
        """
        R1 = join([R0[i] for i in self.__E], '')  # First expand to 48 bits
        R2 = self.xor(R1, keys[m])              # Then add the current subkey
        R3 = join([self.apply_sbox(R2[6 * i: 6 * i + 6], i)
                   for i in range(8)], '')
                                       # Apply the sboxes to each six bits
        R4 = join([R3[i] for i in self.__P], '')  # Then reduce to 32 bits
        return R4

    def encrypt(self, plaintext):
        r"""
        Encrypts a single plaintext with the instance of DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, bin2hex
            sage: k1=8*'01010011'
            sage: k2=8*'10011001'
            sage: D1=DES(k1)
            sage: D2=DES(k2)
            sage: pl=64*'0'
            sage: ct1=D1.encrypt(pl)
            sage: bin2hex(ct1)
            '17d819b45d919a56'
            sage: ct2=D2.encrypt(pl)
            sage: bin2hex(ct2)
            '0f2fcf4aeb6c56d4'
        """
        # Sanity checking of arguments.
        if len(plaintext) != 64:
            raise ValueError("Invalid DES plaintext size. Plaintext must be exactly 64 bits long.")
        pl = join([plaintext[i] for i in self.__IP], '')
            # first permute the bits of the plaintext
        ks = self.keyschedule()  # Create all the subkeys
        L = pl[:32]  # Obtain left and right halves
        R = pl[32:]
        for i in range(16):
            tmp = R
            R = self.xor(self._mixright(tmp, ks, i), L)
            L = tmp
    #       print "i,L+R=",i,hex(Integer('0b'+L)).zfill(8),hex(Integer('0b'+R)).zfill(8)
        out = R+L
        return join([out[i] for i in self.__FP], '')

    def decrypt(self, ciphertext):
        r"""
        Decrypts a single ciphertext with the instance of DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, bin2hex
            sage: k1=8*'01010011'
            sage: k2=8*'10011001'
            sage: D1=DES(k1)
            sage: D2=DES(k2)
            sage: pl=64*'0'
            sage: ct1=D1.encrypt(pl)
            sage: bin2hex(ct1)
            '17d819b45d919a56'
            sage: ct2=D2.encrypt(pl)
            sage: bin2hex(ct2)
            '0f2fcf4aeb6c56d4'
            sage: de1=D1.decrypt(ct1)
            sage: bin2hex(de1)
            '0000000000000000'
            sage: de2=D2.decrypt(ct2)
            sage: bin2hex(de2)
            '0000000000000000'
        """
        # Sanity checking of arguments.
        if len(ciphertext) != 64:
            raise ValueError("Invalid DES ciphertext size. Ciphertext must be exactly 64 bits long.")
        ct = join([ciphertext[i] for i in self.__IP], '')  # first permute the bits of the ciphertext
        ks = self.keyschedule()  # Create all the subkeys
        L = ct[:32]  # Obtain left and right halves
        R = ct[32:]
        for i in range(16):
            tmp = R
            R = self.xor(self._mixright(tmp, ks, 15-i), L)
                # Apply subkeys in opposite order from encryption
            L = tmp
            # print "i,L+R=",i,hex(Integer('0b'+L)).zfill(8),hex(Integer('0b'+R)).zfill(8)
        out = R+L
        return join([out[i] for i in self.__FP], '')

    def getkey(self):  # Returns the key used
        r"""
        Returns the current key.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, bin2hex
            sage: k1=8*'01010011'
            sage: D1=DES(k1)
            sage: D1.getkey()
            '0101001101010011010100110101001101010011010100110101001101010011'
            sage: bin2hex(_)
            '5353535353535353'
        """
        return self.__key

    def setkey(self, newkey):  # Changes the key for this instance of DES
        r"""
        Changes the key for this instance of DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES, bin2hex
            sage: k1=8*'01010011'
            sage: D1=DES(k1)
            sage: D1.setkey(8*'10011001')
            sage: bin2hex(D1.getkey())
            '9999999999999999'
        """
        self.__key = newkey

    def avalanche(self, plaintext1, plaintext2):
        r"""
        Demontrates the avalanche effect for two different plaintexts.

        Given two plaintexts pl1 and pl2, then the call for an instance
        D1 of DES as ``D1.avalanche(pl1,pl2)`` will print 16 rows of
        values, one for each round of DES.  Each row
        contains the round number, the number of differences between the
        current value for each plaintext, and a string representing the
        bit-wise exclusive-ors between the two strings.  So the number
        of differences is just the number of ones in the exclusive-or string.

        The avalanche effect is best shown when the two plaintexts differ
        at a single position.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.des import DES
            sage: k1=8*'01010011'
            sage: D1=DES(k1)
            sage: pl1=8*'10011100'
            sage: pl2=7*'10011100'+'10111011'
            sage: D1.avalanche(pl1, pl2)
            1 13 0000000010000000000000001000000001000010001110011000010110001000
            2 33 0100001000111001100001011000100001001111111001111100011011111101
            3 37 0100111111100111110001101111110111100011100111001100001100011000
            4 28 1110001110011100110000110001100001001100010010000010011010100111
            5 31 0100110001001000001001101010011100001101010101110001111111101100
            6 35 0000110101010111000111111110110001101111110001001100011010001101
            7 34 0110111111000100110001101000110110001111100111111001101100001000
            8 30 1000111110011111100110110000100010100001100010110000010010110110
            9 27 1010000110001011000001001011011000110011111000000110010001101001
            10 27 0011001111100000011001000110100110110100101110001000100100100010
            11 26 1011010010111000100010010010001001111010011001100100100000010100
            12 29 0111101001100110010010000001010000101011010000011011011010111001
            13 29 0010101101000001101101101011100110000001011010000101100101010011
            14 25 1000000101101000010110010101001111100001000110110000001011000001
            15 27 1110000100011011000000101100000101100111001101011000010010110010
            16 33 0110011100110101100001001011001000100001110011011001111111110100
            ['1011110101110101010100111011111101011100110110011010101001011100', '0100100111110011111011101010101101111011001110100011101101000011']
            
        This shows that after the first few rounds half of the bits have
        changed.
        """
        # Checks the avalanche criterion for two different plaintexts
        pl1 = join([plaintext1[i] for i in self.__IP], '')
            # first permute the bits of the plaintext
        pl2 = join([plaintext2[i] for i in self.__IP], '')
            # first permute the bits of the plaintext
        ks = self.keyschedule()  # Create all the subkeys
        L1 = pl1[:32]  # Obtain left and right halves
        R1 = pl1[32:]
        L2 = pl2[:32]  # Obtain left and right halves
        R2 = pl2[32:]
        for i in range(16):
            tmp1 = R1
            tmp2 = R2
            R1 = self.xor(self._mixright(tmp1, ks, i), L1)
            R2 = self.xor(self._mixright(tmp2, ks, i), L2)
            L1 = tmp1
            L2 = tmp2
            X = self.xor(L1+R1, L2+R2)
            wt = sum(Integer(i) for i in X)
            print i+1, wt, X
        out1 = R1+L1
        out2 = R2+L2
        return [join([out1[i] for i in self.__FP], ''), join([out2[i] for i in self.__FP], '')]

    # def avalanche_k(plaintext, key1, key2):
    #     """
    #     Avalanche with two separate keys - not yet properly implemented
    #     """
    #     pl = join([plaintext[i] for i in IP],'') # first permute the bits of the plaintext
    #     ks1=keyschedule(key1) # Create all the subkeys
    #     ks2=keyschedule(key2) # Create all the subkeys
    #     L1=pl[:32] # Obtain left and right halves
    #     R1=pl[32:]
    #     L2=pl[:32]
    #     R2=pl[32:]
    #     for i in range(16):
    #         tmp1=R1
    #         tmp2=R2
    #         R1=xor(f(tmp1,ks1,i),L1)
    #         R2=xor(f(tmp2,ks2,i),L2)
    #         L1=tmp1
    #         L2=tmp2
    #         X=xor(L1+R1,L2+R2)
    #         wt=sum(Integer(i) for i in X)
    #         print i+1,wt,X
    #     out1=R1+L1
    #     out2=R2+L2
    #     return [join([out1[i] for i in FP],''),join([out2[i] for i in FP],'')]
