# -*- coding: utf-8 -*-

from sage.matrix.constructor import matrix, random_matrix
from sage.misc.prandom import shuffle
from sage.coding.linear_code import LinearCode
from sage.coding.channel_constructions import StaticErrorRateChannel

class McElieceCryptosystem:
    r"""
    EXAMPLES::

        sage: C = codes.RandomLinearCode(11,4,GF(2))
        sage: en = C.encoder()
        sage: de = C.decoder()
        sage: crypt = codes.McElieceCryptosystem(C, en, de, 1)
        sage: m = vector([1,0,1,1])
        sage: c = crypt.encrypt(m)
        sage: r = crypt.decrypt(c); r
        (1, 0, 1, 1)
    """

    def __init__(self, code, encoder, decoder, decoding_radius):
        G = encoder.generator_matrix()
        S = random_matrix(code.base_field(), code.dimension(), code.dimension())
        while S.is_singular():
            S = random_matrix(code.base_field(), code.dimension(), code.dimension())
        sigma = range(0, code.length())
        shuffle(sigma)
        P = matrix(code.base_field(), code.length(), code.length())
        for i in sigma:
            P[i,sigma[i]] = 1
        G_pub = S*G*P
        C_pub = LinearCode(G_pub)
        E_pub = C_pub.encoder()  
        self.pubkey = E_pub, decoding_radius
        self.privkey = (S.inverse(), P.inverse(), decoder)
    
    def encrypt(self, m):
        E_pub, t = self.pubkey
        c = E_pub.encode(m)
        ambient = E_pub.code().ambient_space()
        channel = StaticErrorRateChannel(ambient, t)
        y = channel.transmit(c)
        return y
    
    def decrypt(self, y):
        S_inv, P_inv, decoder = self.privkey
        y = y * P_inv
        m = decoder.decode_to_message(y)
        m = m * S_inv
        return m
    
    
    
