# -*- coding: utf-8 -*-
"""
McEliece Cryptosystem
"""
from sage.coding.channel import StaticRankErrorChannel
from sage.crypto.cryptosystem import PublicKeyCryptosystem
from sage.matrix.constructor import matrix, random_matrix
from sage.misc.prandom import randint, sample
from sage.modules.free_module_element import random_vector, vector
from sage.modules.free_module import span

class McElieceCryptosystem(PublicKeyCryptosystem):

    def __init__(self, code, encoder, decoder):
        """
        TESTS::

        This module uses the following experimental feature:
        This test block is here only to trigger the experimental warning so it does not
        interferes with doctests::

            sage: Fqm = GF(2^9)
            sage: Fq = GF(2^3)
            sage: C = codes.GabidulinCode(Fqm, 2, 2, Fq)
            sage: S.<x> = Fqm['x', C.twisting_homomorphism()]
            sage: z9 = Fqm.gen()
            sage: p = (z9^6 + z9^2 + z9 + 1)*x + z9^7 + z9^5 + z9^4 + z9^2
            sage: vector(p.multi_point_evaluation(C.evaluation_points()))
            doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
            See http://trac.sagemath.org/13215 for details.
            (z9^7 + z9^6 + z9^5 + z9^4 + z9 + 1, z9^6 + z9^5 + z9^3 + z9)

        EXAMPLES::
            sage: from sage.crypto.public_key.mceliece import McElieceCryptosystem
            sage: C = codes.GabidulinCode(GF(4096), 6, 2, GF(4))
            sage: E = codes.encoders.GabidulinVectorEvaluationEncoder(C)
            sage: D = codes.decoders.GabidulinGaoDecoder(C)
            sage: M = McElieceCryptosystem(C, E, D)
            sage: m = random_vector(C.base_field(), C.dimension())
            sage: x = M.encrypt(m)
            sage: M.decrypt(x) == m
            True
        """
        G = code.generator_matrix()
        decoding_radius = decoder.decoding_radius()
        S = random_matrix(code.base_field(), code.dimension(), code.dimension())
        while S.is_singular():
            S = random_matrix(code.base_field(), code.dimension(), code.dimension())
        x = randint(1, decoding_radius - 1)
        P = self.scrambler(code, x)
        G_pub = S*G + P
        self.P = P
        self.pubkey = G_pub, decoding_radius - x, encoder
        self.privkey = S.inverse(), decoder

    def scrambler(self, code, x):
        Fqm = code.base_field()
        Fq = code.sub_field()
        m = code.extension_degree()
        k = code.dimension()
        n = code.length()
        km = k * m
        first_col = random_vector(Fq, km)
        while first_col == vector(Fq, [0]*km):
            first_col = random_vector(Fq, km)
        M = [first_col]

        for i in range(x - 1):
            V = span(M)
            W = V.complement()
            a = W.random_element()
            while a == W.zero():
                a = W.random_element()
            b = V.random_element()
            while b == V.zero():
                b = V.random_element()
            M.append(a + b)

        for i in range(n - x):
            V = span(M)
            x = V.random_element()
            M.append(x)

        parts = [e[i:i+m] for e in M for i in range(0, len(e), m)]
        VS, to_big, from_big = Fqm.vector_space(Fq, map=True)
        mat = matrix(Fqm, n, k, [to_big(v) for v in parts])
        return mat.transpose()

    def encrypt(self, m):
        M, t, E = self.pubkey
        c = m * M
        ambient = E.code().ambient_space()
        channel = StaticRankErrorChannel(ambient, E.code().sub_field(), t)
        y = channel.transmit(c)
        return y

    def decrypt(self, y):
        S_inv, decoder = self.privkey
        m = decoder.decode_to_message(y)
        m = vector(m.coefficients()) #Gabidulin Gao decoder returns a polynomial
        m = m * S_inv
        return m
