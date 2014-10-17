from  sage_coding_project.RS import *
from  sage_coding_project.BlockCode import *
import random

def key_gen(C):

    G = C.generator_matrix()
    S = random_matrix(C.F, C.k, C.k)
    while S.is_singular() == True:
        S = random_matrix(C.F, C.k, C.k)

    x = [0..C.n-1]
    random.shuffle(x)
    P=Matrix(C.F, C.n, C.n)
    for i in x:
        P[i,x[i]] = 1
    
    t = floor((C.n-C.k)/2) #I should implement a decoding capacity at some point, when
                  # I'll know how this works in a more specific way

    Gp = S*G*P

    pubkey = (Gp, t)
    privkey = (S.inverse(), G, P.inverse(), C)


    return pubkey, privkey


def encrypt(m, pubkey):

    Gp = pubkey[0]
    c = vector(m)*Gp
    e = vector(Gp.base_ring(), Gp.ncols())
    rand = []
    w = 0
    while w !=pubkey[1]:
        r = randint(0, Gp.ncols()-1)
        if r not in rand:
            w += 1
            rand.append(r)

    for i in range(0, pubkey[1]):
        e[rand[i]] = randint(1, Gp.ncols()-1)
    
    c  = c + e

    return c

def decrypt(c, privkey):
    S,G,P,C = privkey
    c = c*P
    c = C.decode_welch_berlekamp(c)
    m = C.unencode(c)
    m = m*S
    return m

