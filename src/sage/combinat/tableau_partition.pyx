from itertools import permutations, combinations, product
from sage.combinat.partition import Partition
from sage.combinat.tableau import StandardTableaux, Tableau
from sage.misc.misc_c import prod
from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

def matchings(tuple a, tuple b):
    """
    if a = (i1,...,il), b=(j1,...,jm) generate all (k1,...,km) sublists of a
    such that kr<jr for r=1,...,m.
    """
    cdef unsigned short int k, l
    cdef tuple K, L
    l = len(a)
    if len(b)==0:
        yield tuple()
    else:
        for k in range(l):
            if a[k]<b[0]:
                aa = tuple(ii for ii in a if ii!=a[k])
                for K in matchings(aa,b[1:]):
                    L = (a[k],)+K
                    yield L
            else:
                break

def set_partitions(la):
    la = Partition(la)
    for T in StandardTableaux(la):
        for sp in tableau_set_partitions(T):
            yield sp

cdef list extend_set_partition(tuple J0,tuple J1,tuple K, sp):
    cdef int i
    cdef list newsp = []
    for i in J0:
        if i in K:
            newsp.append((i,)+sp[K.index(i)])
        else:
            newsp.append((i,))
    return newsp

def tableau_set_partitions(T):
    T = Tableau(T)
    return tableau_set_partitions_conj(T.conjugate())

def tableau_set_partitions_conj(T):
    if len(T)==1:
        yield([(a,) for a in T[0]])
    else:
        for K in matchings(T[0],T[1]):
            for sp in tableau_set_partitions_conj(T[1:]):
                yield extend_set_partition(T[0],T[1],K,sp)

def non_interlacing_matching(tuple r,tuple rr):
    K = list()
    L = list(r)
    for j in rr:
        m = max(a for a in L if a<j)
        K.append(m)
        L.remove(m)
    return K

def non_interlacing_set_partition(T):
    T = Tableau(T)
    return non_interlacing_set_partition_conj(T.conjugate())

def non_interlacing_set_partition_conj(T):
    """
    Return the non-interlacing set partitions corresponding to ``T``.
    """
    if len(T)==1:
        return [(a,) for a in T[0]]
    else:
        K = non_interlacing_matching(T[0],T[1])
        sp = non_interlacing_set_partition_conj(T[1:])
        newsp = []
        for i in T[0]:
            if i in K:
                newsp.append((i,)+sp[K.index(i)])
            else:
                newsp.append((i,))
        return newsp
    

