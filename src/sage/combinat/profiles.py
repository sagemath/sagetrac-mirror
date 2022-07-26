r"""
Set Partitions and Profiles.

AUTHOR: 

- Amritanshu Prasad (amri@imsc.res.in)

This module provides functions that help verify and understand the manuscript Set Partitions, Tableaux, and Subspace Profiles of Regular Diagonal Operators by Amritanshu Prasad and Samrith Ram.

EXAMPLES:

Computation of `c_q(T)`::

    sage: T = Tableau([[1,2,9],[3,6],[5,8]])
    sage: tableau_stat(T)
    q^3 + 2*q^2 + 2*q + 1

Computation of `b_\lambda(q)`::

    sage: b([3,3,1])
    q^4 + 5*q^3 + 15*q^2 + 28*q + 21

Computation of `b_{(3,3,1)}(q)` using set partitions::

    sage: R.<q> = QQ['q']
    sage: sum(q^interlacing_number(AA) for AA in SetPartitions(7,[3,3,1]))

Computation of qq-Stirling numbers of the second kind using `b_\lambda(q)`::

    sage: R.<q> = QQ['q']
    sage: sum(q^sum(i*(l-1) for i,l in enumerate(la))*b(la) for la in Partitions(6, length=3))
    q^6 + 4*q^5 + 11*q^4 + 19*q^3 + 25*q^2 + 20*q + 10

We can check this against Sage's q-Stirling number function::

    sage: from sage.combinat.q_analogues import q_stirling_number2
    sage: q_stirling_number2(6,3)/q^binomial(3,2) # Sage uses a different normalization for q-Stirling numbers of the second kind
    q^6 + 4*q^5 + 11*q^4 + 19*q^3 + 25*q^2 + 20*q + 10

The number of subspaces of `\mathbf F_q^5` with profile `(3,1)` is::    

    sage: sigma(5,[3,1])
    5*q^3 + 5*q^2 + 5*q - 15

The number of subspaces of `\mathbf F_q^5` with partial profile `(3,1)` is::    

    sage: partial(5,[3,1])
    q^4 + 6*q^3 + 6*q^2 - 4*q - 9
"""
from itertools import permutations, combinations, product
from sage.combinat.partition import Partition
from sage.combinat.tableau import StandardTableaux
from sage.functions.other import binomial
from sage.misc.all import prod
from sage.rings.infinity import infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ

def q_int(n):
    """
    Return the `q`-analog of `n`.
    """
    R= PolynomialRing(QQ,'q')
    q = R.gen()
    return sum(q**i for i in range(n))

def c2(a,b):
    """
    Return `c_q(T)` for a two-row tableau with rows `a` and `b`.
    """
    return prod(q_int(sum(1 for i in a[k:] if i<j)) for k,j in enumerate(b)) 


def tableau_stat(T):
    r"""
    Return the statistic `c_q` associated to the conjugate of tableau ``T``.
    """
    # Easier to work with rows of tableaux
    def tableau_row_stat(T):
        r"""
        Return the statistic `c_q` associated to the conjugate of tableau ``T``.
        """
        if len(T)==1:
            return 1
        else:
            return prod(c2(a,T[i+1]) for i,a in enumerate(T[:-1]))
    return tableau_row_stat(T.conjugate())
    
def b(la):
    r"""
    Return the polynomial `b_\lambda(q)` associated to ``la``.
    """
    return sum((tableau_stat(T) for T in StandardTableaux(la)))

def arcs(A):
    """
    Generate the arcs of a set ``A`` of positive integers.
    """
    for i,a in enumerate(A[:-1]):
        yield (a,A[i+1])
    yield (A[-1],infinity)

def cross(arc1,arc2):
    """
    Decide whether or not arcs cross.
    """
    return (arc1[0]<arc2[0] and arc2[0]<arc1[1] and arc1[1]<arc2[1]) or (arc2[0]<arc1[0] and arc1[0]<arc2[1] and arc2[1]<arc1[1])

def interlacing_sets(A,B):
    """
    Return the number of interlacings of arcs of ``A`` and ``B``.
    """
    return sum(int(cross(*pair)) for pair in zip(arcs(A),arcs(B)))

def interlacing_number(AA):
    """
    Return the interlacing number of the set partition ``AA``.
    """
    AA = sorted([sorted(A) for A in AA]) # Sage sometimes gives set partitions in non-standard notation!
    return sum(interlacing_sets(A,B) for A,B in combinations(AA,2))

def sigma(n,mu):
    r"""
    Return the number of subspaces of `F_q^n` with profile ``mu``.
    """
    mu = Partition(mu)
    R = PolynomialRing(QQ,'q')
    q = R.gen()
    return binomial(n,mu.size())*(q-1)**sum(mu[1:])*q**sum(binomial(m,2) for m in mu[1:])*b(mu.conjugate())
        
