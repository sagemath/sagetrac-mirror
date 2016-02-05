# -*- coding: utf-8 -*-
r"""
S-unit equations over `\QQ`

Let `K` be a number field and `S` a finite set of prime ideals of
`K`. Let `\mathcal O_{K,S}^*:=\{x\in K|v_{\mathfrak p}(x) =0` for all
`\mathfrak p\not\in S\}` be the group of `S`-units of K.  The `S`-unit
equation is the Diophantine equation `x+y=1` where  `x,y \in
\mathcal O_{K,S}^*`.  A classical reult is that this equation has only
finitely many solutions.  An algorithm to solve the equation was
developed by de Weger in his thesis [Weg88]_.

This module contains an implementation of the algorithm in the special
case `K = \QQ`, following [Weg88]_ and the exposition in Smart's book
[Sma98]_ (see also [Sma99]_).

AUTHORS:

- Angelos Koutsianas (2015-2016)
- John Cremona (2016)

EXAMPLES::

    sage: from sage.arith.S_unit_equation import solve_S_unit_equation
    sage: solve_S_unit_equation([])
    []

The function solve_S_unit_equation takes as input a list S of primes,
and returns the finite list of all solutions x; that is, all S-units x
such that y=1-x is also an S-unit::

    sage: solve_S_unit_equation([2,3,5])
    [1/16, 25/16, 1/4, 1/25, 25, 4, 16/25, 16, 125/128, 5/32, 5/8, 1/10, 5/2,
    2/5, 10, 8/5, 32/5, 128/125, -1/80, -5/4, -1/5, -5, -4/5, -80, -1/8, -1/2,
    -25/2, -2/25, -2, -8, 1/81, 3/128, 3/8, 2/27, 6, 32/27, 9/4, 4/9, 27/32,
    1/6, 27/2, 8/3, 128/3, -9/16, -1/9, -9, -16/9, -1/24, -2/3, -1/4, -4, -3/2,
    -24, 1/2, -1, 2, 15/16, -15, 16/15, -1/15, 9/25, 25/9, 3/4, -3, 4/3, -1/3,
    24/25, 25/24, -3/125, -125/3, -27/5, -5/27, -3/5, -5/3, 9/10, 10/9, 3/5, 5/3,
    81/80, 81, 80/81, 9/5, 5/9, 6/5, 5/6, 9/8, 9, 8/9, 1/9, 3/2, 3, 2/3, 1/3, 27/25,
    25/27, 5/4, 5, 4/5, 1/5]

    sage: solve_S_unit_equation([2,3,7])
    [1/64, 1/4, 1/49, 49, 4, 64, 7/16, 1/28, 7/4, 1/7, 7, 4/7, 28, 16/7, -49/32, -1/8,
    -1/2, -2, -8, -32/49, -7/2, -2/7, 81/32, 2/9, 9/16, 1/8, 8, 16/9, 9/2, 32/81, -1/27,
    -27, -6, -1/48, -4/3, -3/4, -48, -1/6, 1/2, -1, 2, 63/64, -63, 64/63, -1/63, 3/4, -3,
    4/3, -1/3, 48/49, 49/48, -9/7, -7/9, 27/28, 28/27, 3/7, 7/3, 6/7, 7/6, 81/49, 49/81,
    9/8, 9, 8/9, 1/9, 3/2, 3, 2/3, 1/3, 9/7, 7/9, 7/8, -7, 8/7, -1/7]

REFERENCES:

..  [Sma98] Nigel P. Smart. The Algorithmic Resolution of Diophantine Equations. Number 41 in Students Texts. London
    Mathematical Society, 1998.

..  [Sma99] SMART, N. , Determine the small solutions to S-unit equations, Mathematics of Computations 68(228):1687-1699,1999

..  [Weg88] B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

"""

#*****************************************************************************
#       Copyright (C) 2013 Angelos Koutsianas <koutsis.jr@hotmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR

def minimal_vector(A,y):
    r"""
    Return a lower bound for the square distance from a vector to a lattice.

    INPUT:

    - ``A`` : an square non-singular integer matrix whose rows generate a lattice `\mathcal L`

    - ``y`` : a row vector with integer coordinates

    OUTPUT:

    A low bound for the square of `\ell (\mathcal L,\vec y)
    =\begin{cases}\displaystyle\min_{\vec x\in\mathcal L} \Vert\vec
    x-\vec y\Vert &, \vec y\not\in\mathcal L. \\
    \displaystyle\min_{0\neq\vec x\in\mathcal L}\Vert\vec x
    \Vert&,\vec y\in\mathcal L.\end{cases}`

    .. NOTE:

        The algorithm is based on V.9 and V.10 of [Sma98]_.

    EXAMPLES::

        sage: from sage.arith.S_unit_equation import minimal_vector
        sage: B = matrix(ZZ,2,[1,1,1,0])
        sage: y = vector(ZZ,[2,1])
        sage: minimal_vector(B,y)
        1/2

        sage: B = random_matrix(ZZ,3)
        sage: y = vector([1,2,100])
        sage: minimal_vector(B,y) # random
        15/28
    """
    if A.is_singular():
        raise ValueError('The matrix A is singular')

    n = len(y)
    c1 = 2**(n-1)
    ALLL = A.LLL()
    ALLLinv = ALLL.inverse()
    ybrace = [(a-a.round()).abs() for a in y * ALLLinv if (a-round(a)) != 0]

    if len(ybrace) == 0:
        return (ALLL.rows()[0].norm())**2 / c1
    else:
        sigma = ybrace[len(ybrace)-1]
        return ((ALLL.rows()[0].norm())**2 * sigma) / c1


def initial_bound(S):
    r"""
    Return an upper bound for the exponents of solutions to the S-unit equation.
    
    INPUT:

    - ``S`` : a list of prime numbers
    
    OUTPUT:

    A (large) upper bound for the absolute value of the exponents of
    the solutions of the `S`-unit equation `x \pm y=1`. This is based
    on the theorem 6.1 of [Weg88]_.
        
    EXAMPLE::
        
        sage: from sage.arith.S_unit_equation import initial_bound
        sage: S = [2,7,11]
        sage: initial_bound(S)
        8.40119849273500e17
    """    
    C1t = [768523,476217,373024,318871,284931,261379,2770008]
    s = len(S)
    t = (2 * s)/3 # rounded down
    from sage.misc.all import prod
    P = prod(S)
    q = 0
    for p in S:
        m = p * (p-1)
        qi = ZZ(3)
        while qi.divides(m):
            qi = qi.next_prime()
        q = max(qi, q)
    q = RR(q)

    e = RR(1).exp()
    if t<8:
        a1 = (56 * e)/15
        c = C1t[t-2]
    else:
        a1 = (8 * e)/3
        c = C1t[6]

    m = max([((qi-1) * (2+1/(qi-1))**t)/((RR(qi).log())**(t+2)) for qi in S])
    # Does this assume that S[s-1] is the largest?
    smax = S[-1]
    smin = S[0]
    log_smax = RR(smax).log()
    log_smin = RR(smin).log()
    U = c * (a1**t) * (t**((t+5)/2)) * (q**(2*t)) * (q-1) * ((RR(t*q).log())**2) * m * (log_smax**t) * ((4*log_smax).log()+log_smax/(8*t))
    C1 = U/(6 * t)
    C2 = U * RR(4).log()
    Omega = 1
    for i in range(s-t,s):
        Vi = Vs_1 = Vs = RR(1)
        Vi = max(Vi, RR(S[i]).log())
        if i == s-2:
            Vs_1 = Vi
        if i == s-1:
            Vs = Vi
        Omega = Omega * Vi
    C3 = 2**(9 * t + 26) * t**(t + 4) * Omega * (1+Vs_1.log())
    P0 = RR(P)/smin
    log_P0 = P0.log()
    C4 = max(RR(7.4),  (C1 * log_P0 + C3)/log_smin)
    C6 = max((C2 * log_P0+C3 * (e * Vs).log()+0.327)/log_smin,
             (C2 * log_P0 + RR(2).log())/log_smin)
    C7 = 2 * (C6 + C4 * C4.log())
    C8 = max([smax, (2 * P0**smax).log()/log_smin, C2 + C1 * C7.log(), C7])
    return C8


def primitive_p_1_root_mod_pn(p,n):
    r"""
    Return a primitive `(p-1)`st root of unity in `\ZZ/p^n`.

    INPUT:

    - ``p`` : a prime number

    - ``n`` : a natural number

    OUTPUT:

    A primitive `(p-1)`-st root of unity `\mod p^n` if it exists, or 1 otherwise.

    EXAMPLES::

        sage: from sage.arith.S_unit_equation import primitive_p_1_root_mod_pn
        sage: primitive_p_1_root_mod_pn(5,1)
            2
        sage: primitive_p_1_root_mod_pn(11,3)
            596
    """
    from sage.rings.finite_rings.integer_mod import mod
    P = p**n
    if p == 2 and n > 2:
        return mod(1,P)

    from sage.arith.all import primitive_root
    ap = mod(primitive_root(p), P) # prim root mod p is enough here!
    for i in range(n-1):
        ap = ap**p

    return ap


def change_basis(v):
    r"""
    Return a unimodular matrix mapping ``v`` to its gcd.

    INPUT:

    - ``v`` : a list of integers

    OUTPUT:

    Let `v=[v_1,...,v_n]` and `g = \gcd(v) = l_1v_1+\cdots +l_nv_n`.
    The output is a unimodular integral whose last row is
    `[l_1,\cdots,l_n]`.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import change_basis
        sage: v = [2,11,4]
        sage: change_basis(v)
            [-11   2   0]
            [ 20  -4   1]
            [ -5   1   0]
    """
    n = len(v)
    from sage.matrix.all import matrix
    v = matrix(v)
    D,U,V = v.smith_form();
    V = V.transpose()
    t = V[0]
    V[0] = V[n-1]
    V[n-1] = t

    return V


def p_adic_approximation_of_a_homogenous_lattice(theta,p,m):
    r"""
    Return technical data needed in solving S-unit equations.

    INPUT:

    - ``theta`` : a list of `p`-adic numbers as defined in section 6.3
      page 121 of [Sma98]_

    - ``p`` : a prime number

    - ``m`` : the precision of the approximation

    OUTPUT:

    - The matrix `B_{\mu}` in the page 68 of [Weg88]_.

    - A copy of ``theta`` such that the last element has the minimal
      valuation.

    - the position of the element in ``theta`` that has the minimal
      valuation and was permuted with the last element of the theta.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import p_adic_approximation_of_a_homogenous_lattice
        sage: R = Qp(5,20,'capped-rel','series')
        sage: theta = [R(-6),R(7),R(14)]
        sage: p_adic_approximation_of_a_homogenous_lattice(theta,5,5)
            (
            [   1    0    0]
            [   0    1    0]
            [1044  522 3125], [4 + 2*5 + O(5^20), 2 + 5 + O(5^20), 4 + 3*5 + 4*5^2 +
            4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11
            + 4*5^12 + 4*5^13 + 4*5^14 + 4*5^15 + 4*5^16 + 4*5^17 + 4*5^18 + 4*5^19
            + O(5^20)],
            0
            )
    """
    n = len(theta)

    #we are going to find the theta_i with the smallest valuation

    ord = min([t.valuation() for t in theta])
    position = [i for i in range(n) if theta[i].valuation() == ord][0]

    # Swap theta[n-1] and theta[position], so that the theta_i
    # with minimal valuation is be last one
    a = theta[position]
    theta[position] = theta[n-1]
    theta[n-1] = a

    # Now the last element of theta has minimal valuation, and we
    # create the matix Bm as defined in De Weger's thesis (page 68)

    from sage.matrix.all import identity_matrix
    from sage.rings.finite_rings.integer_mod import mod

    Bm = copy(identity_matrix(n))
    P = p**m
    Bm[n-1] = [mod(-t/a, P) for t in theta]
    Bm[n-1,n-1] = P

    return Bm,theta,position


def a_base_for_Gmstar(B,A,p,m,m0):
    r"""
    Return a matrix whose columns generate a lattice needed in the algorithm.

    INPUT:

    - ``B`` : the matrix whose columns generate the lattice
      `\Gamma_{\mu}` as it is defined in page 68 of [Sma98]_

    - ``A`` : a list `[a_1,..,a_n]` such that `x=\prod a_i^{x_i}`
      where `a_i\in \QQ_p` and `v_p(a_i)=0`

    - ``p`` : a prime number

    - ``m`` : the precision of the lattice

    - ``m0``: the minimal order of `\log_p(a_i)` for `i=1,\cdots , n`

    OUTPUT:

    A matrix such that its columns generate the lattice `\Gamma_{\mu}`
    as in page 72 of [Weg88]_ when `p>3`

    .. NOTE:

        `v_p(\log_p(a_n))` must be minimal.

    EXAMPLES::

        sage: from sage.arith.S_unit_equation import a_base_for_Gmstar
        sage: B = matrix(ZZ,[[1,0],[1268,15625]])
        sage: Q5 = Qp(5, prec = 200, type = 'capped-rel', print_mode = 'series')
        sage: a_base_for_Gmstar(B,[Q5(3),Q5(2)],5,6,1)
        [    1     6]
        [ 1268 38858]

        sage: B = matrix(ZZ,[[1,0],[757,16384]])
        sage: Q2 = Qp(2, prec = 200, type = 'capped-rel', print_mode = 'series')
        sage: a_base_for_Gmstar(B,[Q2(5),Q2(3)],2,14,2)
        [    1     0]
        [  757 16384]
    """
    if p <= 3:
        return B

    n = len(A)
    zeta = primitive_p_1_root_mod_pn(p,m+m0)
    from sage.misc.all import prod
    xi = [prod([A[j]**B[j,i] for j in range(n)]) for i in range(n)]

    #xi has the values of the products Π ai^xi with x=bi
    #kbi has the values of the exponents k(bi)
    #zeta_powers contains the powers of zeta 

    P = p**(m+m0)
    from sage.rings.finite_rings.integer_mod import mod
    zeta_powers = [mod(zeta**i, P) for i in range(p-1)]
    kbi = [min([k for k in range(p-1) if (mod(xi[i], P)-zeta_powers[k]) == 0]) for i in range(n)]


    #V is the matrix which change the basis of Gamma from the basis b to the basis b'
    V = change_basis(kbi)

    #B2 is the matrix of the Gamma lattice with respect to the basis b'
    B2 = B * (V.inverse())
    from sage.matrix.all import matrix
    kbi = matrix(kbi).transpose()

    #kbi is containing the exponents of the new basis b'
    kbi = V*kbi
    B2 = B2.transpose()
    
    
    #we find bi* for i = 1 up to n-1 
    #Bstar is transpose matrix of the matrix that response to a basis for the Gm* sub-lattice of Gm.
     
    Bstar = matrix(ZZ,n)
    for i in range(n-1):
        a = mod(kbi[i][0] / kbi[n-1][0],(p-1)/2)
        gstar = a.lift_centered()
        Bstar[i] = B2[i]-gstar * B2[n-1]
    
    
    #we find bn*
    from sage.arith.all import lcm
    gstar = lcm(kbi[n-1][0],(p-1)/2)/kbi[n-1][0]
    Bstar[n-1] = gstar * B2[n-1]

    return Bstar.transpose()


def reducing_the_bound(X0,A,p,m):
    r"""
    Return a smaller uppoer bound for the exponent of `p` in a solution to an S-unit equation.

    INPUT:

    - ``X0`` : a big upper bound for the exponents

    - ``A`` : a list `[a_1,..,a_n]` such that `x=\prod a_i^x_i` where `a_i\in \QQ_p` and `v_p(a_i)=0`

    - ``p`` : a prime number

    - ``m`` : the precision of the lattice

    OUTPUT:

    - An new upper bound with respect to the prime ``p``.

    - A boolean variable that is True when the condition of lemma 3.14
      page 68 of [Weg88]_ holds.

    EXAMPLES::

        sage: from sage.arith.S_unit_equation import reducing_the_bound
        sage: Q2 = Qp(2, prec = 200, type = 'capped-rel', print_mode = 'series')
        sage: reducing_the_bound(294667190680076544,[Q2(3),Q2(5)],2,116)
        (294667190680076544, True)

        sage: reducing_the_bound(294667190680076544,[Q2(3),Q2(5)],2,121)
        (122, False)
    """
    n = len(A)
    A_log = [a.log() for a in A]
    Bm = p_adic_approximation_of_a_homogenous_lattice(A_log,p,m)
    A_log = Bm[1]
 
    pos = Bm[2]
    a = A[pos]
    A[pos] = A[n-1]
    A[n-1] = a
    m0 = A_log[n-1].valuation()

    #if p>3 we find a matrix for Gm* lattice. Otherwise Gm=Gm*
    Bmstar = a_base_for_Gmstar(Bm[0],A,p,m,m0)

    #We have to take the transpose of the matrix because of the
    #LLL() function
    Bmstar = Bmstar.transpose()

    #assume that the rows of the matrix generate the lattice
    C = Bmstar.LLL()
    from sage.modules.all import zero_vector
    e = copy(zero_vector(ZZ,n))
    e[0] = 1
    v = e * C
    vnorm = v.norm()**2
    if 2**(1-n) * vnorm > n * X0**2:
        increase_m = False
        X0 = (m-1+m0)
    else:
        increase_m = True

    return RR(X0).floor(),increase_m


def find_the_new_bound_for_all_primes(X0,A,precision):
    r"""
    Return a new exponent bound, given one bound.

    INPUT:

    - ``X0`` : an exponent upper bound for all the primes in A

    - ``A`` :a list of primes

    - ``precision`` : a working precision

    OUTPUT:

    A list with upper bounds for the exponents of each prime in ``A``.

    EXAMPLES::

        sage: from sage.arith.S_unit_equation import find_the_new_bound_for_all_primes
        sage: find_the_new_bound_for_all_primes(1000,[2,3,5],100)
        [24, 15, 10]

        sage: find_the_new_bound_for_all_primes(10000,[2,3,5,7,11,13],250)
        [85, 53, 37, 29, 24, 22]
    """
    B = [1] * len(A)
    for i,p in enumerate(A):
        #for its prime in A we are going to find a new bound

        from sage.rings.padics.all import Qp
        K = Qp(p, prec = precision, type = 'capped-rel', print_mode = 'series')

        #e = a vector with the primes different to p as Qp elements
        e = [K(a) for a in A if a != p]
        m0 = min([a.log().valuation() for a in e])
        m = (2 * RR(X0).log()/RR(p).log()).round()
        newbound = True
        while newbound:
            T = reducing_the_bound(X0,e,p,m)
            newbound = T[1]
            m += 1
            if m + m0 > K.precision_cap():
                # Sieve
                #if m is bigger than the precision we have, we have to increase it an evaluate all the p-adic numbers

                K = Qp(p, prec = 2 * K.precision_cap(), type = 'capped-rel', print_mode = 'series')
                e = [K(A[j]) for j in range(len(A)) if i != j]
        B[i] = T[0]

    return B


def applying_De_Weger_method(A,precision):
    r"""
    Return an upper bound for the exponents of primes in ``A`` for all S-unit solutions.

    INPUT:

    - ``A`` : a list of prime numbers

    - ``precision`` : `p`-adic precision to be used

    OUTPUT:

    An upper bound of the exponents of the primes in ``A``.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import applying_De_Weger_method
        sage: applying_De_Weger_method([2,3,5,11],200)
        32
    """
    X0 = RR(initial_bound(A)).floor()
    Xnew = max(find_the_new_bound_for_all_primes(X0,A,precision))
    while Xnew < X0:
        X0 = Xnew
        M = find_the_new_bound_for_all_primes(Xnew,A,precision);
        Xnew = max(M)
    return Xnew


def simple_loop(S,bounds):
    r"""
    Return a list of solutions `x` to the S-unit equation, given an exponent bound.

    INPUT:

    - ``S`` : a list of primes

    - ``bounds`` : a list of upper bounds of the absolute value of the
      exponents, or a single bound for all primes

    OUTPUT:

    A list of values `x` such that `(x,1-x)` give all solutions of the
    `S`-unit equation `x+y=1`, such that the absolute values of the
    exponents of `x,y` are smaller than ``B``

    .. NOTE:

        Here we use the fact that either `v_p(x)=v_p(y)<0` or
        `v_p(x)>0,v_p(y)=0` or `v_p(x)=0,v_p(y)>0` for all `p\in S`.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import simple_loop
        sage: simple_loop([2,3,5],12)
        [1/16, 15/16, -15, 16, 1/4, 3/4, -3, 4, 1/6, 5/6, -5, 6, 1/10, 9/10, -9,
        10, 1/2, 1/2, -1, 2, 1/81, 80/81, -80, 81, 1/9, 8/9, -8, 9, 1/3, 2/3,
        -2, 3, 1/25, 24/25, -24, 25, 1/5, 4/5, -4, 5]
    """
    solutions = []
    if not isinstance(bounds,list):
        bounds = [bounds]*len(S)
    from sage.misc.all import cartesian_product_iterator
    for v in cartesian_product_iterator([xrange(-b,b+1) if b !=0 else xrange(1) for b in bounds]):
        #for each candidate x we store the potential y in T
        T = [1]
        x = 1
        for pr,exp in zip(S,v):
            x = x * pr**exp
            temp = []
            for y in T:
                if exp < 0:
                    y = y * pr**exp
                    temp = temp + [y]
                elif exp == 0:
                    for j in range(bounds[S.index(pr)]+1):
                        temp = temp + [y * pr**j]
            T = temp
        for y in T:
            if x + y == 1:
                solutions.extend([x,y,-y/x,1/x])
    return solutions


def solve_S_unit_equation(S,precision = 200):
    r"""
    Return all solutions x to the S-unit equation x+y=1.

    INPUT:

    - ``S`` : a list of primes

    - ``precision`` : the precision for the calculations of the
      `p`-adic logarithms

    OUTPUT:

    All the `x` of the pairs of the solutions of the `S`-unit equation `x+y=1`.

    EXAMPLES::

        sage: from sage.arith.S_unit_equation import solve_S_unit_equation
        sage: solve_S_unit_equation([2,3], 100)
        [1/4, 4, -1/8, -1/2, -2, -8, 1/2, -1, 2, 3/4, -3, 4/3, -1/3, 9/8, 9,
        8/9, 1/9, 3/2, 3, 2/3, 1/3]

        sage: solve_S_unit_equation([2,3,5], 100)
        [1/16, 25/16, 1/4, 1/25, 25, 4, 16/25, 16, 125/128, 5/32, 5/8, 1/10,
        5/2, 2/5, 10, 8/5, 32/5, 128/125, -1/80, -5/4, -1/5, -5, -4/5, -80,
        -1/8, -1/2, -25/2, -2/25, -2, -8, 1/81, 3/128, 3/8, 2/27, 6, 32/27,
        9/4, 4/9, 27/32, 1/6, 27/2, 8/3, 128/3, -9/16, -1/9, -9, -16/9, -1/24,
        -2/3, -1/4, -4, -3/2, -24, 1/2, -1, 2, 15/16, -15, 16/15, -1/15, 9/25,
        25/9, 3/4, -3, 4/3, -1/3, 24/25, 25/24, -3/125, -125/3, -27/5, -5/27,
        -3/5, -5/3, 9/10, 10/9, 3/5, 5/3, 81/80, 81, 80/81, 9/5, 5/9, 6/5, 5/6,
        9/8, 9, 8/9, 1/9, 3/2, 3, 2/3, 1/3, 27/25, 25/27, 5/4, 5, 4/5, 1/5]

    TESTS:

    There are no solutions when 2 is not in ``S``::

        sage: solve_S_unit_equation([3,5], 100)
        []
        sage: solve_S_unit_equation([3,5,7], 100)
        []
    """
    if len(S)==0:
        return []
    if len(S)==1:
        p = S[0]
        if p==2:
            return [QQ(1/2), QQ(-1), QQ(2)]
        else:
            return []
    try:
        bad = any([not p.is_prime() for p in S])
    except AttributeError:
        bad = True
    if bad:
        raise ValueError('S contains nonprime values')

    S.sort()

    #we find an upper bound
    B = applying_De_Weger_method(S,precision)
    return sieve_S_unit_equation_over_Q(S,B,precision)


def trivial_Tp_finite_place_over_Q(S,p,bounds,delta,precision):
    r"""Return whether or not a technical inequality has solutions.

    INPUT:

    - ``S``: a list of rational primes

    - ``p``: a rational prime

    - ``bounds``: a list of upper bounds for the exponents of the
      primes in ``S``

    - ``delta``: a real number less than 1

    OUTPUT:

    True, if the inequality `|x-1|_{p}<\delta` does not have
    non-trivial solutions with `x=\prod_{p_i\neq p} p_i^{e_i}` for
    `p_i\in S` and `|e_i|\leq B`. Otherwise, False.

    .. NOTE:

        Here we implement paragraph 3.2 of [Sma99]_.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import trivial_Tp_finite_place_over_Q
        sage: trivial_Tp_finite_place_over_Q([2,3,5],5,[24, 15, 10],1/1000000,200)
        True

    """
    if p in S:
        i = S.index(p)
        S.remove(p)
        bounds.pop(i)

    deltaprime = (-(RR(delta).log()))/RR(p).log()
    if deltaprime < 1:
        return False

    from sage.rings.padics.all import Qp
    K = Qp(p, prec = precision, type = 'capped-rel', print_mode = 'series')
    LogS = [K(x).log() for x in S]

    c6 = min([x.valuation() for x in LogS])
    lam = p**c6
    M = [x/lam for x in LogS]
    c7 = RR(deltaprime - c6).floor()

    from sage.matrix.all import matrix, identity_matrix, zero_matrix, block_matrix
    from sage.modules.all import zero_vector
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    A = copy(identity_matrix(ZZ,len(S)))
    B = copy(zero_matrix(ZZ,len(S),1))

    u = 1
    while u <= c7:
        if u > precision:
            K = Qp(p, prec = precision, type = 'capped-rel', print_mode = 'series')
            LogS = [K(x).log() for x in S]
            M = [x/lam for x in LogS]
        P = p**u
        ZmodP = Zmod(P)
        C = matrix(ZZ,[ZmodP(x) for x in M]+[P])
        L = block_matrix([[block_matrix([[A,B]])],[C]])
        l = minimal_vector(L.transpose(),zero_vector(ZZ,len(S)+1))
        if l**2 > sum([b**2 for b in bounds]):
            return True
        u += 1
    return False


def sieve_S_unit_equation_over_Q(S,B,precision):
    r"""
    Return a list of solutions `x` to the S-unit equation, given an exponent bound.

    INPUT:

    - ``S`` : a list of primes

    - ``B`` : an upper bound of the absolute value of the exponents
      for the solutions

    - ``precision`` : the precision for the calculations of `p`-adic
      logarithms

    OUTPUT:

    A list of `x` of all the solutions of the `S`-unit group `x+y=1`

    .. NOTE:

        The sieve is based on the ideas of [Sma99]_.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import sieve_S_unit_equation_over_Q
        sage: sieve_S_unit_equation_over_Q([2,3,5],24,100)
        [1/16, 25/16, 1/4, 1/25, 25, 4, 16/25, 16, 125/128, 5/32, 5/8, 1/10, 5/2,
        2/5, 10, 8/5, 32/5, 128/125, -1/80, -5/4, -1/5, -5, -4/5, -80, -1/8, -1/2,
        -25/2, -2/25, -2, -8, 1/81, 3/128, 3/8, 2/27, 6, 32/27, 9/4, 4/9, 27/32,
        1/6, 27/2, 8/3, 128/3, -9/16, -1/9, -9, -16/9, -1/24, -2/3, -1/4, -4, -3/2,
        -24, 1/2, -1, 2, 15/16, -15, 16/15, -1/15, 9/25, 25/9, 3/4, -3, 4/3, -1/3,
        24/25, 25/24, -3/125, -125/3, -27/5, -5/27, -3/5, -5/3, 9/10, 10/9, 3/5, 5/3,
        81/80, 81, 80/81, 9/5, 5/9, 6/5, 5/6, 9/8, 9, 8/9, 1/9, 3/2, 3, 2/3, 1/3,
        27/25, 25/27, 5/4, 5, 4/5, 1/5]
    """
    #I define Q as a number field
    bounds = len(S) * [B]

    for i,p in enumerate(S):
        B0 = 1
        B1 = bounds[i]
        Bmid = ((B0+B1)/2).floor()
        Scopy = copy(S)
        boundscopy = copy(bounds)
        while (B0-B1).abs() > 1:
            delta = (p**(-Bmid))
            if trivial_Tp_finite_place_over_Q(Scopy,p,boundscopy,delta,precision):
                B1 = Bmid
                Bmid = ((B0+B1)/2).floor()
            else:
                B0 = Bmid
                Bmid = ((B0+B1)/2).floor()

        bounds[i] = B1
    Sol = []

    #we find solutions which are divible by a `suitable' high power of a prime in S and are positive.
    #We do that only for odd primes

    for p in S:
        if p%2 == 1:
            m = max(1,RR(((bounds[S.index(p)]).log()-RR(2).log()-RR(p-1).log())/RR(p).log()).round())
            Sol += solutions_divible_by_higher_power_of_p(S,p,bounds,m)
            bounds[S.index(p)] = m-1
    Sol += simple_loop(S,bounds)
    for l in Sol:
        temp = [l,1-l,1/l,1-1/l,1/(1-l),l/(l-1)]
        for a in temp:
            if a not in Sol:
                Sol.append(a)
    solutions =[]
    for a in Sol:
        if a not in solutions:
            solutions.append(a)

    return solutions


def solutions_divible_by_higher_power_of_p(S,p,bounds,m):
    r"""
    Return all solutions to the S-unit equation with `p^m\mid y`.

    INPUT:

    - ``S`` : a set of primes which contains ``p``

    - ``p`` : a prime number in ``S``

    - ``bounds`` : a list of natural numbers

    - ``m`` : a natural number

    OUTPUT:

    All `x` of the solutions of the `S`-unit equation `x+y=1` such that `y` is divisible by `p^m`.

    EXAMPLE::

        sage: from sage.arith.S_unit_equation import solutions_divible_by_higher_power_of_p
        sage: solutions_divible_by_higher_power_of_p([2,3,5],3,[24, 15, 10],2)
        [5/32, 10, 25/16, 16/25, 1/10, 32/5, -2/25, -1/80, -4/5, -1/8, -8, -5/4, -80, -25/2]
    """
    if p not in S:
        raise ValueError('p not in S')
    Snp = copy(S)
    Bnp = copy(bounds)
    Snp.remove(p)
    Bnp.remove(Bnp[S.index(p)])

    P = p**m
    from sage.rings.finite_rings.integer_mod_ring import Zmod
    Zpm = Zmod(P)
    order = Zpm.unit_group_order()

    from sage.modules.all import vector
    #we find the powers of primes in Snp that are congruence to 1 mod p^m
    A = vector(ZZ,[Zpm(q).log() for q in [-1]+Snp])
    B = [b//order for b in Bnp]
    X0 = []
    V = []
    from sage.misc.all import cartesian_product_iterator, prod
    for v in cartesian_product_iterator([xrange(2)]+[xrange(order)]*len(Bnp)):
        if (vector(v)*A)%order == 0:
            X0.append(QQ(prod([p**e for e,p in zip(v,[-1]+Snp)])))
            V.append(vector(v))

    #we find the solutions x such that x is congruence to 1 mod p**m
    #applying Smart's ideas using Fincke-Pohst algorithm

    solutions = []
    Snp_power = [q**order for q in Snp]
    for x0 in X0:
        for v in cartesian_product_iterator([xrange(-b-1,b+1) for b in B]):
            x1 = prod([q**e for e,q in zip(v,Snp_power)])
            if (QQ(1-x0*x1)).prime_to_S_part(S).abs() == 1:
                solutions.append(x0*x1)

    return solutions
