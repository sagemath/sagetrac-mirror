r"""
Dimensions of spaces of modular forms

AUTHORS:

- William Stein

- Jordi Quer

ACKNOWLEDGEMENT: The dimension formulas and implementations in this
module grew out of a program that Bruce Kaskel wrote (around 1996)
in PARI, which Kevin Buzzard subsequently extended. I (William
Stein) then implemented it in C++ for Hecke. I also implemented it
in Magma. Also, the functions for dimensions of spaces with
nontrivial character are based on a paper (that has no proofs) by
Cohen and Oesterle (Springer Lecture notes in math, volume 627,
pages 69-78). The formulas for `\Gamma_H(N)` were found
and implemented by Jordi Quer.

The formulas here are more complete than in Hecke or Magma.

Currently the input to each function below is an integer and either a Dirichlet
character `\varepsilon` or a finite index subgroup of `{\rm SL}_2(\ZZ)`.
If the input is a Dirichlet character `\varepsilon`, the dimensions are for
subspaces of `M_k(\Gamma_1(N), \varepsilon)`, where `N` is the modulus of
`\varepsilon`.

These functions mostly call the methods dimension_cusp_forms,
dimension_modular_forms and so on of the corresponding congruence subgroup
classes.
"""

##########################################################################
#       Copyright (C) 2004,2005,2006,2007,2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##########################################################################


from sage.rings.arith import (factor, is_prime,
                              valuation, kronecker_symbol, gcd, euler_phi, lcm, divisors, prime_divisors, jacobi_symbol)

from sage.misc.misc import mul
from sage.rings.all import Mod, Integer, IntegerModRing, ZZ
from sage.rings.rational_field import frac, QQ
from sage.rings.rational import Rational
import dirichlet
Z = ZZ  # useful abbreviation.

from sage.modular.arithgroup.all import Gamma0, Gamma1, is_ArithmeticSubgroup, is_GammaH
from sage.functions.other import sqrt
from sage.modular.dirichlet import DirichletGroup, trivial_character

##########################################################################
# Helper functions for calculating dimensions of spaces of modular forms
##########################################################################

def eisen(p):
    """
    Return the Eisenstein number `n` which is the numerator of
    `(p-1)/12`.

    INPUT:


    -  ``p`` - a prime


    OUTPUT: Integer

    EXAMPLES::

        sage: [(p,sage.modular.dims.eisen(p)) for p in prime_range(24)]
        [(2, 1), (3, 1), (5, 1), (7, 1), (11, 5), (13, 1), (17, 4), (19, 3), (23, 11)]
    """
    if not is_prime(p):
        raise ValueError("p must be prime")
    return frac(p-1,12).numerator()

##########################################################################
# Formula of Cohen-Oesterle for dim S_k(Gamma_1(N),eps).  REF:
# Springer Lecture notes in math, volume 627, pages 69--78.  The
# functions CO_delta and CO_nu, which were first written by Kevin
# Buzzard, are used only by the function CohenOesterle.
##########################################################################

def CO_delta(r,p,N,eps):
    r"""
    This is used as an intermediate value in computations related to
    the paper of Cohen-Oesterle.

    INPUT:


    -  ``r`` - positive integer

    -  ``p`` - a prime

    -  ``N`` - positive integer

    -  ``eps`` - character


    OUTPUT: element of the base ring of the character

    EXAMPLES::

        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CO_delta(1,5,7,eps^3)
        2
    """
    if not is_prime(p):
        raise ValueError("p must be prime")
    K = eps.base_ring()
    if p%4 == 3:
        return K(0)
    if p==2:
        if r==1:
            return K(1)
        return K(0)
    # interesting case: p=1(mod 4).
    # omega is a primitive 4th root of unity mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p-1)//4)
    # this n is within a p-power root of a "local" 4th root of 1 modulo p.
    n = Mod(int(omega.crt(Mod(1,N//(p**r)))),N)
    n = n**(p**(r-1))   # this is correct now
    t = eps(n)
    if t==K(1):
        return K(2)
    if t==K(-1):
        return K(-2)
    return K(0)

def CO_nu(r, p, N, eps):
    r"""
    This is used as an intermediate value in computations related to
    the paper of Cohen-Oesterle.

    INPUT:


    -  ``r`` - positive integer

    -  ``p`` - a prime

    -  ``N`` - positive integer

    -  ``eps`` - character


    OUTPUT: element of the base ring of the character

    EXAMPLES::

        sage: G.<eps> = DirichletGroup(7)
        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CO_nu(1,7,7,eps)
        -1
    """
    K = eps.base_ring()
    if p%3==2:
        return K(0)
    if p==3:
        if r==1:
            return K(1)
        return K(0)
    # interesting case: p=1(mod 3)
    # omega is a cube root of 1 mod p.
    omega = (IntegerModRing(p).unit_gens()[0])**((p-1)//3)
    n = Mod(omega.crt(Mod(1,N//(p**r))), N)  # within a p-power root of a "local" cube root of 1 mod p.
    n = n**(p**(r-1))  # this is right now
    t = eps(n)
    if t==K(1):
        return K(2)
    return K(-1)

def CohenOesterle(eps, k):
    r"""
    Compute the Cohen-Oesterle function associate to eps, `k`.
    This is a summand in the formula for the dimension of the space of
    cusp forms of weight `k` with character `\varepsilon`.

    INPUT:


    -  ``eps`` - Dirichlet character. Its modulus must be divisible by 4
        if k is half-integral.

    -  ``k`` - integer, or half an integer


    OUTPUT: element of the base ring of eps.

    EXAMPLES::

        sage: G.<eps> = DirichletGroup(7)
        sage: sage.modular.dims.CohenOesterle(eps, 2)
        -2/3
        sage: sage.modular.dims.CohenOesterle(eps, 4)
        -1
        sage: G = DirichletGroup(108)
        sage: psi = G.0*G.1^9
        sage: sage.modular.dims.CohenOesterle(psi, 3/2)
        -6
        sage: H = DirichletGroup(576)
        sage: chi = H.0*H.2^3
        sage: sage.modular.dims.CohenOesterle(chi, 1/2)
        -24
    """
    k = QQ(k)
    den = abs(k.denominator())
    if den > 2:
        raise TypeError("The weight must be an integer or half an integer")
    N = eps.modulus()
    f = eps.conductor()
    K = eps.base_ring()
    L = prime_divisors(N)
    r = {}
    s = {}
    for p in L:
        r[p] = valuation(N, p)
        s[p] = valuation(f, p)

    def _lambda(r,s,p):
        """
        Used internally by the CohenOesterle function.

        INPUT:


        -  ``r, s, p`` - integers


        OUTPUT: Integer

        EXAMPLES: (indirect doctest)

        ::

            sage: K = CyclotomicField(3)
            sage: eps = DirichletGroup(7*43,K).0^2
            sage: sage.modular.dims.CohenOesterle(eps,2)
            -4/3
            sage: G = DirichletGroup(108)
            sage: eps = G.0*G.1^9
            sage: sage.modular.dims.CohenOesterle(eps,3/2)
            -6
        """
        if 2*s<=r:
            if r%2==0:
                return p**(r//2) + p**((r//2)-1)
            return 2*p**((r-1)//2)
        return 2*(p**(r-s))
    #end def of lambda

    # We first consider the case of integral weight
    if den == 1:
        gamma_k = 0
        if k%4==2:
            gamma_k = frac(-1,4)
        elif k%4==0:
            gamma_k = frac(1,4)
        mu_k = 0
        if k%3==2:
            mu_k = frac(-1,3)
        elif k%3==0:
            mu_k = frac(1,3)
        return K(frac(-1,2) * mul([_lambda(r[p], s[p], p)     for p in L]) + \
                    gamma_k * mul([CO_delta(r[p], p, N, eps)  for p in L]) + \
                    mu_k    * mul([CO_nu(r[p], p, N, eps)     for p in L]))

    # We now consider the case of half-integral weight
    if den == 2:
        if N%4 != 0:
            raise TypeError("The modulus of the character must be divisible by 4.")
        if r[2] >= 4:
            zeta = _lambda(r[2], s[2], 2)
        if r[2] == 3:
            zeta = 3
        if r[2] == 2:
            C = False
            for p in L:
                if p%4 == 3:
                    if r[p]%2 == 1 or (0 < r[p] < 2*s[p]):
                        C = True
                        break
            if C:
                zeta = 2
            else:
                if 2*k % 4 == 1:
                    if s[2] == 0:
                        zeta = frac(3,2)
                    else:
                        zeta = frac(5,2)
                else:
                    if s[2] == 0:
                        zeta = frac(5,2)
                    else:
                        zeta = frac(3,2)
        return K(frac(-zeta,2) * mul([_lambda(r[p], s[p], p) for p in L if p!= 2]))

def SerreStark(eps, cusp_space=False):
    r"""
    Compute the size of the basis of theta series for the space
    of modular forms of weight 1/2 and character \varepsilon
    given by Serre-Stark.

    INPUT:


    -  ``eps`` - an even Dirichlet character, with conductor
        divisible by 4.

    -  ``cusp_space`` - (optional: default False) if True
        returns the size of a basis for the space of cusp forms.

    OUTPUT: Integer


    EXAMPLES::

        sage: G = DirichletGroup(108)
        sage: sage.modular.dims.SerreStark(G.0*G.1^9)
        2
        sage: H = DirichletGroup(64)
        sage: sage.modular.dims.SerreStark(H.0^2)
        3
        sage: sage.modular.dims.SerreStark(H.1^2)
        0
        sage: sage.modular.dims.SerreStark(H.1^8)
        2
        sage: sage.modular.dims.SerreStark(H.1^8, cusp_space=True)
        0
        sage: H = DirichletGroup(576)
        sage: chi = H.0*H.2^3
        sage: sage.modular.dims.SerreStark(chi, cusp_space=True)
        1
    """
    if eps.is_odd():
        raise TypeError("The character must be even")
    N = eps.modulus()
    if N%4 != 0:
        raise TypeError("The modulus of the character must be divisible by 4")
    d = 0
    for t in divisors(N/4):
        for rr in divisors(N/(4*t)):
            if rr.is_square():
                r = sqrt(rr)
                chars = DirichletGroup(r).list()
                for psi in chars:
                    if cusp_space:
                        if psi.is_totally_even():
                            continue
                    if psi.is_primitive() and psi.is_even():
                        a = 1
                        for n in range(1, N):
                            if gcd(n, N) == 1 and psi(n)*jacobi_symbol(t, n) != eps(n):
                                a = 0
                                break
                        d += a
    return ZZ(d)

####################################################################
# Functions exported to the global namespace.
# These have very flexible inputs.
####################################################################

def dimension_new_cusp_forms(X, k=2, p=0):
    """
    Return the dimension of the new (or `p`-new) subspace of
    cusp forms for the character or group `X`.

    INPUT:


    -  ``X`` - integer, congruence subgroup or Dirichlet
       character

    -  ``k`` - weight (integer)

    -  ``p`` - 0 or a prime


    EXAMPLES::

        sage: dimension_new_cusp_forms(100,2)
        1

    ::

        sage: dimension_new_cusp_forms(Gamma0(100),2)
        1
        sage: dimension_new_cusp_forms(Gamma0(100),4)
        5

    ::

        sage: dimension_new_cusp_forms(Gamma1(100),2)
        141
        sage: dimension_new_cusp_forms(Gamma1(100),4)
        463

    ::

        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,2)
        2
        sage: dimension_new_cusp_forms(DirichletGroup(100).1^2,4)
        8

    ::

        sage: sum(dimension_new_cusp_forms(e,3) for e in DirichletGroup(30))
        12
        sage: dimension_new_cusp_forms(Gamma1(30),3)
        12

    Check that Trac #12640 is fixed::

        sage: dimension_new_cusp_forms(DirichletGroup(1)(1), 12)
        1
        sage: dimension_new_cusp_forms(DirichletGroup(2)(1), 24)
        1
    """
    if is_GammaH(X):
        return X.dimension_new_cusp_forms(k,p=p)
    elif isinstance(X, dirichlet.DirichletCharacter):
        N = X.modulus()
        if N <= 2:
            return Gamma0(N).dimension_new_cusp_forms(k,p=p)
        else:
            # Gamma1(N) for N<=2 just returns Gamma0(N), which has no eps parameter. See Trac #12640.
            return Gamma1(N).dimension_new_cusp_forms(k,eps=X,p=p)
    elif isinstance(X, (int,long,Integer)):
        return Gamma0(X).dimension_new_cusp_forms(k,p=p)
    else:
        raise TypeError("X (=%s) must be an integer, a Dirichlet character or a congruence subgroup of type Gamma0, Gamma1 or GammaH" % X)

def dimension_cusp_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup or Dirichlet character.

    INPUT:


    -  ``X`` - congruence subgroup or Dirichlet character
       or integer

    -  ``k`` - weight (integer, or half an integer)


    EXAMPLES::

        sage: dimension_cusp_forms(5,4)
        1

    ::

        sage: dimension_cusp_forms(DirichletGroup(13).0^2,2)
        1
        sage: dimension_cusp_forms(DirichletGroup(13).0,3)
        1

    ::

        sage: dimension_cusp_forms(Gamma0(11),2)
        1
        sage: dimension_cusp_forms(Gamma0(11),0)
        0
        sage: dimension_cusp_forms(Gamma0(1),12)
        1
        sage: dimension_cusp_forms(Gamma0(1),2)
        0
        sage: dimension_cusp_forms(Gamma0(1),4)
        0

    ::

        sage: dimension_cusp_forms(Gamma0(389),2)
        32
        sage: dimension_cusp_forms(Gamma0(389),4)
        97
        sage: dimension_cusp_forms(Gamma0(2005),2)
        199
        sage: dimension_cusp_forms(Gamma0(11),1)
        0

    ::

        sage: dimension_cusp_forms(Gamma1(11),2)
        1
        sage: dimension_cusp_forms(Gamma1(1),12)
        1
        sage: dimension_cusp_forms(Gamma1(1),2)
        0
        sage: dimension_cusp_forms(Gamma1(1),4)
        0
        sage: dimension_cusp_forms(Gamma1(13),2)
        2

    ::

        sage: dimension_cusp_forms(Gamma1(389),2)
        6112
        sage: dimension_cusp_forms(Gamma1(389),4)
        18721
        sage: dimension_cusp_forms(Gamma1(2005),2)
        159201

    ::

        sage: dimension_cusp_forms(Gamma1(11),1)
        0

    ::

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_cusp_forms(e,2)
        0
        sage: dimension_cusp_forms(e^2,2)
        1

    Check that Trac #12640 is fixed::

        sage: dimension_cusp_forms(DirichletGroup(1)(1), 12)
        1
        sage: dimension_cusp_forms(DirichletGroup(2)(1), 24)
        5

    Examples with half-integral weights: ::

        sage: dimension_cusp_forms(44,3/2)
        2
        sage: dimension_cusp_forms(Gamma0(160),3/2)
        6
        sage: dimension_cusp_forms(Gamma1(44),3/2)
        12
        sage: G = DirichletGroup(108)
        sage: eps = G.0*G.1^9
        sage: dimension_cusp_forms(eps,3/2)
        5
        sage: dimension_cusp_forms(eps,1/2)
        0
        sage: H = DirichletGroup(576)
        sage: chi = H.0*H.2^3
        sage: dimension_cusp_forms(chi,1/2)
        1
        sage: K = DirichletGroup(156)
        sage: dimension_cusp_forms(K.1*K.2, 9/2)
        94
    """

    k = QQ(k)
    den = abs(k.denominator())

    if isinstance(X, dirichlet.DirichletCharacter):
        N = X.modulus()
        if N <= 2 and den == 1:
            return Gamma0(N).dimension_cusp_forms(k)
        else:
            return Gamma1(N).dimension_cusp_forms(k, X)
    elif is_ArithmeticSubgroup(X):
        return X.dimension_cusp_forms(k)
    elif isinstance(X, (Integer,int,long,Rational)):
        # For integral k we use the method given in Diamond--Shurman.
        # For half-integral k, we use the algorithms of Cohen--Oesterle
        # and Serre--Stark.
        if den == 1:
            return Gamma0(X).dimension_cusp_forms(k)
        else:
            return Gamma1(X).dimension_cusp_forms(k, trivial_character(X))
    else:
        raise TypeError("Argument 1 must be a Dirichlet character, an integer or a finite index subgroup of SL2Z")

def dimension_eis(X, k=2):
    """
    The dimension of the space of Eisenstein series for the given
    congruence subgroup.

    INPUT:


    -  ``X`` - congruence subgroup or Dirichlet character
       or integer

    -  ``k`` - weight (integer, or half an integer)


    EXAMPLES::

        sage: dimension_eis(5,4)
        2

    ::

        sage: dimension_eis(Gamma0(11),2)
        1
        sage: dimension_eis(Gamma1(13),2)
        11
        sage: dimension_eis(Gamma1(2006),2)
        3711

    ::

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2

    ::

        sage: e = DirichletGroup(13).0
        sage: e.order()
        12
        sage: dimension_eis(e,2)
        0
        sage: dimension_eis(e^2,2)
        2
        sage: dimension_eis(e,13)
        2

    ::

        sage: G = DirichletGroup(20)
        sage: dimension_eis(G.0,3)
        4
        sage: dimension_eis(G.1,3)
        6
        sage: dimension_eis(G.1^2,2)
        6

    ::

        sage: G = DirichletGroup(200)
        sage: e = prod(G.gens(), G(1))
        sage: e.conductor()
        200
        sage: dimension_eis(e,2)
        4

    ::

        sage: dimension_modular_forms(Gamma1(4), 11)
        6

    Examples with half-integral weights: ::

        sage: dimension_eis(44,3/2)
        3
        sage: dimension_eis(Gamma0(160),3/2)
        14
        sage: dimension_eis(Gamma1(44),3/2)
        38
        sage: G = DirichletGroup(108)
        sage: dimension_eis(G.0*G.1^9, 3/2)
        10
        sage: H = DirichletGroup(64)
        sage: dimension_eis(H.0^2, 1/2)
        3
        sage: dimension_eis(H.1^2, 1/2)
        0
        sage: dimension_eis(H.1^8, 1/2)
        2
        sage: K = DirichletGroup(156)
        sage: dimension_eis(K.1*K.2, 9/2)
        8
    """

    k = QQ(k)
    den = abs(k.denominator())

    if is_ArithmeticSubgroup(X):
        return X.dimension_eis(k)
    elif isinstance(X, dirichlet.DirichletCharacter):
        return Gamma1(X.modulus()).dimension_eis(k, X)
    elif isinstance(X, (int, long, Integer, Rational)):
        # For integral k we use the method given in Diamond--Shurman.
        # For half-integral k, we use the algorithms of Cohen--Oesterle
        # and Serre--Stark.
        if den == 1:
            return Gamma0(X).dimension_eis(k)
        else:
            return Gamma1(X).dimension_eis(k, trivial_character(X))
    else:
        raise TypeError("Argument in dimension_eis must be an integer, a Dirichlet character, or a finite index subgroup of SL2Z (got %s)" % X)

def dimension_modular_forms(X, k=2):
    r"""
    The dimension of the space of cusp forms for the given congruence
    subgroup (either `\Gamma_0(N)`, `\Gamma_1(N)`, or
    `\Gamma_H(N)`) or Dirichlet character.

    INPUT:


    -  ``X`` - congruence subgroup or Dirichlet character

    -  ``k`` - weight (integer, or half an integer)


    EXAMPLES::

        sage: dimension_modular_forms(Gamma0(11),2)
        2
        sage: dimension_modular_forms(Gamma0(11),0)
        1
        sage: dimension_modular_forms(Gamma1(13),2)
        13
        sage: dimension_modular_forms(GammaH(11, [10]), 2)
        10
        sage: dimension_modular_forms(GammaH(11, [10]))
        10
        sage: dimension_modular_forms(GammaH(11, [10]), 4)
        20
        sage: e = DirichletGroup(20).1
        sage: dimension_modular_forms(e,3)
        9
        sage: dimension_cusp_forms(e,3)
        3
        sage: dimension_eis(e,3)
        6
        sage: dimension_modular_forms(11,2)
        2

    Examples with half-integral weights: ::

        sage: dimension_modular_forms(44,3/2)
        5
        sage: dimension_modular_forms(Gamma0(160),3/2)
        20
        sage: dimension_modular_forms(Gamma1(44),3/2)
        50
        sage: G = DirichletGroup(156)
        sage: dimension_modular_forms(G.1*G.2, 9/2)
        102
        sage: H = DirichletGroup(64)
        sage: dimension_modular_forms(H.0^2, 1/2)
        3
        sage: dimension_modular_forms(H.1^2, 1/2)
        0
        sage: dimension_modular_forms(H.1^8, 1/2)
        2
        sage: K = DirichletGroup(108)
        sage: dimension_modular_forms(K.0*K.1^9, 3/2)
        15
    """

    if isinstance(X, (int, long, Integer)):
        if isinstance(k, Integer):
            return Gamma0(X).dimension_modular_forms(k)
        elif isinstance(k, Rational):
            return Gamma1(X).dimension_modular_forms(k, trivial_character(X))
        else: TypeError("Argument 2 must be an integer or half an integer.")
    elif is_ArithmeticSubgroup(X):
        return X.dimension_modular_forms(k)
    elif isinstance(X,dirichlet.DirichletCharacter):
        return Gamma1(X.modulus()).dimension_modular_forms(k, eps=X)
    else:
        raise TypeError("Argument 1 must be an integer, a Dirichlet character or an arithmetic subgroup.")

def sturm_bound(level, weight=2):
    r"""
    Returns the Sturm bound for modular forms with given level and weight. For
    more details, see the documentation for the sturm_bound method of
    sage.modular.arithgroup.CongruenceSubgroup objects.

    INPUT:


    -  ``level`` - an integer (interpreted as a level for Gamma0) or a congruence subgroup

    -  ``weight`` - an integer `\geq 2` (default: 2)

    EXAMPLES::

        sage: sturm_bound(11,2)
        2
        sage: sturm_bound(389,2)
        65
        sage: sturm_bound(1,12)
        1
        sage: sturm_bound(100,2)
        30
        sage: sturm_bound(1,36)
        3
        sage: sturm_bound(11)
        2
    """
    if is_ArithmeticSubgroup(level):
        if level.is_congruence():
            return level.sturm_bound(weight)
        else:
            raise ValueError("No Sturm bound defined for noncongruence subgroups")
    if isinstance(level, (int, long, Integer)):
        return Gamma0(level).sturm_bound(weight)
