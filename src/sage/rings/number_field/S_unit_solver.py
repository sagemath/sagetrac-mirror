r"""
code to solve S-unit equation x + y = 1


REFERENCES::

[TCDF] = N.P. Smart, The Solution of Triangularly Connected Decomposable Form Equations
[Smart] = N.P. Smart, The algorithmic resolution of Diophantine equations

AUTHORS:

- Alejandra Alvarado, Angelos Koutsianas, Beth Malmskog, Christopher Rasmussen, Christelle Vincent, Mckenzie West (2017-01-10): original version

EXAMPLES::


"""


#*****************************************************************************
#       Copyright (C) 2017 Alejandra Alvarado <aalvarado2 at eiu.edu>
#                          Angelos Koutsianas <koutsis.jr at gmail.com>
#                          Beth Malmskog <beth.malmskog at gmail.com>
#                          Christopher Rasmussen <crasmussen at wesleyan.edu>
#                          Christelle Vincent <christelle.vincent at uvm.edu>
#                          Mckenzie West <mckenzierwest at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

import numpy
import math

from sage.rings.ring import Field
#from sage.rings.real_mpfr import RealField

from sage.rings.all import Infinity
from sage.rings.number_field.number_field import NumberField
from sage.rings.number_field.unit_group import UnitGroup
from sage.rings.number_field.number_field_ideal import NumberFieldIdeal
from sage.rings.number_field.number_field_element import NumberFieldElement
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.groups.abelian_gps.values import AbelianGroupWithValues_class
from sage.combinat.multichoose_nk import MultichooseNK
from sage.combinat.combination import Combinations
from sage.arith.all import factorial
from sage.misc.all import prod
from sage.rings.real_mpfr import RealField
from math import pi


def set_R(prec = None):
    r"""
    Return a real field of precision ``prec``

    INPUT:

    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:

    - A real field (either RIF, RDF or RealField(prec))

    EXAMPLES::

        sage: set_R()
        Real Field with 106 bits of precision

        sage: set_R(53)
        Real Field with 53 bits of precision
    """
    if prec == None:
        return RealField(106)
    else:
        return RealField(prec)

def is_real(v, prec = None):
    r"""
    Return ``True`` if `v` is real, ``False`` if `v` is complex

    INPUT:

    - ``v`` -- an infinite place of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:

    - ````True`` or ``False``

    EXAMPLES::

        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: is_real(phi_real)
        True

        sage: is_real(phi_complex)
        False

    It is an error to put in a finite place

    ::

        sage: is_real(v_fin)
        Traceback (most recent call last):
        ...
        AttributeError: 'NumberFieldFractionalIdeal' object has no attribute 'im_gens'

    """
    R = set_R(prec)

    try:
        R(v.im_gens()[0])
        return True
    except TypeError:
        return False

def abs_val(SUK, v, iota, prec = None):
    r"""
    Return the value `|iota|_{v}`.

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of SUK.number_field().places(prec = prec))
    - ``iota`` -- an element of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The absolute value as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: abs_val(SUK,phi_real,xi^2)
        2.080083823051904114530056824358

        sage: abs_val(SUK,phi_complex,xi^2)
        4.326748710922225349406744498992

        sage: abs_val(SUK,v_fin,xi^2)
        0.1111111111111111111111111111111

    It is an error to use a complex embedding rather than a place
    ::
        sage: phi_real = K.complex_embeddings()[2]
        sage: abs_val(SUK,phi_real,xi^2)
        Traceback (most recent call last):
        ...
        AttributeError: 'NumberFieldHomomorphism_im_gens' object has no attribute 'smallest_integer'
    """
    R = set_R(prec)
    K = SUK.number_field()
    primes = SUK.primes()
    num_primes = len(primes)
    rank = SUK.rank()
    if v in K.places(prec = prec):
        if is_real(v):
            return R(v(iota).abs())
        else:
            return R(v(iota).abs()**2)
    else:
        # v is a finite place
        p = v.smallest_integer()
        iota_ideal = K.ideal(K(iota))
        exponent = - v.residue_class_degree() * iota_ideal.valuation(v)
        return R(p**exponent)

def column_Log(SUK, iota, U, prec = None):
    r"""
    Return the "log vector" of iota; i.e., the logs of all the valuations

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``iota`` -- an element of ``K``
    - ``U`` -- a list of places (finite or infinite) of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The log vector as a list of real numbers

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: U = [phi_complex, v_fin]
        sage: column_Log(SUK, xi^2, U)
        [1.464816384890812968648768625966, -2.197224577336219382790490473845]

    REFERENCES:
    .. [TCDF] p. 823
    """
    R = set_R(prec)
    from sage.functions.log import log
    return [ R(log(abs_val(SUK,v,iota,prec))) for v in U]

def c1_func(SUK, prec = None):
    r"""
    Return the constant ``c1`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant c1, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))

        sage: c1_func(SUK)
        1.174298947359884426898831030919

    REFERENCES:
    .. [TCDF] p. 823

    """
    R = set_R(prec)
    all_places = list(SUK.primes()) + SUK.number_field().places(prec = prec)
    Possible_U = Combinations(all_places, SUK.rank())
    c1 = 0
    for U in Possible_U:
        # first, build the matrix C_{i,U}
        columns_of_C = []
        for unit in SUK.fundamental_units():
            columns_of_C.append( column_Log(SUK, unit, U, prec) )
        from sage.matrix.constructor import Matrix
        C = Matrix(SUK.rank(), SUK.rank(), columns_of_C)
        # Is it invertible?
        if abs(C.determinant()) > 10**(-10):
            A = C.inverse().apply_map(abs)
            poss_c1 = max([sum(i) for i in list(A)])
            c1 = max(poss_c1,c1)
    return R(c1)

def c3_func(SUK, prec = None):
    r"""
    Return the constant ``c3`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c3``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))

        sage: c3_func(SUK)
        0.4257859134798034905777797121118

    ..NOTE::
        The numerator should be as close to 1 as possible, especially as the rank of the `S`-units grows large

    REFERENCES:
    .. [TCDF] p. 823

    """

    R = set_R(prec)
    return R(0.9999999/(c1_func(SUK)*SUK.rank()))

def c4_func(SUK,v, A, prec = None):
    r"""
    Return the constant ``c4`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of SUK.number_field().places(prec = prec))
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c4``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c4_func(SUK,phi_real,A)
        1.000000000000000000000000000000

        sage: c4_func(SUK,phi_complex,A)
        1.000000000000000000000000000000

        sage: c4_func(SUK,v_fin,A)
        1.000000000000000000000000000000

    REFERENCES:
    .. [TCDF] p. 824

    """
    R = set_R(prec)
    return max( [abs_val(SUK, v, alpha, prec) for alpha in A])

def c5_func(SUK, v, prec = None):
    r"""
    Return the constant ``c5`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant c5, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: c5_func(SUK,v_fin)
        0.1291890135314859427129849620671

    REFERENCES:
    .. [TCDF] p. 824

    """
    R = set_R(prec)
    from sage.functions.log import log
    return R(c3_func(SUK, prec)/(v.residue_class_degree()*log(v.smallest_integer())*v.ramification_index()))

def c6_func(SUK, v, A, prec = None):
    r"""
    Return the constant ``c6`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c6``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c6_func(SUK,v_fin,A)
        0.0000000000000000000000000000000

    REFERENCES:
    .. [TCDF] p. 824

    """
    R = set_R(prec)
    from sage.functions.log import log
    return R(log(c4_func(SUK, v, A, prec))/(v.residue_class_degree()*log(v.smallest_integer())*v.ramification_index()))

def c7_func(SUK, v, A, prec = None):
    r"""
    Return the constant ``c7`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c7``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c7_func(SUK,v_fin,A)
        0.0000000000000000000000000000000

    REFERENCES:
    .. [TCDF] p. 824

    """
    R = set_R(prec)
    from sage.functions.log import log
    return R(log(c4_func(SUK, v, A, prec))/c3_func(SUK, prec))

def beta_ns(SUK,v):
    r"""
    Return a list of pairs `[beta,val_v(beta)]`, for `beta` in the fundamental units of SUK

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``


    OUTPUT:
    The list of pairs `[beta,val_v(beta)]`, for `beta` in the fundamental units of SUK (elements of K), val_v(beta) is an integer

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: beta_ns(SUK,v_fin)
        [[xi^2 - 2, 0], [xi, 1]]

    ..NOTE::
        n_0 is the valuation of the coefficient alpha_d of the `S`-unit equation such that val_v(alpha_d tau_d) = 0
        We have set n_0 = 0 here since the coefficients are roots of unity

    REFERENCES:
    .. [TCDF] pp. 824-825

    """
    betas = SUK.fundamental_units()
    # n_0 = valuation_v of one of the coefficients of the equation = 0 for x + y = 1 p. 824
    return [[beta,beta.valuation(v)] for beta in betas]

def beta_k(betas_and_ns):
    r"""
    Return a pair `[beta_k,val_v(beta_k)]`, where `beta_k` has the smallest nonzero valuation in absolute value of the list `betas_and_ns`

    INPUT:
    - ``betas_and_ns`` -- a list of pairs ``[beta,val_v(beta)]`` outputted from the function beta_ns

    OUTPUT:
    The pair ``[beta_k,v(beta_k)],`` ``beta_k`` is an element of ``K`` and ``val_v(beta_k)`` is a integer

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: betas = beta_ns(SUK,v_fin)
        sage: beta_k(betas)
        [xi, 1]

    REFERENCES:
    .. [TCDF] pp. 824-825

    """
    for pair in betas_and_ns:
        if pair[1].abs() != 0:
            good_pair = pair
            break
    for pair in betas_and_ns:
        if (pair[1].abs() != 0 and pair[1].abs() < good_pair[1].abs()):
            good_pair = pair
    return good_pair

def mus(SUK,v):
    r"""
    Return a list `[mus]`, for `mus` defined on pp. 824-825 of TCDF

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``


    OUTPUT:
    A list ``[mus]`` where each mu is an element of ``K``

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: mus(SUK,v_fin)
        [xi^2 - 2, 1]

    REFERENCES:
    .. [TCDF] pp. 824-825

    """
    beta_and_ns = beta_ns(SUK,v)
    if all(pair[1]==0 for pair in beta_and_ns):
        return SUK.fundamental_units()
    else:
        good_pair = beta_k(beta_and_ns)
        return [(beta[0]**good_pair[1])*(good_pair[0]**(-beta[1])) for beta in beta_and_ns]

def possible_mu0s(SUK,v):
    r"""
    Return a list `[mu0s]` of all possible `mu0s` defined on pp. 824-825 of TCDF

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``


    OUTPUT:
    A list ``[mu0s]`` where each ``mu0`` is an element of ``K``

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: possible_mu0s(SUK,v_fin)
        [1]

    ..NOTE::
        `n_0` is the valuation of the coefficient alpha_d of the `S`-unit equation such that `val_v(alpha_d tau_d) = 0`
        We have set `n_0 = 0` here since the coefficients are roots of unity
        `alpha_0` is not defined in the paper, we set it to be 1

    REFERENCES:
    .. [TCDF] pp. 824-825, but we modify the definition of ``sigma`` (``sigma_tilde``) to make it easier to code

    """
    beta_and_ns = beta_ns(SUK,v)
    ns = [beta[1] for beta in beta_and_ns]
    betas = [beta[0] for beta in beta_and_ns]
    betak, nk = beta_k(beta_and_ns)
    mu0s = []
    for rs in MultichooseNK(nk.abs(),SUK.rank()):
        # n_0 = valuation_v of one of the coefficients of the equation = 0 for x + y = 1 p. 824
        n_rs = zip(ns,rs)
        sigma_tilde = -(sum([n_r[0]*n_r[1] for n_r in n_rs]))
        if sigma_tilde % nk == 0:
            # alpha0 = 1, but we don't know what alpha0 is
            beta_rs = zip(betas,rs)
            mu0s.append(prod([beta_r[0]**beta_r[1] for beta_r in beta_rs])*betak**(sigma_tilde/nk))
    return mu0s

def modified_height(SUK,v,D,b, prec = None):
    r"""
    Return the modified height at the finite place `v`

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an finite place of ``K``
    - ``D`` -- an auxiliary quantity (an integer)
    - ``b`` -- an element of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The modified height at the place `v`, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: modified_height(SUK,v_fin,2*K.degree(),xi^2)
        0.7324081924454063363683076204325

    REFERENCES:
    .. [Smart] p. 226

    """
    R = set_R(prec)
    d = SUK.number_field().degree()
    f_p = v.residue_class_degree()
    p = v.smallest_integer()
    from sage.functions.log import log
    max_log_b = max([log(phi(b)).abs() for phi in SUK.number_field().places(prec = prec)])
    return R(max([b.global_height(),max_log_b/(2*pi*D),f_p*log(p)/d]))

def local_c3(SUK,v,D, prec = None):
    r"""
    Return a factor of the constant ``c3`` defined on p. 226 of [Smart]

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an finite place of ``K``
    - ``D`` -- an auxiliary quantity (an integer)
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    A factor of the constant ``c3,`` as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: local_c3(SUK,v_fin,2*K.degree())
        0.1128578369553013940126003163566

    REFERENCES:
    .. [Smart] p. 226

    """
    R = set_R(prec)
    mus_prod = prod([modified_height(SUK,v,D,b,prec) for b in mus(SUK,v)])
    return R(max([mus_prod*modified_height(SUK,v,D,mu0,prec) for mu0 in possible_mu0s(SUK,v)]))

def c8_c9_func(SUK, v, A, prec = None):
    r"""
    Return the constants `c8` and `c9` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K`` (a fractional ideal)
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constants ``c8`` and ``c9``, as real numbers

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c8_c9_func(SUK,v_fin,A)
        (3.773087980338467375636381140719e19, 1.352093223400401265533140708466e20)

    REFERENCES:
    .. [TCDF] p. 825
    .. [Smart] p. 226, Theorem A.2 for the local constants

    """
    R = set_R(prec)
    num_mus = len(mus(SUK,v))+1
    p = v.smallest_integer()
    f_p = v.residue_class_degree()
    d = SUK.number_field().degree()
    if p == 2:
        local_c2 = 197142*36**num_mus
    elif p%4 == 1:
        local_c2 = 35009*(45/2)**num_mus
    else:
        local_c2 = 30760*25**num_mus
    x = polygen(SUK.number_field())
    if ( p > 2 and not ((x**2+1).is_irreducible()) ) or ( p==2 and not ((x**2+3).is_irreducible()) ):
        D = d
    else:
        D = 2*d
    from sage.functions.log import log
    l_c3 = (num_mus+1)**(2*num_mus+4)*p**(D * f_p/d)*(f_p*log(p))**(-num_mus-1)*D**(num_mus+2)
    l_c3 *= local_c3(SUK,v,D, prec)
    H = max([modified_height(SUK,v,D,alpha, prec) for alpha in mus(SUK,v)+possible_mu0s(SUK,v)])
    if p == 2:
        local_c4 = log(3*2**10*(num_mus+1)**2*D**2*H)
    else:
        local_c4 = log(2**11*(num_mus+1)**2*D**2*H)
    local_c5 = 2*log(D)
    return R(local_c2*l_c3*local_c4), R(local_c2*l_c3*local_c4*local_c5)

def c10_func(SUK, v, A, prec = None):
    r"""
    Return the constant `c10` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K`` (a fractional ideal)
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
   - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c10,`` as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c10_func(SUK,v_fin,A)
        9.659069151659341551509042989895e21

    REFERENCES:
    .. [TCDF] p. 824

    """
    R = set_R(prec)
    e_h = v.ramification_index()
    c_8, c_9 = c8_c9_func(SUK, v, A, prec)
    from sage.functions.log import log
    return R((2/(e_h*c5_func(SUK,v,prec)))*(e_h*c6_func(SUK, v, A, prec) + c_9 + c_8 * log( c_8/(e_h*c5_func(SUK,v, prec)))))

def c11_func(SUK, v, A, prec = None):
    r"""
    Return the constant `c11` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of SUK.number_field().places(prec = prec))
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c11``, a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c11_func(SUK,phi_real,A)
        3.255848343572896031429547409663

        sage: c11_func(SUK,phi_complex,A)
        6.511696687145792062859094819326

    REFERENCES:
    .. [TCDF] p. 825

    """
    R = set_R(prec)
    assert v in SUK.number_field().places(prec = prec)
    from sage.functions.log import log
    if is_real(v,prec):
        return R(log(4*c4_func(SUK, v, A, prec))/(c3_func(SUK, prec)))
    else:
        from sage.functions.other import sqrt
        return R(2*(log(4*sqrt(c4_func(SUK,v, A, prec))))/(c3_func(SUK, prec)))

def c12_func(SUK, v, A, prec = None):
    r"""
    Return the constant ``c12`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places(prec = prec))
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c12``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c12_func(SUK,phi_real,A)
        2.000000000000000000000000000000
        
        sage: c12_func(SUK,phi_complex,A)
        2.000000000000000000000000000000

    REFERENCES:
    .. [TCDF] p. 825

    """
    R = set_R(prec)
    assert v in SUK.number_field().places(prec = prec)
    if is_real(v,prec):
        return R(2*c4_func(SUK, v, A, prec))
    else:
        from sage.functions.other import sqrt
        return R(2*sqrt(c4_func(SUK,v, A, prec)))

def c13_func(SUK, v, prec = None):
    r"""
    Return the constant ``c13`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places(prec = prec))
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c13``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: c13_func(SUK,phi_real)
        0.4257859134798034905777797121118

        sage: c13_func(SUK,phi_complex)
        0.2128929567399017452888898560559

    REFERENCES:
    .. [TCDF] p. 825

    """
    R = set_R(prec)
    assert v in SUK.number_field().places(prec = prec)
    if is_real(v,prec):
        return c3_func(SUK,prec)
    else:
        return c3_func(SUK,prec)/2

def Baker_C(t,d,prec = None):
    r"""
    Return ``C(t,d)`` from Smart's TCDF paper

    INPUT:
    - `t` -- the rank of the `S`-unit group (an integer)
    - `d` -- the degree of `K` (an integer)
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``C(t,d)`` as a real number

    EXAMPLES::
        sage: Baker_C(3,4)
        3.371389974809293972795830360284e19

    REFERENCES:
    .. [Smart] p. 225, Theorem A.1

    """
    R = set_R(prec)
    from sage.functions.log import log
    return R( 18 * factorial(t+2) * (t+1)**(t+2) * (32*d)**(t + 3) * log( 2*(t+1) * d) )

def hprime(SUK, alpha, v, prec = None):
    r"""
    Return the modified height at the infinite place `v`

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``alpha`` -- an element of ``K``
    - ``v`` -- an infinite place of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The modified height at the place `v`, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: hprime(SUK,xi^2,phi_real)
        0.7324081924454063363683076204325

        sage: hprime(SUK,xi^2,phi_complex)
        0.7395879186929970039443560381187

    REFERENCES:
    .. [Smart] p. 225

    """
    R = set_R(prec)
    assert v in SUK.number_field().places(prec = prec)
    assert alpha in SUK.number_field()
    from sage.functions.log import log
    return R(max(alpha.global_height(), 1/SUK.number_field().degree(), log(v(alpha)).abs()/SUK.number_field().degree()))

def c14_func(SUK,v,A,prec = None):
    r"""
    Return the constant ``c14`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places(prec = prec))
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c14``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c14_func(SUK,phi_real,A)
        2.661434475699529175064652896133e14

        sage: c14_func(SUK,phi_complex,A)
        5.908765797306433904563179346367e14

    REFERENCES:
    .. [TCDF] p. 825
    .. [Smart] p. 225, Theorem A.1

    """
    R = set_R(prec)
    assert v in SUK.number_field().places(prec = prec)
    c_1 = Baker_C(SUK.rank(),SUK.number_field().degree(),prec)
    hproduct = c_1 * prod([hprime(SUK, alpha, v, prec) for alpha in SUK.gens_values()])
    return hproduct

def c15_func(SUK, v, A, prec = None):
    r"""
    Return the constant `c15` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places(prec = prec))
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``c15``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c15_func(SUK,phi_real,A)
        4.396386097852707225469494870675e16

        sage: c15_func(SUK,phi_complex,A)
        2.034870098399844351862009972531e17

    REFERENCES:
    .. [TCDF] p. 825

    """
    R = set_R(prec)
    from sage.functions.log import log
    return R(2*(log(c12_func(SUK,v,A,prec))+c14_func(SUK,v,A,prec)*log((SUK.rank()+1)*c14_func(SUK,v,A,prec)/c13_func(SUK,v,prec)))/c13_func(SUK,v,prec))

def K0_func(SUK, A, prec = None):
    r"""
    Return the constant `K0` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``K0``, a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: K0_func(SUK,A)
        9.659069151659341551509042989895e21

    REFERENCES:
    .. [TCDF] p. 824

    """
    R = set_R(prec)
    return R(max([c10_func(SUK,v, A, prec) for v in SUK.primes()] + [c7_func(SUK,v,A, prec) for v in SUK.primes()]))

def K1_func(SUK, v, A, prec = None):
    r"""
    Return the constant ``K1`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places(prec = prec))
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``
    - ``prec`` (default = None) -- the precision of the real field

    OUTPUT:
    The constant ``K1,`` a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: K1_func(SUK,phi_real,A)
        4.396386097852707225469494870675e16

        sage: K1_func(SUK,phi_complex,A)
        2.034870098399844351862009972531e17

    REFERENCES:
    .. [TCDF] p. 825

    """
    R = set_R(prec)
    return max([c11_func(SUK,v, A, prec), c15_func(SUK,v,A,prec)])
