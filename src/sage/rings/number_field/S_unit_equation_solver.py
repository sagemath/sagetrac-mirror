r"""
code to solve S-unit equation x + y = 1


REFERENCES::

[TCDF] =
[Smart] = N.P. Smart, The algorithmic resolution of Diophantine equations

AUTHORS:

- Alejandra Alvarado, Angelos Koutsianas, Beth Malmskog, Chris Rasmussen, Christelle Vincent, Mckenzie West (2017-01-10): original version

EXAMPLES::


"""


#*****************************************************************************
#       Copyright (C) 2017 Alejandra Alvarado <aalvarado2 at eiu.edu>
#                          Angelos Koutsianas <koutsis.jr at gmail.com>
#                          Beth Malmskog <beth.malmskog at gmail.com>
#                          Chris Rasmussen <crasmussen at wesleyan.edu>
#                          Christelle Vincent <christelle.vincent at uvm.edu>
#                          Mckenzie West <mckenzierwest at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def is_real(v):
    r"""
    Return ``True`` if `v` is real, ``False`` if `v` is complex

    INPUT:
    - ``v`` -- an infinite place of ``K``

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
    try:
        RR(v.im_gens()[0])
        return True
    except TypeError:
        return False

def abs_val(SUK, v, iota):
    r"""
    Return the value `|iota|_{v}`.

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of SUK.number_field().places())
    - ``iota`` -- an element of ``K``

    OUTPUT:
    The absolute value as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: abs_val(SUK,phi_real,xi^2)
        2.08008382305190

        sage: abs_val(SUK,phi_complex,xi^2)
        4.32674871092223

        sage: abs_val(SUK,v_fin,xi^2)
        0.111111111111111

    It is an error to use a complex embedding rather than a place
    ::
        sage: phi_real = K.complex_embeddings()[2]
        sage: abs_val(SUK,phi_real,xi^2)
        Traceback (most recent call last):
        ...
        AttributeError: 'NumberFieldHomomorphism_im_gens' object has no attribute 'smallest_integer'
    """
    K = SUK.number_field()
    primes = SUK.primes()
    num_primes = len(primes)
    rank = SUK.rank()
    if v in K.places():
        if is_real(v):
            return RealField()(v(iota).abs())
        else:
            return RealField()(v(iota).abs()**2)
    else:
        # v is a finite place
        p = v.smallest_integer()
        iota_ideal = ideal(K(iota))
        exponent = - v.residue_class_degree() * iota_ideal.valuation(v)
        return RealField()(p**exponent)

def column_Log(SUK, iota, U):
    r"""
    Return the "log vector" of iota; i.e., the logs of all the valuations

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``iota`` -- an element of ``K``
    - ``U`` -- a list of places (finite or infinite) of ``K``

    OUTPUT:
    The log vector as a list of real numbers

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: U = [phi_complex, v_fin]
        sage: column_Log(SUK, xi^2, U)
        [1.46481638489081, -2.19722457733622]

    REFERENCES:
    .. [TCDF] p. 823
    """
    return [ RealField()(log(abs_val(SUK,v,iota))) for v in U]

def c1(SUK):
    r"""
    Return the constant ``c1`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units

    OUTPUT:
    The constant c1, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))

        sage: c1(SUK)
        1.1742989473598846

    REFERENCES:
    .. [TCDF] p. 823

    """
    all_places = list(SUK.primes()) + SUK.number_field().places()
    Possible_U = Combinations(all_places, SUK.rank())
    c1 = 0
    for U in Possible_U:
        # first, build the matrix C_{i,U}
        columns_of_C = []
        for unit in SUK.fundamental_units():
            columns_of_C.append( column_Log(SUK, unit, U) )
        C = Matrix(SUK.rank(), SUK.rank(), columns_of_C)
        # Is it invertible?
        if abs(C.determinant()) > 10**(-10):
            poss_c1 = C.inverse().norm(Infinity)
            c1 = max(poss_c1,c1)
    return c1

def c3(SUK):
    r"""
    Return the constant ``c3`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units

    OUTPUT:
    The constant ``c3``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))

        sage: c3(SUK)
        0.4257859134798034

    ..NOTE::
        The numerator should be as close to 1 as possible, especially as the rank of the `S`-units grows large

    REFERENCES:
    .. [TCDF] p. 823

    """
    return 0.9999999/(c1(SUK)*SUK.rank())

def c4(SUK,v, A):
    r"""
    Return the constant ``c4`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of SUK.number_field().places())
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c4``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c4(SUK,phi_real,A)
        1.00000000000000

        sage: c4(SUK,phi_real,A)
        1.00000000000000

        sage: c4(SUK,v_fin,A)
        1.00000000000000

    REFERENCES:
    .. [TCDF] p. 824

    """
    return max( [abs_val(SUK, v, alpha) for alpha in A])

def c5(SUK, v):
    r"""
    Return the constant ``c5`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``

    OUTPUT:
    The constant c5, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: c5(SUK,v_fin)
        0.129189013531486

    REFERENCES:
    .. [TCDF] p. 824

    """
    return RR(c3(SUK)/(v.residue_class_degree()*log(v.smallest_integer())*v.ramification_index()))

def c6(SUK, v, A):
    r"""
    Return the constant ``c6`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c6``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c6(SUK,v_fin,A)
        0.000000000000000

    REFERENCES:
    .. [TCDF] p. 824

    """
    return RR(log(c4(SUK, v, A))/(v.residue_class_degree()*log(v.smallest_integer())*v.ramification_index()))

def c7(SUK, v, A):
    r"""
    Return the constant ``c7`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K``
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c7``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c7(SUK,v_fin,A)
        0.000000000000000

    REFERENCES:
    .. [TCDF] p. 824

    """
    return RR(log(c4(SUK, v, A))/c3(SUK))

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

def modified_height(SUK,v,D,b):
    r"""
    Return the modified height at the finite place `v`

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an finite place of ``K``
    - ``D`` -- an auxiliary quantity (an integer)
    - ``b`` -- an element of ``K``

    OUTPUT:
    The modified height at the place `v`, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: modified_height(SUK,v_fin,2*K.degree(),xi^2)
        0.732408192445406

    REFERENCES:
    .. [Smart] p. 226

    """
    d = SUK.number_field().degree()
    f_p = v.residue_class_degree()
    p = v.smallest_integer()
    max_log_b = max([log(phi(b)).abs() for phi in SUK.number_field().places()])
    return RR(max([b.global_height(),max_log_b/(2*pi*D),f_p*log(p)/d]))

def local_c3(SUK,v,D):
    r"""
    Return a factor of the constant ``c3`` defined on p. 226 of [Smart]

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an finite place of ``K``
    - ``D`` -- an auxiliary quantity (an integer)

    OUTPUT:
    A factor of the constant c3, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]

        sage: local_c3(SUK,v_fin,2*K.degree())
        0.112857836955301

    REFERENCES:
    .. [Smart] p. 226

    """
    mus_prod = prod([modified_height(SUK,v,D,b) for b in mus(SUK,v)])
    return RR(max([mus_prod*modified_height(SUK,v,D,mu0) for mu0 in possible_mu0s(SUK,v)]))

def c8_c9(SUK, v, A):
    r"""
    Return the constants `c8` and `c9` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K`` (a fractional ideal)
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constants ``c8`` and ``c9``, as real numbers

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c8_c9(SUK,v_fin,A)
        (3.77308798033847e19, 1.35209322340040e20)

    REFERENCES:
    .. [TCDF] p. 825
    .. [Smart] p. 226, Theorem A.2 for the local constants

    """
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
    l_c3 = (num_mus+1)**(2*num_mus+4)*p**(D * f_p/d)*(f_p*log(p))**(-num_mus-1)*D**(num_mus+2)
    l_c3 *= local_c3(SUK,v,D)
    H = max([modified_height(SUK,v,D,alpha) for alpha in mus(SUK,v)+possible_mu0s(SUK,v)])
    if p == 2:
        local_c4 = log(3*2**10*(num_mus+1)**2*D**2*H)
    else:
        local_c4 = log(2**11*(num_mus+1)**2*D**2*H)
    local_c5 = 2*log(D)
    return RR(local_c2*l_c3*local_c4), RR(local_c2*l_c3*local_c4*local_c5)

def c10(SUK, v, A):
    r"""
    Return the constant `c10` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a finite place of ``K`` (a fractional ideal)
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c10,`` as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: v_fin = tuple(K.primes_above(3))[0]
        sage: A = K.roots_of_unity()

        sage: c10(SUK,v_fin,A)
        9.65906915165934e21

    REFERENCES:
    .. [TCDF] p. 824

    """
    e_h = v.ramification_index()
    c_8, c_9 = c8_c9(SUK, v, A)
    return RR((2/(e_h*c5(SUK,v)))*(e_h*c6(SUK, v, A) + c_9 + c_8 * log( c_8/(e_h*c5(SUK,v)))))

def c11(SUK, v, A):
    r"""
    Return the constant `c11` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- a place of ``K``, finite (a fractional ideal) or infinite (element of SUK.number_field().places())
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c11``, a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c11(SUK,phi_real,A)
        3.25584834357290

        sage: c11(SUK,phi_complex,A)
        6.51169668714579

    REFERENCES:
    .. [TCDF] p. 825

    """
    assert v in SUK.number_field().places()
    if is_real(v):
        return RR(log(4*c4(SUK, v, A))/(c3(SUK)))
    else:
        return RR(2*(log(4*sqrt(c4(SUK,v, A))))/(c3(SUK)))

def c12(SUK, v, A):
    r"""
    Return the constant ``c12`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places())
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c12``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c12(SUK,phi_real,A)
        2.00000000000000

        sage: c12(SUK,phi_complex,A)
        2.00000000000000

    REFERENCES:
    .. [TCDF] p. 825

    """
    assert v in SUK.number_field().places()
    if is_real(v):
        return RealField()(2*c4(SUK, v, A))
    else:
        return RealField()(2*sqrt(c4(SUK,v, A)))

def c13(SUK, v):
    r"""
    Return the constant ``c13`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places())

    OUTPUT:
    The constant ``c13``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: c13(SUK,phi_real)
        0.4257859134798034

        sage: c13(SUK,phi_complex)
        0.2128929567399017

    REFERENCES:
    .. [TCDF] p. 825

    """
    assert v in SUK.number_field().places()
    if is_real(v):
        return c3(SUK)
    else:
        return c3(SUK)/2

def Baker_C(t,d):
    r"""
    Return ``C(t,d)`` from Smart's TCDF paper

    INPUT:
    - `t` -- the rank of the `S`-unit group (an integer)
    - `d` -- the degree of `K` (an integer)

    OUTPUT:
    The constant ``C(t,d)`` as a real number

    EXAMPLES::
        sage: Baker_C(3,4)
        3.37138997480929e19

    REFERENCES:
    .. [Smart] p. 225, Theorem A.1

    """
    return RR( 18 * factorial(t+2) * (t+1)**(t+2) * (32*d)**(t + 3) * log( 2*(t+1) * d) )

def hprime(SUK, alpha, v):
    r"""
    Return the modified height at the infinite place `v`

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``alpha`` -- an element of ``K``
    - ``v`` -- an infinite place of ``K``

    OUTPUT:
    The modified height at the place `v`, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]

        sage: hprime(SUK,xi^2,phi_real)
        0.732408192445406

        sage: hprime(SUK,xi^2,phi_complex)
        0.739587918692997

    REFERENCES:
    .. [Smart] p. 225

    """
    assert v in SUK.number_field().places()
    assert alpha in SUK.number_field()
    return RR(max(alpha.global_height(), 1/K.degree(), log(v(alpha)).abs()/K.degree()))

def c14(SUK,v,A):
    r"""
    Return the constant ``c14`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places())
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c14``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c14(SUK,phi_real,A)
        2.66143447569953e14

        sage: c14(SUK,phi_complex,A)
        5.90876579730643e14

    REFERENCES:
    .. [TCDF] p. 825
    .. [Smart] p. 225, Theorem A.1

    """
    assert v in SUK.number_field().places()
    c_1 = Baker_C(SUK.rank(),SUK.number_field().degree())
    hproduct = c_1 * prod([hprime(SUK, alpha, v) for alpha in SUK.gens_values()])
    return hproduct

def c15(SUK, v, A):
    r"""
    Return the constant `c15` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``v`` -- an infinite place of ``K`` (element of SUK.number_field().places())
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``c15``, as a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: phi_real = K.places()[0]
        sage: phi_complex = K.places()[1]
        sage: A = K.roots_of_unity()

        sage: c15(SUK,phi_real,A)
        4.39638609785271e16

        sage: c15(SUK,phi_complex,A)
        2.03487009839984e17

    REFERENCES:
    .. [TCDF] p. 825

    """
    return RR(2*(log(c12(SUK,v,A))+c14(SUK,v,A)*log((SUK.rank()+1)*c14(SUK,v,A)/c13(SUK,v)))/c13(SUK,v))

def K0(SUK, A):
    r"""
    Return the constant `K0` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``K0``, a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: K0(SUK,A)
        9.65906915165934e21

    REFERENCES:
    .. [TCDF] p. 824

    """
    return RR(max([c10(SUK,v, A) for v in SUK.primes()] + [c7(SUK,v,A) for v in SUK.primes()]))

def K1(SUK, A):
    r"""
    Return the constant ``K1`` from Smart's TCDF paper

    INPUT:
    - ``SUK`` -- a group of ``S``-units
    - ``A`` -- the set of the product of the coefficients of the ``S``-unit equation with each root of unity of ``K``

    OUTPUT:
    The constant ``K1,`` a real number

    EXAMPLES::
        sage: K.<xi> = NumberField(x^3-3)
        sage: SUK = UnitGroup(K,S=tuple(K.primes_above(3)))
        sage: A = K.roots_of_unity()

        sage: K1(SUK,A)
        2.03487009839984e17

    REFERENCES:
    .. [TCDF] p. 825

    """
    return max([c11(SUK,v, A) for v in SUK.number_field().places()]+[c15(SUK,v,A) for v in SUK.number_field().places()])
