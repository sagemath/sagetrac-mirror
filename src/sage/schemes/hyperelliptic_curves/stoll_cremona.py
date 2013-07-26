r"""
Reduction of binary forms and hyperelliptic curve equations

This file contains functions that reduce the coefficients of a binary form or
hyperelliptic curve by integral transformations, leaving the discriminant
intact.

AUTHORS:

- Florian Bouyer (2011 - 09) : covariant and stoll_cremona_reduction
- Marco Streng (2012 - 09) : Hilbert fundamental domains in some cases

"""
#*****************************************************************************
#       Copyright (C) 2011,2012,2013 Florian Bouyer <f.j.s.c.bouyer@gmail.com>
#                                and Marco Streng <marco.streng@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.matrix.all import Matrix
from sage.schemes.plane_conics.constructor import Conic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.arith import gcd, primes_first_n, factor
from sage.rings.all import (ZZ, QQ, CC, ComplexField, RealField,
                            is_NumberField, FiniteField)
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
from sage.misc.misc import prod
from sage.matrix.constructor import identity_matrix
from sage.rings.number_field.arithgroup_nf.all import HilbertModularGroup

def stoll_cremona_reduction(f, degree=None, algorithm='default',
                            precision=53, group='gl', transformation=False):
    r"""
    Reduces the polynomial f using the Stoll-Cremona method
    
    INPUT:
    
    - ``f`` - a univariate polynomial
    
    - ``degree`` - a positive integer (default: None). Interpret ``f`` as a
      homogeneous polynomial of degree ``degree`` in two variables, i.e.,
      as `F(X,Y) = Y^{degree} f(X/Y)`. If None, use the degree of `f`.
    
    - ``algorithm`` - ('default' or 'magma'). If set to 'magma', instead of
      using :func:`covariant_z0` to find the covariant `z_0`, the algorithm
      uses Magma to compute the 'better' covariant `z` (see [SC2002] for more
      details) (requires Magma to be installed)
      
    - ``precision`` - a positive integer (default=53). The precision to use
      for the covariants.

    - ``group`` (default='sl') -- a string: 'sl' for SL_2(O_K), 'gl' for
      GL_2(O_K), 'plus' for the subgroup GL_2(O_K)^+ of elements of totally
      positive determinant in GL_2(O_K).
    
    - ``transformation`` (default: False). Whether to also output a matrix
      (a,b;c,d) such that the output g satisfies
      g(z) = (cz+d)^degree * f((az+b)/(cz+d)).
      
    OUTPUT:
        
    A polynomial whose associated covariant point in the upper half plane
    is in a fundamental domain and that is `SL_2(ZZ)`-equivalent to f
    
    EXAMPLES:
    
    We check with Stoll and Cremona's examples::
        
        sage: from sage.schemes.hyperelliptic_curves.stoll_cremona import stoll_cremona_reduction
        sage: P.<x> = QQ[]
        sage: f = 19*x^8 - 262*x^7 + 1507*x^6 - 4784*x^5 + 9202*x^4 - 10962*x^3 + 7844*x^2 - 3040*x + 475
        sage: stoll_cremona_reduction(f)
        -x^8 + 6*x^7 - 7*x^6 - 12*x^5 + 27*x^4 - 4*x^3 - 19*x^2 + 10*x - 5
        
        sage: Q.<y> = QQ[]
        sage: f = y^6 + 30*y^5 + 371*y^4 + 2422*y^3 + 8813*y^2 + 16968*y + 13524
        sage: stoll_cremona_reduction(f)    
        y^6 - 4*y^4 + 2*y^3 + 8*y^2 - 12*y + 9
        
    An example over a real quadratic number field::
    
        sage: K.<a> = NumberField(x^2-5)
        sage: P.<x> = K[]
        sage: f = (x+a)^6+7; f
        x^6 + 6*a*x^5 + 75*x^4 + 100*a*x^3 + 375*x^2 + 150*a*x + 132
        sage: stoll_cremona_reduction(f)
        x^6 + 7

    The chosen degree is relevant, in the following example, the Igusa-Clebsch
    invariants show that the reduction for degree 5 is not equivalent for
    degree 6::
    
        sage: Q.<x> = QQ[]
        sage: f = 4*x^5 + 26*x^4 + 88*x^3 + 144*x^2 + 112*x + 32
        sage: C = HyperellipticCurve(f)
        sage: i = C.igusa_clebsch_invariants()
        sage: g = stoll_cremona_reduction(f); g
        2*x^5 + 4*x^4 + 4*x^3 + 24*x^2 - 6*x + 4
        
        sage: D = HyperellipticCurve(g)
        sage: j = D.igusa_clebsch_invariants()
        sage: i[0]^5/i[3] == j[0]^5/j[3]
        False
        
        sage: stoll_cremona_reduction(f, degree=6)
        Traceback (most recent call last):
        ...
        NotImplementedError: Please implement covariant_z0 for binary forms with roots at infinity
        
        sage: h = stoll_cremona_reduction(f, degree=6, algorithm='magma'); h # optional - magma
        4*x^5 + 6*x^4 + 24*x^3 - 4*x^2 + 4*x - 2
        sage: D = HyperellipticCurve(h) # optional - magma
        sage: j = D.igusa_clebsch_invariants() # optional - magma
        sage: i[0]^5/i[3] == j[0]^5/j[3]
        True


    With Magma, we can get a better reduction, as the covariant `z` is
    implemented there, rather than just the covariant `z0`::
    
        sage: P.<x> = QQ[]
        sage: f = good example, maybe from [SC2002]
        sage: stoll_cremona_reduction(f)
        
        sage: stoll_cremona_reduction(f, algorithm='magma') # optional - magma
        
    An example with ``transformation=True``::
    
        sage: g = 32*x^6 + 688*x^5 + 6144*x^4 + 29176*x^3 + 77714*x^2 + 110104*x + 64830
        sage: h, T = stoll_cremona_reduction(g, algorithm='magma', transformation=True)
        sage: T
        [ 3  4]
        [-1 -1]
        sage: g((3*x+4)/(-x-1))*(-x-1)^6 == h
        True

    ALGORITHM:
    
    This is the algorithm from [SC2002].
    
    If the base field is `QQ`, then the algorithm is exactly the same. For
    real quadratic fields, we implemented a reduction in the Hilbert upper
    half space HH^2. If the real quadratic field has class number one, then
    this region is the standard fundamental domain. For more general number
    fields, a more general reduction algorithm for the covariant point is
    needed.
    
    REFERENCES
    
    [SC2002] M. Stoll and J. Cremona "On the reduction theory of binary forms"
             J. Reine Angew. Math. 565 (2003) pp. 79-99
    """
#   s = [] 
    K = f.base_ring()
    P = f.parent()
    x = P.gen()
    G = HilbertModularGroup(K,group='GL')
#    h = _HilbertFundDomData(K, embs=precision, group=group)
    h = G.fundamental_domain()
    if degree is None:
        degree = f.degree()
    R = RealField(precision)
    C = ComplexField(precision)
    T = Matrix([[1,0],[0,1]])
    previous_height = R(-2)
    height = R(-1)
    one_minus_epsilon = R(1) - R(2)**(-5)    
    if algorithm == 'default':
        while previous_height < height * one_minus_epsilon:
            previous_height = height
            z = covariant_z0(f, precision=precision, degree=degree)
            M, z = h.step(z)
            [[a,b],[c,d]] = M
            f = P(f((d*x-b)/(-c*x+a)) * (-c*x+a)**degree)
#           s.append((f,z,covariant_z0(f)))
            T = T * Matrix([[d, -b], [-c, a]])
            height = prod([t.imag() for t in covariant_z0(f)]).abs()
            print f
            print height
            print previous_height
        if height < previous_height * one_minus_epsilon:
            raise RuntimeError, "Possibly too low precision in stoll_cremona_reduction"
    elif algorithm == 'magma':
        from sage.interfaces.all import magma
        z = 0; zp = 1
        while z != zp:
            Phi = K.embeddings(R)
            z = []
            for emb in Phi:
                gcoeffs = [emb(f[k]) for k in range(degree+1)]
                P = R['x']
                x = P.gen()
                g = sum([gcoeffs[indx]*x**indx for indx in range(degree+1)])
                zi = magma(magma.Covariant(degree, magma(g))).sage()
                z.append(C(zi))
            M,zp = h.step(z)
            [a,b],[c,d] = M
            P = K['x']
            x = P.gen()
            f = f.parent()(f((d*x-b)/(-c*x+a)) * (-c*x+a)**degree)
            T = T * Matrix([[d, -b], [-c, a]])
    else:
        raise ValueError("Unkown algorithm '%s'" % algorithm)   
    if transformation:
        return f, T
    return f


def covariant_z0(f, precision=53, degree=None):   
    r"""
    Returns the covariant `z_0` of f in H (upper half plane) from [SC2002].
    
    INPUT:
        
     - ``f`` - a polynomial
    
     - ``precision`` - An integer (Default = 53) What precision for real and
       complex numbers should the algorithm use
       
     - ``degree`` - the degree to homogenize with
    
    OUTPUT: 
    
    A point z in the upper half plane associated to f. If f is defined over
    a totally real number field of degree g, then the output is a sequence
    of length g.
    
    EXAMPLES::
    
        sage: z0 = sage.schemes.hyperelliptic_curves.stoll_cremona.covariant_z0
        sage: f = x^6 + 3*x^3 +4*x^2+ 3
        sage: z0(f)
        -0.0478810014556224 + 1.20469810839717*I
        sage: -z0(f)^-1
        
        sage: z0(f.reverse())
        
        sage: z0(f)-1
        
        sage: z0(f(x+1))
        
    Comparing the different answers using different precision ::

        sage: z0 = sage.schemes.hyperelliptic_curves.stoll_cremona.covariant_z0
        sage: P.<x> = QQ[]
        sage: f = x^8 + 24*x^7 + 3*x^5 +4*x^2+ 3
        sage: z0(f)
        -0.00832825511768667 + 0.835323919995269*I
        sage: z0(f,200)
        -0.0083282551176867045891624031810609648523736031152551979546157 + 0.83532391999526880130877046174257268976426032236105109987080*I
    
    If f is not in `\QQ` then ``covariant_z0`` returns a list where every entry
    is a z associated to one embeddings of f in `\RR`::
    
        sage: z0 = sage.schemes.hyperelliptic_curves.stoll_cremona.covariant_z0
        sage: X = var('X')
        sage: k = NumberField(X^2-41,'a')
        sage: P.<x> = k[]
        sage: a = k.an_element()
        sage: f = x^8 + 24*x^7 + 3*x^5 +4*x^2+ 3
        sage: z0(f)
        [-0.00832825511768667 + 0.835323919995269*I, -0.00832825511768667 + 0.835323919995269*I]
        sage: f = x^6 + (a+1)*x^5 + (5*a+2)*x^2 + (a/2+1/2)*x + 5*a
        [0.167978514759694 + 1.68636006242584*I, -0.125921323440448 + 1.64323412027265*I]

    Even though Stoll and Cremona allow roots of higher multiplicity, and roots
    at infinity (i.e. degree > f.degree()), this is currently not fully
    supported::
        
        sage: z0 = sage.schemes.hyperelliptic_curves.stoll_cremona.covariant_z0
        sage: P.<x> = RR[]
        sage: i = (x-3)^2 * (x^2 + 2)^2
        sage: z0(i)
        Traceback (most recent call last):
        ...
        ValueError: Cannot convert NaN or infinity to Pari float
        
        sage: j = x^8 + 3
        sage: z0(j, degree=9)
        Traceback (most recent call last):
        ...
        NotImplementedError
        
    ALGORITHM:
        
    This algorithm is taken straight from [SC2002] page 6
    
    Let n be the degree of the polynomial f. For each embeddings of f into
    `\RR` let `\alpha_j` be the roots of f. Define a quadratic polynomial
    `Q_0` by `Q_0(x) = \sum_{j=1}^{n}
                       (x-\alpha_j)(x-\bar{\alpha_j})/|f'(\alpha_j)|^{2/(n-2)}`
    Finally take `z_0` to be the root with positive imaginary part of `Q_0`
    
    REFERENCE:
    
    [SC2002] M. Stoll and J. Cremona "On the reduction theory of binary forms"
             J. Reine Angew. Math. 565 (2003) pp. 79-99
    """    
    C_out = ComplexField(precision)
 
    internal_precision = 6*precision # TODO: make an informed choice
    C = ComplexField(internal_precision)
    R = RealField(internal_precision)
    

    if degree is None:
        degree = f.degree()
    elif degree != f.degree():
        raise NotImplementedError("Please implement covariant_z0 for "\
                                  "binary forms with roots at infinity")
    
    k = f.base_ring()
    x = R['x'].gen()
    if k is QQ or is_NumberField(k):
        Phi = k.embeddings(C)
    else:
        Phi = [C.convert_map_from(k)]
    I = C.gen()
    
    # Actually, the roots should be taken with multiplicities,
    # but this does not matter, as the relevant coefficients
    # are divided by infinity anyway.
    
    z_list = []
    for emb in Phi:
        gcoeffs = [emb(coeff) for coeff in f.coeffs()]
        g = sum([gcoeffs[indx] * x**indx for indx in range(degree+1)])
        roots = g.roots(C, multiplicities=False)
        gdiv = g.derivative()
        q0 = R(0)
        for ai in roots:
            q0 += ((x - ai) * (x - ai.conjugate()) /
                   (abs(gdiv(ai))**(2 / (degree - 2))))
        z1,z2 = q0.roots(C, multiplicities=False)
        if z1.imag_part() > 0:
            z0 = z1
        else:
            z0 = z2
        z_list.append(z0)

    return [C_out(z) for z in z_list]

