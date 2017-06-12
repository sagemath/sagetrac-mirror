r"""
Gauss-Legendre integration for vector-valued functions

Routine to perform Gauss-Legendre integration for vector-functions.

AUTHORS:

 - Nils Bruin (2017-06-06): initial version

EXAMPLES::

NOTE:

The code here is directly based on mpmath (see http://mpmath.org), but has a highly
optimized routine to compute the nodes.
"""

#*****************************************************************************
#       Copyright (C) 2017 Nils Bruin <nbruin@sfu.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#as it turns out, computing the nodes can easily turn out to be more
#expensive than computing the integrals. So it's worth optimizing this.
#making the function into a cython routine helps a little bit. If we really
#want to we can optimize this further, probably to a point where
#we don't have to bother with node computation routines that have a better order
#than this naive approach (which is quadratic)
from __future__ import absolute_import, division, print_function
from sage.libs.mpfr cimport *
import math
from sage.rings.real_mpfr import RealField
from sage.misc.cachefunc import cached_function
from sage.rings.real_mpfr cimport RealNumber, RealField_class

@cached_function
def nodes(degree,prec):
    r"""
    Compute the integration nodes and weights for the Gauss-Legendre quadrature scheme.

    INPUT:

     - ``degree`` -- integer. The number of nodes. Must be 3 or even.

     - ``prec`` -- integer (minimal value 53). Binary precision with which the nodes and weights are computed.

    OUTPUT:

    A list of (node,weight) pairs.

    EXAMPLE:

    The nodes for the Gauss-Legendre scheme are roots of Legendre polynomials.
    The weights can be computed by a straightforward formula (note that evaluating
    a derivative of a Legendre polynomial isn't particularly numerically stable, so the results
    from this routine are actually more accurate than what the values the closed formula produces)::

        sage: from sage.numerical.gauss_legendre import nodes
        sage: L1=nodes(24,53)
        sage: P=RR['x'](sage.functions.orthogonal_polys.legendre_P(24,x))
        sage: Pdif=P.diff()
        sage: L2=[( (r+1)/2,1/(1-r^2)/Pdif(r)^2) for r,_ in RR['x'](P).roots()]
        sage: all((a[0]-b[0]).abs() < 10^-15 and (a[1]-b[1]).abs() < 10^-9 for a,b in zip(L1,L2))
        True
    """
    cdef long j,j1,n
    cdef RealNumber r,t1,t2,t3,t4,a,w
    cdef mpfr_t u,v
    cdef RealField_class R
    if prec < 53:
        prec = 53
    if degree !=3 and degree % 2 !=0:
        raise ValueError("degree=%s not supported (degree must be 3 or even)"%degree)
    R = RealField(int(prec*3/2))
    Rout = RealField(prec)
    mpfr_init2(u,R.__prec)
    mpfr_init2(v,R.__prec)
    ZERO = R.zero()
    ONE = R.one()
    HALF = ONE/2
    TWO = 2*ONE
    rnd = R.rnd
    epsilon = R(1)>>(prec+8)
    if degree == 1:
        x = (R(3)/5).sqrt()
        w = R(5)/18
        nodes = [((1-x)/2,w),(HALF,R(4)/9),((1+x)/2,w)]
    else:
        nodes = []
        n = degree
        upto = n//2+1
        for j in xrange(1,upto):
            r = R(math.cos(math.pi*(j-0.25)/(n+0.5)))
            while True:
                t1,t2=ONE,ZERO
                for j1 in xrange(1,n+1):
                    mpfr_mul(u,r.value,t1.value,rnd)
                    mpfr_mul_si(u,u,2*j1-1,rnd)
                    mpfr_mul_si(v,t2.value,j1-1,rnd)
                    mpfr_sub(u,u,v,rnd)
                    mpfr_div_si(u,u,j1,rnd)
                    t2=t1
                    t1=R._new()
                    mpfr_set(t1.value,u,rnd)
                t4 = R(n)*(r*t1-t2)/(r**2-ONE)
                a = t1/t4
                r = r-a
                if a.abs()<epsilon:
                    break
            x = r
            w = ONE/((ONE-r**2)*t4**2)
            nodes.append(((ONE+x)/TWO,w))
            nodes.append(((ONE-x)/TWO,w))
    nodes=[(Rout(x),Rout(w)) for x,w in nodes]
    nodes.sort()
    mpfr_clear(u)
    mpfr_clear(v)
    return nodes

def estimate_error(results,prec,epsilon):
    r"""
    Routine to estimate the error in a list of quadrature approximations.
    
    The method used is based on Borwein, Bailey, and Girgensohn. As mentioned in
    mpmath: Although not very conservative, this method seems to be very robust in
    practice.

    The routine takes a list of vector results and, under assumption that these
    vectors approximate a given vector approximately quadratically, gives an estimate
    of the maximum norm of the error in the last approximation.

    INPUT:
    
     - ``results`` -- list. List of approximations to estimate the error from. Should be at least length 2.

     - ``prec`` -- integer. Binary precision at which computations are happening.

     - ``epsilon`` -- multiprecision float. Default error estimate in case of insufficient data.

    OUTPUT:

    An estimate of the error.

    EXAMPLE::

        sage: from sage.numerical.gauss_legendre import estimate_error
        sage: prec=200
        sage: K=RealField(prec)
        sage: V=VectorSpace(K,2)
        sage: a=V([1,-1])
        sage: b=V([1,1/2])
        sage: L=[a+2^(-2^i)*b for i in [0..5]]
        sage: estimate_error(L,prec,K(2^(-prec)))
        2.328235...e-10

    """
    if len(results)==2:
        return max((results[0][i]-results[1][i]).abs() for i in xrange(len(results[0])))
    e = []
    ZERO = 0*epsilon
    for i in xrange(len(results[0])):
        try:
            if results[-1][i] == results[-2][i] == results[-3][i]:
                e.append(0*epsilon)
            D1 = (results[-1][i]-results[-2][i]).abs().log()
            D2 = (results[-1][i]-results[-3][i]).abs().log()
        except ValueError:
            e.append(epsilon)
        #we follow mpmath in clipping the precision
        D4 = min(ZERO,max(D1**2/D2,2*D1,ZERO-prec))
        e.append(D4.exp())
    return max(e)

def integrate_vector(f,prec,epsilon=None):
    r"""
    Integrate a one-argument vector-valued function numerically using Gauss-Legendre.

    This function uses the Gauss-Legendre quadrature scheme to approximate
    the integral of f(t) for t=0..1.

    INPUT:

     - `f` -- callable. Vector-valued integrand.

     - `prec` -- integer. Binary precision to be used.

     - `epsilon` -- Multiprecision float. Target error bound.

    OUTPUT:

    Vector approximating value of the integral.

    EXAMPLE::

        sage: from sage.numerical.gauss_legendre import integrate_vector
        sage: prec=200
        sage: K=RealField(prec)
        sage: V=VectorSpace(K,2)
        sage: epsilon=K(2^(-prec+4))
        sage: f=lambda t:V((1+t^2,1/(1+t^2)))
        sage: I=integrate_vector(f,prec,epsilon)
        sage: J=V((4/3,pi/4))
        sage: max(c.abs() for c in (I-J)) < epsilon
        True

    We can also use complex-valued integrands::

        sage: prec=200
        sage: Kreal=RealField(prec)
        sage: K=ComplexField(prec)
        sage: V=VectorSpace(K,2)
        sage: epsilon=Kreal(2^(-prec+4))
        sage: f=lambda t: V((t,K(exp(2*pi*t*K.0))))
        sage: I=integrate_vector(f,prec,epsilon)
        sage: J=V((1/2,0))
        sage: max(c.abs() for c in (I-J)) < epsilon
        True
    """
    results = []
    cdef long degree = 3
    Rout = RealField(prec)
    if epsilon is None:
        epsilon = Rout(2)**(-prec+3)
    while True:
        nodelist = nodes(degree,prec)
        I = nodelist[0][1]*f(nodelist[0][0])
        for i in xrange(1,len(nodelist)):
            I += nodelist[i][1]*f(nodelist[i][0])
        results.append(I)
        if degree > 3:
            err = estimate_error(results,prec,epsilon)
            if err <= epsilon:
                return I
        #double the degree to double expected precision
        degree *= 2
