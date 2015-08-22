##############################################################################
#       Copyright (C) 2011 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################
"""
Chow-Heegner Points (fast Cython code)

This module implements fast compiled helper code in Cython that is
useful for computing Chow-Heegner points, but may be more generally
useful for fast numerical applications.  This includes code for
arithmetic with the Sage ComplexNumber type that avoids any object
creation, a complex polynomial class optimized for fast evaluation
when the degree is relatively large, and a polynomial class for
polynomials with RDF coefficients that can compute roots and do Newton
iteration via fast code in the GSL (much faster than PARI and numpy).

AUTHORS:

- William Stein (November 2011, January 2012)

"""

##############################################################################
# Implementation Note: one could make a case that most of the
# functions in the file could be moved to other places in the Sage
# library where they might speed up all code.  We do not plan to do
# that for the first version of this patch, since this code could also
# slow down code in low degree, and may have worrisome numerical
# rounding implications for other code.  By keeping it here, we keep
# this code as simple as possible as well, since we don't have to
# optimize for general cases that aren't need for our application.
#
# The right way to do this is to make a completely new complete
# RDF['x'] class that uses GSL.  For example, for an instance of the
# class Polynomial_RDF_gsl of a polynomial of degree 1000, evaluation
# is over 100 times faster than for the same element of RDF['x']!
# So this has to happen.  However, it's a hell of a lot more work
# than implementing Polynomial_RDF_gsl below.
##############################################################################

from sage.libs.gsl.types cimport *
from sage.libs.gsl.complex cimport *
from sage.libs.gsl.poly cimport *

cdef extern from "gsl/gsl_poly.h":
    gsl_complex gsl_poly_complex_eval (double c[], int n, gsl_complex z)

include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'

from sage.rings.all import CDF
from sage.rings.all import RR
from sage.rings.complex_double cimport ComplexDoubleElement
from sage.stats.intlist cimport IntList

from sage.rings.complex_number cimport (ComplexNumber,
     mpfr_t, mpfr_init2, mpfr_mul, mpfr_sub, mpfr_add, mpfr_clear, GMP_RNDN)


def required_series_prec(y, prec):
    """
    Return integer B so that using B terms of the series of a modular
    parameterization the tail end of the sum is less than 2**(-prec)
    in absolute value for any point z with Im(z)>=ymin.

    INPUT:
    
    - ``ymin`` -- positive real number; minimum y coordinate
    - ``prec`` -- positive integer (bits of precision)

    OUTPUT:

    - Python int

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner_fast import required_series_prec
        sage: required_series_prec(1e-4, 53)
        71306
        sage: required_series_prec(1e-5, 53)
        749700

    The following means that for any elliptic curve and any point z in
    the upper half plane with imaginary part at least 1e-4, the tail
    end of the modular parametrization evaluated at z is less than
    2^(-200) in absolute value::
    
        sage: required_series_prec(1e-4, 200)
        233473        
    """
    y = RR(y)
    epsilon = RR(2)**(-(prec+1))
    pi = RR.pi()
    return int((epsilon*(1 - (-2*pi*y).exp())).log() / (-2*pi*y)) + 1

###############################################################
# Cython code to efficiently evalute a polynomial with
# ComplexField(prec) coefficients very efficiently.
# Below we have functions for * and + of Cython ComplexNumber
# objects.  This code is basically same as code in the
# complex_number.pyx file; however, it allows for
# mutating the output argument in place without having
# to create a new ComplexNumber argument, which is something
# the code in complex_number.pyx simply does not allow.
# This results in a huge speedup.   The "right" way to do this
# would probably be to completely refactor complex_number.pyx
# so that all arithmetic functions there call functions
# like the ones below.  However, that would be a lot more work
# for a questionable payback.
###############################################################

cdef int ComplexNumber_mul(ComplexNumber x, ComplexNumber left, ComplexNumber right):
    """
    Multiply two ComplexNumber objects.
    """
    cdef mpfr_t t0, t1
    mpfr_init2(t0, left._prec)
    mpfr_init2(t1, left._prec)
    mpfr_mul(t0, left.__re, right.__re,  GMP_RNDN)
    mpfr_mul(t1, left.__im, right.__im,  GMP_RNDN)
    mpfr_sub(x.__re, t0, t1,  GMP_RNDN)
    mpfr_mul(t0, left.__re, right.__im,  GMP_RNDN)
    mpfr_mul(t1, left.__im, right.__re,  GMP_RNDN)
    mpfr_add(x.__im, t0, t1,  GMP_RNDN)
    mpfr_clear(t0)
    mpfr_clear(t1)

cdef int ComplexNumber_add(ComplexNumber x, ComplexNumber left, ComplexNumber right):
    """
    Add two ComplexNumber objects.
    """
    mpfr_add(x.__re, left.__re, right.__re,  GMP_RNDN)
    mpfr_add(x.__im, left.__im, right.__im,  GMP_RNDN)

cdef class ComplexPolynomial:
    """
    A polynomial with floating point complex number entries that is
    optimized for fast evaluation when the degree is relatively large.

    TESTS::

        sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
        sage: f = ComplexPolynomial(ComplexField(100)['x']([1,2]))
        sage: TestSuite(f).run()
    """
    cdef list a
    cdef readonly object f
    def __init__(self, f):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
            sage: f = ComplexPolynomial(ComplexField(100)['x']([1,2])); f
            2.0000000000000000000000000000*x + 1.0000000000000000000000000000
            sage: type(f)
            <type 'sage.schemes.elliptic_curves.chow_heegner_fast.ComplexPolynomial'>        
        """
        self.a = f.list()
        self.f = f

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
            sage: f = ComplexPolynomial(ComplexField(100)['x']([1,2])); f.__repr__()
            '2.0000000000000000000000000000*x + 1.0000000000000000000000000000'
        """
        return repr(self.f)

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
            sage: f = ComplexPolynomial(ComplexField(100)['x']([1,2])); g = ComplexPolynomial(ComplexField(100)['x']([1,3]))
            sage: cmp(f,f)
            0
            sage: cmp(f,g)
            -1
            sage: cmp(g,f)
            1
        """
        return cmp(self.f, right.f)

    def __reduce__(self):
        """
        Used in pickling.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
            sage: ComplexPolynomial(ComplexField(100)['x']([1,2])).__reduce__()
            (<type 'sage.schemes.elliptic_curves.chow_heegner_fast.ComplexPolynomial'>, (2.0000000000000000000000000000*x + 1.0000000000000000000000000000,))        
        """
        return ComplexPolynomial, (self.f,)

    def __call__(self, z0):
        """
        Evaluate this polynomial at z0.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
            sage: f = ComplexField(100)['x']([1..20]); g = ComplexPolynomial(f)
            sage: f(1)
            210.00000000000000000000000000
            sage: g(1)
            210.00000000000000000000000000        
        """
        cdef list a = self.a
        C = a[0].parent()
        cdef ComplexNumber z = C(z0)
            
        # Use Horner's rule (see http://en.wikipedia.org/wiki/Horner_scheme)
        cdef Py_ssize_t n = len(a)-1
        cdef ComplexNumber t = C(0), b = C(a[n].real(), a[n].imag())  # b must be a new copy -- do not do "b = a[n]" here or you'll get a huge bug!
        
        while n >= 1:
            n -= 1
            # Fast version of: b = a[n] + b*z
            ComplexNumber_mul(t, b, z)      # t = b*z
            ComplexNumber_add(b, a[n], t)   # b = a[n] + t
        return b

    def degree(self):
        """
        OUTPUT:
        
        - integer
            
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import ComplexPolynomial
            sage: ComplexPolynomial(ComplexField(100)['x']([1..10])).degree()
            9
        """
        return self.f.degree()


cdef class Polynomial_RDF_gsl:
    """
    A polynomial with RDF (real double) coefficients, which can
    compute roots and do Newton iteration via fast code in the GSL
    library.

    TESTS::

        sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
        sage: s = Polynomial_RDF_gsl(RDF['x']([1..3]))
        sage: TestSuite(s).run()    
    """
    cdef double *c
    cdef double *c_d
    cdef readonly object f
    cdef int n
    cdef int n_d

    def __cinit__(self, f):
        """
        EXAMPLES::

        Make an object of type Polynomial_RDF_gsl::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..3])); type(s)
            <type 'sage.schemes.elliptic_curves.chow_heegner_fast.Polynomial_RDF_gsl'>
            sage: s.f
            3.0*x^2 + 2.0*x + 1.0

        Check that the f attribute is read only::
        
            sage: s.f = 0
            Traceback (most recent call last):
            ...
            AttributeError: attribute 'f' of 'sage.schemes.elliptic_curves.chow_heegner_fast.Polynomial_RDF_gsl' objects is not writable        
        """
        self.c = NULL
        self.c_d = NULL
        self.f = f

    def __reduce__(self):
        """
        Used in pickling.
        
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..3]))
            sage: s.__reduce__()
            (<type 'sage.schemes.elliptic_curves.chow_heegner_fast.Polynomial_RDF_gsl'>, (3.0*x^2 + 2.0*x + 1.0,))
        """
        return Polynomial_RDF_gsl, (self.f, )

    def __cmp__(self, right):
        """
        Compares underlying polynomials.
        
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..3])); t = Polynomial_RDF_gsl(RDF['x']([3,2,1]))
            sage: s == s
            True
            sage: s == t
            False
            sage: cmp(s,t)
            1
            sage: cmp(t,s)
            -1
            sage: cmp(t.f,s.f)
            -1        
        """
        return cmp(self.f, right.f)
    
    def __init__(self, f):
        """
        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..5])); s
            5.0*x^4 + 4.0*x^3 + 3.0*x^2 + 2.0*x + 1.0
            sage: type(s)
            <type 'sage.schemes.elliptic_curves.chow_heegner_fast.Polynomial_RDF_gsl'>
        """
        self.n = f.degree() + 1
        self.c = <double*> sage_malloc(sizeof(double)*self.n)
        if not self.c:
            self.c = NULL
            raise MemoryError
        cdef list v = f.list()
        cdef Py_ssize_t i        
        for i in range(self.n):
            self.c[i] = v[i]

    def __repr__(self):
        """
        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: Polynomial_RDF_gsl(RDF['x']([-2.5,3,-1])).__repr__()
            '-x^2 + 3.0*x - 2.5'
        """
        return self.f.__repr__()

    def __dealloc__(self):
        if self.c: sage_free(self.c)
        if self.c_d: sage_free(self.c_d)        

    def __call__(self, x):
        """
        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..100]))
            sage: s(1)
            5050.0
            sage: s.f(1)
            5050.0            
        """
        cdef ComplexDoubleElement w, z
        cdef gsl_complex a
        if isinstance(x, ComplexDoubleElement):
            z = x
        else:
            z = CDF(x)
        a = gsl_poly_complex_eval(self.c, self.n, z._complex)
        w = z._new_c(a)
        return w

    def _init_deriv(self):
        """
        Initially derivative attribute, which is used in the __call__ method.
        
        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..5]))
            sage: s._init_deriv()
        """
        if self.c_d: return
        f_d = self.f.derivative()
        self.n_d = f_d.degree()+1
        self.c_d = <double*> sage_malloc(sizeof(double)*self.n_d)
        cdef list v = f_d.list()
        cdef Py_ssize_t i
        for i in range(self.n_d):
            self.c_d[i] = v[i]

    def roots(self):
        """
        OUTPUT:
        
        - list of complex double precision numbers
        
        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..3]))
            sage: s
            3.0*x^2 + 2.0*x + 1.0
            sage: s.roots()  # abs tol 1e-10
            [-0.333333333333 - 0.471404520791*I, -0.333333333333 + 0.471404520791*I]
            sage: s.f.roots(CDF, multiplicities=False)  # abs tol 1e-10
            [-0.333333333333 - 0.471404520791*I, -0.333333333333 + 0.471404520791*I]

        A higher degree example::
        
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..400]))
            sage: v = s.roots()
            sage: w = s.f.roots(CDF, multiplicities=False)
            sage: max([abs(v[i]-w[i]) for i in range(len(v))]) < 1e-12
            True
        """
        return cdf_roots_of_rdf_poly(self.f)

    def degree(self):
        """
        OUTPUT:
        
        - integer
            
        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: Polynomial_RDF_gsl(RDF['x']([1..10])).degree()
            9
        """
        return self.f.degree()

    def newton(self, x, int max_iter=1000, double max_err=1e-14):
        """
        Refine approximate double precision roots x using Newton's
        method.  Return list of triples consisting of the refined
        approximation, number of iterations and the error as output by
        GSL.

        EXAMPLES::
        
            sage: from sage.schemes.elliptic_curves.chow_heegner_fast import Polynomial_RDF_gsl
            sage: s = Polynomial_RDF_gsl(RDF['x']([1..4])); s
            4.0*x^3 + 3.0*x^2 + 2.0*x + 1.0
            sage: s.f.roots(CDF, multiplicities=False)  # abs tol 1e-10
            [-0.605829586188, -0.0720852069059 - 0.638326735148*I, -0.0720852069059 + 0.638326735148*I]
            sage: s.newton(0)  # abs tol 1e-10
            [(-0.605829586188, 7, 0.0)]

        The above used 7 iterations; let's restrict to at most 4 iterations::
        
            sage: s.newton(0, 4)  # abs tol 1e-8
            [(-0.6058300580523304, 4, 0.000552920671073931)]

        A very small max_err::

            sage: s.newton(0, max_err=0.1)  # abs tol 1e-10
            [(-0.606382978723, 3, 0.018617021276595702)]        

        Start at a complex approximate root::
        
            sage: s.newton(-.6*I)  # abs tol 1e-10
            [(-0.0720852069059 - 0.638326735148*I, 6, 1.3877787807814457e-17)]

        Make x a list of approximate roots::
            sage: s.newton([0, -.6*I])  # abs tol 1e-10
            [(-0.605829586188, 7, 0.0), (-0.0720852069059 - 0.638326735148*I, 6, 1.3877787807814457e-17)]
        """
        self._init_deriv()
        cdef int i
        cdef gsl_complex root, last_root
        cdef double err
        if not isinstance(x, list):
            v = [x]
        else:
            v = x

        ans = []
        cdef ComplexDoubleElement z
        for x in v:
            z=CDF(x)
            GSL_SET_COMPLEX(&root, z._complex.dat[0], z._complex.dat[1])
            GSL_SET_COMPLEX(&last_root, z._complex.dat[0], z._complex.dat[1])        
            sig_on()
            for i in range(max_iter):
                # We recode what would be the following simple Python
                # algorithm in GSL instead:
                #    root = root - f(root) / f_prime(root)
                #    if abs(last_root - root) < tiny:
                #       break
                #    else:
                #       last_root = root
                root = gsl_complex_sub(root,
                                gsl_complex_div(gsl_poly_complex_eval(self.c, self.n, root),
                                    gsl_poly_complex_eval(self.c_d, self.n_d, root)))
                err = gsl_complex_abs(gsl_complex_sub(last_root, root))
                if err <= max_err:
                    break
                GSL_SET_COMPLEX(&last_root, root.dat[0], root.dat[1])                    
            sig_off()
            ans.append((z._new_c(root), i+1, err))
        return ans

###############################################################
# Cython code to efficiently find all complex roots of a real
# polynomial using GSL.  This is *much* faster than using numpy.
###############################################################

def cdf_roots_of_rdf_poly(f):
    """
    Return the CDF roots of a polynomial with coefficients in RDF.

    Uses a very fast function in GSL that works by computing the
    eigenvalues of the companion matrix. 
    
    INPUT:

    - f -- polynomial with RDF coefficients
        
    OUTPUT:
    
    - list -- all CDF roots of f

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.chow_heegner_fast import cdf_roots_of_rdf_poly
        sage: f = RDF['x']([1..4])
        sage: cdf_roots_of_rdf_poly(f)  # abs tol 1e-10
        [-0.605829586188, -0.0720852069059 - 0.638326735148*I, -0.0720852069059 + 0.638326735148*I]
        sage: f.roots(CDF, multiplicities=False)  # abs tol 1e-10
        [-0.605829586188, -0.0720852069059 - 0.638326735148*I, -0.0720852069059 + 0.638326735148*I]    
    """
    cdef Py_ssize_t i, n = f.degree() + 1
    cdef double* a = <double*> sage_malloc(sizeof(double)*n)
    if not a:
        raise MemoryError
    cdef double* z = <double*> sage_malloc(sizeof(double)*2*n)
    if not z:
        sage_free(a)
        raise MemoryError
    
    cdef list v = f.list()
    for i in range(n):
        a[i] = v[i]

    cdef gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(n)
    sig_on()
    gsl_poly_complex_solve(a, n, w, z)
    sig_off()
    gsl_poly_complex_workspace_free(w)
    
    rts = [CDF(z[2 * i], z[2 * i + 1]) for i in range(n - 1)]
    rts.sort()
    
    sage_free(a)
    sage_free(z)
    
    return rts
