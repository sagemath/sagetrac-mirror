# -*- coding: utf-8 -*-
r"""
Groebner basis computation in finite fields using F4 algorithm.

AUTHORS:

- Vanessa VITSE (2015)

- Antoine JOUX (2015)

- Titouan COLADON (2015)

EXAMPLES:

We compute a Groebner basis for some given ideal under the degrevlex ordering using F4::

    sage: from sage.rings.polynomial.groebner_basis_f4 import groebner_basis_f4
    sage: R.<x1,x2,x3,x4,x5,x6,x7,x8> = Zmod(65521)[]
    sage: I = sage.rings.ideal.Cyclic(R,8)
    sage: B = groebner_basis_f4(I)
    sage: type(B)
    <class 'sage.rings.polynomial.multi_polynomial_sequence.PolynomialSequence_generic'>

TESTS::

    sage: from sage.rings.polynomial.groebner_basis_f4 import groebner_basis_f4
    sage: R.<a,b,c,d,e,f> = Zmod(65521)[]
    sage: I = sage.rings.ideal.Cyclic(R,6)
    sage: B = groebner_basis_f4(I)

    sage: from sage.rings.polynomial.groebner_basis_f4 import groebner_basis_f4
    sage: F.<t>=GF(2)[]
    sage: K.<t>=GF(2^31, name='t', modulus=t^31+t^3+1)
    sage: R.<x0,x1,x2,x3,x4,x5> = K[]
    sage: I = ideal((t+t^3)*x0+(t+t^3)*x1+(t+t^3)*x2+(t+t^3)*x3+(t+t^3)*x4+(t+t^3)*x5, (t+t^3)*x0*x1+(t+t^3)*x1*x2+(t+t^3)*x2*x3+(t+t^3)*x3*x4+(t+t^3)*x0*x5+(t+t^3)*x4*x5, (t+t^3)*x0*x1*x2+(t+t^3)*x1*x2*x3+(t+t^3)*x2*x3*x4+(t+t^3)*x0*x1*x5+(t+t^3)*x0*x4*x5+(t+t^3)*x3*x4*x5, (t+t^3)*x0*x1*x2*x3+(t+t^3)*x1*x2*x3*x4+(t+t^3)*x0*x1*x2*x5+(t+t^3)*x0*x1*x4*x5+(t+t^3)*x0*x3*x4*x5+(t+t^3)*x2*x3*x4*x5, (t+t^3)*x0*x1*x2*x3*x4+(t+t^3)*x0*x1*x2*x3*x5+(t+t^3)*x0*x1*x2*x4*x5+(t+t^3)*x0*x1*x3*x4*x5+(t+t^3)*x0*x2*x3*x4*x5+(t+t^3)*x1*x2*x3*x4*x5, (t+t^3)*x0*x1*x2*x3*x4*x5-1)
    sage: B = groebner_basis_f4(I)


"""


#*****************************************************************************
#
#                               Sage
#
#       Copyright (C) 2015 Antoine Joux and Vanessa Vitse
#       Copyright (C) 2015 Titouan COLADON
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
#*****************************************************************************


#include <Python.h>

from libcpp.string cimport string
from libcpp.vector cimport vector

from libc.stdint cimport int64_t

include "sage/ext/interrupt.pxi"


cdef extern from "libf4.h":
    cdef vector[string] groebnerBasisF4(int64_t modulus,
                                        int nbVariable,
                                        vector[string] variableName,
                                        vector[string] polynomialList,
                                        int nbThread, int verbose)
    cdef vector[string] groebnerBasisGF2F4(int nbVariable,
                                           vector[string] variableName,
                                           vector[string] polynomialList,
                                           int nbThread, int verbose)
#    cdef vector[string] groebnerBasisGivaroIntegerF4(string modulus, int nbVariable, vector[string] variableName, vector[string] polynomialList, int nbThread, int verbose)
    cdef vector[string] groebnerBasisGF2ExtensionF4(string modulus,
                                                    int nbVariable,
                                                    vector[string] variableName,
                                                    string polyVarName,
                                                    vector[string] polynomialList,
                                                    int nbThread,
                                                    int verbose)

def groebner_basis_f4(self, verbose = 0, nb_thread = 1):
        """
        Computes a Groebner Basis for this ideal using F4

        INPUT:

        EXAMPLES::

            sage: from sage.rings.polynomial.groebner_basis_f4 import groebner_basis_f4
            sage: R.<x1,x2,x3,x4,x5,x6> = Zmod(65521)[]
            sage: I = sage.rings.ideal.Cyclic(R,6)
            sage: B = groebner_basis_f4(I) # optional - f4
            sage: len(B)
            45

        TESTS::

            sage: P = PolynomialRing(GF(next_prime(2^32)), 8, 'x')
            sage: I = sage.rings.ideal.Cyclic(P)
            sage: gb0 = I.groebner_basis('f4') # optional - f4
            Traceback (most recent call last):
            ...
            NotImplementedError: Prime field with characteristic > 2^32 are not handled in Sage for the moment

        """
        R = self.ring()
        polynomial_list = self.gens()
        cdef vector[string] variable_name = [str(v) for v in R.gens()]
        cdef vector[string] polynomial_list_cpp = [str(poly).replace(" ", "") for poly in polynomial_list]
        cdef vector[string] basis;
        cdef int nb_variable = R.ngens()
        cdef int64_t modulus_int
        cdef string modulus_str
        cdef string poly_var_name = str(self.base_ring().gen())

        # Prime field
        if self.base_ring().is_prime_field():
            if R.characteristic() == 2:
                sig_on()
                basis = groebnerBasisGF2F4(nb_variable, variable_name, polynomial_list_cpp, nb_thread, verbose);
                sig_off()
            elif R.characteristic() < 2**32:
                modulus_int = R.characteristic()
                sig_on()
                basis = groebnerBasisF4(modulus_int, nb_variable, variable_name, polynomial_list_cpp, nb_thread, verbose);
                sig_off()
            else:
                # not handled gmp and givaro version in Sage are too old
                raise NotImplementedError("Prime field with characteristic > 2^32 are not handled in Sage for the moment")
#                modulus_str = str(R.characteristic())
#                sig_on()
#                basis = groebnerBasisGivaroIntegerF4(modulus_str, nb_variable, variable_name, polynomial_list_cpp, nb_thread, verbose);
#                sig_off()

        if not self.base_ring().is_prime_field():
            if R.characteristic() != 2:
                raise NotImplementedError("Non prime field with characteristic != 2 are not handled by F4")
            elif self.base_ring().modulus().degree() > 63:
                raise NotImplementedError("GF(2^n) field with n > 63 are not handled by F4")
            else:
                modulus_str =str(self.base_ring().modulus()).replace("x", str(self.base_ring().gen()))
                sig_on()
                basis = groebnerBasisGF2ExtensionF4(modulus_str, nb_variable, variable_name, poly_var_name, polynomial_list_cpp, nb_thread, verbose);
                sig_off()

        from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
        base_ring = R.base_ring()
        gens_dict = R.gens_dict()
        cache = {}
        r = PolynomialSequence([convert_from_f4_string(e, gens_dict, base_ring, cache) for e in basis], R, immutable=True)
        return r


def convert_from_f4_string(f_string, gens_dict, base_ring, cache=None):
    """
    Convert F4 string representation to Sage polynomial.

    INPUT:

    - fs_string - string representation
    - gens_dict - gens dict of a multivariate polynomial ring
    - base_ring - base ring of a multivariate polynomial ring
    - cache - optional variable power cache

    EXAMPLE::

        sage: P.<x,y,z> = PolynomialRing(GF(previous_prime(2^30)))
        sage: fs = '(1*x^1*y^2) + (-32*z^2)'
        sage: from sage.rings.polynomial.groebner_basis_f4 import convert_from_f4_string
        sage: convert_from_f4_string(fs, P.gens_dict(), P.base_ring())
        x*y^2 - 32*z^2

        sage: cache = {}
        sage: convert_from_f4_string(fs, P.gens_dict(), P.base_ring(), cache)
        x*y^2 - 32*z^2

        sage: cache
        {'x^1': x, 'y^2': y^2, 'z^2': z^2}
    """

    M = f_string.split(" + ")
    f = base_ring(0)

    if cache is None:
        cache = {}

    for m in M:
        m = m[1:-1]
        m = m.split("*")
        c, p = m[0], m[1:]
        c = base_ring(c)
        m = c
        for p_ in p:
            if p_ not in cache:
                cache[p_] = eval(p_.replace("^", "**"), gens_dict)
            m *= cache[p_]
        f += m
    return f
