"""
Cython types for elements of path algebras
"""
#*****************************************************************************
#     Copyright (C) 2014 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# This file declares the types. The implementation of the basic on these types
# is in algebra_elements.pxi, the implementation of the Python class is in
# algebra_elements.pyx. The latter file also contains all doctests.

from cpython cimport PyObject
from sage.data_structures.bounded_integer_sequences cimport *
from sage.structure.element cimport RingElement, ModuleElement, Element
from sage.quivers.paths cimport QuiverPath

# Type definitions

cdef struct path_mon_t:
    # if mid==-1, then we encode an element of a path semigroups, i.e., we use
    # it for elements of an ideal.
    #
    # Otherwise, the monomial is of the form "a*idempotent*b" with "a" a path
    # of length "mid".
    int mid
    # In a sub-module of a direct sum, "pos" denotes the direct summand that
    # this monomial belongs to.
    unsigned int pos
    # paths are encoded as lists of integers. We store a*b if the monomial is
    # a*idempotent*b.
    biseq_t path
    # reference counter
    unsigned int ref

cdef struct path_term_t:
    path_mon_t *mon
    # We need to manually take care of the reference count for the
    # coefficient!
    PyObject *coef
    # In a polynomial, the terms are arranged in a pointered list
    path_term_t *nxt

# Type of monomial ordering functions.
ctypedef int (*path_order_t)(path_mon_t*, path_mon_t*)

# Polynomials are decreasingly sorted lists of terms. For convenience, the
# number of terms is directly available.
cdef struct path_poly_t:
    path_term_t *lead
    size_t nterms

# In path_poly_t, the terms need not to have all the same start and end
# points. path_homog_poly_t provides a list of "start and end point
# homogeneous polynomials". They are sorted by increasing (start, end)
cdef struct path_homog_poly_t:
    path_poly_t *poly
    int start, end
    path_homog_poly_t *nxt

cdef class PathAlgebraElement(RingElement):
    cdef path_homog_poly_t *data
    cdef path_order_t cmp_terms
    cdef long _hash
    cpdef ssize_t degree(self) except -2
    cpdef dict monomial_coefficients(self)
    cpdef list coefficients(self)
    cpdef list monomials(self)
    cpdef list support(self)
    cpdef list terms(self)
    cpdef object coefficient(self, QuiverPath P)
    cdef list _sorted_items_for_printing(self)
    cdef inline PathAlgebraElement _new_(self, path_homog_poly_t *h)
