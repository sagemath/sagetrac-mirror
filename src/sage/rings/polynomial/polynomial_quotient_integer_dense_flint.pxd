from sage.libs.flint.fmpz_poly cimport fmpz_poly_t, fmpz_poly_struct
from sage.structure.parent cimport Parent
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint

ctypedef fmpz_poly_struct celement
ctypedef fmpz_poly_struct *cparent

include "polynomial_template_header.pxi"

cdef cparent get_cparent(parent) except? NULL

cdef class PolynomialQuotientRingElement_integer_flint(Polynomial_template):
    cdef Polynomial_template _new(self)
    cdef int _set_list(self, list x) except -1
    cdef int _set_fmpz_poly(self, fmpz_poly_t) except -1
