from bounded_series_ring_element cimport BoundedSeries
from sage.rings.polynomial.polynomial_element cimport Polynomial

from sage.structure.element cimport RingElement
from sage.structure.parent cimport Parent

cdef class BoundedSeries_poly(BoundedSeries):
    cdef Polynomial _f

    cdef _valuation
    cdef int _degree
    cdef _mu
    cdef BoundedSeries_poly _a
    cdef BoundedSeries_poly _b
    cdef BoundedSeries_poly _binv

    cpdef BoundedSeries _new_constant_series(self, RingElement a, Parent P, char check=*)

    cdef _compute_valuation_degree(self)
