from sage.structure.element cimport CommutativeAlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.polydict cimport ETuple, PolyDict
from sage.rings.polynomial.multi_polynomial cimport MPolynomial


cdef class LaurentPolynomial(CommutativeAlgebraElement):
    cdef LaurentPolynomial _new_c(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _floordiv_(self, other)
    cpdef long number_of_terms(self) except -1
    cpdef dict dict(self)

cdef class LaurentPolynomial_univariate(LaurentPolynomial):
    cpdef ModuleElement __u
    cdef long __n
    cpdef __normalize(self)
    cpdef _unsafe_mutate(self, i, value)

cdef class LaurentPolynomial_mpair(LaurentPolynomial):
    cdef ETuple _mon
    cdef MPolynomial _poly
    cdef PolyDict _prod
    cdef _compute_polydict(self)
    cdef _normalize(self, i=*)

