from sage.structure.element cimport CommutativeAlgebraElement, ModuleElement, RingElement, Element
from sage.rings.polynomial.polydict cimport ETuple, PolyDict
from sage.rings.polynomial.multi_polynomial cimport MPolynomial


cdef class PuiseuxPolynomial(CommutativeAlgebraElement):
    cdef PuiseuxPolynomial _new_c(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cpdef _floordiv_(self, other)
    cpdef long number_of_terms(self) except -1
    cpdef dict dict(self)

cdef class PuiseuxPolynomial_univariate(PuiseuxPolynomial):
    cpdef ModuleElement __u
    cdef long __n
    cpdef _unsafe_mutate(self, i, value)
