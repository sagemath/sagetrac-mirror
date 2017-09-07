from sage.libs.flint.types cimport nmod_poly_t

from sage.rings.morphism cimport RingHomomorphism
from sage.categories.morphism cimport Morphism
from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.categories.map cimport Section

cdef class FpTElement(RingElement):
    cdef nmod_poly_t _numer, _denom
    cdef bint initalized
    cdef long p

    cdef FpTElement _new_c(self)
    cpdef _add_(self, other)
    cpdef _mul_(self, other)
    cdef FpTElement _copy_c(self)
    cpdef numerator(self)
    cpdef denominator(self)
    cpdef FpTElement next(self)
    cpdef _sqrt_or_None(self)
    cpdef bint is_square(self)

cdef class FpT_iter:
    cdef parent
    cdef long degree
    cdef FpTElement cur
    cdef nmod_poly_t g
