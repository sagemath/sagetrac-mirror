from sage.libs.flint.types cimport nmod_poly_t

from sage.rings.morphism cimport RingHomomorphism
from sage.categories.morphism cimport Morphism
from sage.structure.element cimport Element, ModuleElement, FieldElement
from sage.categories.map cimport Section

cdef class FpTElement(FieldElement):
    cdef nmod_poly_t _numer, _denom
    cdef bint initialized
    cdef long p

    cdef FpTElement _new_c(self)
    cdef _add_(self, other)
    cdef _mul_(self, other)
    cdef FpTElement _copy_c(self)
    cdef numerator(self)
    cdef denominator(self)
    cdef FpTElement next(self)
    cdef _sqrt_or_None(self)
    cdef bint is_square(self)

cdef class FpT_iter:
    cdef parent
    cdef long degree
    cdef FpTElement cur
    cdef nmod_poly_t g
