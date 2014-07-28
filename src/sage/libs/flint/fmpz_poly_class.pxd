include "sage/ext/stdsage.pxi"

from flint cimport *
from sage.libs.flint.fmpz_poly cimport *

from sage.structure.sage_object cimport SageObject

cdef class Fmpz_poly(SageObject):
    cdef fmpz_poly_t poly

