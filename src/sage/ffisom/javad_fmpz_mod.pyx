include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.finite_rings.finite_field_flint_fq cimport FiniteField_flint_fq

cdef extern from "flint/fq.h":
    ctypedef struct fq_ctx_struct:
        fmpz_mod_poly_t modulus

def isom_javad_fmpz_mod(k1, k2):
    cdef FFIsomorphism *isom

    if not isinstance(k1, FiniteField_flint_fq) \
       or not isinstance(k1, FiniteField_flint_fq):
        raise TypeError

    sig_on() 
    isom = new FFIsomorphism((<fq_ctx_struct*>(<FiniteField_flint_fq>k1)._ctx).modulus, (<fq_ctx_struct*>(<FiniteField_flint_fq>k2)._ctx).modulus)
    sig_off()

    # Do nothing, for now timings only matter
    return 
