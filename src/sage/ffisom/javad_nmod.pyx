include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.finite_rings.finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod

cdef extern from "flint/fq_nmod.h":
    ctypedef struct fq_nmod_ctx_struct:
        nmod_poly_t modulus

def isom_javad_nmod(k1, k2):
    cdef FFIsomorphism *isom

    if not isinstance(k1, FiniteField_flint_fq_nmod) \
       or not isinstance(k1, FiniteField_flint_fq_nmod):
        raise TypeError

    sig_on()
    isom = new FFIsomorphism((<fq_nmod_ctx_struct*>(<FiniteField_flint_fq_nmod>k1)._ctx).modulus, (<fq_nmod_ctx_struct*>(<FiniteField_flint_fq_nmod>k2)._ctx).modulus)
    sig_off()

    # Do nothing, for now timings only matter
    return 
