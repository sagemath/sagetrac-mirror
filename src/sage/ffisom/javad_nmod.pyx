include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.finite_rings.finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod
from sage.rings.finite_rings.element_flint_fq_nmod cimport FiniteFieldElement_flint_fq_nmod
from sage.libs.flint.fq_nmod cimport *

def isom_javad_nmod(k1, k2):
    cdef FFIsomorphism *isom
    cdef FiniteFieldElement_flint_fq_nmod res

    if not isinstance(k1, FiniteField_flint_fq_nmod) \
       or not isinstance(k1, FiniteField_flint_fq_nmod):
        raise TypeError

    sig_on()
    isom = new FFIsomorphism((<fq_nmod_ctx_struct*>(<FiniteField_flint_fq_nmod>k1)._ctx).modulus, (<fq_nmod_ctx_struct*>(<FiniteField_flint_fq_nmod>k2)._ctx).modulus)
    sig_off()

    res = k2(0)
    res.set_from_fq_nmod(<fq_nmod_struct *>(isom.x_image))    
    return res
