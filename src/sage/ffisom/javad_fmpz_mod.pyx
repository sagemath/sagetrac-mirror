include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.finite_rings.finite_field_flint_fq cimport FiniteField_flint_fq
from sage.rings.finite_rings.element_flint_fq cimport FiniteFieldElement_flint_fq
from sage.libs.flint.fq cimport *

def isom_javad_fmpz_mod(k1, k2):
    cdef FFIsomorphism *isom
    cdef FiniteFieldElement_flint_fq res

    if not isinstance(k1, FiniteField_flint_fq) \
       or not isinstance(k1, FiniteField_flint_fq):
        raise TypeError

    sig_on() 
    isom = new FFIsomorphism((<fq_ctx_struct*>(<FiniteField_flint_fq>k1)._ctx).modulus, (<fq_ctx_struct*>(<FiniteField_flint_fq>k2)._ctx).modulus)
    sig_off()

    res = k2(0)
    res.set_from_fq(<fq_struct *>(isom.x_image))
    return res
