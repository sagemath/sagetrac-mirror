include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.finite_rings.finite_field_flint_fq_nmod cimport FiniteField_flint_fq_nmod
from sage.rings.finite_rings.element_flint_fq_nmod cimport FiniteFieldElement_flint_fq_nmod
from sage.libs.flint.fq_nmod cimport *

def isom_javad_nmod(k1, k2):
    cdef FFIsomorphism *isom
    cdef FiniteFieldElement_flint_fq_nmod res
    cdef nmod_poly_t g1, g2

    if not isinstance(k1, FiniteField_flint_fq_nmod) \
       or not isinstance(k2, FiniteField_flint_fq_nmod):
        raise TypeError

    nmod_poly_init(g1, ((<FiniteField_flint_fq_nmod>k1)._ctx).modulus.mod.n)
    nmod_poly_init(g2, ((<FiniteField_flint_fq_nmod>k2)._ctx).modulus.mod.n)

    sig_on()
    isom = new FFIsomorphism((<FiniteField_flint_fq_nmod>k1)._ctx.modulus, (<FiniteField_flint_fq_nmod>k2)._ctx.modulus)
    isom.compute_generators(g1, g2)
    isom.build_isomorphism(g1, g2)
    sig_off()

    res = k2(0)
    res.set_from_fq_nmod(<fq_nmod_t>(isom.x_image))

    del isom
    nmod_poly_clear(g1)
    nmod_poly_clear(g2)

    return res

def change_basis_javad_nmod(FiniteFieldElement_flint_fq_nmod x, FiniteFieldElement_flint_fq_nmod u):
    cdef FFIsomBaseChange *bc
    cdef fq_nmod_t x_image
    cdef list l
    bc = new FFIsomBaseChange()
    fq_nmod_init(x_image, x._cparent)
    sig_on()
    bc.change_basis(<nmod_poly_struct*> (x_image), <nmod_poly_struct*> (u.val), <nmod_poly_struct*> (x.val), <nmod_poly_struct*> (x._cparent.modulus))
    sig_off()
    n = x.parent().degree()
    l = []
    for i in xrange(n):
        l.append(nmod_poly_get_coeff_ui(<nmod_poly_struct*> x_image, i))
    del bc
    fq_nmod_clear(x_image, x._cparent)
    return tuple(l)

