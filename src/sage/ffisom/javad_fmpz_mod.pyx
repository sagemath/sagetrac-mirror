include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.rings.integer cimport Integer
from sage.rings.finite_rings.finite_field_flint_fq cimport FiniteField_flint_fq
from sage.rings.finite_rings.element_flint_fq cimport FiniteFieldElement_flint_fq
from sage.libs.flint.fmpz cimport *
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
    del isom
    return res

def change_basis_javad_fmpz_mod(FiniteFieldElement_flint_fq x, FiniteFieldElement_flint_fq u):
    cdef FFIsomBaseChange *bc
    cdef fmpz_mod_poly_t x_image
    cdef mpz_t cgmp
    cdef fmpz_t cflint
    cdef Integer cint
    cdef list l
    bc = new FFIsomBaseChange()
    fmpz_mod_poly_init(x_image, fq_ctx_prime(x._cparent))
    sig_on()
    bc.change_basis(<fmpz_mod_poly_struct*> (x_image), <fmpz_mod_poly_struct*> (u.val), <fmpz_mod_poly_struct*> (x.val), <fmpz_mod_poly_struct*> (x._cparent.modulus))
    sig_off()
    n = x.parent().degree()
    l = []
    fmpz_init(cflint)
    cint = Integer.__new__(Integer)
    for i in xrange(n):
        fmpz_mod_poly_get_coeff_fmpz(cflint, x_image, i)
        flint_mpz_init_set_readonly(cgmp, cflint)
        cint.set_from_mpz(cgmp)
        flint_mpz_clear_readonly(cgmp)
        l.append(Integer(cint))
    fmpz_clear(cflint)
    fmpz_mod_poly_clear(x_image)
    del bc
    return tuple(l)

