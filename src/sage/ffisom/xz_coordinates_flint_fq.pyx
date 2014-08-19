include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fq cimport *
from sage.rings.integer cimport Integer
from sage.rings.finite_rings.finite_field_flint_fq cimport FiniteField_flint_fq
from sage.rings.finite_rings.element_flint_fq cimport FiniteFieldElement_flint_fq

cpdef ladder(tuple P, Integer m, FiniteFieldElement_flint_fq a, FiniteFieldElement_flint_fq b):
    cdef fq_weierstrass_xz_t E
    cdef fq_ctx_struct *K
    cdef FiniteFieldElement_flint_fq x, z
    cdef fmpz_t m_fmpz 

    K = a._cparent
    x = a._new()
    z = a._new()
    fq_weierstrass_xz_init(E, K)
    fq_weierstrass_xz_set_fq(E, a.val, b.val, K)
    fmpz_init_set_readonly(m_fmpz, m.value)

    sig_on()
    fq_weierstrass_xz_mul_ltr(x.val, z.val, (<FiniteFieldElement_flint_fq>(P[0])).val, (<FiniteFieldElement_flint_fq>(P[1])).val, m_fmpz, E, K)
    sig_off()

    fmpz_clear_readonly(m_fmpz)
    fq_weierstrass_xz_clear(E, K)

    return (x, z)

cpdef find_ordm(object E, object m):
    cdef FiniteField_flint_fq K
    K = <FiniteField_flint_fq>(E.base_ring())
    cofactor = E.cardinality()//m
    coprime = m.prime_divisors()

    while True:
        P = ladder((E.random_point()[0], K(1)), cofactor, E.a4(), E.a6())
        if P[1] == 0:
            continue
        for a in coprime:
            m_a = m//a
            if ladder((P[0], K(1)), m_a, E.a4(), E.a6())[1] == 0:
                break
        else:
            return P

