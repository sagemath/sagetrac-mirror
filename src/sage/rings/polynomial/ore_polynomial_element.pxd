from sage.structure.element cimport AlgebraElement
from sage.structure.parent cimport Parent
from sage.rings.morphism cimport Morphism
from sage.structure.element cimport RingElement
from sage.rings.polynomial.polynomial_element cimport Polynomial_generic_dense

cdef class OrePolynomial(AlgebraElement):
    cdef _is_gen

    cdef long _hash_c(self)
    cdef OrePolynomial _new_c(self, list coeffs, Parent P, char check=*)
    cdef OrePolynomial _new_constant_poly(self, RingElement a, Parent P, char check=*)
    cdef _neg_(self)
    cdef _floordiv_(self, right)
    cdef _mod_(self, right)

    cdef bint is_zero(self)
    cdef bint is_one(self)
 
    cdef _left_quo_rem(self, OrePolynomial other)
    cdef _right_quo_rem(self, OrePolynomial other)
    cdef OrePolynomial _left_lcm_cofactor(self, OrePolynomial other)
    cdef OrePolynomial _right_lcm_cofactor(self, OrePolynomial other)

    # Abstract methods
    cdef int degree(self)
    cdef list coefficients(self, sparse=*)


cdef void lmul_gen(list A, Morphism m, d)

cdef class OrePolynomial_generic_dense(OrePolynomial):
    cdef list _coeffs

    cdef void __normalize(self)
    cdef _add_(self, other)
    cdef list _mul_list(self, list A)
    cdef _mul_(self, other)

    cdef dict dict(self)
    cdef list list(self, bint copy=*)


cdef class OrePolynomialBaseringInjection(Morphism):
    cdef RingElement _an_element
    cdef object _new_constant_poly_
