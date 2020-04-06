from sage.rings.polynomial.skew_polynomial_finite_order cimport SkewPolynomial_finite_order_dense
from sage.rings.polynomial.skew_polynomial_element cimport CenterSkewPolynomial_generic_dense
from sage.matrix.matrix_dense cimport Matrix_dense

cdef class SkewPolynomial_finite_field_dense (SkewPolynomial_finite_order_dense):
    cdef _norm_factor
    cdef dict _rdivisors
    cdef dict _types
    cdef _factorization

    cdef SkewPolynomial_finite_field_dense _rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rrem(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_rfloordiv(self, SkewPolynomial_finite_field_dense other)
    cdef void _inplace_lmonic(self)
    cdef void _inplace_rmonic(self)
    cdef void _inplace_rgcd(self,SkewPolynomial_finite_field_dense other)
    cdef Py_ssize_t _val_inplace_unit(self)
    cdef SkewPolynomial_finite_field_dense _rquo_inplace_rem(self, SkewPolynomial_finite_field_dense other)

    # Finding divisors
    cdef SkewPolynomial_finite_field_dense _rdivisor_c(P, CenterSkewPolynomial_generic_dense N)

    # Finding factorizations
    cdef _factor_c(self)
    cdef _factor_uniform_c(self)
